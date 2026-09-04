// SPDX-License-Identifier: LGPL-3.0-or-later

//! The advised modular inverse, against the in-circuit one it stands in for.
//!
//! [`MontgomeryWordRef::inv_advised`] does not compute an inverse; it takes the one the prover
//! computed on the host and asserts that it is one. So there are two things to check: that an
//! honest prover's answer agrees with [`MontgomeryWordRef::inv`] everywhere, including at zero;
//! and that a dishonest one is caught. The second is checked by computing the advice with a field
//! whose native inversion lies — which is exactly where a dishonest prover's lie would live, since
//! the advice is an ordinary input the prover fills in.

use zkboo::{
    backend::{Backend, Frontend, WordRef},
    circuit::{Assertions, Circuit},
    executor::{OwnedFlexibleWordPool, exec},
    word::{CompositeWord, Words},
};
use zkboo_modular::field::FieldRep;
use zkboo_modular::montgomery::{MontgomeryFrontendIO, MontgomeryMod, MontgomeryWordRef};
use zkboo_profiling::profile;
use zkboo::executor::ExecOptions;

type WP = OwnedFlexibleWordPool<usize>;
type Word2 = CompositeWord<u128, 2>;

/// secp256k1 base field modulus, as two u128 limbs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Secp256k1;

impl MontgomeryMod<u128, 2> for Secp256k1 {
    fn n(&self) -> Word2 {
        CompositeWord::from_be_words([
            0xffffffffffffffffffffffffffffffff,
            0xfffffffffffffffffffffffefffffc2f,
        ])
    }
    fn rr_mod_n(&self) -> Word2 {
        CompositeWord::from_be_words([0x0, 0x1000007a2000e90a1])
    }
    fn n_neg_inv(&self) -> Word2 {
        CompositeWord::from_be_words([
            0xc9bd1905155383999c46c2c295f2b761,
            0xbcb223fedc24a059d838091dd2253531,
        ])
    }
}

/// The same field, with a native inversion that lies — a prover computing bad advice.
///
/// It cannot be a [`MontgomeryMod`], since those get their [`FieldRep`] from a blanket
/// implementation and so cannot override the inversion; it delegates to one instead.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Lying {
    /// Whether to return zero rather than the true inverse of a nonzero argument.
    lie_about_nonzero: bool,
    /// Whether to return one rather than zero for the inverse of zero.
    lie_about_zero: bool,
}

impl FieldRep<u128, 2> for Lying {
    fn modulus(&self) -> Word2 {
        return FieldRep::modulus(&Secp256k1);
    }
    fn encode_const(&self, value: Word2) -> Word2 {
        return Secp256k1.encode_const(value);
    }
    fn decode_const(&self, internal: Word2) -> Word2 {
        return Secp256k1.decode_const(internal);
    }
    fn mul_reduce_const(&self, lo: Word2, hi: Word2) -> Word2 {
        return Secp256k1.mul_reduce_const(lo, hi);
    }
    fn encode<B: Backend>(&self, value: WordRef<B, u128, 2>) -> WordRef<B, u128, 2> {
        return Secp256k1.encode(value);
    }
    fn decode<B: Backend>(&self, internal: WordRef<B, u128, 2>) -> WordRef<B, u128, 2> {
        return Secp256k1.decode(internal);
    }
    fn mul_reduce<B: Backend>(
        &self,
        lo: WordRef<B, u128, 2>,
        hi: WordRef<B, u128, 2>,
    ) -> WordRef<B, u128, 2> {
        return Secp256k1.mul_reduce(lo, hi);
    }
    fn invert_const(&self, internal: Word2) -> Word2 {
        if internal.is_nonzero() {
            if self.lie_about_nonzero {
                // Zero is the inverse of nothing, so this is a lie whatever the argument was.
                return Word2::ZERO;
            }
            return Secp256k1.invert_const(internal);
        }
        if self.lie_about_zero {
            // Zero has no inverse, so any nonzero answer is a lie — and the product assertion
            // cannot see it.
            return self.encode_const(Word2::ONE);
        }
        return Word2::ZERO;
    }
}

/// The prover's host-side work: the advice it will pass in, in internal representation.
fn advice_for<F: FieldRep<u128, 2>>(field: F, a: Word2) -> Word2 {
    return field.invert_const(field.encode_const(a));
}

/// Outputs the advised inverse and the in-circuit one for the same value, in that order.
struct BothInverses {
    a: Word2,
    advice: Word2,
}

impl Circuit for BothInverses {
    fn exec<B: Backend>(&self, fe: &Frontend<B>) {
        let mut asserts = Assertions::new();
        let a = MontgomeryWordRef::new(fe.input(self.a), Secp256k1);
        let advice = MontgomeryWordRef::from_inner(fe.input(self.advice), Secp256k1);
        fe.montgomery_output(a.clone().inv_advised(advice, &mut asserts));
        fe.montgomery_output(a.inv());
        asserts.output(fe);
    }
}

/// Inverts a value over a field whose advice may be a lie, outputting only the assertion flag.
struct AdvisedInverse {
    a: Word2,
    advice: Word2,
    field: Lying,
}

impl AdvisedInverse {
    fn new(a: Word2, field: Lying) -> Self {
        return Self {
            a,
            advice: advice_for(field, a),
            field,
        };
    }
}

impl Circuit for AdvisedInverse {
    fn exec<B: Backend>(&self, fe: &Frontend<B>) {
        let mut asserts = Assertions::new();
        let a = MontgomeryWordRef::new(fe.input(self.a), self.field);
        let advice = MontgomeryWordRef::from_inner(fe.input(self.advice), self.field);
        let _ = a.inv_advised(advice, &mut asserts);
        asserts.output(fe);
    }
}

/// The structured values where inversion bookkeeping is most delicate, plus zero.
fn probe_values() -> Vec<Word2> {
    let p = Secp256k1.n();
    let one = Word2::ONE;
    let two = CompositeWord::from_le_words([2, 0]);
    return vec![
        Word2::ZERO,
        one,
        two,
        CompositeWord::from_le_words([7, 0]),
        p.wrapping_sub(one),
        p.wrapping_sub(two),
        CompositeWord::from_le_words([0, 1]),
    ];
}

fn flag(words: &Words) -> u8 {
    return *words
        .u8
        .last()
        .expect("the circuit outputs its assertion flag");
}

#[test]
fn the_advised_inverse_agrees_with_the_in_circuit_one() {
    for a in probe_values() {
        let out = exec::<_, WP, _>(&BothInverses {
            a,
            advice: advice_for(Secp256k1, a),
        }, ExecOptions::new());
        assert_eq!(flag(&out), 1, "assertions failed for a={a:?}");
        let limbs = &out.u128;
        assert_eq!(limbs.len(), 4, "two inverses of two limbs each");
        assert_eq!(
            limbs[..2],
            limbs[2..],
            "the advised inverse differs from the in-circuit one for a={a:?}"
        );
    }
}

#[test]
fn a_lie_about_an_invertible_value_is_caught() {
    let field = Lying {
        lie_about_nonzero: true,
        lie_about_zero: false,
    };
    for a in probe_values().into_iter().filter(|a| a.is_nonzero()) {
        let out = exec::<_, WP, _>(&AdvisedInverse::new(a, field), ExecOptions::new());
        assert_eq!(
            flag(&out),
            0,
            "a wrong inverse of a={a:?} passed the assertions"
        );
    }
}

#[test]
fn a_lie_about_the_inverse_of_zero_is_caught() {
    // Zero has no inverse, so the product assertion is vacuous there and says nothing about the
    // advice. Without the second assertion a prover could return anything, and everything computed
    // from it — the affine coordinates of a point at infinity, say — would be its to choose.
    let field = Lying {
        lie_about_nonzero: false,
        lie_about_zero: true,
    };
    let out = exec::<_, WP, _>(&AdvisedInverse::new(Word2::ZERO, field), ExecOptions::new());
    assert_eq!(flag(&out), 0, "a nonzero inverse of zero passed the assertions");
}

#[test]
fn an_honest_prover_over_the_delegating_field_is_accepted() {
    // The lying field only lies when told to: this pins that the two tests above fail for the
    // reason claimed and not because the delegation is broken.
    let field = Lying {
        lie_about_nonzero: false,
        lie_about_zero: false,
    };
    for a in probe_values() {
        let out = exec::<_, WP, _>(&AdvisedInverse::new(a, field), ExecOptions::new());
        assert_eq!(flag(&out), 1, "honest advice rejected for a={a:?}");
    }
}

/// Outputs its input unchanged, with an empty assertion flag: the cost both circuits below pay
/// before either inverts anything.
struct NoInverse {
    a: Word2,
}

impl Circuit for NoInverse {
    fn exec<B: Backend>(&self, fe: &Frontend<B>) {
        let a = MontgomeryWordRef::new(fe.input(self.a), Secp256k1);
        fe.montgomery_output(a);
        Assertions::new().output(fe);
    }
}

/// Outputs the advised inverse alone.
struct AdvisedOnly {
    a: Word2,
    advice: Word2,
}

impl Circuit for AdvisedOnly {
    fn exec<B: Backend>(&self, fe: &Frontend<B>) {
        let mut asserts = Assertions::new();
        let a = MontgomeryWordRef::new(fe.input(self.a), Secp256k1);
        let advice = MontgomeryWordRef::from_inner(fe.input(self.advice), Secp256k1);
        fe.montgomery_output(a.inv_advised(advice, &mut asserts));
        asserts.output(fe);
    }
}

/// Outputs the in-circuit inverse alone.
struct InCircuitOnly {
    a: Word2,
}

impl Circuit for InCircuitOnly {
    fn exec<B: Backend>(&self, fe: &Frontend<B>) {
        let a = MontgomeryWordRef::new(fe.input(self.a), Secp256k1);
        fe.montgomery_output(a.inv());
        Assertions::new().output(fe);
    }
}

#[test]
fn asking_for_an_inverse_is_far_cheaper_than_computing_one() {
    let a = Word2::ONE;
    let baseline = profile(&NoInverse { a }).and_msg_size().sum();
    let advised = profile(&AdvisedOnly {
        a,
        advice: advice_for(Secp256k1, a),
    })
    .and_msg_size()
    .sum()
        - baseline;
    let in_circuit = profile(&InCircuitOnly { a }).and_msg_size().sum() - baseline;
    // One field multiplication, two zero tests, a two-way constant select and a comparison,
    // against a full run of divsteps. These are the same figures the hint-gate form cost: advice
    // as an ordinary input is byte-for-byte as cheap, and the circuit-side accumulator conjoins
    // the two assertions at the same single AND the backend one did.
    assert_eq!(advised, 1_515, "the advised inverse changed cost");
    assert_eq!(in_circuit, 60_649, "the in-circuit inverse changed cost");
    assert!(
        in_circuit / advised >= 40,
        "the advised inverse is not forty times cheaper: {advised} against {in_circuit}"
    );
}
