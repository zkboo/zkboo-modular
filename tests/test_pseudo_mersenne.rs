// SPDX-License-Identifier: LGPL-3.0-or-later

//! Validates pseudo-Mersenne (Solinas) reduction and field multiplication against an arbitrary-
//! precision reference, for the secp256k1 modulus p = 2^256 - (2^32 + 977).

use dashu_int::UBig;
use zkboo::{
    backend::{Backend, Frontend},
    circuit::Circuit,
    executor::{exec, OwnedFlexibleWordPool},
    word::{CompositeWord, Words},
};
use zkboo_modular::pseudo_mersenne::{self, PseudoMersenneMod};
use zkboo::executor::ExecOptions;

type WP = OwnedFlexibleWordPool<usize>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Secp256k1;

impl PseudoMersenneMod<u64, 4> for Secp256k1 {
    fn p(&self) -> CompositeWord<u64, 4> {
        CompositeWord::from_be_words([
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xfffffffefffffc2f,
        ])
    }
    fn c(&self) -> CompositeWord<u64, 4> {
        // 2^256 - p = 2^32 + 977
        CompositeWord::from_be_words([0, 0, 0, 0x00000001_000003d1])
    }
}

// Verify the concrete FieldRep impl coheres alongside the MontgomeryMod blanket impl.
zkboo_modular::impl_pseudo_mersenne_field_rep!(Secp256k1, u64, 4);

fn to_ubig(w: CompositeWord<u64, 4>) -> UBig {
    let mut v = UBig::ZERO;
    for &limb in w.to_le_words().iter().rev() {
        v = (v << 64) + UBig::from(limb);
    }
    v
}
fn p_ubig() -> UBig {
    to_ubig(Secp256k1.p())
}

/// Circuit: `a · b mod p` via pseudo-Mersenne reduction.
struct PmMul {
    a: CompositeWord<u64, 4>,
    b: CompositeWord<u64, 4>,
}
impl Circuit for PmMul {
    fn exec<B: Backend>(&self, fe: &Frontend<B>) {
        let a = fe.input(self.a);
        let b = fe.input(self.b);
        fe.output(pseudo_mersenne::mul(&Secp256k1, a, b));
    }
}

fn check_mul(a: CompositeWord<u64, 4>, b: CompositeWord<u64, 4>) {
    let outputs = exec::<_, WP, _>(&PmMul { a, b }, ExecOptions::new());
    let want = (to_ubig(a) * to_ubig(b)) % p_ubig();
    // Encode the expected residue back into 4 little-endian u64 words.
    let mut words = [0u64; 4];
    let mask = UBig::from(u64::MAX);
    let mut t = want.clone();
    for w in words.iter_mut() {
        *w = u64::try_from(&t & &mask).unwrap();
        t >>= 64;
    }
    let mut expected = Words::new();
    expected
        .as_vec_mut::<u64>()
        .extend(CompositeWord::<u64, 4>::from_le_words(words).to_le_words());
    assert_eq!(outputs, expected, "pseudo-Mersenne mul wrong for a={a:?} b={b:?}");
}

#[test]
fn pseudo_mersenne_mul_edges() {
    let p = Secp256k1.p();
    let one = CompositeWord::<u64, 4>::ONE;
    let zero = CompositeWord::<u64, 4>::ZERO;
    let pm1 = p.wrapping_sub(one);
    check_mul(zero, pm1);
    check_mul(one, pm1);
    check_mul(pm1, pm1); // (p-1)^2 — stresses the high half of the product
    check_mul(p.wrapping_sub(CompositeWord::from_le_words([2, 0, 0, 0])), pm1);
    check_mul(CompositeWord::from_le_words([0, 0, 0, 0xffffffffffffffff]), pm1);
    check_mul(CompositeWord::from_le_words([u64::MAX, u64::MAX, u64::MAX, u64::MAX]), one);
}

#[test]
fn pseudo_mersenne_mul_random() {
    // deterministic pseudo-random pairs reduced below p
    let p = p_ubig();
    let mut lo: u64 = 0x9e3779b97f4a7c15;
    let mut hi: u64 = 0xc2b2ae3d27d4eb4f;
    let mut next = || {
        lo = lo.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        hi = hi.wrapping_mul(6364136223846793005).wrapping_add(lo);
        let w = CompositeWord::<u64, 4>::from_le_words([lo, hi, lo ^ hi, hi.wrapping_add(lo)]);
        // reduce below p
        let r = to_ubig(w) % &p;
        let mut words = [0u64; 4];
        let mask = UBig::from(u64::MAX);
        let mut t = r;
        for x in words.iter_mut() {
            *x = u64::try_from(&t & &mask).unwrap();
            t >>= 64;
        }
        CompositeWord::<u64, 4>::from_le_words(words)
    };
    for _ in 0..64 {
        check_mul(next(), next());
    }
}
