// SPDX-License-Identifier: LGPL-3.0-or-later

//! Edge-case validation for the constant-time safegcd modular inverse.
//!
//! The random-sample inverse test (`test_secp256k1::test_inv`) never lands on the structured edge
//! values where safegcd's divstep/coefficient bookkeeping is most delicate. Here we pin those
//! exactly at the 256-bit secp256k1 width — `1`, `2`, `p-1` (≡ -1), `p-2`, a tiny value, and `0`
//! (non-invertible) — by checking the defining property `a · a⁻¹ ≡ 1 (mod p)` directly in-circuit
//! (and `a⁻¹ = 0` for `a = 0`, matching the Fermat convention), so no reference inverse is needed.

use zkboo::{
    backend::{Backend, Frontend},
    circuit::Circuit,
    executor::{exec, OwnedFlexibleWordPool},
    word::{CompositeWord, Words},
};
use zkboo_modular::montgomery::{MontgomeryFrontendIO, MontgomeryMod, MontgomeryWordRef};

type WP = OwnedFlexibleWordPool<usize>;

/// secp256k1 base field modulus, as two u128 limbs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Secp256k1;

impl MontgomeryMod<u128, 2> for Secp256k1 {
    fn n(&self) -> CompositeWord<u128, 2> {
        CompositeWord::from_be_words([
            0xffffffffffffffffffffffffffffffff,
            0xfffffffffffffffffffffffefffffc2f,
        ])
    }
    fn rr_mod_n(&self) -> CompositeWord<u128, 2> {
        CompositeWord::from_be_words([0x0, 0x1000007a2000e90a1])
    }
    fn n_neg_inv(&self) -> CompositeWord<u128, 2> {
        CompositeWord::from_be_words([
            0xc9bd1905155383999c46c2c295f2b761,
            0xbcb223fedc24a059d838091dd2253531,
        ])
    }
}

/// Circuit: input `a`, output `a · a⁻¹` in canonical form.
struct InvCheck {
    a: CompositeWord<u128, 2>,
}

impl Circuit for InvCheck {
    fn exec<B: Backend>(&self, fe: &Frontend<B>) {
        let a = MontgomeryWordRef::new(fe.input(self.a), Secp256k1);
        let prod = a.clone() * a.inv();
        fe.montgomery_output(prod);
    }
}

fn check(a: CompositeWord<u128, 2>, expect_one: bool) {
    let outputs = exec::<_, WP>(&InvCheck { a });
    let mut expected = Words::new();
    let want = if expect_one {
        CompositeWord::<u128, 2>::ONE
    } else {
        CompositeWord::<u128, 2>::ZERO
    };
    expected.as_vec_mut::<u128>().extend(want.to_le_words());
    assert_eq!(outputs, expected, "a·a⁻¹ wrong for a={a:?}");
}

#[test]
fn safegcd_edges_secp256k1() {
    let p = Secp256k1.n();
    let one = CompositeWord::<u128, 2>::ONE;
    let two = CompositeWord::<u128, 2>::from_le_words([2, 0]);
    // invertible edges: a · a⁻¹ == 1
    check(one, true);
    check(two, true);
    check(p.wrapping_sub(one), true); // p - 1 ≡ -1
    check(p.wrapping_sub(two), true); // p - 2
    check(CompositeWord::from_le_words([0xdeadbeef, 0]), true);
    check(CompositeWord::from_le_words([3, 0]), true);
    check(CompositeWord::from_le_words([0, 1]), true); // 2^128
    // non-invertible: 0⁻¹ defined as 0, so 0·0⁻¹ == 0
    check(CompositeWord::<u128, 2>::ZERO, false);
}

#[test]
fn safegcd_random_secp256k1() {
    // Deterministic 256-bit pseudo-random spread, reduced below p, checked via a·a⁻¹ ≡ 1.
    let p = Secp256k1.n();
    let mut lo: u128 = 0x9e3779b97f4a7c15f39cc0605cedc834;
    let mut hi: u128 = 0xc2b2ae3d27d4eb4f165667b19e3779b9;
    for _ in 0..32 {
        // xorshift-ish mixing
        lo = lo.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        hi = hi.wrapping_mul(6364136223846793005).wrapping_add(lo);
        let mut a = CompositeWord::<u128, 2>::from_le_words([lo, hi]);
        // ensure a < p and a != 0
        if a.ge(p) {
            a = a.wrapping_sub(p);
        }
        if a.is_zero() {
            a = CompositeWord::from_le_words([1, 0]);
        }
        check(a, true);
    }
}
