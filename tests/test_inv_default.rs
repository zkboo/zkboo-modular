// SPDX-License-Identifier: LGPL-3.0-or-later

//! Validates the `MontgomeryMod::inv_exp` trait DEFAULT (`n - 2`, Fermat's little theorem for a
//! prime modulus). `PrimeMod` below deliberately does NOT override `inv_exp`, so `inv()` works out
//! of the box for prime fields without the caller supplying an exponent — the ergonomic fix.

mod common;

use crate::common::utils::modinv_i128;
use zkboo::{
    executor::{ExecutionBackend, OwnedFlexibleWordPool},
    word::CompositeWord,
};
use zkboo_modular::montgomery::{MontgomeryMod, MontgomeryWordRef};

type WP = OwnedFlexibleWordPool<usize>;

/// A prime modulus over `u32`. Crucially, it does NOT override [MontgomeryMod::inv_exp]: it relies
/// on the trait default of `n - 2`, which is the correct inverse exponent for a prime modulus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PrimeMod {
    n: u32,
}

impl MontgomeryMod<u32, 1> for PrimeMod {
    #[inline]
    fn n(&self) -> CompositeWord<u32, 1> {
        return self.n.into();
    }

    // NB: inv_exp is intentionally NOT implemented here — exercising the trait default (n - 2).

    fn rr_mod_n(&self) -> CompositeWord<u32, 1> {
        let n = self.n as u128;
        let r = 1u128 << 32;
        let r_mod_n = r % n;
        return (((r_mod_n * r_mod_n) % n) as u32).into();
    }

    #[inline]
    fn n_neg_inv(&self) -> CompositeWord<u32, 1> {
        let n = self.n as u128;
        let r = 1u128 << 32;
        let n_inv = modinv_i128(n as i128, r as i128).unwrap() as u128;
        assert!((n * n_inv) % r == 1, "n_inv is not correct");
        return (r.wrapping_sub(n_inv) as u32).into();
    }
}

/// `a * a^{-1} == 1 (mod p)`, computed through the circuit, using the DEFAULT inv_exp.
fn check_inv(p: u32, a: u32) {
    let m = PrimeMod { n: p };
    let executor = ExecutionBackend::<WP>::new().into_executor();
    let aw = MontgomeryWordRef::new(executor.input(a), m);
    let inv = aw.clone().inv();
    let prod = aw * inv;
    executor.output(prod.value());
    let out: u32 = executor.finalize().as_vec()[0];
    assert_eq!(out, 1, "a * a^-1 should be 1 mod {p} (a = {a})");
}

#[test]
fn test_default_inv_exp_value() {
    // The default inv_exp must equal n - 2.
    for &p in &[97u32, 65537, 2147483647] {
        let m = PrimeMod { n: p };
        let exp: u32 = m.inv_exp().expect("default inv_exp must be Some").into();
        assert_eq!(exp, p - 2, "default inv_exp should be n - 2");
    }
}

#[test]
fn test_inv_small_prime() {
    let p = 97u32;
    for a in 1..p {
        check_inv(p, a);
    }
}

#[test]
fn test_inv_fermat_prime() {
    let p = 65537u32; // F4
    for a in [1, 2, 3, 100, 256, 65535, 65536] {
        check_inv(p, a);
    }
}

#[test]
fn test_inv_mersenne_prime() {
    let p = 2147483647u32; // M31 = 2^31 - 1
    for a in [1u32, 2, 3, 7, 12345, 1 << 20, p - 1, p / 2] {
        check_inv(p, a);
    }
}
