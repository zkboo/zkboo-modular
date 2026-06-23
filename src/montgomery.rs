// SPDX-License-Identifier: LGPL-3.0-or-later

//! Montgomery modular arithmetic for the [zkboo] crate.

use core::fmt::Debug;
use core::ops::{Add, Mul, Neg, Sub};
use zkboo::backend::{BooleanWordRef, Frontend};
use zkboo::{
    backend::{Backend, WordRef},
    word::{CompositeWord, Word, WordLike},
};

/// Trait for an odd modulus, implementing Montgomery reduction for efficient modular arithmetic.
pub trait MontgomeryMod<W: Word, const N: usize>: Clone + Copy + Debug + PartialEq + Eq {
    /// The value of the modulus.
    ///
    /// This must be an odd number, i.e. its least significant bit must be 1.
    fn n(&self) -> CompositeWord<W, N>;

    /// Modular inversion exponent. Optional, enables inversion by repeated squaring.
    ///
    /// The default assumes a **prime** modulus and returns `n - 2`: by Fermat's little theorem,
    /// `a^(n-2) ≡ a^(-1) (mod n)` for any `a` coprime to a prime `n`. This is correct for the
    /// prime fields used by the ECC curves (e.g. secp256k1). Composite moduli (where Fermat does
    /// not apply) must override this to return `None` — disabling [MontgomeryWord::inv] — or to a
    /// Carmichael-based exponent if one is known.
    #[inline]
    fn inv_exp(&self) -> Option<CompositeWord<W, N>> {
        let one = CompositeWord::<W, N>::ONE;
        return Some(self.n().wrapping_sub(one).wrapping_sub(one));
    }

    /// The residue by [Self::n] for the square of the Montgomery radix `r`.
    ///
    /// Specifically, working with arbitrary precision unsigned integers and letting the
    /// Montgomery radix be `r := 2**Self::W::WIDTH`, we have `rr_mod_n := r**2 % n`.
    fn rr_mod_n(&self) -> CompositeWord<W, N>;

    /// The negative of the multiplicative inverse of [Self::n] modulo the Montgomery radix `r`.
    ///
    /// Specifically, working with arbitrary precision unsigned integers and letting the Montgomery
    /// radix be `r := 2**Self::W::WIDTH`, we have `n * n_neg_inv % r = r - 1`.
    fn n_neg_inv(&self) -> CompositeWord<W, N>;

    /// Creates a [MontgomeryWord] representing the given value in Montgomery form.
    #[inline]
    fn const_word<U: WordLike<W, N>>(&self, value: U) -> MontgomeryWord<W, N, Self> {
        return MontgomeryWord::new(value, *self);
    }

    /// Creates a [MontgomeryWordRef] representing the zero value in Montgomery form.
    #[inline]
    fn zero_word(&self) -> MontgomeryWord<W, N, Self> {
        return MontgomeryWord::from_inner(CompositeWord::<W, N>::ZERO, *self);
    }

    /// Creates a [MontgomeryWordRef] representing the one value in Montgomery form.
    #[inline]
    fn one_word(&self) -> MontgomeryWord<W, N, Self> {
        return MontgomeryWord::new(CompositeWord::<W, N>::ONE, *self);
    }

    /// Performs Montgomery reduction on the given value.
    /// Returns the reduced value in Montgomery form.
    ///
    /// Equivalent to [Self::redc_wide] with `t_hi` set to `Self::W::ZERO`.
    fn redc<B: Backend>(&self, lo: WordRef<B, W, N>) -> WordRef<B, W, N> {
        let n = self.n();
        // Montgomery's `m = lo * n_neg_inv mod r` only needs the low `r`-residue, so a low-half
        // (wrapping) multiply suffices — the discarded high half of a full `wide_mul_const` is
        // pure wasted work.
        let m = lo.clone().wrapping_mul_const(self.n_neg_inv());
        let (mn_lo, mn_hi) = m.wide_mul_const(n);
        let (_, t_lo_carry) = lo.overflowing_add(mn_lo);
        let (mut t, t_hi_lo_carry) = mn_hi.overflowing_add(WordRef::from_bool(t_lo_carry));
        t = t_hi_lo_carry.select(t.clone() + n.wrapping_neg(), t);
        t = t.clone().ge_const(n).select(t.clone() - n, t);
        return t;
    }

    /// Performs Montgomery reduction on the given wide value, specified as two limbs `(lo, hi)`.
    /// Returns the reduced value in Montgomery form.
    ///
    /// ⚠️ Safety: The caller must ensure that the value represented by `(lo, hi)` is guaranteed to
    /// be less than the wide value `N.wide_mul(R)`, where `N` is the modulus and `R` is the
    /// Montgomery radix. Failure to do so may result in incorrect behaviour.
    fn redc_wide<B: Backend>(
        &self,
        lo: WordRef<B, W, N>,
        hi: WordRef<B, W, N>,
    ) -> WordRef<B, W, N> {
        let n = self.n();
        // Montgomery's `m = lo * n_neg_inv mod r` only needs the low `r`-residue, so a low-half
        // (wrapping) multiply suffices — the discarded high half of a full `wide_mul_const` is
        // pure wasted work.
        let m = lo.clone().wrapping_mul_const(self.n_neg_inv());
        let (mn_lo, mn_hi) = m.wide_mul_const(n);
        let (_, t_lo_carry) = lo.overflowing_add(mn_lo);
        let (t, t_hi_carry) = hi.overflowing_add(mn_hi);
        let (mut t, t_hi_lo_carry) = t.overflowing_add(WordRef::from_bool(t_lo_carry));
        t = t_hi_carry.select(t.clone() + n.wrapping_neg(), t);
        t = t_hi_lo_carry.select(t.clone() + n.wrapping_neg(), t);
        t = t.clone().ge_const(n).select(t.clone() - n, t);
        return t;
    }

    /// Converts a given value to Montgomery form.
    fn to_montgomery<B: Backend>(&self, value: WordRef<B, W, N>) -> WordRef<B, W, N> {
        let (lo, hi) = value.wide_mul_const(self.rr_mod_n());
        return self.redc_wide(lo, hi);
    }

    /// Converts a given value from Montgomery form to canonical form.
    ///
    /// ⚠️ Safety: The caller must ensure that the value is guaranteed to be in Montgomery form.
    /// Failure to do so may result in incorrect behaviour.
    fn from_montgomery<B: Backend>(&self, value: WordRef<B, W, N>) -> WordRef<B, W, N> {
        return self.redc(value);
    }

    /// Reduces the given value modulo the modulus.
    fn reduce<B: Backend>(&self, value: WordRef<B, W, N>) -> WordRef<B, W, N> {
        return self.from_montgomery(self.to_montgomery(value));
    }

    /// Performs Montgomery reduction on the given value.
    /// Returns the reduced value in Montgomery form.
    ///
    /// Equivalent to [Self::redc_wide] with `t_hi` set to `Self::W::ZERO`.
    fn redc_const(&self, lo: CompositeWord<W, N>) -> CompositeWord<W, N> {
        let n = self.n();
        let (m, _) = lo.wide_mul(self.n_neg_inv());
        let (mn_lo, mn_hi) = m.wide_mul(n);
        let (_, t_lo_carry) = lo.overflowing_add(mn_lo);
        let (mut t, t_hi_lo_carry) = mn_hi.overflowing_add(CompositeWord::from_bool(t_lo_carry));
        if t_hi_lo_carry {
            t = t.wrapping_add(n.wrapping_neg());
        }
        if t.ge(n) {
            t = t.wrapping_sub(n);
        }
        return t;
    }

    /// Performs Montgomery reduction on the given wide value, specified as two limbs `(lo, hi)`.
    /// Returns the reduced value in Montgomery form.
    ///
    /// ⚠️ Safety: The caller must ensure that the value represented by `(lo, hi)` is guaranteed to
    /// be less than the wide value `N.wide_mul(R)`, where `N` is the modulus and `R` is the
    /// Montgomery radix. Failure to do so may result in incorrect behaviour.
    fn redc_wide_const(
        &self,
        lo: CompositeWord<W, N>,
        hi: CompositeWord<W, N>,
    ) -> CompositeWord<W, N> {
        let n = self.n();
        let (m, _) = lo.clone().wide_mul(self.n_neg_inv());
        let (mn_lo, mn_hi) = m.wide_mul(n);
        let (_, t_lo_carry) = lo.overflowing_add(mn_lo);
        let (t, t_hi_carry) = hi.overflowing_add(mn_hi);
        let (mut t, t_hi_lo_carry) = t.overflowing_add(CompositeWord::from_bool(t_lo_carry));
        if t_hi_carry {
            t = t.wrapping_add(n.wrapping_neg());
        }
        if t_hi_lo_carry {
            t = t.wrapping_add(n.wrapping_neg());
        }
        if t.ge(n) {
            t = t.wrapping_sub(n);
        }
        return t;
    }

    /// Converts a given value to Montgomery form.
    fn to_montgomery_const(&self, value: CompositeWord<W, N>) -> CompositeWord<W, N> {
        let (lo, hi) = value.wide_mul(self.rr_mod_n());
        return self.redc_wide_const(lo, hi);
    }

    /// Converts a given value from Montgomery form to canonical form.
    ///
    /// ⚠️ Safety: The caller must ensure that the value is guaranteed to be in Montgomery form.
    /// Failure to do so may result in incorrect behaviour.
    fn from_montgomery_const(&self, value: CompositeWord<W, N>) -> CompositeWord<W, N> {
        return self.redc_const(value);
    }

    /// Reduces the given constant value modulo the modulus.
    fn reduce_const(&self, value: CompositeWord<W, N>) -> CompositeWord<W, N> {
        return self.from_montgomery_const(self.to_montgomery_const(value));
    }
}

/// Modular word value, with generic modulus.
/// It is a wrapper around a [WordRef] enforcing the guarantee that the value is in Montgomery form.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct MontgomeryWord<W: Word, const N: usize, M: MontgomeryMod<W, N>> {
    montgomery_val: CompositeWord<W, N>,
    modulus: M,
}

impl<W: Word, const N: usize, M: MontgomeryMod<W, N>> MontgomeryWord<W, N, M> {
    /// Converts the given value to Montgomery form.
    #[inline]
    pub fn new<U: WordLike<W, N>>(value: U, modulus: M) -> Self {
        return MontgomeryWord {
            montgomery_val: modulus.to_montgomery_const(value.to_word()),
            modulus,
        };
    }

    /// Wraps a value in Montgomery form, without converting it.
    ///
    /// ⚠️ Safety: The caller must ensure that the value is in Montgomery form.
    /// Failure to do so may result in incorrect behaviour.
    pub fn from_inner(montgomery_val: CompositeWord<W, N>, modulus: M) -> Self {
        return Self {
            montgomery_val,
            modulus,
        };
    }

    /// Returns the modulus associated with this Montgomery word.
    pub fn modulus(&self) -> M {
        return self.modulus;
    }

    /// Consumes this Montgomery word and returns its modulus.
    pub fn into_modulus(self) -> M {
        return self.modulus;
    }

    /// Converts this value out of Montgomery form.
    pub fn value(self) -> CompositeWord<W, N> {
        return self.modulus.from_montgomery_const(self.montgomery_val);
    }

    /// Returns the inner Montgomery value.
    pub fn inner(&self) -> &CompositeWord<W, N> {
        return &self.montgomery_val;
    }

    /// Consumes this Montgomery word and returns the inner Montgomery value.
    pub fn into_inner(self) -> CompositeWord<W, N> {
        return self.montgomery_val;
    }

    /// Consumes this Montgomery word and returns its inner Montgomery value and modulus.
    pub fn destructure(self) -> (CompositeWord<W, N>, M) {
        return (self.montgomery_val, self.modulus);
    }

    /// Returns true if this Montgomery word is zero, false otherwise.
    pub fn is_zero(self) -> bool {
        return self.montgomery_val.is_zero();
    }

    /// Returns true if this Montgomery word is nonzero, false otherwise.
    pub fn is_nonzero(self) -> bool {
        return self.montgomery_val.is_nonzero();
    }

    /// Modular inverse of this Montgomery word.
    pub fn inv(self) -> Self {
        let e = self.modulus.inv_exp();
        return match e {
            Some(e) => self.inv_by_rep_squaring(e),
            None => unimplemented!("Cannot yet compute inverse without explicit inverse exponent."),
        };
    }

    /// Helper function implementing modular inversion by repeated squaring.
    fn inv_by_rep_squaring(self, mut e: CompositeWord<W, N>) -> Self {
        let modulus = self.modulus;
        let mut res = MontgomeryWord {
            montgomery_val: modulus.to_montgomery_const(CompositeWord::<W, N>::ONE),
            modulus: modulus,
        };
        let mut base = self;
        while e.is_nonzero() {
            if e.lsb() {
                res = res * base.clone();
            }
            base = base.clone() * base;
            e = e >> 1;
        }
        return res;
    }
}

impl<W: Word, const N: usize, M: MontgomeryMod<W, N>> Neg for MontgomeryWord<W, N, M> {
    type Output = Self;

    /// Modular negation.
    fn neg(self) -> Self::Output {
        if self.montgomery_val == CompositeWord::<W, N>::ZERO {
            return self;
        } else {
            return MontgomeryWord {
                montgomery_val: self.modulus.n().wrapping_sub(self.montgomery_val),
                modulus: self.modulus,
            };
        }
    }
}

impl<W: Word, const N: usize, M: MontgomeryMod<W, N>> Add<Self> for MontgomeryWord<W, N, M> {
    type Output = Self;

    /// Modular addition.
    fn add(self, rhs: Self) -> Self::Output {
        if self.modulus != rhs.modulus {
            panic!("Cannot add modular words with different moduli");
        }
        let (res, carry) = self.montgomery_val.overflowing_add(rhs.montgomery_val);
        if carry || res.ge(self.modulus.n()) {
            return MontgomeryWord {
                montgomery_val: res.wrapping_sub(self.modulus.n()),
                modulus: self.modulus,
            };
        } else {
            return MontgomeryWord {
                montgomery_val: res,
                modulus: self.modulus,
            };
        }
    }
}

impl<W: Word, const N: usize, M: MontgomeryMod<W, N>> Sub<Self> for MontgomeryWord<W, N, M> {
    type Output = Self;
    /// Modular subtraction.
    fn sub(self, rhs: Self) -> Self::Output {
        if self.modulus != rhs.modulus {
            panic!("Cannot subtract modular words with different moduli");
        }
        let (res, borrow) = self.montgomery_val.overflowing_sub(rhs.montgomery_val);
        if borrow {
            return MontgomeryWord {
                montgomery_val: res.wrapping_add(self.modulus.n()),
                modulus: self.modulus,
            };
        } else {
            return MontgomeryWord {
                montgomery_val: res,
                modulus: self.modulus,
            };
        }
    }
}

impl<W: Word, const N: usize, M: MontgomeryMod<W, N>> Mul<Self> for MontgomeryWord<W, N, M> {
    type Output = Self;
    /// Modular multiplication.
    fn mul(self, rhs: Self) -> Self::Output {
        if self.modulus != rhs.modulus {
            panic!("Cannot multiply modular words with different moduli");
        }
        let (res_lo, res_hi) = self.montgomery_val.wide_mul(rhs.montgomery_val);
        return MontgomeryWord {
            montgomery_val: self.modulus.redc_wide_const(res_lo, res_hi),
            modulus: self.modulus,
        };
    }
}

/// Word reference to a modular value, with generic modulus.
/// It is a wrapper around a [WordRef] enforcing the guarantee that the value is in Montgomery form.
#[derive(Debug)]
pub struct MontgomeryWordRef<B: Backend, W: Word, const N: usize, M: MontgomeryMod<W, N>> {
    montgomery_val: WordRef<B, W, N>,
    modulus: M,
}

impl<B: Backend, M: MontgomeryMod<W, N>, W: Word, const N: usize> Clone
    for MontgomeryWordRef<B, W, N, M>
{
    fn clone(&self) -> Self {
        return Self {
            montgomery_val: self.montgomery_val.clone(),
            modulus: self.modulus,
        };
    }
}

impl<B: Backend, M: MontgomeryMod<W, N>, W: Word, const N: usize> MontgomeryWordRef<B, W, N, M> {
    /// Converts the given value to Montgomery form.
    pub fn new(value: WordRef<B, W, N>, modulus: M) -> Self {
        return MontgomeryWordRef {
            montgomery_val: modulus.to_montgomery(value),
            modulus,
        };
    }

    /// Wraps a value in Montgomery form, without converting it.
    ///
    /// ⚠️ Safety: The caller must ensure that the value is in Montgomery form.
    /// Failure to do so may result in incorrect behaviour.
    pub fn from_inner(montgomery_val: WordRef<B, W, N>, modulus: M) -> Self {
        return Self {
            montgomery_val,
            modulus,
        };
    }

    /// Returns the modulus associated with this Montgomery word.
    pub fn modulus(&self) -> M {
        return self.modulus;
    }

    /// Consumes this Montgomery word and returns its modulus.
    pub fn into_modulus(self) -> M {
        return self.modulus;
    }

    /// Converts this value out of Montgomery form.
    pub fn value(self) -> WordRef<B, W, N> {
        return self.modulus.from_montgomery(self.montgomery_val);
    }

    /// Returns the inner Montgomery value.
    pub fn inner(&self) -> &WordRef<B, W, N> {
        return &self.montgomery_val;
    }

    /// Consumes this Montgomery word and returns the inner Montgomery value.
    pub fn into_inner(self) -> WordRef<B, W, N> {
        return self.montgomery_val;
    }

    /// Consumes this Montgomery word and returns its inner Montgomery value and modulus.
    pub fn destructure(self) -> (WordRef<B, W, N>, M) {
        return (self.montgomery_val, self.modulus);
    }

    /// Returns a boolean reference that is true if this Montgomery word is zero, false otherwise.
    pub fn is_zero(self) -> BooleanWordRef<B> {
        return self.montgomery_val.is_zero();
    }

    /// Returns a boolean reference that is true if this Montgomery word is nonzero, false otherwise.
    pub fn is_nonzero(self) -> BooleanWordRef<B> {
        return self.montgomery_val.is_nonzero();
    }

    /// Compares this Montgomery word with another for equality, returning a boolean reference.
    pub fn eq(self, other: Self) -> BooleanWordRef<B> {
        if self.modulus != other.modulus {
            panic!("Cannot compare modular words with different moduli");
        }
        return self.montgomery_val.eq(other.montgomery_val);
    }

    /// Compares this Montgomery word with another for inequality, returning a boolean reference.
    pub fn ne(self, other: Self) -> BooleanWordRef<B> {
        if self.modulus != other.modulus {
            panic!("Cannot compare modular words with different moduli");
        }
        return self.montgomery_val.ne(other.montgomery_val);
    }

    /// Compares this Montgomery word with a constant for equality, returning a boolean reference.
    pub fn eq_const(self, other: MontgomeryWord<W, N, M>) -> BooleanWordRef<B> {
        return self.montgomery_val.eq_const(other.into_inner());
    }

    /// Compares this Montgomery word with a constant for inequality, returning a boolean reference.
    pub fn ne_const(self, other: MontgomeryWord<W, N, M>) -> BooleanWordRef<B> {
        return self.montgomery_val.ne_const(other.into_inner());
    }

    /// Modular addition with a constant.
    pub fn add_const(self, rhs: CompositeWord<W, N>) -> Self {
        let (res, carry) = self
            .montgomery_val
            .overflowing_add_const(self.modulus.to_montgomery_const(rhs));
        return Self {
            montgomery_val: (carry | res.clone().ge_const(self.modulus.n()))
                .select(res.clone() - self.modulus.n(), res),
            modulus: self.modulus,
        };
    }

    /// Modular multiplication by a constant.
    pub fn mul_const(self, rhs: CompositeWord<W, N>) -> Self {
        let (res_lo, res_hi) = self
            .montgomery_val
            .wide_mul_const(self.modulus.to_montgomery_const(rhs));
        return Self {
            montgomery_val: self.modulus.redc_wide(res_lo, res_hi),
            modulus: self.modulus,
        };
    }

    /// Reduces this Montgomery word modulo the modulus,
    /// returning a new Montgomery word in Montgomery form.
    pub fn into_zero(self) -> Self {
        let modulus = self.modulus;
        let montgomery_val = self.montgomery_val;
        return Self {
            montgomery_val: montgomery_val.into_zero(),
            modulus,
        };
    }

    /// Consumes this Montgomery word and returns a new Montgomery word in Montgomery form
    /// representing the given constant value.
    pub fn into_const(self, word: CompositeWord<W, N>) -> Self {
        let modulus = self.modulus;
        let montgomery_val = self.montgomery_val;
        let word = modulus.to_montgomery_const(word);
        return Self {
            montgomery_val: montgomery_val.into_const_same_width(word),
            modulus,
        };
    }

    /// Consumes this Montgomery word and returns a new Montgomery word in Montgomery form
    /// representing the given constant value already in Montgomery form.
    ///
    /// ⚠️ Safety: The caller must ensure that the given constant value is in Montgomery form.
    /// Failure to do so may result in incorrect behaviour.
    pub fn into_montgomery_const(self, word: MontgomeryWord<W, N, M>) -> Self {
        let modulus = self.modulus;
        let montgomery_val = self.montgomery_val;
        return Self {
            montgomery_val: montgomery_val.into_const_same_width(word.into_inner()),
            modulus,
        };
    }

    /// Modular inverse of this Montgomery word.
    ///
    /// Uses the constant-time Bernstein–Yang safegcd algorithm (see [crate::safegcd]), which needs
    /// only that the modulus be odd and `self` be coprime to it — no [MontgomeryMod::inv_exp] is
    /// required. For an invertible value it returns the inverse in Montgomery form; for a
    /// non-invertible value (e.g. zero) it returns zero, matching the Fermat convention.
    pub fn inv(self) -> Self {
        let modulus = self.modulus;
        // safegcd works on the canonical integer residue; convert in, invert, convert back.
        let a_int = modulus.from_montgomery(self.montgomery_val);
        let inv_int = crate::safegcd::safegcd_invert(a_int, modulus.n());
        return Self {
            montgomery_val: modulus.to_montgomery(inv_int),
            modulus,
        };
    }
}

impl<B: Backend, W: Word, const N: usize, M: MontgomeryMod<W, N>> Neg
    for MontgomeryWordRef<B, W, N, M>
{
    type Output = Self;

    /// Modular negation.
    fn neg(self) -> Self::Output {
        let value = self.montgomery_val;
        return Self {
            montgomery_val: value
                .clone()
                .is_zero()
                .select(value.clone(), -value + self.modulus.n()),
            modulus: self.modulus,
        };
    }
}

impl<B: Backend, W: Word, const N: usize, M: MontgomeryMod<W, N>> Add<Self>
    for MontgomeryWordRef<B, W, N, M>
{
    type Output = Self;

    /// Modular addition.
    fn add(self, rhs: Self) -> Self::Output {
        if self.modulus != rhs.modulus {
            panic!("Cannot add modular words with different moduli");
        }
        let (res, carry) = self.montgomery_val.overflowing_add(rhs.montgomery_val);
        return Self {
            montgomery_val: (carry | res.clone().ge_const(self.modulus.n()))
                .select(res.clone() - self.modulus.n(), res),
            modulus: self.modulus,
        };
    }
}

impl<B: Backend, W: Word, const N: usize, M: MontgomeryMod<W, N>> Add<MontgomeryWord<W, N, M>>
    for MontgomeryWordRef<B, W, N, M>
{
    type Output = Self;

    /// Modular addition.
    fn add(self, rhs: MontgomeryWord<W, N, M>) -> Self::Output {
        if self.modulus != rhs.modulus {
            panic!("Cannot add modular words with different moduli");
        }
        let (res, carry) = self
            .montgomery_val
            .overflowing_add_const(rhs.montgomery_val);
        return Self {
            montgomery_val: (carry | res.clone().ge_const(self.modulus.n()))
                .select(res.clone() - self.modulus.n(), res),
            modulus: self.modulus,
        };
    }
}

impl<B: Backend, W: Word, const N: usize, M: MontgomeryMod<W, N>> Sub<Self>
    for MontgomeryWordRef<B, W, N, M>
{
    type Output = Self;
    /// Modular subtraction.
    fn sub(self, rhs: Self) -> Self::Output {
        if self.modulus != rhs.modulus {
            panic!("Cannot subtract modular words with different moduli");
        }

        let (res, borrow) = self.montgomery_val.overflowing_sub(rhs.montgomery_val);
        return Self {
            montgomery_val: borrow.select(res.clone() + self.modulus.n(), res),
            modulus: self.modulus,
        };
    }
}

impl<B: Backend, W: Word, const N: usize, M: MontgomeryMod<W, N>> Sub<MontgomeryWord<W, N, M>>
    for MontgomeryWordRef<B, W, N, M>
{
    type Output = Self;
    /// Modular subtraction.
    fn sub(self, rhs: MontgomeryWord<W, N, M>) -> Self::Output {
        if self.modulus != rhs.modulus {
            panic!("Cannot subtract modular words with different moduli");
        }

        let (res, borrow) = self
            .montgomery_val
            .overflowing_sub_const(rhs.montgomery_val);
        return Self {
            montgomery_val: borrow.select(res.clone() + self.modulus.n(), res),
            modulus: self.modulus,
        };
    }
}

impl<B: Backend, W: Word, const N: usize, M: MontgomeryMod<W, N>> Mul<Self>
    for MontgomeryWordRef<B, W, N, M>
{
    type Output = Self;
    /// Modular multiplication.
    fn mul(self, rhs: Self) -> Self::Output {
        if self.modulus != rhs.modulus {
            panic!("Cannot multiply modular words with different moduli");
        }
        let (res_lo, res_hi) = self.montgomery_val.wide_mul(rhs.montgomery_val);
        return Self {
            montgomery_val: self.modulus.redc_wide(res_lo, res_hi),
            modulus: self.modulus,
        };
    }
}

impl<B: Backend, W: Word, const N: usize, M: MontgomeryMod<W, N>> Mul<MontgomeryWord<W, N, M>>
    for MontgomeryWordRef<B, W, N, M>
{
    type Output = Self;
    /// Modular multiplication.
    fn mul(self, rhs: MontgomeryWord<W, N, M>) -> Self::Output {
        if self.modulus != rhs.modulus {
            panic!("Cannot multiply modular words with different moduli");
        }
        let (res_lo, res_hi) = self.montgomery_val.wide_mul_const(rhs.montgomery_val);
        return Self {
            montgomery_val: self.modulus.redc_wide(res_lo, res_hi),
            modulus: self.modulus,
        };
    }
}

/// Helper implementing Montgomery word selection by a [BooleanWordRef].
pub trait MontgomeryBooleanWordRefSelector<
    B: Backend,
    W: Word,
    const N: usize,
    M: MontgomeryMod<W, N>,
>
{
    /// Selects between two constant Montgomery words.
    fn montgomery_select_const_const(
        self,
        then: MontgomeryWord<W, N, M>,
        else_: MontgomeryWord<W, N, M>,
    ) -> MontgomeryWordRef<B, W, N, M>;

    /// Selects between a constant Montgomery word and a variable Montgomery word.
    fn montgomery_select_const_var(
        self,
        then: MontgomeryWord<W, N, M>,
        else_: MontgomeryWordRef<B, W, N, M>,
    ) -> MontgomeryWordRef<B, W, N, M>;

    /// Selects between a variable Montgomery word and a constant Montgomery word.
    fn montgomery_select_var_const(
        self,
        then: MontgomeryWordRef<B, W, N, M>,
        else_: MontgomeryWord<W, N, M>,
    ) -> MontgomeryWordRef<B, W, N, M>;

    /// Selects between two variable Montgomery words.
    fn montgomery_select(
        self,
        then: MontgomeryWordRef<B, W, N, M>,
        else_: MontgomeryWordRef<B, W, N, M>,
    ) -> MontgomeryWordRef<B, W, N, M>;
}

impl<B: Backend, W: Word, const N: usize, M: MontgomeryMod<W, N>>
    MontgomeryBooleanWordRefSelector<B, W, N, M> for BooleanWordRef<B>
{
    fn montgomery_select_const_const(
        self,
        then: MontgomeryWord<W, N, M>,
        else_: MontgomeryWord<W, N, M>,
    ) -> MontgomeryWordRef<B, W, N, M> {
        if then.modulus != else_.modulus {
            panic!("Cannot select between words with different moduli");
        }
        let (then, modulus) = then.destructure();
        let else_ = else_.into_inner();
        return MontgomeryWordRef {
            montgomery_val: self.select_const_const(then, else_),
            modulus,
        };
    }

    fn montgomery_select_const_var(
        self,
        then: MontgomeryWord<W, N, M>,
        else_: MontgomeryWordRef<B, W, N, M>,
    ) -> MontgomeryWordRef<B, W, N, M> {
        if then.modulus != else_.modulus {
            panic!("Cannot select between words with different moduli");
        }
        let (then, modulus) = then.destructure();
        let else_ = else_.into_inner();
        return MontgomeryWordRef {
            montgomery_val: self.select_const_var(then, else_),
            modulus,
        };
    }

    fn montgomery_select_var_const(
        self,
        then: MontgomeryWordRef<B, W, N, M>,
        else_: MontgomeryWord<W, N, M>,
    ) -> MontgomeryWordRef<B, W, N, M> {
        if then.modulus != else_.modulus {
            panic!("Cannot select between words with different moduli");
        }
        let (then, modulus) = then.destructure();
        let else_ = else_.into_inner();
        return MontgomeryWordRef {
            montgomery_val: self.select_var_const(then, else_),
            modulus,
        };
    }

    fn montgomery_select(
        self,
        then: MontgomeryWordRef<B, W, N, M>,
        else_: MontgomeryWordRef<B, W, N, M>,
    ) -> MontgomeryWordRef<B, W, N, M> {
        if then.modulus != else_.modulus {
            panic!("Cannot select between words with different moduli");
        }
        let (then, modulus) = then.destructure();
        let else_ = else_.into_inner();
        return MontgomeryWordRef {
            montgomery_val: self.select(then, else_),
            modulus,
        };
    }
}

/// Helper trait implementing Montgomery word allocation from a [WordRef].
pub trait MontgomeryWordRefAllocator<B: Backend, W: Word, const N: usize, M: MontgomeryMod<W, N>> {
    /// Variant of [WordRef::alloc_constant] for Montgomery words.
    fn alloc_montgomery_constant(
        &self,
        word: MontgomeryWord<W, N, M>,
    ) -> MontgomeryWordRef<B, W, N, M>;
}

impl<_B: Backend, W: Word, const N: usize, M: MontgomeryMod<W, N>, _W: Word, const _N: usize>
    MontgomeryWordRefAllocator<_B, W, N, M> for WordRef<_B, _W, _N>
{
    fn alloc_montgomery_constant(
        &self,
        word: MontgomeryWord<W, N, M>,
    ) -> MontgomeryWordRef<_B, W, N, M> {
        if word.is_zero() {
            return MontgomeryWordRef::new(self.alloc_new_zero(), word.modulus());
        }
        return MontgomeryWordRef::new(self.alloc_new_word(word.value()), word.modulus());
    }
}

/// Helper trait implementing Montgomery word allocation for a [Frontend].
pub trait MontgomeryFrontendIO<B: Backend, W: Word, const N: usize, M: MontgomeryMod<W, N>> {
    /// Variant of [Frontend::input] for Montgomery words.
    fn montgomery_input(&self, in_: MontgomeryWord<W, N, M>) -> MontgomeryWordRef<B, W, N, M>;

    /// Variant of [Frontend::alloc] for Montgomery words.
    fn montgomery_alloc(&self, in_: MontgomeryWord<W, N, M>) -> MontgomeryWordRef<B, W, N, M>;

    /// Variant of [Frontend::output] for Montgomery words.
    fn montgomery_output(&self, out: MontgomeryWordRef<B, W, N, M>);
}

impl<B: Backend, W: Word, const N: usize, M: MontgomeryMod<W, N>> MontgomeryFrontendIO<B, W, N, M>
    for Frontend<B>
{
    fn montgomery_input(&self, in_: MontgomeryWord<W, N, M>) -> MontgomeryWordRef<B, W, N, M> {
        let (in_, modulus) = in_.destructure();
        let in_ = self.input(in_);
        return MontgomeryWordRef::from_inner(in_, modulus);
    }

    fn montgomery_alloc(&self, in_: MontgomeryWord<W, N, M>) -> MontgomeryWordRef<B, W, N, M> {
        let (in_, modulus) = in_.destructure();
        let in_ = self.alloc(in_);
        return MontgomeryWordRef::from_inner(in_, modulus);
    }

    fn montgomery_output(&self, out: MontgomeryWordRef<B, W, N, M>) {
        self.output(out.value());
    }
}
