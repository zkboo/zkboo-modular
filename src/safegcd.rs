// SPDX-License-Identifier: LGPL-3.0-or-later

//! Constant-time modular inversion via the Bernstein–Yang "safegcd" algorithm.

use zkboo::backend::{Backend, BooleanWordRef, WordRef};
use zkboo::word::{CompositeWord, Word};

/// Number of divsteps that guarantees `gcd` completion for `width`-bit inputs, with a small margin.
pub fn divstep_count(width: usize) -> usize {
    let n = if width < 46 { 46 } else { width };
    return (49 * n + 80) / 17 + 2;
}

/// A two's-complement integer one machine word wider than the field, held as `lo + ext·2^WIDTH`.
struct Wide<B: Backend, W: Word, const N: usize> {
    lo: WordRef<B, W, N>,
    ext: WordRef<B, W, 1>,
}

impl<B: Backend, W: Word, const N: usize> Wide<B, W, N> {
    const WIDTH: usize = W::WIDTH * N;

    /// Zero-extends a non-negative `N`-word value (the sign word is zero, so the value is positive
    /// regardless of the top bit of `lo` — this is exactly why the extra word is needed).
    fn from_unsigned(lo: WordRef<B, W, N>) -> Self {
        let ext = lo.alloc_new_zero::<W, 1>();
        return Self { lo, ext };
    }

    fn clone_of(&self) -> Self {
        return Self {
            lo: self.lo.clone(),
            ext: self.ext.clone(),
        };
    }

    /// Two's-complement addition across both limbs.
    fn add(self, rhs: Self) -> Self {
        let (lo, carry) = self.lo.overflowing_add(rhs.lo);
        let ext = self
            .ext
            .wrapping_add(rhs.ext)
            .wrapping_add(WordRef::<B, W, 1>::from_bool(carry));
        return Self { lo, ext };
    }

    /// Two's-complement negation: bitwise-not then add one.
    fn neg(self) -> Self {
        let one = Self {
            lo: self.lo.alloc_new_word(CompositeWord::<W, N>::ONE),
            ext: self.ext.alloc_new_zero::<W, 1>(),
        };
        let not = Self {
            lo: self.lo.not(),
            ext: self.ext.not(),
        };
        return not.add(one);
    }

    /// Arithmetic shift right by one (sign-preserving exact halving of an even value).
    fn ashr1(self) -> Self {
        let width = Self::WIDTH;
        // Bit `WIDTH-1` of the new low word is bit 0 of `ext`.
        let ext_lsb_into_lo = WordRef::<B, W, N>::from_bool(self.ext.clone().lsb()) << (width - 1);
        let lo = (self.lo >> 1).bitxor(ext_lsb_into_lo);
        // The sign bit of `ext` is replicated (arithmetic shift of the top word).
        let ext_msb_into_top = WordRef::<B, W, 1>::from_bool(self.ext.clone().msb()) << (W::WIDTH - 1);
        let ext = (self.ext >> 1).bitxor(ext_msb_into_top);
        return Self { lo, ext };
    }

    /// `cond ? a : b`, selecting both limbs.
    fn select(cond: BooleanWordRef<B>, a: Self, b: Self) -> Self {
        let lo = cond.clone().select(a.lo, b.lo);
        let ext = cond.select(a.ext, b.ext);
        return Self { lo, ext };
    }

    /// True iff the value is negative (its sign bit, the top bit of `ext`).
    fn is_neg(&self) -> BooleanWordRef<B> {
        return self.ext.clone().msb();
    }
}

/// `(-x) mod p` for `x ∈ [0, p)`: `0` if `x == 0`, else `p − x`.
fn neg_mod_p<B: Backend, W: Word, const N: usize>(
    x: WordRef<B, W, N>,
    p: CompositeWord<W, N>,
) -> WordRef<B, W, N> {
    let is_zero = x.clone().is_zero();
    let (p_minus_x, _) = x.overflowing_sub_from_const(p);
    return is_zero.select_const_var(CompositeWord::<W, N>::ZERO, p_minus_x);
}

/// `(x / 2) mod p` for `x ∈ [0, p)`, with `p` odd: `x/2` if even, else `(x + p)/2`.
fn halve_mod_p<B: Backend, W: Word, const N: usize>(
    x: WordRef<B, W, N>,
    p: CompositeWord<W, N>,
) -> WordRef<B, W, N> {
    let width = W::WIDTH * N;
    let odd = x.clone().lsb();
    // t = x + (odd ? p : 0); this can carry out of the top (x + p < 2p < 2^{WIDTH+1}).
    let addend = odd.clone().select_const_const(p, CompositeWord::<W, N>::ZERO);
    let (t, carry) = x.overflowing_add(addend);
    // (carry·2^WIDTH + t) is even; divide by two, shifting the carry bit into the top.
    let carry_into_top = WordRef::<B, W, N>::from_bool(carry) << (width - 1);
    return (t >> 1).bitxor(carry_into_top);
}

/// Computes `a^{-1} mod p` for an `N`-word value `a ∈ [0, p)` and odd modulus `p`, in constant
/// time.
pub fn safegcd_invert<B: Backend, W: Word, const N: usize>(
    a: WordRef<B, W, N>,
    p: CompositeWord<W, N>,
) -> WordRef<B, W, N> {
    let width = W::WIDTH * N;
    assert!(
        W::WIDTH >= 16 || width + 2 < (1usize << (W::WIDTH - 1)),
        "safegcd: field width too large for the signed divstep counter in a {}-bit word",
        W::WIDTH,
    );
    let iters = divstep_count(width);

    // State.
    let mut delta = a.alloc_new_word::<W, 1, _>(CompositeWord::<W, 1>::ONE);
    let mut f = Wide::from_unsigned(a.alloc_new_word(p));
    let mut g = Wide::from_unsigned(a);
    let mut d = f.lo.alloc_new_zero::<W, N>();
    let mut e = f.lo.alloc_new_word(CompositeWord::<W, N>::ONE);

    let one_w1 = CompositeWord::<W, 1>::ONE;

    for _ in 0..iters {
        let g_odd = g.lo.clone().lsb();
        // δ > 0  ⇔  δ is non-negative and nonzero.
        let delta_pos = (!delta.clone().msb()) & delta.clone().is_nonzero();
        let cond = delta_pos & g_odd.clone();

        // s·f for the g update: cond ? −f : (g odd ? f : 0).
        let sf = Wide::select(
            cond.clone(),
            f.clone_of().neg(),
            Wide::select(
                g_odd.clone(),
                f.clone_of(),
                Wide {
                    lo: f.lo.alloc_new_zero::<W, N>(),
                    ext: f.ext.alloc_new_zero::<W, 1>(),
                },
            ),
        );
        let new_g = g.clone_of().add(sf).ashr1();
        let new_f = Wide::select(cond.clone(), g.clone_of(), f.clone_of());

        // Coefficient update: s·d as a residue mod p (cond ? −d : (g odd ? d : 0)).
        let neg_d = neg_mod_p(d.clone(), p);
        let sd = cond.clone().select(
            neg_d,
            g_odd.select_var_const(d.clone(), CompositeWord::<W, N>::ZERO),
        );
        // e + s·d, reduced to [0, p), then halved mod p.
        let (e_plus, carry) = e.clone().overflowing_add(sd);
        let ge = carry | e_plus.clone().ge_const(p);
        let e_reduced = ge.select(e_plus.clone() - p, e_plus);
        let new_e = halve_mod_p(e_reduced, p);
        let new_d = cond.clone().select(e, d);

        // δ update: cond ? 1 − δ : 1 + δ.
        let one_plus = delta.clone().wrapping_add_const(one_w1);
        let one_minus = delta.wrapping_neg().wrapping_add_const(one_w1);
        delta = cond.select(one_minus, one_plus);

        f = new_f;
        g = new_g;
        d = new_d;
        e = new_e;
    }

    // f == ±1 now (for invertible a); a^{-1} = d if f > 0 else −d mod p.
    let f_neg = f.is_neg();
    return f_neg.select(neg_mod_p(d.clone(), p), d);
}
