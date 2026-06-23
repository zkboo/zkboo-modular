// SPDX-License-Identifier: LGPL-3.0-or-later

//! Pseudo-Mersenne (Solinas) reduction for moduli of the form `p = 2^WIDTH - c` with a small,
//! low-Hamming-weight `c`, as an alternative to Montgomery reduction.

use zkboo::backend::{Backend, WordRef};
use zkboo::word::{CompositeWord, Word};

/// An odd modulus of pseudo-Mersenne form `p = 2^(W::WIDTH·N) - c`, with `c` small.
pub trait PseudoMersenneMod<W: Word, const N: usize>: Clone + Copy {
    /// The modulus `p`.
    fn p(&self) -> CompositeWord<W, N>;

    /// The folding constant `c = 2^WIDTH - p`.
    fn c(&self) -> CompositeWord<W, N>;
}

/// Reduces a double-width value `(lo, hi)` (representing `hi·2^WIDTH + lo`) modulo `p = 2^WIDTH−c`,
/// returning a canonical residue in `[0, p)`.
pub fn reduce_wide<B: Backend, W: Word, const N: usize, M: PseudoMersenneMod<W, N>>(
    modulus: &M,
    lo: WordRef<B, W, N>,
    hi: WordRef<B, W, N>,
) -> WordRef<B, W, N> {
    let c = modulus.c();
    let p = modulus.p();

    // Fold 1: value ≡ lo + hi·c. `hi·c` has a low word `h1_lo` and a small high word `h1_hi`.
    let (h1_lo, h1_hi) = hi.wide_mul_const(c);
    let (s1, carry1) = lo.overflowing_add(h1_lo);
    // The bits above WIDTH are `h1_hi + carry1`; fold them in too (still ≡ ·c).
    let top1 = h1_hi.wrapping_add(WordRef::<B, W, N>::from_bool(carry1));

    // Fold 2: value ≡ s1 + top1·c. `top1` is small, so `top1·c` fits in WIDTH bits (low half only).
    let t2 = top1.wrapping_mul_const(c);
    let (s2, carry2) = s1.overflowing_add(t2);

    // Fold 3: at most one extra multiple of `2^WIDTH ≡ c` remains.
    let add_c = carry2.select_const_const(c, CompositeWord::<W, N>::ZERO);
    let (s3, carry3) = s2.overflowing_add(add_c);
    let add_c2 = carry3.select_const_const(c, CompositeWord::<W, N>::ZERO);
    let s4 = s3.wrapping_add(add_c2);

    // s4 is now < 2^WIDTH + (small); bring it into [0, p) with up to two conditional subtractions.
    let ge1 = s4.clone().ge_const(p);
    let s5 = ge1.select(s4.clone() - p, s4);
    let ge2 = s5.clone().ge_const(p);
    return ge2.select(s5.clone() - p, s5);
}

/// Field multiplication `a · b mod p` in canonical (non-Montgomery) form.
pub fn mul<B: Backend, W: Word, const N: usize, M: PseudoMersenneMod<W, N>>(
    modulus: &M,
    a: WordRef<B, W, N>,
    b: WordRef<B, W, N>,
) -> WordRef<B, W, N> {
    let (lo, hi) = a.wide_mul(b);
    return reduce_wide(modulus, lo, hi);
}
