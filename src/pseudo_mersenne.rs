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

/// Builds-time (native) analogue of [reduce_wide].
pub fn reduce_wide_const<W: Word, const N: usize, M: PseudoMersenneMod<W, N>>(
    modulus: &M,
    lo: CompositeWord<W, N>,
    hi: CompositeWord<W, N>,
) -> CompositeWord<W, N> {
    let c = modulus.c();
    let p = modulus.p();
    let (h1_lo, h1_hi) = hi.wide_mul(c);
    let (s1, carry1) = lo.overflowing_add(h1_lo);
    let top1 = h1_hi.wrapping_add(CompositeWord::<W, N>::from_bool(carry1));
    let t2 = top1.wrapping_mul(c);
    let (s2, carry2) = s1.overflowing_add(t2);
    let (s3, carry3) = s2.overflowing_add(if carry2 { c } else { CompositeWord::ZERO });
    let s4 = if carry3 { s3.wrapping_add(c) } else { s3 };
    let s5 = if s4.ge(p) { s4.wrapping_sub(p) } else { s4 };
    return if s5.ge(p) { s5.wrapping_sub(p) } else { s5 };
}

/// Implements [crate::field::FieldRep] for a concrete [PseudoMersenneMod] type, so the shared field
/// element type ([crate::montgomery::MontgomeryWord] / [crate::montgomery::MontgomeryWordRef]) can
/// be backed by pseudo-Mersenne reduction.
#[macro_export]
macro_rules! impl_pseudo_mersenne_field_rep {
    ($t:ty, $w:ty, $n:literal) => {
        impl $crate::field::FieldRep<$w, $n> for $t {
            fn modulus(&self) -> ::zkboo::word::CompositeWord<$w, $n> {
                return $crate::pseudo_mersenne::PseudoMersenneMod::p(self);
            }
            fn encode_const(
                &self,
                value: ::zkboo::word::CompositeWord<$w, $n>,
            ) -> ::zkboo::word::CompositeWord<$w, $n> {
                let p = $crate::pseudo_mersenne::PseudoMersenneMod::p(self);
                return if value.ge(p) { value.wrapping_sub(p) } else { value };
            }
            fn decode_const(
                &self,
                internal: ::zkboo::word::CompositeWord<$w, $n>,
            ) -> ::zkboo::word::CompositeWord<$w, $n> {
                return internal;
            }
            fn mul_reduce_const(
                &self,
                lo: ::zkboo::word::CompositeWord<$w, $n>,
                hi: ::zkboo::word::CompositeWord<$w, $n>,
            ) -> ::zkboo::word::CompositeWord<$w, $n> {
                return $crate::pseudo_mersenne::reduce_wide_const(self, lo, hi);
            }
            fn encode<B: ::zkboo::backend::Backend>(
                &self,
                value: ::zkboo::backend::WordRef<B, $w, $n>,
            ) -> ::zkboo::backend::WordRef<B, $w, $n> {
                let p = $crate::pseudo_mersenne::PseudoMersenneMod::p(self);
                return value.clone().ge_const(p).select(value.clone() - p, value);
            }
            fn decode<B: ::zkboo::backend::Backend>(
                &self,
                internal: ::zkboo::backend::WordRef<B, $w, $n>,
            ) -> ::zkboo::backend::WordRef<B, $w, $n> {
                return internal;
            }
            fn mul_reduce<B: ::zkboo::backend::Backend>(
                &self,
                lo: ::zkboo::backend::WordRef<B, $w, $n>,
                hi: ::zkboo::backend::WordRef<B, $w, $n>,
            ) -> ::zkboo::backend::WordRef<B, $w, $n> {
                return $crate::pseudo_mersenne::reduce_wide(self, lo, hi);
            }
        }
    };
}
