// SPDX-License-Identifier: LGPL-3.0-or-later

//! The [FieldRep] trait abstracts a prime field's *internal representation* and reduction strategy,
//! so that a single field-element type ([crate::montgomery::MontgomeryWord] /
//! [crate::montgomery::MontgomeryWordRef]) can be backed by either Montgomery reduction
//! ([crate::montgomery::MontgomeryMod]) or pseudo-Mersenne reduction
//! ([crate::pseudo_mersenne::PseudoMersenneMod]).

use crate::montgomery::MontgomeryWord;
use zkboo::backend::{Backend, WordRef};
use zkboo::word::{CompositeWord, Word, WordLike};
use core::fmt::Debug;

/// Reduction strategy / internal representation for a prime field of width `W::WIDTH · N`.
pub trait FieldRep<W: Word, const N: usize>: Clone + Copy + Debug + PartialEq + Eq {
    /// The field modulus `p`.
    fn modulus(&self) -> CompositeWord<W, N>;

    /// Encodes a canonical value `[0, p)` into the internal representation (build-time).
    fn encode_const(&self, value: CompositeWord<W, N>) -> CompositeWord<W, N>;

    /// Decodes an internal value back to its canonical residue `[0, p)` (build-time).
    fn decode_const(&self, internal: CompositeWord<W, N>) -> CompositeWord<W, N>;

    /// Reduces a double-width internal product `(lo, hi)` to a single internal value (build-time).
    fn mul_reduce_const(
        &self,
        lo: CompositeWord<W, N>,
        hi: CompositeWord<W, N>,
    ) -> CompositeWord<W, N>;

    /// Encodes a canonical value into the internal representation (in-circuit).
    fn encode<B: Backend>(&self, value: WordRef<B, W, N>) -> WordRef<B, W, N>;

    /// Decodes an internal value back to its canonical residue (in-circuit).
    fn decode<B: Backend>(&self, internal: WordRef<B, W, N>) -> WordRef<B, W, N>;

    /// Reduces a double-width internal product `(lo, hi)` to a single internal value (in-circuit).
    fn mul_reduce<B: Backend>(
        &self,
        lo: WordRef<B, W, N>,
        hi: WordRef<B, W, N>,
    ) -> WordRef<B, W, N>;

    /// Native modular inverse of an internal value, in the internal representation.
    #[inline]
    fn invert_const(&self, internal: CompositeWord<W, N>) -> CompositeWord<W, N>
    where
        Self: Sized,
    {
        return MontgomeryWord::from_inner(internal, *self)
            .fermat_inv()
            .into_inner();
    }

    /// Creates a constant field element representing the given canonical value.
    #[inline]
    fn const_word<U: WordLike<W, N>>(&self, value: U) -> MontgomeryWord<W, N, Self>
    where
        Self: Sized,
    {
        return MontgomeryWord::new(value, *self);
    }

    /// Creates the constant field element zero.
    #[inline]
    fn zero_word(&self) -> MontgomeryWord<W, N, Self>
    where
        Self: Sized,
    {
        return MontgomeryWord::from_inner(CompositeWord::<W, N>::ZERO, *self);
    }

    /// Creates the constant field element one.
    #[inline]
    fn one_word(&self) -> MontgomeryWord<W, N, Self>
    where
        Self: Sized,
    {
        return MontgomeryWord::new(CompositeWord::<W, N>::ONE, *self);
    }
}
