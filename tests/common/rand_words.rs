// SPDX-License-Identifier: LGPL-3.0-or-later

use rand::{Rng, RngExt, SeedableRng, rngs::StdRng};
use zkboo::word::{CompositeWord, Word, WordLike};

pub fn rand_word<W: Word, const N: usize, U: WordLike<W, N>>(rng: &mut impl Rng) -> U {
    let special_cases: [U; _] = [
        U::from_word(CompositeWord::<W, N>::ZERO),
        U::from_word(CompositeWord::<W, N>::ONE),
        U::from_word(CompositeWord::<W, N>::ONE << 1),
        U::from_word(CompositeWord::<W, N>::ONE << CompositeWord::<W, N>::WIDTH / 2),
        U::from_word(CompositeWord::<W, N>::ONE << CompositeWord::<W, N>::WIDTH - 1),
        U::from_word(CompositeWord::<W, N>::MAX),
    ];
    let choice = rng.random_bool(0.5);
    if choice {
        return special_cases[rng.random_range(0..special_cases.len())];
    }
    let mut buf = [W::Bytes::default(); N];
    for b in buf.iter_mut() {
        rng.fill_bytes(b.as_mut());
    }
    return U::from_word(CompositeWord::<W, N>::from_le_bytes(buf));
}

pub fn test_vec<W: Word, const N: usize, U: WordLike<W, N>>(
    num_samples: usize,
    seed: u64,
) -> Vec<U> {
    let mut samples = Vec::new();
    let mut rng = StdRng::seed_from_u64(seed);
    samples.extend((0..num_samples).map(|_| rand_word::<W, N, U>(&mut rng)));
    return samples;
}
