# Changelog

All notable changes to this crate are documented in this file, starting at 1.2.0.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- `MontgomeryWord::value` and `MontgomeryWordRef::value` are `canonical`.
  The old name read like the value one would hand back to `from_inner`, when it is the one that does not round-trip through it.
- The tests build their assertion accumulators with `Assertions::scope`, following the core crate.
  The advised-inverse API is unchanged: it still takes the accumulator its caller supplies.

### Added

- `MontgomeryFrontendIO::montgomery_output_inner`, emitting a field element in the representation it is stored in, for a cleartext pass whose output is read back as advice.
  Its sibling `montgomery_output` emits the canonical residue, and the two look equally reasonable at the point where one has to be chosen.
- `Montgomery`, a wrapper for a field element's value in the representation its modulus stores it in.
  `MontgomeryWord::inner`, `into_inner` and `from_inner` carry it, so a stored value and a canonical residue can no longer be passed for one another: both are a `CompositeWord`, both lie in `[0, p)`, and under a pseudo-Mersenne modulus they are the same word, so nothing but a type can tell them apart.
- `Zeroize` for `MontgomeryWord`, which erases the value and leaves the modulus alone.
  A field element can hold witness-derived material, and a caller holding one had no way to erase it.

## [1.2.0] — 2026-09-04

### Changed

- Adapted to the single entry points for proving, verifying and executing.
