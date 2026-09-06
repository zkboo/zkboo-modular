# Changelog

All notable changes to this crate are documented in this file, starting at 1.2.0.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- The tests build their assertion accumulators with `Assertions::scope`, following the core crate.
  The advised-inverse API is unchanged: it still takes the accumulator its caller supplies.

### Added

- `Zeroize` for `MontgomeryWord`, which erases the value and leaves the modulus alone.
  A field element can hold witness-derived material, and a caller holding one had no way to erase it.

## [1.2.0] — 2026-09-04

### Changed

- Adapted to the single entry points for proving, verifying and executing.
