# Changelog

All notable changes to this project are documented in this file.

## [1.0.0] - Unreleased

### Breaking changes

- Constructors now return `Result<Self, QuantileError>` and validate zero sizes and non-positive, NaN, or infinite epsilon values.
- `query` now takes `&self`, returns `Result<T, QuantileError>`, and validates empty summaries and ranks outside `0.0..=1.0`.
- `FixedSizeEpsilonSummary::update` now panics if more than its declared `n` elements are inserted, preventing silent loss of rank accuracy.
- The crate now uses Rust 2021 and declares Rust 1.65 as its library MSRV.

### Added

- Optional `serde` serialization support, with query caches skipped and rebuilt lazily.
- `Clone` support for both summary types, `Debug` implementations, and `#[must_use]` annotations on result-bearing APIs.
- Rank-error contract tests for fixed and unbounded summaries, including streams beyond two million elements.
- Expanded API, algorithm, architecture, error-handling, and benchmark documentation.

### Fixed

- Prevented fixed-size summaries from silently discarding a compressed block when their configured levels were full.
- Sized unbounded substream summaries to the next Zhang-Wang boundary and corrected fixed-summary level allocation for binary block carries.
- Fixed the unbounded boundary index overflow for large streams.
- Fixed empty-summary underflow and the inverted scan termination condition in quantile lookup.
- Made `RankInfo` ordering internally consistent and removed an unnecessary `Debug` bound.
- Corrected malformed rustdoc and a misnamed fixed-summary test.

### Performance and engineering

- Reused owned input blocks during compression and used unstable sorting where stability is irrelevant.
- Updated development dependencies and modernized CI, including all-target clippy, MSRV, formatting, test, coverage, and lower-cost benchmark jobs.
- Added black-boxed warm-cache and cold-cache Criterion query benchmarks and refreshed the README results.

## [0.2.3] - 2021-12-30

- Last release before the 1.0 API and correctness modernization.
