# Changelog

All notable changes to this project are documented in this file.

## [1.0.1] - 2026-09-06

### Fixed

- `FixedSizeEpsilonSummary` no longer panics or answers mid-range ranks with the maximum element when `epsilon * n` is small (roughly `1 < epsilon * n < 5`). The Zhang-Wang block size degenerates there, so such streams are now stored exactly (`epsilon * n <= 8`). For example `new(13, 0.1)` used to answer every rank with the largest value, and `new(101, 0.01)` answered the median with the maximum.
- `UnboundEpsilonSummary` inherited that bug through its early sub-stream summaries, which broke most `epsilon` values on short streams (`new(0.03)` panicked after 233 elements of a zigzag sequence, `new(0.3)` after 23). Random-stream sweeps over 99 `epsilon` values are now clean.
- Merging now uses the same one-based rank convention as the sorted blocks it receives. Every merge level previously under-counted ranks by one, which consumed part of the error budget and caused off-by-one contract violations for large `epsilon`.
- A query whose landing element overshoots the error window now also considers the elements on the low side of the window and, failing that, returns the element whose rank interval is closest to the target instead of silently returning the maximum.
- `UnboundEpsilonSummary` no longer panics after `(2^31 - 1) / epsilon` elements. Sub-stream boundaries are computed on demand instead of read from a 32-entry table whose exhaustion branch was unreachable.
- Both constructors reject `epsilon` above `1.0` with `QuantileError::InvalidEpsilon`; `UnboundEpsilonSummary::new` previously reported `InvalidN` for such values and `FixedSizeEpsilonSummary::new` accepted them.

### Changed

- The `serde` layout of both summaries changed: `FixedSizeEpsilonSummary` no longer stores a level count and `UnboundEpsilonSummary` stores a boundary counter instead of a boundary table. Summaries serialized by 1.0.0 do not load in 1.0.1 and, because of the rank fix above, should be rebuilt from their source stream anyway. The layout is documented as unstable.
- The `Display` text of `QuantileError::InvalidEpsilon` is now `epsilon must be in (0.0, 1.0]`.
- The compression precision condition and the level-capacity invariant are checked with `debug_assert!` instead of `assert!`.

### Performance

- Level merging and compression move elements instead of cloning them, merging reads both inputs through a small cursor abstraction, and `UnboundEpsilonSummary::update` compares one integer instead of binary-searching a boundary table. On the bundled Criterion benchmark (5,000 `usize` values, `epsilon = 0.01`) inserts are 15-22% faster, warm-cache queries 8-23% faster and the fixed-size cold-cache query 19% faster than 1.0.0; inserting 200,000 `String` values is about 30% faster.
- The unbounded cold-cache query in that benchmark is about 11% slower, because short early sub-streams are now kept exactly and contribute more elements to the merged summary. The effect shrinks as the stream grows.

### Documentation

- Every public item now has rustdoc with `# Errors` and `# Panics` sections; the crate enables `#![warn(missing_docs)]` and `#![forbid(unsafe_code)]`, and CI builds the docs with warnings denied.
- Corrected the algorithm and architecture documents (level count, sub-stream sizing, boundary limit, a memory-tracking test that did not exist) and documented exact mode, the `epsilon` range and the serde caveat.
- Added short-stream contract tests over a grid of `(n, epsilon)` values, regression tests for the cases above, and boundary-tracking tests.
- Planning documents under `docs/plans/` and the CI configuration are no longer packaged into the crate.

## [1.0.0] - 2026-07-12

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
