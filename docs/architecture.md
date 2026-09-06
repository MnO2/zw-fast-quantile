# Architecture Overview

This document describes the code organization and design decisions in `zw-fast-quantile`.

## File Structure

```
src/
  lib.rs          — entire library implementation (single-file crate)
benches/
  zw_benchmark.rs — criterion benchmarks comparing ZW vs GK01
docs/
  algorithm.md    — algorithm internals and theory
  architecture.md — this file
```

## Public API

The library exports three items:

| Type | Description |
|------|-------------|
| `FixedSizeEpsilonSummary<T>` | Quantile summary for streams of known size |
| `UnboundEpsilonSummary<T>` | Quantile summary for streams of unknown size |
| `QuantileError` | Error type for construction and query failures |

Both summary types share the same interface:

- `new(...)` — constructor returning `Result<Self, QuantileError>`
- `update(e: T)` — insert a single element
- `query(r: f64)` — query a quantile rank, returns `Result<T, QuantileError>`
- `size() -> usize` — number of elements inserted

Every public item carries rustdoc with `# Errors` and `# Panics` sections where applicable; the crate enables `#![warn(missing_docs)]` and `#![forbid(unsafe_code)]`.

### Generic Constraints

Both types are generic over `T: Clone + Ord`. This means:

- Integer types (`i32`, `u64`, etc.) work directly.
- Floating-point types require a wrapper like `ordered_float::NotNan<f64>` since `f64` does not implement `Ord`.
- Custom types work as long as they implement `Clone + Ord`.

## Internal Components

### RankInfo\<T\>

```rust
struct RankInfo<T> {
    val: T,       // the element value
    rmin: i64,    // minimum possible rank (one-based)
    rmax: i64,    // maximum possible rank (one-based)
}
```

The core data unit stored in summaries. Implements `Ord` (delegating to `val`) for sorting during merge and compress operations. Ranks are one-based everywhere, matching the paper's merge formulas.

### Helper Functions

| Function | Signature | Purpose |
|----------|-----------|---------|
| `assign_exact_ranks` | `(&mut [RankInfo<T>])` | Gives a sorted block exact ranks `1..=len` |
| `merge` | `(impl MergeSource<T>, impl MergeSource<T>) -> Vec<RankInfo<T>>` | Merges two sorted summaries with rank bound recalculation. `MergeSource` is a two-method cursor implemented for `vec::IntoIter` (moves elements) and `&[RankInfo<T>]` (clones them), so the update cascade never clones and the query path never copies a level up front |
| `merged_bounds` | `((i64, i64), Option<(i64, i64)>, Option<(i64, i64)>) -> (i64, i64)` | Rank bounds of one element after merging, from its neighbours in the other summary |
| `compress` | `(Vec<RankInfo<T>>, usize, f64) -> Vec<RankInfo<T>>` | Reduces summary size in place (swap, not clone) while maintaining the epsilon guarantee |
| `query_rank` | `(&[RankInfo<T>], f64, usize, f64) -> Option<T>` | Turns a rank in `0.0..=1.0` into a one-based target rank and error window, shared by both summary types |
| `find_idx` | `(&[RankInfo<T>], i64, i64) -> Option<T>` | Binary search plus forward and backward scans for query resolution, with a nearest-interval fallback |
| `boundary` | `(u32, f64) -> usize` | Position of sub-stream boundary `x`, `floor((2^x - 1) / epsilon)`, saturating at `usize::MAX` |

## Design Decisions

### Single-File Crate

The entire implementation is in `src/lib.rs`. At roughly 1,200 lines (including tests), the algorithm is cohesive enough that splitting into multiple files would add indirection without improving clarity.

### Exact Mode for Short Streams

`FixedSizeEpsilonSummary::new` only uses the Zhang-Wang block size when `epsilon * n > 8` (`EXACT_MODE_THRESHOLD`). Below that it keeps every element in a single level. The formula `floor(log2(epsilon * n) / epsilon)` collapses to a handful of elements when `epsilon * n` is close to 1, which breaks both the compression precision condition and the rank-error guarantee; exact storage costs at most `8 / epsilon` elements, the same order as one block, so nothing is lost. `UnboundEpsilonSummary` inherits this for its short early sub-streams.

### Interior Mutability for Query Caching

Both summary types use `RefCell<Option<Vec<RankInfo<T>>>>` to cache the merged summary. This allows `query()` to take `&self` rather than `&mut self`, which is more ergonomic for callers (e.g., querying multiple quantiles without needing a mutable reference).

The cache is invalidated on every `update()` call via `self.cached_s_m.get_mut().take()`.

### Error Handling with QuantileError

Rather than panicking on invalid inputs, constructors and queries return `Result<_, QuantileError>`. The error type is `Copy + Clone + PartialEq + Eq` for easy matching, and implements both `Display` and `std::error::Error`.

Error variants:
- `InvalidN` — `n` must be > 0 (fixed-size only)
- `InvalidEpsilon` — epsilon must be in `(0.0, 1.0]`, so NaN and infinity are rejected
- `InvalidRank` — rank must be in `[0.0, 1.0]` and not NaN
- `EmptySummary` — cannot query before inserting values

The only panic in normal use is `FixedSizeEpsilonSummary::update` when more than the declared `n` elements are inserted; it is documented under `# Panics`. Internal invariants (the compression precision condition, the level-capacity invariant) are `debug_assert!`s.

### Serde Support

Serde is feature-gated behind `features = ["serde"]`. Cache fields (`cached_s_m`) are skipped during serialization with `#[serde(skip)]` since they are lazily rebuilt on the next query. This keeps serialized output minimal and avoids issues with `RefCell` serialization.

The serialized layout mirrors the private fields and is not a stable API: it changed in 1.0.1 and may change again. Summaries should be rebuilt from the source stream rather than migrated across crate versions.

### Sub-Stream Boundaries Computed on Demand

`UnboundEpsilonSummary` stores the index and position of the boundary that ends the current sub-stream (`next_boundary_index`, `next_boundary`) instead of a precomputed table. Each `update` compares the count with one integer; when they match, the next boundary is computed from `boundary(x + 1, epsilon)`. The computation saturates at `usize::MAX`, so the summary keeps working for as many elements as the platform can count. Earlier versions used a 32-entry table and panicked once the stream passed `(2^31 - 1) / epsilon` elements.

## Testing Strategy

Tests in `src/lib.rs` cover:

1. **Unit tests** — merge/compress correctness with hand-computed rank bounds, and a check that two exact blocks merge into exact ranks
2. **Constructor validation** — error cases for invalid `n`, negative/NaN/infinite epsilon, and epsilon above 1
3. **Small-n correctness** — exact rank recovery on tiny streams, including the parameters that used to derive a degenerate block size
4. **Rank-error contract** — random `u64` streams checked against the `epsilon * n` guarantee at 101 ranks, for large streams (10k, 100k and 2M elements) and for a grid of short streams that steps through every `epsilon * n` regime for both summary types
5. **Statistical distribution tests** — normal and Pareto distributions with 1M elements, verifying that estimated quantiles are within epsilon of exact values
6. **Sub-stream boundaries** — boundary positions, saturation, and the number of frozen summaries after each insertion
7. **Edge cases** — empty summaries, boundary ranks (0.0, 1.0), immutability of query, the query fallback path
8. **Serde round-trip** — serialization/deserialization preserves query results and boundary tracking (feature-gated)

## Benchmarks

`benches/zw_benchmark.rs` uses [criterion](https://github.com/bheisler/criterion.rs) to benchmark:

- **Update**: inserting 5000 sequential values with `epsilon = 0.01`
- **Query**: querying 10 quantiles (0.1 through 1.0) on a pre-built summary, with a warm cache and with a cold cache

Both fixed-size and unbound ZW variants are benchmarked alongside GK01 from the [quantiles](https://github.com/postmates/quantiles) crate.

Run with:
```bash
cargo bench
```
