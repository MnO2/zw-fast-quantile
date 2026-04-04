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

### Generic Constraints

Both types are generic over `T: Clone + Ord`. This means:

- Integer types (`i32`, `u64`, etc.) work directly.
- Floating-point types require a wrapper like `ordered_float::NotNan<f64>` since `f64` does not implement `Ord`.
- Custom types work as long as they implement `Clone + Ord`.

## Internal Components

### RankInfo\<T\>

```rust
struct RankInfo<T: Clone> {
    val: T,       // the element value
    rmin: i64,    // minimum possible rank
    rmax: i64,    // maximum possible rank
}
```

The core data unit stored in summaries. Implements `Ord` (delegating to `val`) for sorting during merge and compress operations.

### Helper Functions

| Function | Signature | Purpose |
|----------|-----------|---------|
| `merge` | `(Vec<RankInfo<T>>, &[RankInfo<T>]) -> Vec<RankInfo<T>>` | Merges two sorted summaries with rank bound recalculation |
| `compress` | `(Vec<RankInfo<T>>, usize, f64) -> Vec<RankInfo<T>>` | Reduces summary size while maintaining epsilon guarantee |
| `find_idx` | `(&[RankInfo<T>], i64, i64) -> Option<T>` | Binary search + scan for query resolution |
| `is_boundary` | `(usize, &[usize; 32]) -> Option<usize>` | Binary search over pre-computed boundary array (used by `UnboundEpsilonSummary`) |

## Design Decisions

### Single-File Crate

The entire implementation is in `src/lib.rs`. At ~950 lines (including tests), the algorithm is cohesive enough that splitting into multiple files would add indirection without improving clarity.

### Interior Mutability for Query Caching

Both summary types use `RefCell<Option<Vec<RankInfo<T>>>>` to cache the merged summary. This allows `query()` to take `&self` rather than `&mut self`, which is more ergonomic for callers (e.g., querying multiple quantiles without needing a mutable reference).

The cache is invalidated on every `update()` call via `self.cached_s_m.get_mut().take()`.

### Error Handling with QuantileError

Rather than panicking on invalid inputs, constructors and queries return `Result<_, QuantileError>`. The error type is `Copy + Clone + PartialEq + Eq` for easy matching, and implements both `Display` and `std::error::Error`.

Error variants:
- `InvalidN` — `n` must be > 0 (fixed-size only)
- `InvalidEpsilon` — epsilon must be positive and finite
- `InvalidRank` — rank must be in `[0.0, 1.0]` and not NaN
- `EmptySummary` — cannot query before inserting values

### Serde Support

Serde is feature-gated behind `features = ["serde"]`. Cache fields (`cached_s_m`) are skipped during serialization with `#[serde(skip)]` since they are lazily rebuilt on the next query. This keeps serialized output minimal and avoids issues with `RefCell` serialization.

### Fixed-Size Boundary Array

`UnboundEpsilonSummary` uses `[usize; 32]` for its boundary table rather than a `Vec`. Since boundaries grow as `(2^i - 1) / epsilon`, 32 entries cover streams up to `2^32 / epsilon` elements — far beyond practical limits. The fixed array avoids heap allocation and enables efficient binary search via `is_boundary()`.

## Testing Strategy

Tests in `src/lib.rs` cover:

1. **Unit tests** — merge/compress correctness with hand-computed rank bounds
2. **Constructor validation** — error cases for invalid `n`, negative/NaN/infinite epsilon
3. **Small-n correctness** — exact rank recovery on tiny streams
4. **Statistical distribution tests** — normal and Pareto distributions with 1M elements, verifying that estimated quantiles are within epsilon of exact values
5. **Edge cases** — empty summaries, boundary ranks (0.0, 1.0), immutability of query
6. **Serde round-trip** — serialization/deserialization preserves query results (feature-gated)
7. **Memory tracking** — uses `stats_alloc` to measure allocations during insertion

## Benchmarks

`benches/zw_benchmark.rs` uses [criterion](https://github.com/bheisler/criterion.rs) to benchmark:

- **Update**: inserting 5000 sequential values with `epsilon = 0.01`
- **Query**: querying 10 quantiles (0.1 through 1.0) on a pre-built summary

Both fixed-size and unbound ZW variants are benchmarked alongside GK01 from the [quantiles](https://github.com/postmates/quantiles) crate.

Run with:
```bash
cargo bench
```
