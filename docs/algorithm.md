# Algorithm Internals

This document explains how the Zhang-Wang (ZW) fast approximate quantile algorithm works, as implemented in this crate.

## Background

Computing exact quantiles over a data stream requires storing all elements and sorting them — `O(n)` space and `O(n log n)` time. For large or unbounded streams, this is impractical. Approximate quantile algorithms trade a bounded error tolerance (`epsilon`) for dramatically reduced memory usage.

The ZW algorithm, described in [Zhang & Wang, SSDBM 2007](http://web.cs.ucla.edu/~weiwang/paper/SSDBM07_2.pdf), achieves `O((1/epsilon) * log(epsilon * n))` space with `O(log(1/epsilon))` amortized update time.

## Core Concepts

### Rank and Epsilon Guarantee

For a stream of `n` elements, querying rank `r` (where `0.0 <= r <= 1.0`) should return the element at position `floor(r * n)` in the sorted order. The epsilon guarantee means the returned element's true rank is within `[r - epsilon, r + epsilon]`, i.e., within `epsilon * n` positions of the exact answer.

### RankInfo

Each element in a summary is stored as a `RankInfo<T>` tuple:

```
RankInfo { val: T, rmin: i64, rmax: i64 }
```

- `val` — the data element
- `rmin` — the minimum possible rank of this element
- `rmax` — the maximum possible rank of this element

The interval `[rmin, rmax]` captures the uncertainty introduced by compression. As summaries are merged and compressed, these bounds widen but remain within the epsilon guarantee.

## Operations

### Update (Insertion)

New elements are inserted into a level-0 buffer (`s[0]`). When this buffer reaches the block size `b`:

1. **Sort** the buffer and assign exact ranks (`rmin = rmax = position`).
2. **Compress** the sorted buffer to roughly `b/2` entries by selecting a subset that preserves the rank coverage.
3. **Propagate** the compressed summary upward through levels: if level `k` is empty, store there; otherwise, merge with level `k`, compress the result, clear level `k`, and continue to level `k+1`.

This cascading merge is analogous to binary counting — level `k` is "full" roughly every `2^k` blocks, giving logarithmic levels.

### Merge

Merging two sorted summaries `S_a` and `S_b` produces a combined summary `S_m` where each element's rank bounds are updated to account for elements from the other summary:

- **rmin**: the element's own `rmin` plus the `rmin` of the predecessor element in the other summary (or 0 if it comes first).
- **rmax**: the element's own `rmax` plus the `rmax` of the next element in the other summary, minus 1 (to avoid double-counting).

The merge proceeds like a standard sorted-list merge, advancing through both inputs in order.

### Compress

Compression reduces a summary of size `n` to at most `block_size` entries while preserving the epsilon guarantee. It works by:

1. Computing `s0_range`, the maximum `rmax` in the summary.
2. Dividing the rank space `[0, s0_range]` into `block_size` equal intervals.
3. For each interval boundary, selecting the first element whose `rmax` reaches or exceeds that boundary.

This ensures uniform coverage of the rank space. The precision condition `2 * epsilon * range >= max(rmax - rmin)` is asserted to verify the error bound is maintained.

### Query

Querying rank `r`:

1. Compute the target rank: `rank = floor(r * n)`.
2. Compute the error window: `epsilon_n = floor(epsilon * n)`.
3. Build a merged summary `S_m` from all levels (cached after first build, invalidated on update).
4. Binary search for the first element with `rmin >= rank`.
5. Scan forward to find an element with `rmax <= rank + epsilon_n` — this element's true rank is guaranteed to be within the tolerance.
6. If no such element is found, return the last element.

## FixedSizeEpsilonSummary

When the stream size `n` is known upfront:

- **Number of levels** = `ceil(log2(epsilon * n))` (or 1 when `epsilon * n <= 1`)
- **Block size** `b` = `floor(log2(epsilon * n) / epsilon)` (or `n + 1` when only 1 level)

Knowing `n` allows pre-allocating the exact number of levels and optimal block size.

## UnboundEpsilonSummary

When the stream size is unknown, the algorithm uses a doubling strategy:

1. Maintain a "current" `FixedSizeEpsilonSummary` (`s_c`) sized for the next power-of-two boundary.
2. Track **boundaries** at positions `floor((2^i - 1) / epsilon)` for `i = 0, 1, 2, ...`
3. When the count hits a boundary at index `x`:
   - Finalize `s_c` (compress all its levels into a single summary).
   - Push `s_c` onto a stack of completed summaries.
   - Create a new `s_c` sized for the interval up to boundary `2x`.
4. At query time, merge all completed summaries with the current `s_c`.

This ensures each completed summary covers a geometrically growing interval, maintaining the epsilon/2 guarantee per summary and epsilon overall after merging.

## Query Caching

Both summary types cache the merged summary (`cached_s_m`) using `RefCell<Option<...>>`. This allows `query()` to take `&self` (immutable reference) while lazily computing and caching the merged result. The cache is invalidated on every `update()` call.

## Space Complexity

- **FixedSizeEpsilonSummary**: `O((1/epsilon) * log(epsilon * n))` elements stored across all levels.
- **UnboundEpsilonSummary**: `O((1/epsilon) * log^2(epsilon * n))` due to the stack of completed summaries.

In practice, for `epsilon = 0.01` and `n = 1,000,000`, the summaries store a few thousand entries rather than the full million.
