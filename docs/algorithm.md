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

Ranks are one-based throughout the crate, matching the paper: a freshly sorted block of `b` elements gets `rmin = rmax = 1, 2, ..., b`. The interval `[rmin, rmax]` captures the uncertainty introduced by compression. As summaries are merged and compressed, these bounds widen but remain within the epsilon guarantee.

## Operations

### Update (Insertion)

New elements are inserted into a level-0 buffer (`s[0]`). When this buffer reaches the block size `b`:

1. **Sort** the buffer and assign exact one-based ranks (`rmin = rmax = position`).
2. **Compress** the sorted buffer to roughly `b/2` entries by selecting a subset that preserves the rank coverage.
3. **Propagate** the compressed summary upward through levels: if level `k` is empty, store there; otherwise, take the summary stored at level `k`, merge it with the incoming one, compress the result, and continue to level `k+1`.

This cascading merge is analogous to binary counting — level `k` is "full" roughly every `2^k` blocks, giving logarithmic levels.

### Merge

Merging two sorted summaries `S_a` and `S_b` produces a combined summary `S_m` where each element's rank bounds are updated to account for elements from the other summary. For an element `e` from `S_a`, let `pred` be the largest element of `S_b` below `e` and `succ` the smallest element of `S_b` above it:

- **rmin**: `rmin_a(e) + rmin_b(pred)`, or just `rmin_a(e)` when there is no predecessor.
- **rmax**: `rmax_a(e) + rmax_b(succ) - 1`, or `rmax_a(e) + rmax_b(pred)` when there is no successor.

These are the one-based formulas from the paper; they only hold because sorted blocks are ranked from 1. The merge proceeds like a standard sorted-list merge, advancing through both inputs in order and moving elements rather than cloning them whenever the input is owned.

### Compress

Compression reduces a summary of size `n` to at most `block_size + 1` entries while preserving the epsilon guarantee. It works by:

1. Computing `s0_range`, the maximum `rmax` in the summary.
2. Dividing the rank space `[0, s0_range]` into `block_size` equal intervals.
3. For each interval boundary, selecting the first element whose `rmax` reaches or exceeds that boundary.

This ensures uniform coverage of the rank space. Every input to compress must already satisfy the precision condition `2 * epsilon * range >= max(rmax - rmin)`; the constructor's block sizing guarantees it, and debug builds assert it.

### Query

Querying rank `r` on a summary of `n` elements:

1. Compute the one-based target rank: `rank = min(floor(r * n) + 1, n)`.
2. Compute the error window: `epsilon_n = floor(epsilon * n)`.
3. Build a merged summary `S_m` from all levels (cached after the first build, invalidated on update).
4. Binary search for the first element with `rmin >= rank`, then scan forward while `rmin <= rank + epsilon_n` for an element with `rmax <= rank + epsilon_n`.
5. If none is found, scan backward over the elements with `rmin >= rank - epsilon_n`; any of them with `rmax <= rank + epsilon_n` also satisfies the guarantee.
6. As a last resort, return the element whose rank interval is centred closest to `rank`. This never happens for summaries built by the crate, but it keeps a query from silently returning the maximum element.

## FixedSizeEpsilonSummary

When the stream size `n` is known upfront and `epsilon * n > 8`:

- **Block size** `b` = `floor(log2(epsilon * n) / epsilon)`
- **Number of levels** = bit width of `floor(n / b)`, plus one for level 0. Full blocks propagate through the levels like binary carries, so the level count is the number of bits needed to count them.

When `epsilon * n <= 8` the summary switches to **exact mode**: a single level with `b = n + 1` that never fills, so every element is kept and queries are exact. The ZW block size is only meaningful once `log2(epsilon * n)` is comfortably above 1; below that the derived blocks become too small (down to a single element) and neither the precision condition nor the rank-error guarantee can hold. Exact mode costs at most `8 / epsilon` stored elements, the same order as a single block.

## UnboundEpsilonSummary

When the stream size is unknown, the algorithm uses a doubling strategy:

1. Maintain a "current" `FixedSizeEpsilonSummary` (`s_c`) sized for the current sub-stream.
2. Sub-stream boundaries lie at positions `floor((2^x - 1) / epsilon)` for `x = 1, 2, ...`; the next boundary is computed on demand and tracked with a counter, so there is no table to exhaust.
3. When the count reaches boundary `x`:
   - Finalize `s_c` (compress all its levels into a single summary).
   - Push `s_c` onto a list of completed summaries.
   - Create a new `s_c` sized for the interval up to boundary `x + 1`.
4. At query time, merge all completed summaries with the current `s_c`.

Each sub-stream summary is built with `epsilon / 2`, so the merged result stays within `epsilon` overall. Short sub-streams (`epsilon * n_i <= 8`) use exact mode like any other fixed-size summary.

## Query Caching

Both summary types cache the merged summary (`cached_s_m`) using `RefCell<Option<...>>`. This allows `query()` to take `&self` (immutable reference) while lazily computing and caching the merged result. The cache is invalidated on every `update()` call.

## Space Complexity

- **FixedSizeEpsilonSummary**: `O((1/epsilon) * log(epsilon * n))` elements stored across all levels.
- **UnboundEpsilonSummary**: `O((1/epsilon) * log^2(epsilon * n))` due to the list of completed summaries.

In practice, for `epsilon = 0.01` and `n = 1,000,000`, the summaries store a few thousand entries rather than the full million.
