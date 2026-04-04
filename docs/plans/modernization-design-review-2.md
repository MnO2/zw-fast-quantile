VERDICT: NEEDS_REVISION

## Summary Assessment

The v2 design addresses all four critical issues from the first review and promotes the `find_idx` termination bug to a proper fix. However, the `RefCell` interior mutability approach for `query(&self)` is fundamentally incomplete: `query` also mutates `self.s[0]` (sorting and assigning rank values) and `self.s_c`, which are not covered by the `RefCell` wrapper on `cached_s_m`. This blocks the core API improvement. One additional minor issue with the proposed `find_idx` fix is noted.

## Critical Issues (must fix)

### C1. RefCell on `cached_s_m` alone is insufficient for `query(&self)`

The design proposes wrapping only `cached_s_m` in `RefCell` to enable `query(&self)`. However, the current `query` implementation in `FixedSizeEpsilonSummary` (lines 171-184) also:

1. Sorts `self.s[0]` in place (line 172: `self.s[0].sort()`)
2. Mutates `rmin` and `rmax` of every element in `self.s[0]` (lines 173-176)

These mutations are to `self.s`, not to `cached_s_m`. With only `cached_s_m` in a `RefCell`, the code would still require `&mut self` to sort and assign ranks to `s[0]`.

The same problem exists in `UnboundEpsilonSummary::query` (line 444), which calls `self.s_c.calc_s_m(self.epsilon / 2.0)`. The `calc_s_m` method (line 195) takes `&mut self` because it sorts and mutates `self.s[0]`.

**Fix options:**
1. Move the `s[0]` sorting and rank assignment into the cache-population path and wrap the entire computed result (including the processed `s[0]` data) into the `RefCell` cache. The approach would be: clone `s[0]`, sort the clone, assign ranks to the clone, then merge -- all producing a `Vec<RankInfo<T>>` that goes into the `RefCell` cache without mutating `self.s`. This requires modifying the cache population logic, not just wrapping a field.
2. Wrap `s` (or the entire mutable state) in `RefCell` -- but this is heavy-handed and defeats the purpose.
3. Use `Cell<Option<Vec<RankInfo<T>>>>` for the cache (move semantics, no runtime borrow tracking), combined with option (1) above to avoid mutating `s`.

Option 1 is the cleanest. The design should explicitly state that the cache population path must be restructured to avoid mutating `self.s[0]` in place, by sorting a clone instead. This is a small performance cost (one extra clone of `s[0]` during cache population) but it only happens once per update cycle due to the cache.

### C2. The `find_idx` proposed fix has incorrect search direction for small-rank queries

The proposed fix (Section 1.6) changes the linear scan condition to:
```rust
while l < s_m.len() && s_m[l].rmin <= rank + epsilon_n {
```

Consider the case where `rank = 0` and `epsilon_n = 5`. The binary search converges to `l = 0` (the first element where `rmin >= 0`). The fix correctly limits the forward scan to elements with `rmin <= 5`.

However, there is a pre-existing issue that the design's analysis does not address: the binary search finds the leftmost element where `rmin >= rank`, but elements *before* this position (with `rmin < rank`) could be valid epsilon-approximate answers if their `rmax <= rank + epsilon_n`. The binary search skips them because it only searches forward from `l`. For a correct epsilon-approximate query, the scan might need to start slightly before the binary search result.

This is technically a pre-existing correctness issue in the algorithm, not introduced by the v2 fix. However, since Section 1.6 is described as both a "correctness and performance fix" and the design is performing a major version bump to 1.0, this is the right time to also address the binary search starting point. At minimum, the design should acknowledge this limitation or verify that the Zhang-Wang algorithm guarantees the binary search landing point is always a valid starting point for the epsilon-window scan.

**Severity note:** This is borderline between critical and suggestion. I mark it critical because the design explicitly claims Section 1.6 is a "correctness fix" but the analysis is incomplete. If the Zhang-Wang paper guarantees the binary search always lands within the epsilon window (which I believe it does, due to the compress invariant), then documenting this reasoning would resolve the issue.

## Suggestions (nice to have)

### S1. RefCell clone semantics should be documented

`RefCell<T>` implements `Clone` by calling `borrow().clone()`. If the `RefCell` is currently mutably borrowed (which cannot happen in safe code from a single thread unless the borrow is held in a re-entrant call), `clone()` will panic at runtime. Since both structs `#[derive(Clone)]`, the design should note that the `Clone` impl remains safe because `clone()` requires `&self`, which is incompatible with an active `borrow_mut()` on the same instance in safe single-threaded code. This is a non-issue in practice but worth a brief note.

### S2. `update` should use `self.cached_s_m.get_mut()` instead of `borrow_mut()`

The design mentions both `*self.cached_s_m.borrow_mut() = None` and `self.cached_s_m.get_mut().take()` for cache invalidation in `update(&mut self)`. The `get_mut()` approach is strictly better: since `update` has `&mut self`, it can use `get_mut()` which avoids runtime borrow checking entirely. The design should commit to `get_mut()` in `update` and note the reasoning.

### S3. Consider documenting `!Sync` implication more prominently

The design correctly notes that `RefCell` makes the types `!Sync`. However, the current types are also not `Sync` in practice (since `query` takes `&mut self`, you cannot call it from multiple threads without external synchronization). The design should note that this is a change in *kind* but not in *practical impact*: before, the types were `Sync` but `query` was unusable without `&mut` access; after, the types are `!Sync` but this is actually more honest about their thread-safety story.

### S4. Serde `#[serde(skip)]` needs `Default`-compatible initialization

The design says `cached_s_m` will use `#[serde(skip)]` and will be initialized to `RefCell::new(None)` on deserialization. This works automatically only if the field type implements `Default`. `RefCell<Option<Vec<RankInfo<T>>>>` does implement `Default` (producing `RefCell::new(None)`), so this is fine. But the design should explicitly confirm this, since a reader might wonder whether a custom `Deserialize` impl is needed.

### S5. The `QuantileError` should consider implementing `Copy`

The `QuantileError` enum has no heap-allocated fields (all variants are unit variants). It could derive `Copy` in addition to `Clone`, which is more ergonomic for callers. This is a minor API quality improvement.

### S6. Edition 2021 requires Rust 1.56+, not 1.65

The design sets `edition = "2021"` and `rust-version = "1.65"`. Edition 2021 requires Rust 1.56+, so the edition is compatible with the MSRV. However, the design doesn't explain *why* 1.65 specifically was chosen (as opposed to, say, 1.60 or 1.70). It would help to document the rationale -- is 1.65 chosen for a specific language feature or stability milestone? If there is no specific feature requiring 1.65, a lower MSRV would be more inclusive.

## Verified Claims (things I confirmed are correct)

### V1. C1 from review 1 (MSRV split-CI) -- ADEQUATELY ADDRESSED

The split-CI approach is sound: `cargo check` and `cargo build --all-features` at MSRV 1.65 will not compile dev-dependencies (criterion, rand, etc.), so they do not affect the MSRV verification. Tests and benchmarks run only on stable. This is a well-established pattern in the Rust ecosystem.

One important detail: `cargo build --all-features` will compile the `serde` feature, and serde 1.x supports Rust 1.31+, which is within MSRV 1.65. This is correct.

### V2. C2 from review 1 (cache regression) -- PARTIALLY ADDRESSED

The decision to keep the cache is correct and addresses the performance concern. However, the implementation approach (RefCell on `cached_s_m` alone) is incomplete, as described in C1 above.

### V3. C3 from review 1 (missing `InvalidRank`) -- ADEQUATELY ADDRESSED

The `InvalidRank` variant is added, and the rank validation logic `!(r >= 0.0 && r <= 1.0)` correctly rejects NaN (since all NaN comparisons return false, the inner expression evaluates to false, and the negation makes it true, triggering the error).

### V4. C4 from review 1 (NaN epsilon) -- ADEQUATELY ADDRESSED

The epsilon validation `!(epsilon > 0.0 && epsilon.is_finite())` correctly rejects NaN (since `NaN > 0.0` is false), negative values, zero, and infinity. This is correct.

### V5. S1 promoted to 1.6 (`find_idx` termination) -- PARTIALLY VERIFIED

The change from `s_m[l].rmin >= rank - epsilon_n` to `s_m[l].rmin <= rank + epsilon_n` correctly fixes the non-termination problem. With non-decreasing `rmin` values, the old condition was always true past the binary search point, causing a full linear scan. The new condition provides a proper upper bound and the loop will terminate once `rmin` exceeds `rank + epsilon_n`. The correctness of the overall epsilon-approximate query depends on the algorithm's compress invariant, which should be documented (see C2 above).

### V6. S4 from review 1 (serde cache skip) -- ADEQUATELY ADDRESSED

`#[serde(skip)]` on `cached_s_m` with `RefCell<Option<...>>` will deserialize to `Default::default()` which is `RefCell::new(None)`. The cache is lazily repopulated on first `query()` call. The note about serialization format instability is a good addition.

### V7. S7 from review 1 (quantiles dependency) -- ADEQUATELY ADDRESSED

The `quantiles` crate is now included in the dependency update table in Section 3.1.

### V8. The `query` return type change (S2) -- ADEQUATELY ADDRESSED

Using `Result<T, QuantileError>` consistently for `query` is a reasonable choice. Both `EmptySummary` and `InvalidRank` are error conditions that callers should handle, and `Result` communicates this clearly. The design chose a consistent approach rather than mixing `Option` and `Result`.
