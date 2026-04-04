VERDICT: NEEDS_REVISION

## Summary Assessment

The design correctly identifies several real bugs and proposes sensible modernization steps, but it has a critical MSRV compatibility gap with the proposed dependencies, underspecifies the cache removal impact on benchmarks, and misses additional bugs and edge cases in the existing code.

## Critical Issues (must fix)

### C1. MSRV 1.65 is incompatible with proposed dependency versions

The design proposes `rust-version = "1.65"` alongside `rand 0.9.x` and `criterion 0.5.x`. Based on known MSRV requirements:
- `rand 0.9` requires Rust 1.71+
- `criterion 0.5` requires Rust 1.67+ (edition 2021 + dependencies)
- `ordered-float 4.x` requires Rust 1.60+

Since these are dev-dependencies, they technically don't affect consumers, but the CI MSRV job proposed in Section 3.5 would run `cargo test --all-features --all --benches`, which compiles dev-dependencies. The MSRV job would fail immediately.

**Fix options:**
1. Raise MSRV to at least 1.71 (to match `rand 0.9`).
2. Keep MSRV 1.65 but split CI: run `cargo check`/`cargo build` at MSRV 1.65 and run tests only on stable. The MSRV then only covers the library itself (no dev-deps). This is a common pattern but should be stated explicitly.
3. Don't upgrade `rand`/`criterion` and keep current versions if a low MSRV is important.

The design must pick an approach and document it. Currently the MSRV claim and the dependency table contradict each other.

### C2. Cache removal causes ~10x regression in the multi-query benchmark with no mitigation plan

The benchmark (`benches/zw_benchmark.rs` lines 52-62) queries 10 different ranks on the same summary without intervening updates. With the current cache, calls 2-10 are essentially free (just a binary search over the already-merged summary). Without the cache, each call recomputes the full merge: sort `s[0]`, clone it, then iteratively merge with each higher level.

For `FixedSizeEpsilonSummary` with n=5000 and epsilon=0.01, the merge involves multiple levels. The 10-query benchmark would go from "1 merge + 9 binary searches" to "10 merges," roughly a 10x slowdown for the query benchmark.

The design says "Users who need repeated queries without intervening updates can store the result themselves," but:
- The merged summary type (`Vec<RankInfo<T>>`) is private and internal. Users cannot store it.
- The design doesn't propose exposing any mechanism for users to avoid repeated merges.
- The benchmark itself would need to be rewritten to not look like a regression.

**Fix:** Either (a) keep an internal cache but make `query(&self)` by using interior mutability (`OnceCell`/`Cell`), (b) provide a public `prepare()` method that returns a queryable snapshot, or (c) explicitly acknowledge the regression and update the benchmark to reflect realistic single-query usage. Option (a) achieves the stated goal (non-mutating reads) while preserving performance.

### C3. `QuantileError` is missing an `InvalidRank` variant

The design proposes `query` returns `Result<T, QuantileError>` but the error enum only has `InvalidN`, `InvalidEpsilon`, and `EmptySummary`. There is no validation or error variant for the rank parameter `r`. Currently, passing `r = -1.0`, `r = 2.0`, or `r = NaN` produces garbage results silently. Since this is a 1.0 release with a proper error type, rank validation should be included.

**Fix:** Add `InvalidRank` to `QuantileError` and validate `0.0 <= r <= 1.0` and `!r.is_nan()` in `query()`.

### C4. NaN epsilon passes constructor validation

Both constructors check `epsilon <= 0.0`, but `f64::NAN <= 0.0` is `false`, so `NaN` passes the guard. This stores `NaN` in `self.epsilon`, causing all subsequent arithmetic to produce `NaN` and all comparisons to behave unpredictably. The design proposes `Result<Self, QuantileError>` for constructors but doesn't mention NaN.

**Fix:** Change the epsilon check to `!(epsilon > 0.0 && epsilon.is_finite())` or equivalent positive assertion, and return `Err(QuantileError::InvalidEpsilon)`.

## Suggestions (nice to have)

### S1. `find_idx` linear scan has a missing early-termination condition

At `lib.rs:373`, the loop condition is:
```rust
while l < s_m.len() && s_m[l].rmin >= rank - epsilon_n {
```

Since the summary is sorted and `rmin` values are non-decreasing, once `l` passes the binary search point, `rmin` will only increase. The condition `rmin >= rank - epsilon_n` will remain true for all subsequent elements, so the loop degenerates to a linear scan of the entire remaining array. It should also check `s_m[l].rmin <= rank + epsilon_n` as a termination condition:

```rust
while l < s_m.len() && s_m[l].rmin <= rank + epsilon_n {
```

This is both a correctness improvement (the function can currently return a value far outside the epsilon-approximate range) and a performance improvement. The design does not mention this.

### S2. Consider `Option<T>` instead of `Result<T, QuantileError>` for `query`

The design proposes `query` returns `Result<T, QuantileError>` (or `Option<T>` for empty-summary). The parenthetical suggests uncertainty. For an empty summary, `Option<T>` is more idiomatic (there is no error -- there is simply no answer). For invalid rank, `Result` is appropriate. Consider splitting: validate rank (returning `Result`), then return `Option<T>` for the actual lookup. Or use `Result` consistently but be clear about it.

### S3. The `compress` function panics on `block_size = 0`

In `compress` at line 315, `block_size` is used as a divisor: `(s0_range as f64) / (block_size as f64)`. If `block_size = 0`, this is a division by zero (producing `Inf` or `NaN` for f64, not a panic, but it produces nonsensical results). While `block_size = 0` is unlikely in normal usage, the design's move to `Result`-based error handling should cover this edge case, or at least document the precondition.

### S4. Serde on `RankInfo` exposes internal representation

The design notes `RankInfo<T>` is internal but needs `Serialize`/`Deserialize` for parent types. This means the serialized format of `FixedSizeEpsilonSummary` and `UnboundEpsilonSummary` will include field names like `rmin`, `rmax`, `cached_s_m`, etc. Since this is a 1.0 release, changing internal field names later becomes a breaking change for anyone who stored serialized data. Consider:
- Adding `#[serde(rename = "...")]` attributes to stabilize field names.
- Skipping the cache field: `#[serde(skip)]` on `cached_s_m`.
- Documenting that the serialized format is not stable, or making it stable.

### S5. `stats_alloc` crate usage in tests may conflict with serde tests

The `#[global_allocator]` attribute in the test module (line 476) sets a custom global allocator. This is process-wide and affects all tests. If future serde round-trip tests are added, they would run under this allocator, which is fine but worth noting. Also, `stats_alloc` is flagged for "check if maintained" -- if it is dropped, the test code needs restructuring.

### S6. Update benchmarks to also cover the "query after update" interleaved pattern

The current benchmarks separate update and query phases entirely. A benchmark that interleaves updates and queries would better reflect real-world streaming usage and would be unaffected by cache removal.

### S7. The design should mention the `quantiles` dev-dependency

The `quantiles` crate (used in benchmarks for GK comparison) is not mentioned in the dependency update table (Section 3.1). It is at version `0.7.1`. Verify it is still maintained and compatible.

## Verified Claims (things you confirmed are correct)

### V1. Boundary array index-out-of-bounds bug (Section 1.1) -- CONFIRMED

The `is_boundary` function can return values up to 30 (the check is `l < 31`). When it returns `x`, the code accesses `self.boundaries[x + x]` (i.e., `self.boundaries[2*x]`). For `x >= 16`, this means index >= 32, which is out of bounds for a `[usize; 32]` array. With `epsilon = 0.1`, `boundaries[16] = floor((2^16 - 1) / 0.1) = 655350`, so inserting 655,350 elements triggers the panic. The design's analysis is correct.

### V2. `find_idx` underflow on empty slice (Section 1.2) -- CONFIRMED

At `lib.rs:361`, `s_m.len() - 1` on an empty slice underflows `usize` to `usize::MAX`. This would cause the binary search to access out-of-bounds memory (panic in debug, UB in release). The fix is correct.

### V3. `PartialOrd`/`Ord` inconsistency (Section 1.3) -- CONFIRMED, but low severity

The implementations at lines 65-75 use different delegation paths (`val.cmp()` vs `val.partial_cmp()`). Since `T: Ord` is required, `partial_cmp` and `cmp` are contractually consistent, so this cannot produce disagreement in practice. The fix to delegate `PartialOrd` through `Some(self.cmp(other))` is good practice but is a code quality issue, not a functional bug.

### V4. Malformed doc comment (Section 1.4) -- CONFIRMED

Line 29 has `//! //!` which breaks the second example block in rustdoc. The first `//!` is the doc comment prefix; the second `//!` becomes literal text, preventing rustdoc from recognizing the code fence.

### V5. Unnecessary `T: Debug` bound (Section 1.5) -- CONFIRMED

`UnboundEpsilonSummary` impl at line 398 requires `T: Clone + Ord + std::fmt::Debug` while `FixedSizeEpsilonSummary` at line 99 only requires `T: Clone + Ord`. The `Debug` bound is not used in any method body and is unnecessary.

### V6. Redundant `epsilon` parameter in `calc_s_m`/`finalize` (Section 2.4) -- CONFIRMED

`UnboundEpsilonSummary` calls `self.s_c.finalize(self.epsilon / 2.0)` and `self.s_c.calc_s_m(self.epsilon / 2.0)`. Since `s_c` was constructed with `FixedSizeEpsilonSummary::new(n, epsilon / 2.0)`, the passed value always equals `self.s_c.epsilon`. The parameter is redundant and can be replaced with `self.epsilon`.

### V7. CI action versions are outdated (Section 3.4) -- CONFIRMED

All GitHub Actions in `.github/workflows/CI.yml` use outdated versions: `actions/checkout@v2`, `actions-rs/toolchain@v1` (archived), `actions-rs/cargo@v1` (archived), `actions-rs/tarpaulin@v0.1`, `codecov/codecov-action@v1`, `actions/upload-artifact@v1`. The `actions-rs` org is no longer maintained. The proposed replacements are appropriate.
