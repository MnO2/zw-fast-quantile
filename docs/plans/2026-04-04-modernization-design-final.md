# zw-fast-quantile Modernization Design (Final)

**Date:** 2026-04-04
**Version:** 0.2.3 -> 1.0.0 (major version bump, breaking changes)
**Previous revisions:** `2026-04-04-modernization-design.md`, `2026-04-04-modernization-design-v2.md`
**Reviews:** `modernization-design-review-1.md`, `modernization-design-review-2.md`

## Changes from v2

This final revision addresses the two remaining issues from review 2:

- **C1 (RefCell insufficient alone):** The cache population path is restructured to clone and sort `s[0]` instead of sorting it in place. This means `query(&self)` never mutates `self.s` — all mutable state is confined to the `RefCell`-wrapped cache. Detailed in Section 2.1.
- **C2 (find_idx correctness documentation):** Added explicit reasoning in Section 1.6 about why the binary search landing point is always a valid starting point, based on the compress invariant. Also added `Copy` to `QuantileError` (reviewer suggestion S5).

---

## 1. Bug Fixes & Correctness

### 1.1 Index out-of-bounds in `UnboundEpsilonSummary::update`
- **Location:** `src/lib.rs:431` — `self.boundaries[x + x]`
- **Problem:** `is_boundary` can return values up to 30, so `x + x` can reach 60, far past the array length of 32. For `epsilon=0.1`, this triggers at ~655K elements.
- **Fix:** Guard against `x + x >= 32` or restructure the boundary computation to avoid fixed-size arrays.

### 1.2 Underflow panic in `find_idx`
- **Location:** `src/lib.rs:362` — `s_m.len() - 1` underflows on empty slice.
- **Fix:** Early return `None` when `s_m.is_empty()`.

### 1.3 `PartialOrd`/`Ord` inconsistency on `RankInfo`
- **Location:** Lines 65-75.
- **Problem:** `PartialOrd` delegates to `val.partial_cmp()`, `Ord` delegates to `val.cmp()`. These can disagree.
- **Fix:** `PartialOrd` should return `Some(self.cmp(other))`.

### 1.4 Malformed doc comment
- **Location:** Line 29 — `//! //!` produces broken rustdoc.
- **Fix:** Remove the doubled prefix so the second example renders correctly.

### 1.5 Unnecessary `T: Debug` bound
- **Location:** `UnboundEpsilonSummary` impl (line 399) requires `T: Debug` but `FixedSizeEpsilonSummary` does not.
- **Fix:** Remove the `Debug` bound from `UnboundEpsilonSummary`.

### 1.6 `find_idx` linear scan has inverted termination condition
- **Location:** `src/lib.rs:373`
- **Problem:** The loop condition is `s_m[l].rmin >= rank - epsilon_n`. Since the summary is sorted and `rmin` values are non-decreasing, once `l` passes the binary search point, `rmin` only increases. The condition `rmin >= rank - epsilon_n` remains true for all subsequent elements, so the loop degenerates to a linear scan of the entire remaining array. Additionally, the function can return a value far outside the epsilon-approximate range.
- **Fix:** Change the termination condition to `s_m[l].rmin <= rank + epsilon_n`:
  ```rust
  while l < s_m.len() && s_m[l].rmin <= rank + epsilon_n {
      if s_m[l].rmax <= rank + epsilon_n {
          return Some(s_m[l].val.clone());
      }
      l += 1;
  }
  ```
  This is both a correctness fix (prevents returning values outside the epsilon-approximate range) and a performance improvement (early termination instead of linear scan).
- **Correctness note:** The binary search finds the leftmost element where `rmin >= rank`. Elements *before* this position have `rmin < rank`, meaning their minimum possible rank is below the target. The compress invariant guarantees that the spacing between consecutive summary elements is at most `2 * epsilon * N`, so there is always an element within the epsilon window of any valid rank. Combined with the binary search targeting `rmin >= rank`, the first candidate at or after the landing point is guaranteed to be within the epsilon-approximate range. Therefore, starting the linear scan at the binary search result is correct and no backward scan is needed.

## 2. API Improvements

### 2.1 `query(&mut self)` -> `query(&self)` via interior mutability
- Keep the internal cache (`cached_s_m`) to preserve multi-query performance.
- Change the cache field type from `Option<Vec<RankInfo<T>>>` to `RefCell<Option<Vec<RankInfo<T>>>>`.
- **Critical implementation detail:** The current `query` implementation sorts `self.s[0]` in place and assigns `rmin`/`rmax` values to elements in `s[0]` (lines 172-176). This mutates `self.s`, which would still require `&mut self`. **The fix:** The cache population path must clone `s[0]`, sort the clone, assign ranks to the clone, then merge — producing the `Vec<RankInfo<T>>` that goes into the `RefCell` cache without ever mutating `self.s`. This is a small additional cost (one clone of `s[0]` per cache miss), but it only happens once per update cycle due to the cache.
  ```rust
  // In query(&self):
  if self.cached_s_m.borrow().is_none() {
      let mut s0_clone = self.s[0].clone();
      s0_clone.sort();
      for (i, r) in s0_clone.iter_mut().enumerate() {
          r.rmin = i as i64;
          r.rmax = i as i64;
      }
      let mut s_m = s0_clone;
      for i in 1..self.level {
          s_m = merge(s_m, &self.s[i]);
      }
      *self.cached_s_m.borrow_mut() = Some(s_m);
  }
  ```
- The same restructuring applies to `UnboundEpsilonSummary::query`, which calls `self.s_c.calc_s_m()`. The `calc_s_m` method must be refactored to not take `&mut self` — it should clone and sort internally rather than mutating `self.s[0]`.
- On `update(&mut self)`, clear the cache using `self.cached_s_m.get_mut().take()` (avoids runtime borrow tracking since `update` has `&mut self`).
- **Thread safety:** `RefCell` makes the types `!Sync`. This matches the current practical behavior — currently the types are technically `Sync` but `query` requires `&mut self`, making concurrent access impossible anyway. The `RefCell` approach is more honest about the thread-safety story. Users who need `Sync` should wrap in their own `Mutex`.

### 2.2 Return `Result` instead of panicking
- Constructors return `Result<Self, QuantileError>`.
- `query` returns `Result<T, QuantileError>`.
- Epsilon validation uses `!(epsilon > 0.0 && epsilon.is_finite())` to correctly reject NaN, negative values, zero, and infinity.
- Rank validation in `query` checks `!(r >= 0.0 && r <= 1.0)` (correctly rejects NaN since all NaN comparisons are false).
- Error type:
  ```rust
  #[derive(Debug, Clone, Copy, PartialEq, Eq)]
  pub enum QuantileError {
      InvalidN,
      InvalidEpsilon,
      InvalidRank,
      EmptySummary,
  }
  impl std::fmt::Display for QuantileError { ... }
  impl std::error::Error for QuantileError {}
  ```

### 2.3 `#[must_use]` attributes
- Add to `query()`, `size()`, and constructors.

### 2.4 Clean up internal methods
- Remove redundant `epsilon` parameter from `calc_s_m` and `finalize` — use `self.epsilon`.
- Refactor `calc_s_m` to not require `&mut self` (clone and sort instead of in-place mutation — same approach as Section 2.1).

### 2.5 Derive `Debug` on public types
- Both `FixedSizeEpsilonSummary` and `UnboundEpsilonSummary` should derive `Debug`.

## 3. Dependency & CI Modernization

### 3.1 Dev dependency updates
| Crate | Current | Target | Notes |
|-------|---------|--------|-------|
| criterion | 0.3.5 | 0.5.x | Requires Rust 1.67+; only used at test/bench time |
| rand | 0.8.4 | 0.9.x | Requires Rust 1.71+; only used at test time |
| rand_distr | 0.4.2 | 0.5.x | Only used at test time |
| rand_pcg | 0.3.1 | check/remove | May be merged into rand 0.9; only used at test time |
| ordered-float | 2.8.0 | 4.x | Requires Rust 1.60+; only used at test time |
| quantiles | 0.7.1 | check latest | Used in benchmarks for GK comparison; verify maintained |
| stats_alloc | 0.1.8 | check if maintained | Used for memory tracking in tests |

**Important:** All items in this table are dev-dependencies. They do not affect the library's MSRV. They are only compiled when running tests and benchmarks, which run on stable Rust in CI (see Section 3.5).

### 3.2 New dependency
- `serde = { version = "1", features = ["derive"], optional = true }` in `[dependencies]`
- This is the only non-dev dependency. serde 1.x supports Rust 1.31+, well within our MSRV 1.65.

### 3.3 Cargo.toml metadata
- `version = "1.0.0"`
- `edition = "2021"`
- `rust-version = "1.65"`

The `rust-version` field declares the MSRV for the library crate itself (i.e., `[dependencies]` only). Dev-dependencies may require a higher Rust version; this is standard practice and does not affect downstream consumers. See Section 3.5 for how CI enforces this split.

### 3.4 CI updates
| Component | Current | Target |
|-----------|---------|--------|
| actions/checkout | v2 | v4 |
| actions-rs/toolchain | v1 | dtolnay/rust-toolchain@stable |
| actions-rs/cargo | v1 | direct `run: cargo ...` |
| actions-rs/tarpaulin | v0.1 | cargo install or maintained action |
| actions/upload-artifact | v1 | v4 |
| codecov/codecov-action | v1 | v4 |

### 3.5 New CI jobs

- **Clippy** lint job (runs on stable).
- **MSRV** job:
  - Uses `dtolnay/rust-toolchain` with `toolchain: "1.65"`.
  - Runs `cargo check` and `cargo build --all-features` only (no `--all-targets`, no `--benches`, no `cargo test`).
  - This verifies the library compiles at MSRV 1.65 without pulling in dev-dependencies, which may require a higher Rust version.
  - **Does not** run tests or benchmarks at MSRV. Tests and benchmarks use dev-dependencies (criterion 0.5, rand 0.9, etc.) that require Rust 1.71+.
- **Test** job (runs on stable):
  - Runs `cargo test --all-features`.
  - Runs `cargo bench --no-run` (compile-check only for benchmarks).

This split-CI approach is a common pattern in the Rust ecosystem (e.g., used by serde, tokio). The MSRV guarantee applies to consumers of the library; contributors and CI use stable for the full test suite.

## 4. Serde Support

- Add `[features] serde = ["dep:serde"]` to Cargo.toml.
- Derive `Serialize`/`Deserialize` behind `#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]` on:
  - `RankInfo<T>` (internal, but needed for parent types)
  - `FixedSizeEpsilonSummary<T>`
  - `UnboundEpsilonSummary<T>`
- Conditional serde bounds on `T` when feature is active.
- **Cache field handling:** The `cached_s_m` field (type `RefCell<Option<Vec<RankInfo<T>>>>`) must be annotated with `#[serde(skip)]`. On deserialization, the cache will be initialized to `RefCell::new(None)` (via `Default`) and lazily repopulated on the first `query()` call. This works automatically because `RefCell<Option<Vec<RankInfo<T>>>>` implements `Default`.
- **Serialization format stability:** The serialized representation exposes internal field names (`rmin`, `rmax`, `s`, `b`, `level`, etc.). The 1.0 release does **not** guarantee serialization format stability. This will be documented in the crate-level docs:
  > The serialized format of summary types (when the `serde` feature is enabled) is not considered part of the stable API. Field names and structure may change in future minor versions. Do not rely on the serialized format for long-term storage without pinning to a specific crate version.

## 5. README & Documentation

- Update examples to use `Result`-returning API.
- Update installation to `zw-fast-quantile = "1.0"`.
- Mention the `serde` feature.
- Remove the Rust 2015 `extern crate` note.
- Update or remove stale benchmark numbers.
- Fix doc comment in lib.rs (Section 1.4).
- Document MSRV policy: "The library crate supports Rust 1.65+. Dev-dependencies (used for tests and benchmarks) may require a newer Rust version."
