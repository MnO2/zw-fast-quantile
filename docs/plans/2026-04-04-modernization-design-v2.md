# zw-fast-quantile Modernization Design (v2)

**Date:** 2026-04-04
**Version:** 0.2.3 -> 1.0.0 (major version bump, breaking changes)
**Previous revision:** `2026-04-04-modernization-design.md`

## Changes from previous version

This revision addresses all critical issues and several suggestions from `modernization-design-review-1.md`:

- **C1 (MSRV incompatibility):** Adopted split-CI approach. MSRV 1.65 applies to the library crate only. CI runs `cargo check`/`cargo build` (no dev-deps) at MSRV 1.65, and full tests/benches on stable only. This is now documented explicitly in Section 3.3 and 3.5.
- **C2 (Cache removal regression):** Reversed the cache removal decision. The cache is kept, but the `cached_s_m` field is now wrapped in `RefCell<Option<Vec<RankInfo<T>>>>` to enable `query(&self)` via interior mutability. This preserves multi-query performance while making reads non-mutating from the caller's perspective.
- **C3 (Missing `InvalidRank`):** Added `InvalidRank` variant to `QuantileError`. `query()` now validates `0.0 <= r <= 1.0` and `!r.is_nan()`.
- **C4 (NaN epsilon):** Changed epsilon validation to `!(epsilon > 0.0 && epsilon.is_finite())`, which correctly rejects NaN and infinity.
- **S1 (find_idx termination) -- promoted to bug fix:** Added as Section 1.6. The `find_idx` linear scan condition is inverted: it should be `s_m[l].rmin <= rank + epsilon_n` (upper bound termination), not `s_m[l].rmin >= rank - epsilon_n` (which never terminates early). This is both a correctness and performance fix.
- **S2 (query return type):** `query` now returns `Result<T, QuantileError>` consistently. Both `EmptySummary` and `InvalidRank` are error conditions conveyed through `Result`.
- **S4 (Serde cache skip):** Added `#[serde(skip)]` on `cached_s_m` fields. Added note that serialization format is not considered stable in 1.0.
- **S7 (quantiles dependency):** Added `quantiles` to the dependency update table in Section 3.1.

---

## 1. Bug Fixes & Correctness

### 1.1 Index out-of-bounds in `UnboundEpsilonSummary::update`
- **Location:** `src/lib.rs:431` -- `self.boundaries[x + x]`
- **Problem:** `is_boundary` can return values up to 30, so `x + x` can reach 60, far past the array length of 32. For `epsilon=0.1`, this triggers at ~655K elements.
- **Fix:** Guard against `x + x >= 32` or restructure the boundary computation to avoid fixed-size arrays.

### 1.2 Underflow panic in `find_idx`
- **Location:** `src/lib.rs:362` -- `s_m.len() - 1` underflows on empty slice.
- **Fix:** Early return `None` when `s_m.is_empty()`.

### 1.3 `PartialOrd`/`Ord` inconsistency on `RankInfo`
- **Location:** Lines 65-75.
- **Problem:** `PartialOrd` delegates to `val.partial_cmp()`, `Ord` delegates to `val.cmp()`. These can disagree.
- **Fix:** `PartialOrd` should return `Some(self.cmp(other))`.

### 1.4 Malformed doc comment
- **Location:** Line 29 -- `//! //!` produces broken rustdoc.
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

## 2. API Improvements

### 2.1 `query(&mut self)` -> `query(&self)` via interior mutability
- Keep the internal cache (`cached_s_m`) to preserve multi-query performance.
- Change the cache field type from `Option<Vec<RankInfo<T>>>` to `RefCell<Option<Vec<RankInfo<T>>>>` (using `std::cell::RefCell`, available since Rust 1.0).
- On `query(&self)`, check `self.cached_s_m.borrow().is_none()` and lazily populate via `self.cached_s_m.borrow_mut()`.
- On `update(&mut self)`, clear the cache: `*self.cached_s_m.borrow_mut() = None;` (or equivalently `self.cached_s_m.get_mut().take()` since `update` has `&mut self`).
- This makes the API unsurprising (reads are non-mutating from the caller's perspective) while preserving the "1 merge + 9 binary searches" performance for the multi-query benchmark.
- **Note:** `RefCell` is not `Sync`, so `FixedSizeEpsilonSummary` and `UnboundEpsilonSummary` will not be `Sync`. This matches the current behavior (they are not `Sync` today since `query` requires `&mut self`). Users who need thread-safe access should wrap in their own `Mutex`.

### 2.2 Return `Result` instead of panicking
- Constructors return `Result<Self, QuantileError>`.
- `query` returns `Result<T, QuantileError>`.
- Epsilon validation uses `!(epsilon > 0.0 && epsilon.is_finite())` to correctly reject NaN, negative values, zero, and infinity.
- Rank validation in `query` checks `!(r >= 0.0 && r <= 1.0)` (this correctly rejects NaN since all comparisons with NaN are false).
- Error type:
  ```rust
  #[derive(Debug, Clone, PartialEq)]
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
- Remove redundant `epsilon` parameter from `calc_s_m` and `finalize` -- use `self.epsilon`.

### 2.5 Derive `Debug` on public types
- Both `FixedSizeEpsilonSummary` and `UnboundEpsilonSummary` should derive `Debug`.
- The `RefCell` cache field should be excluded from `Debug` output (or shown as opaque) via a manual `Debug` impl or by using `#[derive(Debug)]` which will print the `RefCell` contents -- acceptable since `RankInfo` already derives `Debug`.

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
- **Cache field handling:** The `cached_s_m` field (type `RefCell<Option<Vec<RankInfo<T>>>>`) must be annotated with `#[serde(skip)]`. On deserialization, the cache will be initialized to `RefCell::new(None)` and lazily repopulated on the first `query()` call. This is correct because the cache is a derived value.
- **Serialization format stability:** The serialized representation exposes internal field names (`rmin`, `rmax`, `s`, `b`, `level`, etc.). The 1.0 release does **not** guarantee serialization format stability. This will be documented in the crate-level docs:
  > The serialized format of summary types (when the `serde` feature is enabled) is not considered part of the stable API. Field names and structure may change in future minor versions. Do not rely on the serialized format for long-term storage without pinning to a specific crate version.

  If format stability becomes needed in the future, we can add `#[serde(rename = "...")]` attributes in a subsequent release.

## 5. README & Documentation

- Update examples to use `Result`-returning API.
- Update installation to `zw-fast-quantile = "1.0"`.
- Mention the `serde` feature.
- Remove the Rust 2015 `extern crate` note.
- Update or remove stale benchmark numbers.
- Fix doc comment in lib.rs (Section 1.4).
- Document MSRV policy: "The library crate supports Rust 1.65+. Dev-dependencies (used for tests and benchmarks) may require a newer Rust version."
