# zw-fast-quantile Modernization Design

**Date:** 2026-04-04
**Version:** 0.2.3 → 1.0.0 (major version bump, breaking changes)

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

## 2. API Improvements

### 2.1 `query(&mut self)` → `query(&self)`
- Remove the internal cache (`cached_s_m`) and recompute the merged summary on each `query()` call.
- Users who need repeated queries without intervening updates can store the result themselves.
- This makes the API unsurprising — reads are non-mutating.

### 2.2 Return `Result` instead of panicking
- Constructors return `Result<Self, QuantileError>`.
- `query` returns `Result<T, QuantileError>` (or `Option<T>` for empty-summary case).
- Error type:
  ```rust
  #[derive(Debug, Clone, PartialEq)]
  pub enum QuantileError {
      InvalidN,
      InvalidEpsilon,
      EmptySummary,
  }
  impl std::fmt::Display for QuantileError { ... }
  impl std::error::Error for QuantileError {}
  ```

### 2.3 `#[must_use]` attributes
- Add to `query()`, `size()`, and constructors.

### 2.4 Clean up internal methods
- Remove redundant `epsilon` parameter from `calc_s_m` and `finalize` — use `self.epsilon`.

### 2.5 Derive `Debug` on public types
- Both `FixedSizeEpsilonSummary` and `UnboundEpsilonSummary` should derive `Debug`.

## 3. Dependency & CI Modernization

### 3.1 Dev dependency updates
| Crate | Current | Target |
|-------|---------|--------|
| criterion | 0.3.5 | 0.5.x |
| rand | 0.8.4 | 0.9.x |
| rand_distr | 0.4.2 | 0.5.x |
| rand_pcg | 0.3.1 | check/remove (may be merged into rand 0.9) |
| ordered-float | 2.8.0 | 4.x |
| stats_alloc | 0.1.8 | check if maintained |

### 3.2 New dependency
- `serde = { version = "1", features = ["derive"], optional = true }` in `[dependencies]`

### 3.3 Cargo.toml metadata
- `version = "1.0.0"`
- `edition = "2021"`
- `rust-version = "1.65"`

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
- **Clippy** lint job
- **MSRV** job testing against Rust 1.65

## 4. Serde Support

- Add `[features] serde = ["dep:serde"]` to Cargo.toml.
- Derive `Serialize`/`Deserialize` behind `#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]` on:
  - `RankInfo<T>` (internal, but needed for parent types)
  - `FixedSizeEpsilonSummary<T>`
  - `UnboundEpsilonSummary<T>`
- Conditional serde bounds on `T` when feature is active.

## 5. README & Documentation

- Update examples to use `Result`-returning API.
- Update installation to `zw-fast-quantile = "1.0"`.
- Mention the `serde` feature.
- Remove the Rust 2015 `extern crate` note.
- Update or remove stale benchmark numbers.
- Fix doc comment in lib.rs.
