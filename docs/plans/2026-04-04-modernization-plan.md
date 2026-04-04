# Plan: Modernize zw-fast-quantile to 1.0.0

**Goal**: Fix all known bugs, modernize the API (Result-based errors, immutable query, serde support), update dependencies and CI, and release as 1.0.0.
**Architecture**: Single-file Rust library (`src/lib.rs`) with `FixedSizeEpsilonSummary<T>` and `UnboundEpsilonSummary<T>`. Zero runtime deps (serde optional). Interior mutability via `RefCell` for query cache.
**Tech Stack**: Rust (edition 2021, MSRV 1.65), GitHub Actions CI, criterion benchmarks, serde (optional)
**Design Doc**: `docs/plans/2026-04-04-modernization-design-final.md`

## Task Dependencies

| Group | Steps | Can Parallelize | Description |
|-------|-------|-----------------|-------------|
| 1 | Steps 1-2 | Yes | Error type + doc comment fix |
| 2 | Steps 3-5 | Yes | Bug fixes (find_idx, PartialOrd, Debug bound) |
| 3 | Step 6 | No (depends on 1) | Constructors return Result |
| 4 | Step 7 | No (depends on 2,3) | RefCell cache + query(&self) returning Result |
| 5 | Step 8 | No (depends on 4) | Boundary array OOB fix |
| 6 | Step 9 | No (depends on 5) | Internal method cleanup |
| 7 | Step 10 | No (depends on 6) | #[must_use] + Debug derives |
| 8 | Step 11 | No (depends on 7) | Serde support |
| 9 | Step 12 | No (depends on 8) | Update dev-dependencies |
| 10 | Step 13 | No (depends on 9) | Update CI |
| 11 | Step 14 | No (depends on 10) | Update README + doc comments |
| 12 | Step 15 | No (depends on 11) | Version bump + final verification |

---

## Step 1: Add `QuantileError` type

**File**: `src/lib.rs`

### 1a. Write failing test
```rust
// At the top of mod tests:
#[test]
fn test_error_display_and_debug() {
    let e = QuantileError::InvalidN;
    assert_eq!(format!("{}", e), "n must be greater than 0");
    let e = QuantileError::InvalidEpsilon;
    assert_eq!(format!("{}", e), "epsilon must be positive and finite");
    let e = QuantileError::InvalidRank;
    assert_eq!(format!("{}", e), "rank must be between 0.0 and 1.0");
    let e = QuantileError::EmptySummary;
    assert_eq!(format!("{}", e), "cannot query an empty summary");

    // Verify Copy + Clone + PartialEq + Eq
    let e2 = e;
    assert_eq!(e, e2);
}
```

### 1b. Run test to verify it fails
```bash
cargo test test_error_display_and_debug
```

### 1c. Write implementation
```rust
use std::fmt;

/// Errors that can occur when constructing or querying a quantile summary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QuantileError {
    /// `n` must be greater than 0.
    InvalidN,
    /// `epsilon` must be positive and finite.
    InvalidEpsilon,
    /// Rank must be between 0.0 and 1.0 (inclusive), and not NaN.
    InvalidRank,
    /// Cannot query an empty summary (no values have been inserted).
    EmptySummary,
}

impl fmt::Display for QuantileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            QuantileError::InvalidN => write!(f, "n must be greater than 0"),
            QuantileError::InvalidEpsilon => write!(f, "epsilon must be positive and finite"),
            QuantileError::InvalidRank => write!(f, "rank must be between 0.0 and 1.0"),
            QuantileError::EmptySummary => write!(f, "cannot query an empty summary"),
        }
    }
}

impl std::error::Error for QuantileError {}
```

### 1d. Run test to verify it passes
```bash
cargo test test_error_display_and_debug
```

### 1e. Commit
```bash
git commit -m "Add QuantileError type with Display, Error, Copy, Clone, PartialEq, Eq"
```

---

## Step 2: Fix malformed doc comment (line 29)

**File**: `src/lib.rs`

### 2a. No test needed — this is a doc fix

### 2b. Write fix
Change line 29 from:
```rust
//! //! ```rust
```
to:
```rust
//!
//! ```rust
```

### 2c. Verify doc-tests still pass
```bash
cargo test --doc
```

### 2d. Commit
```bash
git commit -m "Fix malformed doc comment that broke rustdoc rendering"
```

---

## Step 3: Fix `find_idx` underflow on empty slice + termination condition

**File**: `src/lib.rs`

### 3a. Write failing test
```rust
#[test]
fn test_find_idx_empty_slice() {
    let empty: Vec<RankInfo<i32>> = vec![];
    assert!(find_idx(&empty, 0, 1).is_none());
}
```

### 3b. Run test to verify it fails
```bash
cargo test test_find_idx_empty_slice
```

### 3c. Write implementation
Replace the `find_idx` function:
```rust
fn find_idx<T: Clone + Ord>(s_m: &[RankInfo<T>], rank: i64, epsilon_n: i64) -> Option<T> {
    if s_m.is_empty() {
        return None;
    }

    let mut l = 0usize;
    let mut r = s_m.len() - 1;
    while l < r {
        let m = l + (r - l) / 2;
        if s_m[m].rmin < rank {
            l = m + 1;
        } else {
            r = m;
        }
    }

    // Scan forward within the epsilon window.
    // The compress invariant guarantees that if a valid element exists,
    // it is at or after the binary search landing point.
    while l < s_m.len() && s_m[l].rmin <= rank + epsilon_n {
        if s_m[l].rmax <= rank + epsilon_n {
            return Some(s_m[l].val.clone());
        }
        l += 1;
    }

    None
}
```

### 3d. Run all tests to verify no regression
```bash
cargo test
```

### 3e. Commit
```bash
git commit -m "Fix find_idx: guard empty slice underflow, fix linear scan termination"
```

---

## Step 4: Fix `PartialOrd`/`Ord` inconsistency on `RankInfo`

**File**: `src/lib.rs`

### 4a. No new test needed — existing tests cover sorting behavior

### 4b. Write fix
Replace the `PartialOrd` impl:
```rust
impl<T: Clone + Ord> PartialOrd for RankInfo<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
```

### 4c. Run all tests
```bash
cargo test
```

### 4d. Commit
```bash
git commit -m "Fix PartialOrd to delegate to Ord for RankInfo"
```

---

## Step 5: Remove unnecessary `T: Debug` bound from `UnboundEpsilonSummary`

**File**: `src/lib.rs`

### 5a. Write test to verify non-Debug types work
```rust
#[test]
fn test_unbound_summary_with_non_debug_type() {
    // A type that is Clone + Ord but NOT Debug
    #[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
    struct Opaque(i32);

    let mut s = UnboundEpsilonSummary::new(0.1).unwrap(); // will use Result after Step 6
    s.update(Opaque(1));
    s.update(Opaque(2));
    // Just verify it compiles and runs
}
```
Note: this test will need adjustment after Step 6 changes the constructor to return Result. For now, use the current panicking constructor.

### 5b. Write fix
Change line 398-399 from:
```rust
impl<T> UnboundEpsilonSummary<T>
where
    T: Clone + Ord + std::fmt::Debug,
```
to:
```rust
impl<T> UnboundEpsilonSummary<T>
where
    T: Clone + Ord,
```

### 5c. Run all tests
```bash
cargo test
```

### 5d. Commit
```bash
git commit -m "Remove unnecessary T: Debug bound from UnboundEpsilonSummary"
```

---

## Step 6: Constructors return `Result<Self, QuantileError>`

**File**: `src/lib.rs`

### 6a. Write failing tests
```rust
#[test]
fn test_fixedsize_constructor_returns_error_on_zero_n() {
    assert_eq!(
        FixedSizeEpsilonSummary::<usize>::new(0, 0.1),
        Err(QuantileError::InvalidN)
    );
}

#[test]
fn test_fixedsize_constructor_returns_error_on_negative_epsilon() {
    assert_eq!(
        FixedSizeEpsilonSummary::<usize>::new(100, -0.1),
        Err(QuantileError::InvalidEpsilon)
    );
}

#[test]
fn test_fixedsize_constructor_returns_error_on_nan_epsilon() {
    assert_eq!(
        FixedSizeEpsilonSummary::<usize>::new(100, f64::NAN),
        Err(QuantileError::InvalidEpsilon)
    );
}

#[test]
fn test_fixedsize_constructor_returns_error_on_inf_epsilon() {
    assert_eq!(
        FixedSizeEpsilonSummary::<usize>::new(100, f64::INFINITY),
        Err(QuantileError::InvalidEpsilon)
    );
}

#[test]
fn test_unbound_constructor_returns_error_on_negative_epsilon() {
    assert_eq!(
        UnboundEpsilonSummary::<usize>::new(-0.1),
        Err(QuantileError::InvalidEpsilon)
    );
}

#[test]
fn test_unbound_constructor_returns_error_on_nan_epsilon() {
    assert_eq!(
        UnboundEpsilonSummary::<usize>::new(f64::NAN),
        Err(QuantileError::InvalidEpsilon)
    );
}
```

### 6b. Run tests to verify they fail
```bash
cargo test test_fixedsize_constructor_returns_error
cargo test test_unbound_constructor_returns_error
```

### 6c. Write implementation
Change `FixedSizeEpsilonSummary::new`:
```rust
#[must_use]
pub fn new(n: usize, epsilon: f64) -> Result<Self, QuantileError> {
    if n == 0 {
        return Err(QuantileError::InvalidN);
    }
    if !(epsilon > 0.0 && epsilon.is_finite()) {
        return Err(QuantileError::InvalidEpsilon);
    }
    // ... rest unchanged, but wrap return in Ok(...)
}
```

Change `UnboundEpsilonSummary::new`:
```rust
#[must_use]
pub fn new(epsilon: f64) -> Result<Self, QuantileError> {
    if !(epsilon > 0.0 && epsilon.is_finite()) {
        return Err(QuantileError::InvalidEpsilon);
    }
    // ... rest unchanged, but wrap return in Ok(...)
    // Internal call: FixedSizeEpsilonSummary::new(n, epsilon / 2.0)?
}
```

Update ALL existing tests and the `UnboundEpsilonSummary::update` method to use `.unwrap()` on constructor calls. Remove the old `#[should_panic]` tests.

### 6d. Run all tests
```bash
cargo test
```

### 6e. Commit
```bash
git commit -m "Change constructors to return Result, reject NaN/Inf epsilon"
```

---

## Step 7: `query(&self)` via RefCell + return `Result<T, QuantileError>`

**File**: `src/lib.rs`

### 7a. Write failing tests
```rust
#[test]
fn test_query_returns_error_on_empty_summary() {
    let s = FixedSizeEpsilonSummary::<usize>::new(10, 0.1).unwrap();
    assert_eq!(s.query(0.5), Err(QuantileError::EmptySummary));
}

#[test]
fn test_query_returns_error_on_invalid_rank() {
    let mut s = FixedSizeEpsilonSummary::new(10, 0.1).unwrap();
    s.update(1);
    assert_eq!(s.query(-0.1), Err(QuantileError::InvalidRank));
    assert_eq!(s.query(1.1), Err(QuantileError::InvalidRank));
    assert_eq!(s.query(f64::NAN), Err(QuantileError::InvalidRank));
}

#[test]
fn test_query_is_immutable() {
    let mut s = FixedSizeEpsilonSummary::new(10, 0.1).unwrap();
    for i in 1..=10 {
        s.update(i);
    }
    // query takes &self, not &mut self
    let s_ref = &s;
    let _ = s_ref.query(0.5);
    let _ = s_ref.query(0.9);
}
```

### 7b. Run tests to verify they fail
```bash
cargo test test_query_returns_error
cargo test test_query_is_immutable
```

### 7c. Write implementation

Add `use std::cell::RefCell;` to imports.

Change `cached_s_m` field in `FixedSizeEpsilonSummary`:
```rust
cached_s_m: RefCell<Option<Vec<RankInfo<T>>>>,
```

Change `FixedSizeEpsilonSummary::new` to initialize:
```rust
cached_s_m: RefCell::new(None),
```

Change `update` to invalidate cache:
```rust
pub fn update(&mut self, e: T) {
    self.cached_s_m.get_mut().take(); // get_mut avoids runtime borrow check
    // ... rest unchanged
}
```

Change `query` signature and implementation:
```rust
pub fn query(&self, r: f64) -> Result<T, QuantileError> {
    if !(r >= 0.0 && r <= 1.0) {
        return Err(QuantileError::InvalidRank);
    }
    if self.cnt == 0 {
        return Err(QuantileError::EmptySummary);
    }

    if self.cached_s_m.borrow().is_none() {
        // Clone s[0] and sort — do NOT mutate self.s
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

    let rank: i64 = ((self.cnt as f64) * r).floor() as i64;
    let epsilon_n: i64 = ((self.cnt as f64) * self.epsilon).floor() as i64;
    let cache = self.cached_s_m.borrow();
    let s_m = cache.as_ref().unwrap();
    if let Some(e) = find_idx(s_m, rank, epsilon_n) {
        Ok(e)
    } else {
        Ok(s_m.last().unwrap().val.clone())
    }
}
```

Apply same changes to `UnboundEpsilonSummary`:
- `cached_s_m: RefCell<Option<Vec<RankInfo<T>>>>`
- `update`: `self.cached_s_m.get_mut().take();`
- `query(&self, r: f64) -> Result<T, QuantileError>` with clone-and-sort pattern
- Refactor `calc_s_m` to take `&self` (clone and sort `s[0]` instead of in-place mutation)

Update ALL existing tests to call `.unwrap()` on query results.

### 7d. Run all tests
```bash
cargo test
```

### 7e. Commit
```bash
git commit -m "Make query(&self) via RefCell cache, return Result, validate rank"
```

---

## Step 8: Fix boundary array index out-of-bounds

**File**: `src/lib.rs`

### 8a. Write failing test
```rust
#[test]
fn test_unbound_summary_large_n_no_panic() {
    // With epsilon=0.1, boundaries[16] = 655350
    // Inserting 655351 elements should not panic
    let epsilon = 0.1;
    let mut s = UnboundEpsilonSummary::new(epsilon).unwrap();
    for i in 0..=655_351 {
        s.update(i);
    }
    let ans = s.query(0.5).unwrap();
    // Just verify it doesn't panic; exact value depends on algorithm
    assert!(ans > 0);
}
```

### 8b. Run test to verify it fails (should panic with index out of bounds)
```bash
cargo test test_unbound_summary_large_n_no_panic -- --nocapture
```

### 8c. Write fix
In `UnboundEpsilonSummary::update`, guard the boundary access:
```rust
pub fn update(&mut self, e: T) {
    self.cached_s_m.get_mut().take();
    self.s_c.update(e);

    if let Some(x) = is_boundary(self.cnt + 1, &self.boundaries) {
        if x + x >= self.boundaries.len() {
            // Beyond the pre-computed boundary range; skip restructuring
            self.cnt += 1;
            return;
        }
        self.s_c.finalize();
        let upper_bound = self.boundaries[x + x];
        let n = upper_bound - self.cnt - 1;
        if n > 0 {
            let mut summary = FixedSizeEpsilonSummary::new(n, self.epsilon / 2.0).unwrap();
            std::mem::swap(&mut self.s_c, &mut summary);
            self.s.push(summary);
        }
    }
    self.cnt += 1;
}
```

### 8d. Run all tests
```bash
cargo test
```
Note: this test may be slow (~5-10 seconds) due to inserting 655K elements. That's acceptable.

### 8e. Commit
```bash
git commit -m "Fix index out-of-bounds in UnboundEpsilonSummary for large inputs"
```

---

## Step 9: Clean up internal methods (`calc_s_m`, `finalize`)

**File**: `src/lib.rs`

### 9a. No new tests — existing tests cover behavior

### 9b. Write changes
Remove the `epsilon` parameter from `calc_s_m` and `finalize` on `FixedSizeEpsilonSummary`:
```rust
fn calc_s_m(&self) -> Vec<RankInfo<T>> {
    let mut s0_clone = self.s[0].clone();
    s0_clone.sort();
    for (i, r) in s0_clone.iter_mut().enumerate() {
        r.rmin = i as i64;
        r.rmax = i as i64;
    }
    let mut s_m = s0_clone;
    for i in 1..self.level {
        s_m = merge(s_m, &self.s[i])
    }
    compress(s_m, self.b, self.epsilon)
}

fn finalize(&mut self) {
    let s_m = self.calc_s_m();
    self.s.clear();
    self.s.push(s_m);
}
```

Update call sites in `UnboundEpsilonSummary`:
- `self.s_c.finalize()` (no arg)
- `self.s_c.calc_s_m()` (no arg)

### 9c. Run all tests
```bash
cargo test
```

### 9d. Commit
```bash
git commit -m "Remove redundant epsilon parameter from internal methods"
```

---

## Step 10: Add `#[must_use]`, `Debug` derives

**File**: `src/lib.rs`

### 10a. Write test
```rust
#[test]
fn test_debug_impls() {
    let s = FixedSizeEpsilonSummary::new(10, 0.1).unwrap();
    let debug_str = format!("{:?}", s);
    assert!(debug_str.contains("FixedSizeEpsilonSummary"));

    let s = UnboundEpsilonSummary::<usize>::new(0.1).unwrap();
    let debug_str = format!("{:?}", s);
    assert!(debug_str.contains("UnboundEpsilonSummary"));
}
```

### 10b. Write implementation
Add `Debug` derive to `FixedSizeEpsilonSummary`. Since `RefCell` field makes auto-derive tricky, implement manually:
```rust
impl<T: Clone + Ord + std::fmt::Debug> std::fmt::Debug for FixedSizeEpsilonSummary<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FixedSizeEpsilonSummary")
            .field("epsilon", &self.epsilon)
            .field("b", &self.b)
            .field("level", &self.level)
            .field("cnt", &self.cnt)
            .finish()
    }
}
```

Same for `UnboundEpsilonSummary`.

Add `#[must_use]` to:
- `FixedSizeEpsilonSummary::new`
- `FixedSizeEpsilonSummary::query`
- `FixedSizeEpsilonSummary::size`
- `UnboundEpsilonSummary::new`
- `UnboundEpsilonSummary::query`
- `UnboundEpsilonSummary::size`

### 10c. Run all tests
```bash
cargo test
```

### 10d. Commit
```bash
git commit -m "Add Debug impls and #[must_use] to public API"
```

---

## Step 11: Add serde support (feature-gated)

**File**: `Cargo.toml`, `src/lib.rs`

### 11a. Update Cargo.toml
Add to `[dependencies]`:
```toml
serde = { version = "1", features = ["derive"], optional = true }
```

Add features section:
```toml
[features]
serde = ["dep:serde"]
```

### 11b. Write test
```rust
#[cfg(feature = "serde")]
#[test]
fn test_serde_roundtrip_fixedsize() {
    let mut s = FixedSizeEpsilonSummary::new(100, 0.1).unwrap();
    for i in 1..=100 {
        s.update(i);
    }
    let serialized = serde_json::to_string(&s).unwrap();
    let deserialized: FixedSizeEpsilonSummary<usize> = serde_json::from_str(&serialized).unwrap();
    assert_eq!(s.query(0.5).unwrap(), deserialized.query(0.5).unwrap());
}

#[cfg(feature = "serde")]
#[test]
fn test_serde_roundtrip_unbound() {
    let mut s = UnboundEpsilonSummary::new(0.1).unwrap();
    for i in 1..=100usize {
        s.update(i);
    }
    let serialized = serde_json::to_string(&s).unwrap();
    let deserialized: UnboundEpsilonSummary<usize> = serde_json::from_str(&serialized).unwrap();
    assert_eq!(s.query(0.5).unwrap(), deserialized.query(0.5).unwrap());
}
```

Add `serde_json` to dev-dependencies:
```toml
[dev-dependencies]
serde_json = "1"
```

### 11c. Write implementation
Add conditional derives to `RankInfo`:
```rust
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
struct RankInfo<T>
where
    T: Clone,
{
    val: T,
    rmin: i64,
    rmax: i64,
}
```

Add to `FixedSizeEpsilonSummary`:
```rust
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FixedSizeEpsilonSummary<T>
where
    T: Clone + Ord,
{
    epsilon: f64,
    b: usize,
    level: usize,
    cnt: usize,
    s: Vec<Vec<RankInfo<T>>>,
    #[cfg_attr(feature = "serde", serde(skip))]
    cached_s_m: RefCell<Option<Vec<RankInfo<T>>>>,
}
```

Same pattern for `UnboundEpsilonSummary` (with `#[serde(skip)]` on `cached_s_m`).

Note: serde bounds on `T` are needed. When `serde` feature is active, the derive macros will auto-add `T: Serialize`/`T: Deserialize` bounds.

### 11d. Run tests
```bash
cargo test --all-features
```

### 11e. Commit
```bash
git commit -m "Add feature-gated serde support with skip on cache fields"
```

---

## Step 12: Update dev-dependencies

**File**: `Cargo.toml`, `benches/zw_benchmark.rs`, `src/lib.rs` (test module)

### 12a. Update Cargo.toml dev-dependencies
```toml
[dev-dependencies]
criterion = "0.5"
rand = "0.9"
rand_distr = "0.5"
ordered-float = "4"
stats_alloc = "0.1.8"
quantiles = "0.7.1"
serde_json = "1"
```

Remove `rand_pcg` if it's merged into `rand 0.9`, or update to compatible version.

### 12b. Fix compilation errors from API changes

**`benches/zw_benchmark.rs`**: Remove `#[macro_use] extern crate criterion;` and `extern crate` lines (not needed in edition 2021). Update benchmark function calls if API changed (query now returns Result).

**`src/lib.rs` test module**: Update `rand_pcg` usage for `rand 0.9` API (may need `rand::SeedableRng` instead). Update `ordered-float` imports for v4. Check `stats_alloc` compatibility.

### 12c. Run all tests and benchmarks compile
```bash
cargo test --all-features
cargo bench --no-run
```

### 12d. Commit
```bash
git commit -m "Update dev-dependencies to latest versions"
```

---

## Step 13: Modernize CI

**File**: `.github/workflows/CI.yml`

### 13a. Rewrite CI workflow
```yaml
on: [push, pull_request]

name: CI

jobs:
  check:
    name: Check
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - run: cargo check --all-features

  test:
    name: Test Suite
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest, windows-latest]
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - run: cargo build
      - run: cargo test --all-features

  msrv:
    name: MSRV (1.65)
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@1.65
      - run: cargo check
      - run: cargo check --all-features

  clippy:
    name: Clippy
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
        with:
          components: clippy
      - run: cargo clippy --all-features -- -D warnings

  fmt:
    name: Rustfmt
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
        with:
          components: rustfmt
      - run: cargo fmt --all -- --check

  codecov:
    name: Code Coverage
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - run: cargo install cargo-tarpaulin
      - run: cargo tarpaulin --all-features --out xml
      - uses: codecov/codecov-action@v4
        with:
          token: ${{ secrets.CODECOV_TOKEN }}
          file: cobertura.xml
```

### 13b. Verify locally (what we can)
```bash
cargo clippy --all-features -- -D warnings
cargo fmt --all -- --check
```

### 13c. Commit
```bash
git commit -m "Modernize CI: update actions, add clippy + MSRV jobs"
```

---

## Step 14: Update README and doc comments

**Files**: `README.md`, `src/lib.rs`

### 14a. Update `src/lib.rs` doc comment
Fix the crate-level doc comment to show the new API:
```rust
//! ## Example
//!
//! ```rust
//! use zw_fast_quantile::FixedSizeEpsilonSummary;
//!
//! let epsilon = 0.1;
//! let n = 10;
//! let mut s = FixedSizeEpsilonSummary::new(n, epsilon).unwrap();
//! for i in 1..=n {
//!     s.update(i);
//! }
//!
//! let ans = s.query(0.0).unwrap();
//! let expected = 1;
//! assert!(expected == ans);
//! ```
//!
//! ```rust
//! use zw_fast_quantile::UnboundEpsilonSummary;
//!
//! let epsilon = 0.1;
//! let n = 10;
//! let mut s = UnboundEpsilonSummary::new(epsilon).unwrap();
//! for i in 1..=n {
//!     s.update(i);
//! }
//!
//! let ans = s.query(0.0).unwrap();
//! let expected = 1;
//! assert!(expected == ans);
//! ```
```

Also update version in doc comment from `"0.1"` to `"1.0"`.

### 14b. Update `README.md`
- Change installation to `zw-fast-quantile = "1.0"`
- Update examples to use `.unwrap()` on `new()` and `query()`
- Remove the Rust 2015 `extern crate` note
- Add `serde` feature mention
- Add MSRV note
- Remove or update stale benchmark numbers

### 14c. Run doc-tests
```bash
cargo test --doc
```

### 14d. Commit
```bash
git commit -m "Update README and doc comments for 1.0 API"
```

---

## Step 15: Version bump + final end-to-end verification

**File**: `Cargo.toml`

### 15a. Bump version
Change `version = "0.2.3"` to `version = "1.0.0"`.

### 15b. Run full verification
```bash
# All tests pass
cargo test --all-features

# Doc-tests pass
cargo test --doc

# Clippy clean
cargo clippy --all-features -- -D warnings

# Formatting clean
cargo fmt --all -- --check

# Benchmarks compile
cargo bench --no-run

# No serde feature still compiles
cargo check

# serde feature compiles
cargo check --all-features

# Docs build without warnings
cargo doc --all-features --no-deps 2>&1 | grep -i warning || echo "No doc warnings"
```

### 15c. Commit
```bash
git commit -m "Bump version to 1.0.0"
```

---

## CHANGELOG (to be maintained during implementation)

Track progress, completed tasks, failed approaches, and known limitations in this section as work proceeds.

### Status: Not started

### Completed:
(none yet)

### Failed approaches:
(none yet)

### Known limitations:
(none yet)
