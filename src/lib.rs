//! Zhang Wang Fast Approximate Quantiles Algorithm in Rust
//!
//! ## Installation
//!
//! Add this to your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! zw-fast-quantile = "1.0"
//! ```
//!
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
//!
use std::cell::RefCell;
use std::cmp::Ordering;
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

impl<T> RankInfo<T>
where
    T: Clone,
{
    fn new(val: T, rmin: i64, rmax: i64) -> Self {
        RankInfo { val, rmin, rmax }
    }
}

impl<T: Clone + Ord> Ord for RankInfo<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.val.cmp(&other.val)
    }
}

impl<T: Clone + Ord> PartialOrd for RankInfo<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T: Clone + PartialEq> PartialEq for RankInfo<T> {
    fn eq(&self, other: &Self) -> bool {
        self.val == other.val
    }
}

impl<T: Clone + PartialEq> Eq for RankInfo<T> {}

#[derive(Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
/// An epsilon-approximate quantile summary for a stream with a known size.
///
/// The summary accepts at most the `n` elements declared in [`Self::new`].
pub struct FixedSizeEpsilonSummary<T>
where
    T: Clone + Ord,
{
    epsilon: f64,
    capacity: usize,
    b: usize,
    level: usize,
    cnt: usize,
    s: Vec<Vec<RankInfo<T>>>,
    #[cfg_attr(feature = "serde", serde(skip))]
    cached_s_m: RefCell<Option<Vec<RankInfo<T>>>>,
}

impl<T> FixedSizeEpsilonSummary<T>
where
    T: Clone + Ord,
{
    pub fn new(n: usize, epsilon: f64) -> Result<Self, QuantileError> {
        if n == 0 {
            return Err(QuantileError::InvalidN);
        }

        if !(epsilon > 0.0 && epsilon.is_finite()) {
            return Err(QuantileError::InvalidEpsilon);
        }

        // number_of_leves = ceil(log2(n))
        // block_size = floor(log2(epsilon * N) / epsilon)
        let epsilon_n: f64 = (n as f64) * epsilon;
        let number_of_levels = if epsilon_n > 1.0 {
            epsilon_n.log2().ceil() as usize
        } else {
            1
        };
        let block_size = if number_of_levels > 1 {
            (epsilon_n.log2() / epsilon).floor() as usize
        } else {
            n + 1
        };

        let mut s = vec![vec![]; number_of_levels];
        s[0].reserve_exact(block_size);

        Ok(FixedSizeEpsilonSummary {
            epsilon,
            capacity: n,
            b: block_size,
            level: number_of_levels,
            cnt: 0,
            s,
            cached_s_m: RefCell::new(None),
        })
    }

    /// Adds an element to the summary.
    ///
    /// # Panics
    ///
    /// Panics if more than the `n` elements declared in [`Self::new`] are added.
    pub fn update(&mut self, e: T) {
        assert!(
            self.cnt < self.capacity,
            "FixedSizeEpsilonSummary capacity exceeded: constructed for n={} elements",
            self.capacity
        );
        self.cached_s_m.get_mut().take();
        let rank_info = RankInfo::new(e, 0, 0);
        self.s[0].push(rank_info);

        self.cnt += 1;
        if self.s[0].len() < self.b {
            return;
        }

        self.s[0].sort();
        for (i, r) in self.s[0].iter_mut().enumerate() {
            r.rmin = i as i64;
            r.rmax = i as i64;
        }

        let compressed_size = self.b / 2;
        let mut s_c = compress(self.s[0].clone(), compressed_size, self.epsilon);

        self.s[0].clear();
        let mut stored = false;
        for k in 1..self.level {
            if self.s[k].is_empty() {
                self.s[k] = s_c;
                stored = true;
                break;
            } else {
                let t = merge(s_c, &self.s[k]);
                s_c = compress(t, compressed_size, self.epsilon);
                self.s[k].clear();
            }
        }
        debug_assert!(stored, "capacity invariant guarantees an empty summary level");
    }

    pub fn query(&self, r: f64) -> Result<T, QuantileError> {
        if !(0.0..=1.0).contains(&r) {
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

    #[inline]
    #[must_use]
    pub fn size(&self) -> usize {
        self.cnt
    }
}

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

fn merge<T: Clone + Ord>(s_a: Vec<RankInfo<T>>, s_b: &[RankInfo<T>]) -> Vec<RankInfo<T>> {
    if s_a.is_empty() {
        return s_b.to_vec();
    }

    if s_b.is_empty() {
        return s_a;
    }

    let mut s_m = Vec::with_capacity(s_a.len() + s_b.len());

    let mut i1 = 0;
    let mut i2 = 0;
    let mut from;

    while i1 < s_a.len() || i2 < s_b.len() {
        let val;
        let rmin;
        let rmax;

        if i1 < s_a.len() && i2 < s_b.len() {
            if s_a[i1].val < s_b[i2].val {
                val = s_a[i1].val.clone();
                from = 1;
            } else {
                val = s_b[i2].val.clone();
                from = 2;
            }
        } else if i1 < s_a.len() && i2 >= s_b.len() {
            val = s_a[i1].val.clone();
            from = 1;
        } else {
            val = s_b[i2].val.clone();
            from = 2;
        }

        if from == 1 {
            if 0 < i2 && i2 < s_b.len() {
                rmin = s_a[i1].rmin + s_b[i2 - 1].rmin;
                rmax = s_a[i1].rmax + s_b[i2].rmax - 1;
            } else if i2 == 0 {
                rmin = s_a[i1].rmin;
                rmax = s_a[i1].rmax + s_b[i2].rmax - 1;
            } else {
                rmin = s_a[i1].rmin + s_b[i2 - 1].rmin;
                rmax = s_a[i1].rmax + s_b[i2 - 1].rmax;
            }

            i1 += 1;
        } else {
            if 0 < i1 && i1 < s_a.len() {
                rmin = s_a[i1 - 1].rmin + s_b[i2].rmin;
                rmax = s_a[i1].rmax + s_b[i2].rmax - 1;
            } else if i1 == 0 {
                rmin = s_b[i2].rmin;
                rmax = s_a[i1].rmax + s_b[i2].rmax - 1;
            } else {
                rmin = s_a[i1 - 1].rmin + s_b[i2].rmin;
                rmax = s_a[i1 - 1].rmax + s_b[i2].rmax;
            }

            i2 += 1;
        }

        let rank_info = RankInfo::new(val, rmin, rmax);
        s_m.push(rank_info);
    }

    s_m
}

fn compress<T: Clone>(mut s0: Vec<RankInfo<T>>, block_size: usize, epsilon: f64) -> Vec<RankInfo<T>> {
    let mut s0_range = 0;
    let mut e: f64 = 0.0;

    for r in &s0 {
        if s0_range < r.rmax {
            s0_range = r.rmax;
        }

        if (r.rmax - r.rmin) as f64 > e {
            e = (r.rmax - r.rmin) as f64;
        }
    }

    let epsilon_n: f64 = epsilon * (s0_range as f64);
    assert!(2.0 * epsilon_n >= e, "precision condition violated.");

    let mut i = 0;
    let mut j = 0;
    let mut k = 0;
    let n = s0.len();
    while i <= block_size && j < n {
        let r = ((i as f64) * (s0_range as f64) / (block_size as f64)).floor() as i64;

        while j < n {
            if s0[j].rmax >= r {
                break;
            }

            j += 1;
        }

        assert!(j < n, "unable to find the summary with precision given.");
        s0[k] = s0[j].clone();
        k += 1;
        j += 1;
        i += 1;
    }

    s0.truncate(k);
    s0
}

#[inline]
#[allow(clippy::many_single_char_names)]
#[allow(clippy::comparison_chain)]
fn is_boundary(x: usize, boundaries: &[usize; 32]) -> Option<usize> {
    let mut l = 0;
    let mut r = 31;

    while l < r {
        let m = l + (r - l) / 2;
        if boundaries[m] < x {
            l = m + 1;
        } else {
            r = m;
        }
    }

    if l < 31 && boundaries[l] == x {
        return Some(l);
    }

    None
}

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

#[derive(Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct UnboundEpsilonSummary<T>
where
    T: Clone + Ord,
{
    epsilon: f64,
    cnt: usize,
    s: Vec<FixedSizeEpsilonSummary<T>>,
    s_c: FixedSizeEpsilonSummary<T>,
    boundaries: [usize; 32],
    #[cfg_attr(feature = "serde", serde(skip))]
    cached_s_m: RefCell<Option<Vec<RankInfo<T>>>>,
}

impl<T> UnboundEpsilonSummary<T>
where
    T: Clone + Ord,
{
    pub fn new(epsilon: f64) -> Result<Self, QuantileError> {
        if !(epsilon > 0.0 && epsilon.is_finite()) {
            return Err(QuantileError::InvalidEpsilon);
        }

        let s = vec![];

        let n = (1.0_f64 / epsilon).floor() as usize;
        let s_c = FixedSizeEpsilonSummary::new(n, epsilon / 2.0)?;
        let mut boundaries: [usize; 32] = [0; 32];
        for (i, b) in boundaries.iter_mut().enumerate() {
            *b = (((usize::pow(2, i as u32) - 1) as f64) / epsilon).floor() as usize;
        }

        Ok(UnboundEpsilonSummary {
            epsilon,
            cnt: 0,
            s,
            s_c,
            boundaries,
            cached_s_m: RefCell::new(None),
        })
    }

    pub fn update(&mut self, e: T) {
        self.cached_s_m.get_mut().take();
        self.s_c.update(e);

        if let Some(x) = is_boundary(self.cnt + 1, &self.boundaries) {
            if x + x >= self.boundaries.len() {
                // The summary created at boundary 15 is sized through boundary 30,
                // so it has headroom after boundary rotation stops at x = 16.
                self.cnt += 1;
                return;
            }
            self.s_c.finalize();

            let upper_bound = self.boundaries[x + x];
            let n = upper_bound - self.cnt - 1;
            let mut summary = FixedSizeEpsilonSummary::new(n, self.epsilon / 2.0).unwrap();
            std::mem::swap(&mut self.s_c, &mut summary);

            self.s.push(summary);
        }
        self.cnt += 1;
    }

    pub fn query(&self, r: f64) -> Result<T, QuantileError> {
        if !(0.0..=1.0).contains(&r) {
            return Err(QuantileError::InvalidRank);
        }
        if self.cnt == 0 {
            return Err(QuantileError::EmptySummary);
        }

        if self.cached_s_m.borrow().is_none() {
            let mut s_m = self.s_c.calc_s_m();
            for i in 0..self.s.len() {
                for j in 0..self.s[i].s.len() {
                    s_m = merge(s_m, &self.s[i].s[j])
                }
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

    #[inline]
    #[must_use]
    pub fn size(&self) -> usize {
        self.cnt
    }
}

impl<T: Clone + Ord + std::fmt::Debug> std::fmt::Debug for UnboundEpsilonSummary<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("UnboundEpsilonSummary")
            .field("epsilon", &self.epsilon)
            .field("cnt", &self.cnt)
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::Rng;
    use rand::SeedableRng;
    use rand_distr::Distribution;

    #[test]
    fn test_merge_and_compress() {
        let mut s0 = Vec::with_capacity(4);
        let mut s1 = Vec::with_capacity(4);

        s0.push(RankInfo::new(2, 1, 1));
        s0.push(RankInfo::new(4, 3, 4));
        s0.push(RankInfo::new(8, 5, 6));
        s0.push(RankInfo::new(17, 8, 8));

        s1.push(RankInfo::new(1, 1, 1));
        s1.push(RankInfo::new(7, 3, 3));
        s1.push(RankInfo::new(12, 5, 6));
        s1.push(RankInfo::new(15, 8, 8));

        let merged = merge(s0, &s1);

        assert_eq!(merged.len(), 8);
        let merged_vals: Vec<i32> = merged.iter().map(|x| x.val).collect();
        let merged_rmins: Vec<i64> = merged.iter().map(|x| x.rmin).collect();
        let merged_rmaxs: Vec<i64> = merged.iter().map(|x| x.rmax).collect();
        assert_eq!(merged_vals, vec![1, 2, 4, 7, 8, 12, 15, 17]);
        assert_eq!(merged_rmins, vec![1, 2, 4, 6, 8, 10, 13, 16]);
        assert_eq!(merged_rmaxs, vec![1, 3, 6, 8, 11, 13, 15, 16]);

        let epsilon: f64 = 0.2;
        let compressed = compress(merged, 4, epsilon);
        let compressed_vals: Vec<i32> = compressed.iter().map(|x| x.val).collect();
        assert_eq!(compressed_vals, vec![1, 4, 7, 12, 17]);
    }

    #[test]
    fn test_fixedsize_constructor_returns_error_on_zero_n() {
        assert!(matches!(
            FixedSizeEpsilonSummary::<usize>::new(0, 0.1),
            Err(QuantileError::InvalidN)
        ));
    }

    #[test]
    fn test_fixedsize_constructor_returns_error_on_negative_epsilon() {
        assert!(matches!(
            FixedSizeEpsilonSummary::<usize>::new(100, -0.1),
            Err(QuantileError::InvalidEpsilon)
        ));
    }

    #[test]
    fn test_fixedsize_constructor_returns_error_on_nan_epsilon() {
        assert!(matches!(
            FixedSizeEpsilonSummary::<usize>::new(100, f64::NAN),
            Err(QuantileError::InvalidEpsilon)
        ));
    }

    #[test]
    fn test_fixedsize_constructor_returns_error_on_inf_epsilon() {
        assert!(matches!(
            FixedSizeEpsilonSummary::<usize>::new(100, f64::INFINITY),
            Err(QuantileError::InvalidEpsilon)
        ));
    }

    #[test]
    fn test_unbound_constructor_returns_error_on_negative_epsilon() {
        assert!(matches!(
            UnboundEpsilonSummary::<usize>::new(-0.1),
            Err(QuantileError::InvalidEpsilon)
        ));
    }

    #[test]
    fn test_unbound_constructor_returns_error_on_nan_epsilon() {
        assert!(matches!(
            UnboundEpsilonSummary::<usize>::new(f64::NAN),
            Err(QuantileError::InvalidEpsilon)
        ));
    }

    #[test]
    fn test_query_fixed_summary_with_insufficient_values() {
        let epsilon = 0.1;
        let n = 3;
        let mut s = FixedSizeEpsilonSummary::new(n, epsilon).unwrap();
        for i in 1..=n {
            s.update(i);
        }

        let rank: f64 = 1.0;
        let ans = s.query(rank).unwrap();
        assert!(n == ans);
    }

    #[test]
    fn test_query_with_small_n_on_fixedsize_summary() {
        let epsilon = 0.1;
        let n = 10;
        let mut s = FixedSizeEpsilonSummary::new(n, epsilon).unwrap();
        for i in 1..=n {
            s.update(i);
        }

        for i in 1..=n {
            let rank: f64 = ((i - 1) as f64) / (n as f64);
            let ans = s.query(rank).unwrap();
            assert!(i == ans);
        }
    }

    #[test]
    #[should_panic(expected = "FixedSizeEpsilonSummary capacity exceeded")]
    fn test_fixedsize_update_panics_when_capacity_is_exceeded() {
        let mut s = FixedSizeEpsilonSummary::new(10, 0.1).unwrap();
        for i in 0..=10 {
            s.update(i);
        }
    }

    fn assert_rank_error<F>(n: usize, epsilon: f64, mut update: F, query: impl Fn(f64) -> u64)
    where
        F: FnMut(u64),
    {
        let mut rng = StdRng::seed_from_u64(42);
        let mut records = Vec::with_capacity(n);
        for _ in 0..n {
            let value = rng.random::<u64>();
            records.push(value);
            update(value);
        }
        records.sort_unstable();

        let allowed_error = (epsilon * n as f64).ceil() as usize;
        for step in 0..=100 {
            let rank = step as f64 / 100.0;
            let target = (rank * n as f64).floor() as usize;
            let value = query(rank);
            let low = records.partition_point(|x| *x < value);
            let high = records.partition_point(|x| *x <= value);

            assert!(
                low <= target.saturating_add(allowed_error) && target <= high.saturating_add(allowed_error),
                "rank {rank}: value rank interval [{low}, {high}] misses target {target} by more than {allowed_error}"
            );
        }
    }

    #[test]
    fn test_fixedsize_summary_respects_rank_error_contract() {
        let n = 10_000;
        let epsilon = 0.01;
        let summary = RefCell::new(FixedSizeEpsilonSummary::new(n, epsilon).unwrap());
        assert_rank_error(
            n,
            epsilon,
            |value| summary.borrow_mut().update(value),
            |rank| summary.borrow().query(rank).unwrap(),
        );
    }

    #[test]
    fn test_unbound_summary_respects_rank_error_contract() {
        let n = 100_000;
        let epsilon = 0.01;
        let summary = RefCell::new(UnboundEpsilonSummary::new(epsilon).unwrap());
        assert_rank_error(
            n,
            epsilon,
            |value| summary.borrow_mut().update(value),
            |rank| summary.borrow().query(rank).unwrap(),
        );
    }

    #[test]
    fn test_query_with_small_n_on_unbound_summary() {
        let epsilon = 0.1;
        let n = 10;
        let mut s = UnboundEpsilonSummary::new(epsilon).unwrap();
        for i in 1..=n {
            s.update(i);
        }

        for i in 1..=n {
            let rank: f64 = ((i - 1) as f64) / (n as f64);
            let ans = s.query(rank).unwrap();
            assert!(i == ans);
        }
    }

    trait TestSummary<T> {
        fn update(&mut self, value: T);
        fn query(&self, rank: f64) -> T;
    }

    impl<T: Clone + Ord> TestSummary<T> for FixedSizeEpsilonSummary<T> {
        fn update(&mut self, value: T) {
            self.update(value);
        }

        fn query(&self, rank: f64) -> T {
            self.query(rank).unwrap()
        }
    }

    impl<T: Clone + Ord> TestSummary<T> for UnboundEpsilonSummary<T> {
        fn update(&mut self, value: T) {
            self.update(value);
        }

        fn query(&self, rank: f64) -> T {
            self.query(rank).unwrap()
        }
    }

    fn assert_distribution_queries<D, S>(mut summary: S, distribution: D, n: usize, tolerances: [f64; 4])
    where
        D: Distribution<f64>,
        S: TestSummary<ordered_float::NotNan<f64>>,
    {
        let mut rng = StdRng::seed_from_u64(42);
        let mut records = Vec::with_capacity(n);
        for _ in 0..n {
            records.push(ordered_float::NotNan::new(distribution.sample(&mut rng)).unwrap());
        }
        records.sort_unstable();
        for &value in &records {
            summary.update(value);
        }

        for (rank, index, tolerance) in [
            (0.5, n / 2, tolerances[0]),
            (0.0, 0, tolerances[1]),
            (0.99, n * 99 / 100, tolerances[2]),
            (1.0, n - 1, tolerances[3]),
        ] {
            assert!((summary.query(rank) - records[index]).abs() < tolerance);
        }
    }

    #[test]
    fn test_normal_distribution_generated_seq_on_fixed_summary() {
        let n = 1_000_000;
        assert_distribution_queries(
            FixedSizeEpsilonSummary::new(n, 0.01).unwrap(),
            rand_distr::Normal::new(0.5, 0.2).unwrap(),
            n,
            [0.01, 0.1, 0.01, 0.01],
        );
    }

    #[test]
    fn test_pareto_distribution_generated_seq_on_fixed_summary() {
        let n = 1_000_000;
        assert_distribution_queries(
            FixedSizeEpsilonSummary::new(n, 0.001).unwrap(),
            rand_distr::Pareto::new(5.0, 10.0).unwrap(),
            n,
            [0.01; 4],
        );
    }

    #[test]
    fn test_normal_distribution_generated_seq_on_unbound_summary() {
        assert_distribution_queries(
            UnboundEpsilonSummary::new(0.01).unwrap(),
            rand_distr::Normal::new(0.5, 0.2).unwrap(),
            1_000_000,
            [0.01; 4],
        );
    }

    #[test]
    fn test_pareto_distribution_generated_seq_on_unbound_summary() {
        assert_distribution_queries(
            UnboundEpsilonSummary::new(0.001).unwrap(),
            rand_distr::Pareto::new(5.0, 10.0).unwrap(),
            1_000_000,
            [0.01; 4],
        );
    }

    #[test]
    fn test_unbound_summary_clone_preserves_queries() {
        let mut summary = UnboundEpsilonSummary::new(0.1).unwrap();
        summary.update(1);
        assert_eq!(summary.query(0.5), summary.clone().query(0.5));
    }

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

    #[test]
    fn test_find_idx_empty_slice() {
        let empty: Vec<RankInfo<i32>> = vec![];
        assert!(find_idx(&empty, 0, 1).is_none());
    }

    #[test]
    fn test_query_returns_error_on_empty_summary() {
        let s = FixedSizeEpsilonSummary::<usize>::new(10, 0.1).unwrap();
        assert!(matches!(s.query(0.5), Err(QuantileError::EmptySummary)));
    }

    #[test]
    fn test_query_returns_error_on_invalid_rank() {
        let mut s = FixedSizeEpsilonSummary::new(10, 0.1).unwrap();
        s.update(1);
        assert!(matches!(s.query(-0.1), Err(QuantileError::InvalidRank)));
        assert!(matches!(s.query(1.1), Err(QuantileError::InvalidRank)));
        assert!(matches!(s.query(f64::NAN), Err(QuantileError::InvalidRank)));
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

    #[test]
    fn test_debug_impls() {
        let s = FixedSizeEpsilonSummary::<usize>::new(10, 0.1).unwrap();
        let debug_str = format!("{:?}", s);
        assert!(debug_str.contains("FixedSizeEpsilonSummary"));

        let s = UnboundEpsilonSummary::<usize>::new(0.1).unwrap();
        let debug_str = format!("{:?}", s);
        assert!(debug_str.contains("UnboundEpsilonSummary"));
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_serde_roundtrip_fixedsize() {
        let mut s = FixedSizeEpsilonSummary::new(100, 0.1).unwrap();
        for i in 1..=100usize {
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
}
