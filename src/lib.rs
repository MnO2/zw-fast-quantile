//! Zhang-Wang fast approximate quantiles algorithm in Rust.
//!
//! This crate implements the multi-level summary described by Zhang and Wang in
//! [*An Efficient Quantile Computation Technique for Approximate Query Processing*](http://web.cs.ucla.edu/~weiwang/paper/SSDBM07_2.pdf)
//! (SSDBM 2007). A summary built with error bound `epsilon` over `n` elements
//! answers a rank query `r` with an inserted element whose position in the sorted
//! stream is within `epsilon * n` of `floor(r * n)`, while storing far fewer than
//! `n` elements.
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
//! ## Choosing a summary
//!
//! - [`FixedSizeEpsilonSummary`] needs the stream length up front and uses the
//!   least memory. Inserting more than the declared number of elements panics.
//! - [`UnboundEpsilonSummary`] accepts any number of elements and grows its
//!   summaries as the stream continues.
//!
//! `epsilon` must be in `(0.0, 1.0]`. Very short streams, where
//! `epsilon * n` is at most `8.0`, are stored exactly rather than summarized.
//!
//! ## Serialization
//!
//! With the `serde` feature both summaries implement `Serialize` and
//! `Deserialize`. Query caches are skipped and rebuilt on the next query. The
//! serialized layout is an implementation detail: it can change between
//! releases, and a summary written by one crate version is not guaranteed to
//! load in another.
#![forbid(unsafe_code)]
#![warn(missing_docs)]

use std::cell::RefCell;
use std::cmp::Ordering;
use std::fmt;

/// Streams with `epsilon * n` at or below this value are stored exactly in a
/// single level instead of being summarized.
///
/// The Zhang-Wang block size `floor(log2(epsilon * n) / epsilon)` is only
/// meaningful once `epsilon * n` is well above 1. Below this threshold the
/// derived blocks are too small to keep the rank-error guarantee, so the
/// summary keeps every element instead. That costs at most
/// `EXACT_MODE_THRESHOLD / epsilon` stored elements, the same order as a
/// single summary block.
const EXACT_MODE_THRESHOLD: f64 = 8.0;

/// Errors that can occur when constructing or querying a quantile summary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QuantileError {
    /// `n` must be greater than 0.
    InvalidN,
    /// `epsilon` must be in `(0.0, 1.0]`.
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
            QuantileError::InvalidEpsilon => write!(f, "epsilon must be in (0.0, 1.0]"),
            QuantileError::InvalidRank => write!(f, "rank must be between 0.0 and 1.0"),
            QuantileError::EmptySummary => write!(f, "cannot query an empty summary"),
        }
    }
}

impl std::error::Error for QuantileError {}

/// An element together with the range of one-based ranks it may occupy in the
/// stream summarized so far.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
struct RankInfo<T> {
    val: T,
    rmin: i64,
    rmax: i64,
}

impl<T> RankInfo<T> {
    fn new(val: T, rmin: i64, rmax: i64) -> Self {
        RankInfo { val, rmin, rmax }
    }
}

impl<T: Ord> Ord for RankInfo<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.val.cmp(&other.val)
    }
}

impl<T: Ord> PartialOrd for RankInfo<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T: PartialEq> PartialEq for RankInfo<T> {
    fn eq(&self, other: &Self) -> bool {
        self.val == other.val
    }
}

impl<T: Eq> Eq for RankInfo<T> {}

/// Assigns exact one-based ranks to a sorted block.
fn assign_exact_ranks<T>(block: &mut [RankInfo<T>]) {
    for (i, r) in block.iter_mut().enumerate() {
        let rank = i as i64 + 1;
        r.rmin = rank;
        r.rmax = rank;
    }
}

/// An epsilon-approximate quantile summary for a stream with a known size.
///
/// The summary accepts at most the `n` elements declared in [`Self::new`].
/// Elements are buffered in blocks; each full block is sorted, compressed and
/// merged upwards through a small number of levels, so the memory used grows
/// with `log2(epsilon * n) / epsilon` rather than with `n`.
#[derive(Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FixedSizeEpsilonSummary<T>
where
    T: Clone + Ord,
{
    epsilon: f64,
    capacity: usize,
    b: usize,
    cnt: usize,
    s: Vec<Vec<RankInfo<T>>>,
    #[cfg_attr(feature = "serde", serde(skip))]
    cached_s_m: RefCell<Option<Vec<RankInfo<T>>>>,
}

impl<T> FixedSizeEpsilonSummary<T>
where
    T: Clone + Ord,
{
    /// Creates a summary for a stream of exactly `n` elements with error bound
    /// `epsilon`.
    ///
    /// A query for rank `r` returns an element whose position in the sorted
    /// stream is within `epsilon * n` of `floor(r * n)`. Smaller `epsilon`
    /// means more accuracy and more memory: the summary stores roughly
    /// `log2(epsilon * n) / epsilon` elements. When `epsilon * n` is at most
    /// `8.0` the summary keeps every element instead, because the Zhang-Wang
    /// block size is not meaningful for such short streams.
    ///
    /// # Errors
    ///
    /// Returns [`QuantileError::InvalidN`] if `n` is zero and
    /// [`QuantileError::InvalidEpsilon`] if `epsilon` is not in `(0.0, 1.0]`
    /// (NaN and infinite values are rejected as well).
    pub fn new(n: usize, epsilon: f64) -> Result<Self, QuantileError> {
        if n == 0 {
            return Err(QuantileError::InvalidN);
        }

        if !(epsilon > 0.0 && epsilon <= 1.0) {
            return Err(QuantileError::InvalidEpsilon);
        }

        let epsilon_n: f64 = (n as f64) * epsilon;
        let (block_size, number_of_levels) = if epsilon_n > EXACT_MODE_THRESHOLD {
            // block_size = floor(log2(epsilon * N) / epsilon); with
            // epsilon * N > 8 this is at least 3 / epsilon.
            let block_size = ((epsilon_n.log2() / epsilon).floor() as usize).max(2);
            // Full blocks propagate through the merge levels like binary carries.
            // The bit width of their maximum count gives the required merge
            // levels; add one for level 0.
            let full_blocks = n / block_size;
            let levels = (usize::BITS - full_blocks.leading_zeros()) as usize + 1;
            (block_size, levels)
        } else {
            // Exact mode: a single level that never fills up.
            (n + 1, 1)
        };

        let mut s = vec![Vec::new(); number_of_levels];
        s[0].reserve_exact(block_size.min(n));

        Ok(FixedSizeEpsilonSummary {
            epsilon,
            capacity: n,
            b: block_size,
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
        self.s[0].push(RankInfo::new(e, 0, 0));

        self.cnt += 1;
        if self.s[0].len() <= self.b {
            return;
        }

        let mut block = std::mem::replace(&mut self.s[0], Vec::with_capacity(self.b));
        block.sort_unstable();
        assign_exact_ranks(&mut block);

        let compressed_size = self.b / 2;
        let mut s_c = compress(block, compressed_size, self.epsilon);
        let mut stored = false;
        for level in self.s.iter_mut().skip(1) {
            if level.is_empty() {
                *level = s_c;
                stored = true;
                break;
            }
            let occupied = std::mem::take(level);
            s_c = compress(
                merge(s_c.into_iter(), occupied.into_iter()),
                compressed_size,
                self.epsilon,
            );
        }
        debug_assert!(
            stored,
            "capacity invariant failed: capacity={}, count={}, block_size={}, levels={}",
            self.capacity,
            self.cnt,
            self.b,
            self.s.len()
        );
    }

    /// Returns an approximate quantile.
    ///
    /// `r` is a rank in `0.0..=1.0`: `0.0` is the minimum, `0.5` the median
    /// and `1.0` the maximum. The result is an inserted element whose position
    /// in the sorted stream is within `epsilon * n` of `floor(r * n)`, where
    /// `n` is the number of elements inserted so far.
    ///
    /// The first query after an [`update`](Self::update) rebuilds the merged
    /// summary; later queries reuse it until the next update.
    ///
    /// # Errors
    ///
    /// Returns [`QuantileError::InvalidRank`] if `r` is outside `0.0..=1.0` or
    /// NaN, and [`QuantileError::EmptySummary`] if nothing has been inserted.
    pub fn query(&self, r: f64) -> Result<T, QuantileError> {
        if !(0.0..=1.0).contains(&r) {
            return Err(QuantileError::InvalidRank);
        }
        if self.cnt == 0 {
            return Err(QuantileError::EmptySummary);
        }

        let mut cache = self.cached_s_m.borrow_mut();
        let s_m = cache.get_or_insert_with(|| self.merged_levels());
        query_rank(s_m, r, self.cnt, self.epsilon).ok_or(QuantileError::EmptySummary)
    }

    /// Merges every level into one summary without compressing it. Level 0 is
    /// cloned and sorted so the summary itself is left untouched.
    fn merged_levels(&self) -> Vec<RankInfo<T>> {
        let mut s_m = self.s[0].clone();
        s_m.sort_unstable();
        assign_exact_ranks(&mut s_m);
        for level in self.s[1..].iter().filter(|level| !level.is_empty()) {
            s_m = merge(s_m.into_iter(), level.as_slice());
        }
        s_m
    }

    /// Merges every level and compresses the result back to one block.
    fn calc_s_m(&self) -> Vec<RankInfo<T>> {
        compress(self.merged_levels(), self.b, self.epsilon)
    }

    /// Collapses the summary into a single compressed level. Used by
    /// [`UnboundEpsilonSummary`] once a sub-stream is complete.
    fn finalize(&mut self) {
        let s_m = self.calc_s_m();
        self.s = vec![s_m];
    }

    /// Returns the number of elements inserted so far.
    #[inline]
    #[must_use]
    pub fn size(&self) -> usize {
        self.cnt
    }
}

impl<T: Clone + Ord + fmt::Debug> fmt::Debug for FixedSizeEpsilonSummary<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FixedSizeEpsilonSummary")
            .field("epsilon", &self.epsilon)
            .field("capacity", &self.capacity)
            .field("b", &self.b)
            .field("levels", &self.s.len())
            .field("cnt", &self.cnt)
            .finish()
    }
}

/// Resolves a rank query against a merged summary of `cnt` elements. Returns
/// `None` only when the summary is empty.
fn query_rank<T: Clone + Ord>(s_m: &[RankInfo<T>], r: f64, cnt: usize, epsilon: f64) -> Option<T> {
    // One-based target rank: the element at position floor(r * n) of the
    // sorted stream, clamped to the maximum.
    let rank: i64 = (((cnt as f64) * r).floor() as i64 + 1).min(cnt as i64);
    let epsilon_n: i64 = ((cnt as f64) * epsilon).floor() as i64;
    find_idx(s_m, rank, epsilon_n)
}

/// One side of a merge: a cursor that can look at the next element and hand
/// it over. Owned input moves its elements; borrowed input clones them.
trait MergeSource<T> {
    fn peek(&self) -> Option<&RankInfo<T>>;
    fn take(&mut self) -> Option<RankInfo<T>>;
    fn remaining(&self) -> usize;
}

impl<T> MergeSource<T> for std::vec::IntoIter<RankInfo<T>> {
    #[inline]
    fn peek(&self) -> Option<&RankInfo<T>> {
        self.as_slice().first()
    }

    #[inline]
    fn take(&mut self) -> Option<RankInfo<T>> {
        self.next()
    }

    #[inline]
    fn remaining(&self) -> usize {
        self.len()
    }
}

impl<T: Clone> MergeSource<T> for &[RankInfo<T>] {
    #[inline]
    fn peek(&self) -> Option<&RankInfo<T>> {
        self.first()
    }

    #[inline]
    fn take(&mut self) -> Option<RankInfo<T>> {
        let (head, tail) = self.split_first()?;
        *self = tail;
        Some(head.clone())
    }

    #[inline]
    fn remaining(&self) -> usize {
        self.len()
    }
}

/// Merges two sorted summaries into one, combining each element's rank bounds
/// with the bounds of its neighbours in the other summary (Zhang and Wang,
/// SSDBM 2007, section 3.1). Ranks are one-based throughout.
fn merge<T: Ord>(mut a: impl MergeSource<T>, mut b: impl MergeSource<T>) -> Vec<RankInfo<T>> {
    let mut s_m = Vec::with_capacity(a.remaining() + b.remaining());
    // Original bounds of the most recently consumed element on each side.
    let mut prev_a: Option<(i64, i64)> = None;
    let mut prev_b: Option<(i64, i64)> = None;

    loop {
        let take_a = match (a.peek(), b.peek()) {
            (Some(x), Some(y)) => x.val < y.val,
            (Some(_), None) => true,
            (None, Some(_)) => false,
            (None, None) => break,
        };

        if take_a {
            if let Some(x) = a.take() {
                let succ = b.peek().map(|y| (y.rmin, y.rmax));
                let (rmin, rmax) = merged_bounds((x.rmin, x.rmax), prev_b, succ);
                prev_a = Some((x.rmin, x.rmax));
                s_m.push(RankInfo::new(x.val, rmin, rmax));
            }
        } else if let Some(y) = b.take() {
            let succ = a.peek().map(|x| (x.rmin, x.rmax));
            let (rmin, rmax) = merged_bounds((y.rmin, y.rmax), prev_a, succ);
            prev_b = Some((y.rmin, y.rmax));
            s_m.push(RankInfo::new(y.val, rmin, rmax));
        }
    }

    s_m
}

/// Rank bounds of an element after merging, given the original bounds of its
/// predecessor and successor in the other summary.
fn merged_bounds(own: (i64, i64), pred: Option<(i64, i64)>, succ: Option<(i64, i64)>) -> (i64, i64) {
    let (rmin, rmax) = own;
    match (pred, succ) {
        (None, None) => (rmin, rmax),
        (None, Some((_, succ_max))) => (rmin, rmax + succ_max - 1),
        (Some((pred_min, _)), Some((_, succ_max))) => (rmin + pred_min, rmax + succ_max - 1),
        (Some((pred_min, pred_max)), None) => (rmin + pred_min, rmax + pred_max),
    }
}

/// Thins a merged summary down to about `block_size` elements spread evenly
/// over its rank range, keeping the first element that reaches each rank
/// boundary.
fn compress<T>(mut s0: Vec<RankInfo<T>>, block_size: usize, epsilon: f64) -> Vec<RankInfo<T>> {
    let mut s0_range = 0;
    let mut max_width = 0;
    for r in &s0 {
        s0_range = s0_range.max(r.rmax);
        max_width = max_width.max(r.rmax - r.rmin);
    }

    // Every input to compress must already be a 2 * epsilon-approximate
    // summary of its own range; the constructor's block sizing guarantees it.
    debug_assert!(
        2.0 * epsilon * (s0_range as f64) >= max_width as f64,
        "precision condition violated: range={s0_range}, max width={max_width}, epsilon={epsilon}"
    );

    let n = s0.len();
    let mut j = 0;
    let mut k = 0;
    for i in 0..=block_size {
        let r = ((i as f64) * (s0_range as f64) / (block_size as f64)).floor() as i64;

        while j < n && s0[j].rmax < r {
            j += 1;
        }
        if j >= n {
            break;
        }

        s0.swap(k, j);
        k += 1;
        j += 1;
    }

    s0.truncate(k);
    s0
}

/// Finds an element whose rank interval fits inside `[rank - epsilon_n, rank + epsilon_n]`.
///
/// The search starts at the first element whose `rmin` reaches `rank`, scans
/// forward, then scans backward over the elements whose `rmin` is still inside
/// the window. If no element satisfies both bounds, the element whose rank
/// interval is centred closest to `rank` is returned instead. Returns `None`
/// only for an empty summary.
fn find_idx<T: Clone + Ord>(s_m: &[RankInfo<T>], rank: i64, epsilon_n: i64) -> Option<T> {
    if s_m.is_empty() {
        return None;
    }

    let hi = rank + epsilon_n;
    let lo = rank - epsilon_n;
    let landing = s_m.partition_point(|e| e.rmin < rank);

    let mut i = landing;
    while i < s_m.len() && s_m[i].rmin <= hi {
        if s_m[i].rmax <= hi {
            return Some(s_m[i].val.clone());
        }
        i += 1;
    }

    let mut i = landing;
    while i > 0 && s_m[i - 1].rmin >= lo {
        i -= 1;
        if s_m[i].rmax <= hi {
            return Some(s_m[i].val.clone());
        }
    }

    s_m.iter()
        .min_by_key(|e| ((e.rmin + e.rmax) - 2 * rank).abs())
        .map(|e| e.val.clone())
}

/// Position of the `x`-th sub-stream boundary, `floor((2^x - 1) / epsilon)`,
/// saturating at `usize::MAX` once it no longer fits.
fn boundary(x: u32, epsilon: f64) -> usize {
    ((2f64.powi(x as i32) - 1.0) / epsilon).floor() as usize
}

/// An epsilon-approximate quantile summary for a stream of unknown size.
///
/// The stream is cut into sub-streams at positions `floor((2^x - 1) / epsilon)`
/// for `x = 1, 2, ...`, so each sub-stream is about twice as long as the
/// previous one. Every completed sub-stream is frozen into a compressed
/// [`FixedSizeEpsilonSummary`] built with error bound `epsilon / 2`, and a
/// query merges all of them with the summary of the current sub-stream. The
/// combined answer stays within `epsilon` of the requested rank.
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
    /// Index `x` of the boundary that ends the current sub-stream.
    next_boundary_index: u32,
    /// Stream position of that boundary, `floor((2^x - 1) / epsilon)`.
    next_boundary: usize,
    #[cfg_attr(feature = "serde", serde(skip))]
    cached_s_m: RefCell<Option<Vec<RankInfo<T>>>>,
}

impl<T> UnboundEpsilonSummary<T>
where
    T: Clone + Ord,
{
    /// Creates a summary with error bound `epsilon` for a stream of unknown
    /// length.
    ///
    /// A query for rank `r` after `n` insertions returns an element whose
    /// position in the sorted stream is within `epsilon * n` of `floor(r * n)`.
    ///
    /// # Errors
    ///
    /// Returns [`QuantileError::InvalidEpsilon`] if `epsilon` is not in
    /// `(0.0, 1.0]` (NaN and infinite values are rejected as well).
    pub fn new(epsilon: f64) -> Result<Self, QuantileError> {
        if !(epsilon > 0.0 && epsilon <= 1.0) {
            return Err(QuantileError::InvalidEpsilon);
        }

        // The first sub-stream ends at boundary 1, floor(1 / epsilon) >= 1.
        let next_boundary = boundary(1, epsilon);
        let s_c = FixedSizeEpsilonSummary::new(next_boundary, epsilon / 2.0)?;

        Ok(UnboundEpsilonSummary {
            epsilon,
            cnt: 0,
            s: Vec::new(),
            s_c,
            next_boundary_index: 1,
            next_boundary,
            cached_s_m: RefCell::new(None),
        })
    }

    /// Adds an element to the summary.
    ///
    /// When the element completes a sub-stream, that sub-stream's summary is
    /// compressed and frozen, and a new summary is started for the next one.
    pub fn update(&mut self, e: T) {
        self.cached_s_m.get_mut().take();
        self.s_c.update(e);
        self.cnt += 1;

        if self.cnt == self.next_boundary {
            // Zhang and Wang, SSDBM 2007, section 3.2: sub-stream P_x spans the
            // elements between boundaries x and x + 1.
            self.s_c.finalize();

            let upper = boundary(self.next_boundary_index + 1, self.epsilon);
            let n = upper - self.cnt;
            let next = FixedSizeEpsilonSummary::new(n, self.epsilon / 2.0)
                .expect("sub-stream length and epsilon are valid by construction");
            let finished = std::mem::replace(&mut self.s_c, next);
            self.s.push(finished);

            self.next_boundary_index += 1;
            self.next_boundary = upper;
        }
    }

    /// Returns an approximate quantile.
    ///
    /// `r` is a rank in `0.0..=1.0`: `0.0` is the minimum, `0.5` the median
    /// and `1.0` the maximum. The result is an inserted element whose position
    /// in the sorted stream is within `epsilon * n` of `floor(r * n)`, where
    /// `n` is the number of elements inserted so far.
    ///
    /// The first query after an [`update`](Self::update) rebuilds the merged
    /// summary; later queries reuse it until the next update.
    ///
    /// # Errors
    ///
    /// Returns [`QuantileError::InvalidRank`] if `r` is outside `0.0..=1.0` or
    /// NaN, and [`QuantileError::EmptySummary`] if nothing has been inserted.
    pub fn query(&self, r: f64) -> Result<T, QuantileError> {
        if !(0.0..=1.0).contains(&r) {
            return Err(QuantileError::InvalidRank);
        }
        if self.cnt == 0 {
            return Err(QuantileError::EmptySummary);
        }

        let mut cache = self.cached_s_m.borrow_mut();
        let s_m = cache.get_or_insert_with(|| self.merged_summaries());
        query_rank(s_m, r, self.cnt, self.epsilon).ok_or(QuantileError::EmptySummary)
    }

    /// Merges the current sub-stream summary with every frozen one.
    fn merged_summaries(&self) -> Vec<RankInfo<T>> {
        let mut s_m = self.s_c.calc_s_m();
        for summary in &self.s {
            for level in summary.s.iter().filter(|level| !level.is_empty()) {
                s_m = merge(s_m.into_iter(), level.as_slice());
            }
        }
        s_m
    }

    /// Returns the number of elements inserted so far.
    #[inline]
    #[must_use]
    pub fn size(&self) -> usize {
        self.cnt
    }
}

impl<T: Clone + Ord + fmt::Debug> fmt::Debug for UnboundEpsilonSummary<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("UnboundEpsilonSummary")
            .field("epsilon", &self.epsilon)
            .field("cnt", &self.cnt)
            .field("sub_streams", &(self.s.len() + 1))
            .field("next_boundary", &self.next_boundary)
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

        let merged = merge(s0.into_iter(), s1.as_slice());

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
    fn test_merge_of_exact_blocks_keeps_exact_ranks() {
        // Two sorted blocks with the one-based ranks update() assigns must
        // merge into exact ranks 1..=6; a 0/1-indexing mismatch shifts them.
        let mut a: Vec<RankInfo<i32>> = [1, 3, 5].iter().map(|&v| RankInfo::new(v, 0, 0)).collect();
        let mut b: Vec<RankInfo<i32>> = [2, 4, 6].iter().map(|&v| RankInfo::new(v, 0, 0)).collect();
        assign_exact_ranks(&mut a);
        assign_exact_ranks(&mut b);

        let merged = merge(a.into_iter(), b.as_slice());
        let ranks: Vec<(i32, i64, i64)> = merged.iter().map(|r| (r.val, r.rmin, r.rmax)).collect();
        assert_eq!(
            ranks,
            vec![(1, 1, 1), (2, 2, 2), (3, 3, 3), (4, 4, 4), (5, 5, 5), (6, 6, 6)]
        );

        let merged = merge(Vec::new().into_iter(), b.into_iter());
        assert_eq!(merged.iter().map(|r| r.rmin).collect::<Vec<_>>(), vec![1, 2, 3]);
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
    fn test_constructors_reject_epsilon_above_one() {
        assert!(matches!(
            FixedSizeEpsilonSummary::<usize>::new(100, 1.5),
            Err(QuantileError::InvalidEpsilon)
        ));
        assert!(matches!(
            UnboundEpsilonSummary::<usize>::new(1.5),
            Err(QuantileError::InvalidEpsilon)
        ));
        assert!(FixedSizeEpsilonSummary::<usize>::new(100, 1.0).is_ok());
        assert!(UnboundEpsilonSummary::<usize>::new(1.0).is_ok());
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
    fn test_fixedsize_short_streams_are_exact() {
        // Regression: these parameters used to derive a block size of 1 or 2
        // and answered rank 0.0 or 0.5 with the maximum element.
        for (n, epsilon) in [(13usize, 0.1f64), (101, 0.01), (24, 0.05), (16, 0.3)] {
            let mut s = FixedSizeEpsilonSummary::new(n, epsilon).unwrap();
            for i in 1..=n {
                s.update(i);
            }
            assert_eq!(s.query(0.0).unwrap(), 1, "n={n} epsilon={epsilon}");
            assert_eq!(s.query(0.5).unwrap(), n / 2 + 1, "n={n} epsilon={epsilon}");
            assert_eq!(s.query(1.0).unwrap(), n, "n={n} epsilon={epsilon}");
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
        assert_rank_error_against(&mut records, n, epsilon, query);
    }

    fn assert_rank_error_against(records: &mut [u64], n: usize, epsilon: f64, query: impl Fn(f64) -> u64) {
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
                "n={n} epsilon={epsilon} rank {rank}: value rank interval [{low}, {high}] misses target {target} by more than {allowed_error}"
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
    fn test_unbound_summary_respects_rank_error_contract_past_boundary_16() {
        let n = 2_000_000;
        let epsilon = 0.1;
        let summary = RefCell::new(UnboundEpsilonSummary::new(epsilon).unwrap());
        assert_rank_error(
            n,
            epsilon,
            |value| summary.borrow_mut().update(value),
            |rank| summary.borrow().query(rank).unwrap(),
        );
    }

    /// Sizes that step through every short-stream regime, including the
    /// `1 < epsilon * n < 5` window where the block-size formula degenerates.
    fn contract_grid() -> Vec<(usize, f64)> {
        let mut grid = Vec::new();
        for epsilon in [0.5, 0.3, 0.15, 0.085, 0.06, 0.03, 0.015] {
            let max_n = (60.0 / epsilon) as usize;
            let mut sizes: Vec<usize> = (1..=40).collect();
            sizes.extend((41..max_n).step_by((max_n / 25).max(1)));
            sizes.push(max_n);
            grid.extend(sizes.into_iter().map(|n| (n, epsilon)));
        }
        grid
    }

    #[test]
    fn test_fixedsize_summary_respects_rank_error_contract_on_short_streams() {
        for (n, epsilon) in contract_grid() {
            let summary = RefCell::new(FixedSizeEpsilonSummary::new(n, epsilon).unwrap());
            assert_rank_error(
                n,
                epsilon,
                |value| summary.borrow_mut().update(value),
                |rank| summary.borrow().query(rank).unwrap(),
            );
        }
    }

    #[test]
    fn test_unbound_summary_respects_rank_error_contract_on_short_streams() {
        for (n, epsilon) in contract_grid() {
            let summary = RefCell::new(UnboundEpsilonSummary::new(epsilon).unwrap());
            assert_rank_error(
                n,
                epsilon,
                |value| summary.borrow_mut().update(value),
                |rank| summary.borrow().query(rank).unwrap(),
            );
        }
    }

    #[test]
    fn test_unbound_summary_survives_zigzag_stream() {
        // Regression: this stream used to trip the compress precision
        // assertion at element 233 for epsilon = 0.03.
        let n = 2000usize;
        let epsilon = 0.03;
        let mut records: Vec<u64> = (0..n as u64)
            .map(|i| if i % 2 == 0 { i / 2 } else { n as u64 - i / 2 })
            .collect();
        let mut summary = UnboundEpsilonSummary::new(epsilon).unwrap();
        for &value in &records {
            summary.update(value);
        }
        assert_rank_error_against(&mut records, n, epsilon, |rank| summary.query(rank).unwrap());
    }

    #[test]
    fn test_boundary_positions_and_saturation() {
        assert_eq!(boundary(1, 0.01), 100);
        assert_eq!(boundary(2, 0.01), 300);
        assert_eq!(boundary(3, 0.01), 700);
        assert_eq!(boundary(1, 1.0), 1);
        assert_eq!(boundary(200, 0.5), usize::MAX);
    }

    #[test]
    fn test_unbound_summary_tracks_sub_stream_boundaries() {
        // epsilon = 1.0 places boundaries at 1, 3, 7, ..., 2^x - 1.
        let mut summary = UnboundEpsilonSummary::new(1.0).unwrap();
        for i in 1..=2000u64 {
            summary.update(i);
            let expected_frozen = (usize::BITS - (i as usize + 1).leading_zeros()) as usize - 1;
            assert_eq!(summary.s.len(), expected_frozen, "after {i} elements");
        }
        assert_eq!(summary.next_boundary_index, 11);
        assert_eq!(summary.next_boundary, 2047);
        assert_eq!(summary.size(), 2000);
        let median = summary.query(0.5).unwrap();
        assert!((1..=2000).contains(&median));
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
        assert_eq!(format!("{}", e), "epsilon must be in (0.0, 1.0]");
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
    fn test_find_idx_falls_back_to_nearest_interval() {
        // No interval fits inside [rank - 1, rank + 1]; the closest one wins
        // instead of the last element.
        let s_m = vec![
            RankInfo::new(10, 1, 4),
            RankInfo::new(20, 5, 9),
            RankInfo::new(30, 12, 16),
        ];
        assert_eq!(find_idx(&s_m, 7, 1), Some(20));
        assert_eq!(find_idx(&s_m, 3, 1), Some(10));
        // An element left of the landing point that fits the window is used.
        let s_m = vec![RankInfo::new(10, 4, 5), RankInfo::new(20, 7, 12)];
        assert_eq!(find_idx(&s_m, 6, 2), Some(10));
    }

    #[test]
    fn test_query_returns_error_on_empty_summary() {
        let s = FixedSizeEpsilonSummary::<usize>::new(10, 0.1).unwrap();
        assert!(matches!(s.query(0.5), Err(QuantileError::EmptySummary)));
        let s = UnboundEpsilonSummary::<usize>::new(0.1).unwrap();
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
        let mut deserialized: UnboundEpsilonSummary<usize> = serde_json::from_str(&serialized).unwrap();
        assert_eq!(s.query(0.5).unwrap(), deserialized.query(0.5).unwrap());
        // The restored summary keeps tracking sub-stream boundaries.
        for i in 101..=1000usize {
            s.update(i);
            deserialized.update(i);
        }
        assert_eq!(s.s.len(), deserialized.s.len());
        assert_eq!(s.query(0.5).unwrap(), deserialized.query(0.5).unwrap());
    }

    #[test]
    fn test_update_exactly_block_size() {
        let n = 400;
        let epsilon = 0.01;
        let mut s = FixedSizeEpsilonSummary::new(n, epsilon);
        let block_size = s.b;
        for i in 0..block_size {
            s.update(i);
        }

        assert_eq!(s.query(0.0), 0);
        assert_eq!(s.query(0.5), block_size / 2);
        assert_eq!(s.query(1.0), block_size - 1);
    }
}
