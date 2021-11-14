//! Zhang Wang Fast Approximate Quantiles Algorithm in Rust
//!
//! ## Installation
//!
//! Add this to your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! zw-fast-quantile = "0.1"
//! ```
//!
//! ## Example
//!
//! ```rust
//! use zw_fast_quantile::FixedSizeEpsilonSummary;
//!
//! let epsilon = 0.2;
//! let n = 10;
//! let mut s = FixedSizeEpsilonSummary::new(n, epsilon).unwrap();
//! for i in 1..=n {
//!     s.update(i);
//! }
//!
//! let ans = s.query(1);
//! let expected: f64 = 0.0;
//!
//! let error: f64 = (expected - ans).abs();
//! assert!(error < epsilon);
//! ```
//!
use std::cmp::Ordering;
use std::result::Result;

#[derive(Debug, Clone)]
pub enum FixedSizeEpsilonSummaryError {
    EpsilonTooSmall,
}

#[derive(Debug, Clone)]
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
        self.val.partial_cmp(&other.val)
    }
}

impl<T: Clone + PartialEq> PartialEq for RankInfo<T> {
    fn eq(&self, other: &Self) -> bool {
        self.val == other.val
    }
}

impl<T: Clone + PartialEq> Eq for RankInfo<T> {}

#[derive(Clone)]
pub struct FixedSizeEpsilonSummary<T>
where
    T: Clone + Ord,
{
    epsilon: f64,
    b: usize,
    level: usize,
    cnt: usize,
    s: Vec<Vec<RankInfo<T>>>,
}

impl<T> FixedSizeEpsilonSummary<T>
where
    T: Clone + Ord,
{
    pub fn new(n: usize, epsilon: f64) -> Result<Self, FixedSizeEpsilonSummaryError> {
        // number_of_leves = ceil(log2(n))
        // block_size = floor(log2(epsilon * N) / epsilon)

        let epsilon_n: f64 = (n as f64) * epsilon;
        if epsilon_n <= 1.0 {
            return Err(FixedSizeEpsilonSummaryError::EpsilonTooSmall);
        }

        let number_of_levels = (epsilon_n as f64).log2().ceil() as usize;
        let block_size = if number_of_levels > 1 {
            (epsilon_n.log2() / epsilon).floor() as usize
        } else {
            n + 1
        };

        let s = vec![vec![]; number_of_levels];

        Ok(FixedSizeEpsilonSummary {
            epsilon,
            b: block_size,
            level: number_of_levels,
            cnt: 0,
            s,
        })
    }

    pub fn update(&mut self, e: T) {
        let rank_info = RankInfo::new(e, 0, 0);
        self.s[0].push(rank_info);

        self.cnt += 1;
        let n = self.s[0].len();
        if n < self.b {
            return;
        }

        self.s[0].sort();
        for (i, r) in self.s[0].iter_mut().enumerate() {
            r.rmin = i as i64;
            r.rmax = i as i64;
        }

        let compressed_size = self.b / 2;
        let mut s_c = compress(&self.s[0], compressed_size, self.epsilon);

        self.s[0].clear();
        for k in 1..self.level {
            if self.s[k].is_empty() {
                self.s[k] = s_c.clone();
                break;
            } else {
                let t = merge(&self.s[k], &s_c);
                s_c = compress(&t, compressed_size, self.epsilon);
                self.s[k].clear();
            }
        }
    }

    pub fn query(&mut self, e: T) -> f64 {
        self.s[0].sort();
        for (i, r) in self.s[0].iter_mut().enumerate() {
            r.rmin = i as i64;
            r.rmax = i as i64;
        }

        let mut s_m = self.s[0].clone();
        for i in 1..self.level {
            s_m = merge(&s_m, &self.s[i])
        }

        let mut i = 0;
        while i < s_m.len() {
            if s_m[i].val >= e {
                break;
            }

            i += 1;
        }

        let quantile: f64 = ((s_m[i].rmin + s_m[i].rmax) as f64) / (2.0_f64 * self.cnt as f64);
        if quantile < 0.0_f64 {
            0.0_f64
        } else {
            quantile
        }
    }

    fn calc_s_m(&mut self, epsilon: f64) -> Vec<RankInfo<T>> {
        self.s[0].sort();
        for (i, r) in self.s[0].iter_mut().enumerate() {
            r.rmin = i as i64;
            r.rmax = i as i64;
        }

        let mut s_m = self.s[0].clone();
        for i in 1..self.level {
            s_m = merge(&s_m, &self.s[i])
        }

        compress(&s_m, self.b, epsilon)
    }

    fn finalize(&mut self, epsilon: f64) {
        let s_m = self.calc_s_m(epsilon);
        self.s.clear();
        self.s.push(s_m);
    }

    #[inline]
    pub fn size(&self) -> usize {
        self.cnt
    }
}

fn merge<T: Clone + Ord>(s_a: &[RankInfo<T>], s_b: &[RankInfo<T>]) -> Vec<RankInfo<T>> {
    if s_a.is_empty() {
        return s_b.to_vec();
    }

    if s_b.is_empty() {
        return s_a.to_vec();
    }

    let mut s_m = Vec::new();

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

fn compress<T: Clone>(s0: &[RankInfo<T>], block_size: usize, epsilon: f64) -> Vec<RankInfo<T>> {
    let mut s_c = Vec::new();

    let mut s0_range = 0;
    let mut e: f64 = 0.0;

    for r in s0 {
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
    while i <= block_size && j < s0.len() {
        let r = ((i as f64) * (s0_range as f64) / (block_size as f64)).floor() as i64;

        while j < s0.len() {
            if s0[j].rmax >= r {
                break;
            }

            j += 1;
        }

        assert!(j < s0.len(), "unable to find the summary with precision given.");
        s_c.push(s0[j].clone());
        j += 1;
        i += 1;
    }

    s_c
}

fn is_power_of_two(x: usize) -> bool {
    (x & (x - 1)) == 0
}

#[derive(Debug, Clone)]
pub enum UnboundEpsilonSummaryError {
    EpsilonTooSmall,
}

pub struct UnboundEpsilonSummary<T>
where
    T: Clone + Ord,
{
    epsilon: f64,
    cnt: usize,
    s: Vec<FixedSizeEpsilonSummary<T>>,
    s_c: FixedSizeEpsilonSummary<T>,
}

impl<T> UnboundEpsilonSummary<T>
where
    T: Clone + Ord,
{
    pub fn new(epsilon: f64) -> Result<Self, UnboundEpsilonSummaryError> {
        let s = vec![];

        let n = (1.0_f64 / epsilon).floor() as usize;
        let s_c = FixedSizeEpsilonSummary::new(n, epsilon / 2.0).unwrap();
        Ok(UnboundEpsilonSummary {
            epsilon,
            cnt: 0,
            s,
            s_c,
        })
    }

    pub fn update(&mut self, e: T) {
        let x = (((self.cnt + 1) as f64) * self.epsilon + 1.0).ceil() as usize;

        if is_power_of_two(x) {
            self.s_c.finalize(self.epsilon / 2.0);
            self.s.push(self.s_c.clone());

            let upper_bound = (((x + x - 1) as f64) / self.epsilon).floor() as usize;
            let n = upper_bound - self.cnt - 1;
            let summary = FixedSizeEpsilonSummary::new(n, self.epsilon / 2.0).unwrap();
            self.s_c = summary;
        } else {
            self.s_c.update(e);
        }
    }

    pub fn query(&mut self, e: T) -> f64 {
        let mut s_m = self.s_c.calc_s_m(self.epsilon / 2.0);

        for i in 0..self.s.len() {
            for j in 0..self.s[i].s.len() {
                s_m = merge(&s_m, &self.s[i].s[j])
            }
        }

        let mut i = 0;
        while i < s_m.len() {
            if s_m[i].val >= e {
                break;
            }

            i += 1;
        }

        let quantile: f64 = ((s_m[i].rmin + s_m[i].rmax) as f64) / (2.0_f64 * self.cnt as f64);
        if quantile < 0.0_f64 {
            0.0_f64
        } else {
            quantile
        }
    }

    #[inline]
    pub fn size(&self) -> usize {
        self.cnt
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    #[test]
    fn test_merge_and_compress() {
        let mut s0 = Vec::new();
        let mut s1 = Vec::new();

        s0.push(RankInfo::new(2, 1, 1));
        s0.push(RankInfo::new(4, 3, 4));
        s0.push(RankInfo::new(8, 5, 6));
        s0.push(RankInfo::new(17, 8, 8));

        s1.push(RankInfo::new(1, 1, 1));
        s1.push(RankInfo::new(7, 3, 3));
        s1.push(RankInfo::new(12, 5, 6));
        s1.push(RankInfo::new(15, 8, 8));

        let merged = merge(&s0, &s1);

        assert_eq!(merged.len(), 8);
        let merged_vals: Vec<i32> = merged.iter().map(|x| x.val).collect();
        let merged_rmins: Vec<i64> = merged.iter().map(|x| x.rmin).collect();
        let merged_rmaxs: Vec<i64> = merged.iter().map(|x| x.rmax).collect();
        assert_eq!(merged_vals, vec![1, 2, 4, 7, 8, 12, 15, 17]);
        assert_eq!(merged_rmins, vec![1, 2, 4, 6, 8, 10, 13, 16]);
        assert_eq!(merged_rmaxs, vec![1, 3, 6, 8, 11, 13, 15, 16]);

        let epsilon: f64 = 0.2;
        let compressed = compress(&merged, 4, epsilon);
        let compressed_vals: Vec<i32> = compressed.iter().map(|x| x.val).collect();
        assert_eq!(compressed_vals, vec![1, 4, 7, 12, 17]);
    }

    #[test]
    fn test_randomly_generated_seq_on_fixedsize_summary() {
        let mut rng = rand::thread_rng();
        let n = rng.gen_range(100..10000);
        let epsilon: f64 = rng.gen_range(0.01..0.2);

        let mut s = FixedSizeEpsilonSummary::new(n, epsilon).unwrap();

        let mut records = Vec::with_capacity(n);
        let mut quantile_ans: Vec<f64> = Vec::with_capacity(n);
        for i in 0..n {
            let x = rand::random::<u32>();
            records.push(x);

            let mut real_rank = 0;
            for j in 0..i {
                if records[j] <= records[i] {
                    real_rank += 1;
                }
            }

            quantile_ans.push((real_rank as f64) / ((i + 1) as f64))
        }

        for i in 0..n {
            s.update(records[i]);

            let quantile_estimated = s.query(records[i]);
            assert!((quantile_ans[i] - quantile_estimated).abs() < epsilon);
        }
    }

    #[test]
    fn test_query_with_small_n_on_fixedsize_summary() {
        let epsilon = 0.2;
        let n = 10;
        let mut s = FixedSizeEpsilonSummary::new(n, epsilon).unwrap();
        for i in 1..=n {
            s.update(i);
        }

        for i in 1..=n {
            let ans = s.query(i);
            let expected: f64 = (i as f64) / (n as f64);
            let error: f64 = (expected - ans).abs();
            assert!(error < epsilon);
        }
    }

    #[test]
    fn test_query_on_unbound_summary() {
        let epsilon = 0.2;
        let n = 100;
        let mut s = UnboundEpsilonSummary::new(epsilon).unwrap();
        for i in 1..=n {
            s.update(i);
        }

        for i in 1..=n {
            let ans = s.query(i);
            let expected: f64 = (i as f64) / (n as f64);
            let error: f64 = (expected - ans).abs();
            assert!(error < epsilon);
        }
    }
}
