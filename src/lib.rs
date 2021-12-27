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
//! let epsilon = 0.1;
//! let n = 10;
//! let mut s = FixedSizeEpsilonSummary::new(n, epsilon);
//! for i in 1..=n {
//!     s.update(i);
//! }
//!
//! let ans = s.query(0.0);
//! let expected = 1;
//! assert!(expected == ans);
//! ```
//!
//! //! ```rust
//! use zw_fast_quantile::UnboundEpsilonSummary;
//!
//! let epsilon = 0.1;
//! let n = 10;
//! let mut s = UnboundEpsilonSummary::new(epsilon);
//! for i in 1..=n {
//!     s.update(i);
//! }
//!
//! let ans = s.query(0.0);
//! let expected = 1;
//! assert!(expected == ans);
//! ```
//!
use std::cmp::Ordering;

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
    cached_s_m: Option<Vec<RankInfo<T>>>,
}

impl<T> FixedSizeEpsilonSummary<T>
where
    T: Clone + Ord,
{
    pub fn new(n: usize, epsilon: f64) -> Self {
        // number_of_leves = ceil(log2(n))
        // block_size = floor(log2(epsilon * N) / epsilon)

        let epsilon_n: f64 = (n as f64) * epsilon;
        let number_of_levels = if epsilon_n > 1.0 {
            (epsilon_n as f64).log2().ceil() as usize
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

        FixedSizeEpsilonSummary {
            epsilon,
            b: block_size,
            level: number_of_levels,
            cnt: 0,
            s,
            cached_s_m: None,
        }
    }

    pub fn update(&mut self, e: T) {
        self.cached_s_m = None;
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
        for k in 1..self.level {
            if self.s[k].is_empty() {
                self.s[k] = s_c;
                break;
            } else {
                let t = merge(s_c, &self.s[k]);
                s_c = compress(t, compressed_size, self.epsilon);
                self.s[k].clear();
            }
        }
    }

    pub fn query(&mut self, r: f64) -> T {
        if self.cached_s_m.is_none() {
            self.s[0].sort();
            for (i, r) in self.s[0].iter_mut().enumerate() {
                r.rmin = i as i64;
                r.rmax = i as i64;
            }

            let mut s_m = self.s[0].clone();
            for i in 1..self.level {
                s_m = merge(s_m, &self.s[i])
            }

            self.cached_s_m = Some(s_m);
        }

        let rank: i64 = ((self.cnt as f64) * r).floor() as i64;
        let epsilon_n: i64 = ((self.cnt as f64) * self.epsilon).floor() as i64;
        let e = find_idx(&self.cached_s_m.as_ref().unwrap(), rank, epsilon_n).unwrap();
        return e;
    }

    fn calc_s_m(&mut self, epsilon: f64) -> Vec<RankInfo<T>> {
        self.s[0].sort();
        for (i, r) in self.s[0].iter_mut().enumerate() {
            r.rmin = i as i64;
            r.rmax = i as i64;
        }

        let mut s_m = self.s[0].clone();
        for i in 1..self.level {
            s_m = merge(s_m, &self.s[i])
        }

        compress(s_m, self.b, epsilon)
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
    let mut l = 0usize;
    let mut r = s_m.len() - 1;
    while l < r {
        let m = l + (r - l) / 2;
        if s_m[m].rmin < rank {
            l = m + 1;
        } else if s_m[m].rmax > rank {
            r = m;
        } else {
            r = m;
        }
    }

    while l < s_m.len() && s_m[l].rmin >= rank - epsilon_n {
        if s_m[l].rmax <= rank + epsilon_n {
            return Some(s_m[l].val.clone());
        }

        l += 1;
    }

    None
}

pub struct UnboundEpsilonSummary<T>
where
    T: Clone + Ord,
{
    epsilon: f64,
    cnt: usize,
    s: Vec<FixedSizeEpsilonSummary<T>>,
    s_c: FixedSizeEpsilonSummary<T>,
    boundaries: [usize; 32],
    cached_s_m: Option<Vec<RankInfo<T>>>,
}

impl<T> UnboundEpsilonSummary<T>
where
    T: Clone + Ord + std::fmt::Debug,
{
    pub fn new(epsilon: f64) -> Self {
        let s = vec![];

        let n = (1.0_f64 / epsilon).floor() as usize;
        let s_c = FixedSizeEpsilonSummary::new(n, epsilon / 2.0);
        let mut boundaries: [usize; 32] = [0; 32];
        for i in 0..32usize {
            let boundary = (((usize::pow(2, i as u32) - 1) as f64) / epsilon).floor() as usize;
            boundaries[i] = boundary
        }

        UnboundEpsilonSummary {
            epsilon,
            cnt: 0,
            s,
            s_c,
            boundaries,
            cached_s_m: None,
        }
    }

    pub fn update(&mut self, e: T) {
        self.cached_s_m = None;
        self.s_c.update(e);

        if let Some(x) = is_boundary(self.cnt + 1, &self.boundaries) {
            self.s_c.finalize(self.epsilon / 2.0);

            let upper_bound = (((usize::pow(2, (x + x) as u32) - 1) as f64) / self.epsilon).floor() as usize;
            let n = upper_bound - self.cnt - 1;
            let mut summary = FixedSizeEpsilonSummary::new(n, self.epsilon / 2.0);
            std::mem::swap(&mut self.s_c, &mut summary);

            self.s.push(summary);
        }
        self.cnt += 1;
    }

    pub fn query(&mut self, r: f64) -> T {
        if self.cached_s_m.is_none() {
            let mut s_m = self.s_c.calc_s_m(self.epsilon / 2.0);
            for i in 0..self.s.len() {
                for j in 0..self.s[i].s.len() {
                    s_m = merge(s_m, &self.s[i].s[j])
                }
            }

            self.cached_s_m = Some(s_m);
        }

        let rank: i64 = ((self.cnt as f64) * r).floor() as i64;
        let epsilon_n: i64 = ((self.cnt as f64) * self.epsilon).floor() as i64;
        let e = find_idx(&self.cached_s_m.as_ref().unwrap(), rank, epsilon_n).unwrap();
        return e;
    }

    #[inline]
    pub fn size(&self) -> usize {
        self.cnt
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_distr::Distribution;

    // INSTRUMENTED_SYSTEM is an instrumented instance of the system allocator
    #[global_allocator]
    static GLOBAL: &stats_alloc::StatsAlloc<std::alloc::System> = &stats_alloc::INSTRUMENTED_SYSTEM;

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
    fn test_query_with_small_n_on_fixedsize_summary() {
        let epsilon = 0.1;
        let n = 10;
        let mut s = FixedSizeEpsilonSummary::new(n, epsilon);
        for i in 1..=n {
            s.update(i);
        }

        for i in 1..=n {
            let rank: f64 = ((i - 1) as f64) / (n as f64);
            let ans = s.query(rank);
            assert!(i == ans);
        }
    }

    #[test]
    fn test_query_with_small_n_on_unbound_summary() {
        let epsilon = 0.1;
        let n = 10;
        let mut s = UnboundEpsilonSummary::new(epsilon);
        for i in 1..=n {
            s.update(i);
        }

        for i in 1..=n {
            let rank: f64 = ((i - 1) as f64) / (n as f64);
            let ans = s.query(rank);
            assert!(i == ans);
        }
    }

    #[test]
    fn test_normal_distribution_generated_seq_on_fixed_summary() {
        let n = 1000000;
        let epsilon: f64 = 0.01;
        let mut s = FixedSizeEpsilonSummary::new(n, epsilon);

        let dn = rand_distr::Normal::new(0.5f64, 0.2f64).unwrap();
        let mut normal_rng = rand_pcg::Pcg64::new(0xcafef00dd15ea5e5, 0xa02bdbf7bb3c0a7ac28fa16a64abf96);

        let startmem = stats_alloc::Region::new(&GLOBAL);
        let mut records = Vec::with_capacity(n);
        for _ in 0..n {
            let x = dn.sample(&mut normal_rng);
            records.push(unsafe { ordered_float::NotNan::new_unchecked(x) });
        }
        let mem = startmem.change();
        let bytes_change = mem.bytes_allocated as isize - mem.bytes_deallocated as isize + mem.bytes_reallocated;

        println!("{}k", bytes_change / 1024);

        records.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let startmem = stats_alloc::Region::new(&GLOBAL);
        for i in 0..n {
            s.update(records[i]);
        }
        let mem = startmem.change();
        let bytes_change = mem.bytes_allocated as isize - mem.bytes_deallocated as isize + mem.bytes_reallocated;

        println!("{}k", bytes_change / 1024);

        let quantile_estimated = s.query(0.5);
        assert!((quantile_estimated - records[n / 2]).abs() < 0.01);

        let quantile_estimated = s.query(0.0);
        assert!((quantile_estimated - records[0]).abs() < 0.1);

        let quantile_estimated = s.query(0.99);
        assert!((quantile_estimated - records[n * 99 / 100]).abs() < 0.01);

        let quantile_estimated = s.query(1.0);
        assert!((quantile_estimated - records[n - 1]).abs() < 0.01);
    }

    #[test]
    fn test_pareto_distribution_generated_seq_on_fixed_summary() {
        let n = 1000000;
        let epsilon: f64 = 0.001;
        let mut s = FixedSizeEpsilonSummary::new(n, epsilon);

        let dn = rand_distr::Pareto::new(5f64, 10f64).unwrap();
        let mut normal_rng = rand_pcg::Pcg64::new(0xcafef00dd15ea5e5, 0xa02bdbf7bb3c0a7ac28fa16a64abf96);

        let startmem = stats_alloc::Region::new(&GLOBAL);
        let mut records = Vec::with_capacity(n);
        for _ in 0..n {
            let x = dn.sample(&mut normal_rng);
            records.push(unsafe { ordered_float::NotNan::new_unchecked(x) });
        }
        let mem = startmem.change();
        let bytes_change = mem.bytes_allocated as isize - mem.bytes_deallocated as isize + mem.bytes_reallocated;

        println!("{}k", bytes_change / 1024);

        records.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let startmem = stats_alloc::Region::new(&GLOBAL);
        for i in 0..n {
            s.update(records[i]);
        }
        let mem = startmem.change();
        let bytes_change = mem.bytes_allocated as isize - mem.bytes_deallocated as isize + mem.bytes_reallocated;

        println!("{}k", bytes_change / 1024);

        let quantile_estimated = s.query(0.5);
        assert!((quantile_estimated - records[n / 2]).abs() < 0.01);

        let quantile_estimated = s.query(0.0);
        assert!((quantile_estimated - records[0]).abs() < 0.01);

        let quantile_estimated = s.query(0.99);
        assert!((quantile_estimated - records[n * 99 / 100]).abs() < 0.01);

        let quantile_estimated = s.query(1.0);
        assert!((quantile_estimated - records[n - 1]).abs() < 0.01);
    }

    #[test]
    fn test_normal_distribution_generated_seq_on_unbound_summary() {
        let n = 1000000;
        let epsilon: f64 = 0.01;
        let mut s = UnboundEpsilonSummary::new(epsilon);

        let dn = rand_distr::Normal::new(0.5f64, 0.2f64).unwrap();
        let mut normal_rng = rand_pcg::Pcg64::new(0xcafef00dd15ea5e5, 0xa02bdbf7bb3c0a7ac28fa16a64abf96);

        let startmem = stats_alloc::Region::new(&GLOBAL);
        let mut records = Vec::with_capacity(n);
        for _ in 0..n {
            let x = dn.sample(&mut normal_rng);
            records.push(unsafe { ordered_float::NotNan::new_unchecked(x) });
        }
        let mem = startmem.change();
        let bytes_change = mem.bytes_allocated as isize - mem.bytes_deallocated as isize + mem.bytes_reallocated;

        println!("{}k", bytes_change / 1024);

        records.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let startmem = stats_alloc::Region::new(&GLOBAL);
        for i in 0..n {
            s.update(records[i]);
        }
        let mem = startmem.change();
        let bytes_change = mem.bytes_allocated as isize - mem.bytes_deallocated as isize + mem.bytes_reallocated;

        println!("{}k", bytes_change / 1024);

        let quantile_estimated = s.query(0.5);
        assert!((quantile_estimated - records[n / 2]).abs() < 0.01);

        let quantile_estimated = s.query(0.0);
        assert!((quantile_estimated - records[0]).abs() < 0.01);

        let quantile_estimated = s.query(0.99);
        assert!((quantile_estimated - records[n * 99 / 100]).abs() < 0.01);

        let quantile_estimated = s.query(1.0);
        assert!((quantile_estimated - records[n - 1]).abs() < 0.01);
    }

    #[test]
    fn test_pareto_distribution_generated_seq_on_unbound_summary() {
        let n = 1000000;
        let epsilon: f64 = 0.001;
        let mut s = UnboundEpsilonSummary::new(epsilon);

        let dn = rand_distr::Pareto::new(5f64, 10f64).unwrap();
        let mut normal_rng = rand_pcg::Pcg64::new(0xcafef00dd15ea5e5, 0xa02bdbf7bb3c0a7ac28fa16a64abf96);

        let startmem = stats_alloc::Region::new(&GLOBAL);
        let mut records = Vec::with_capacity(n);
        for _ in 0..n {
            let x = dn.sample(&mut normal_rng);
            records.push(unsafe { ordered_float::NotNan::new_unchecked(x) });
        }
        let mem = startmem.change();
        let bytes_change = mem.bytes_allocated as isize - mem.bytes_deallocated as isize + mem.bytes_reallocated;

        println!("{}k", bytes_change / 1024);

        records.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let startmem = stats_alloc::Region::new(&GLOBAL);
        for i in 0..n {
            s.update(records[i]);
        }
        let mem = startmem.change();
        let bytes_change = mem.bytes_allocated as isize - mem.bytes_deallocated as isize + mem.bytes_reallocated;

        println!("{}k", bytes_change / 1024);

        let quantile_estimated = s.query(0.5);
        assert!((quantile_estimated - records[n / 2]).abs() < 0.01);

        let quantile_estimated = s.query(0.0);
        assert!((quantile_estimated - records[0]).abs() < 0.01);

        let quantile_estimated = s.query(0.99);
        assert!((quantile_estimated - records[n * 99 / 100]).abs() < 0.01);

        let quantile_estimated = s.query(1.0);
        assert!((quantile_estimated - records[n - 1]).abs() < 0.01);
    }
}
