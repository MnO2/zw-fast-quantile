use std::cmp::Ordering;
use superslice::*;

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
        return RankInfo { val, rmin, rmax };
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

pub struct EpsilonSummary<T>
where
    T: Clone + Ord,
{
    epsilon: f64,
    b: usize,
    L: usize,
    cnt: usize,
    S: Vec<Vec<RankInfo<T>>>,
}

impl<T> EpsilonSummary<T>
where
    T: Clone + Ord + std::fmt::Debug,
{
    pub fn new(level: usize, epsilon: f64) -> Self {
        let b = ((level as f64) / epsilon + epsilon.ln() / (epsilon * 2.0_f64.ln())).floor() as usize;
        let S = vec![vec![]; level];

        EpsilonSummary {
            epsilon: epsilon,
            b: b,
            L: level,
            cnt: 0,
            S: S,
        }
    }

    pub fn update(&mut self, e: T) {
        let idx = self.S[0].upper_bound(&RankInfo::new(e.clone(), 0, 0));
        let rank_info = RankInfo::new(e, idx as i64, idx as i64);
        self.S[0].insert(idx, rank_info);

        self.cnt += 1;
        let n = self.S[0].len();
        if n < self.b {
            return;
        }

        for (i, r) in self.S[0].iter_mut().enumerate() {
            r.rmin = i as i64;
            r.rmax = i as i64;
        }

        let compressed_size = self.b / 2;
        let mut Sc = compress(&self.S[0], compressed_size, self.epsilon);

        self.S[0].clear();
        let L = self.L;
        for k in 1..L {
            if self.S[k].is_empty() {
                self.S[k] = Sc.clone();
            } else {
                let t = merge(&self.S[k], &Sc);
                Sc = compress(&t, compressed_size, self.epsilon);
                self.S[k].clear();
            }
        }
    }

    pub fn query(&mut self, e: T) -> f64 {
        for (i, r) in self.S[0].iter_mut().enumerate() {
            r.rmin = i as i64;
            r.rmax = i as i64;
        }

        let mut Sm = self.S[0].clone();
        for i in 1..self.L {
            Sm = merge(&Sm, &self.S[i])
        }

        let mut i = 0;
        while i < Sm.len() {
            if Sm[i].val >= e {
                break;
            }

            i += 1;
        }

        let quantile: f64 = ((Sm[i].rmin + Sm[i].rmax) as f64) / (2.0_f64 * self.cnt as f64);
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

fn merge<T: Clone + Ord>(Sa: &Vec<RankInfo<T>>, Sb: &Vec<RankInfo<T>>) -> Vec<RankInfo<T>> {
    if Sa.len() == 0 {
        return Sb.clone();
    }

    if Sb.len() == 0 {
        return Sa.clone();
    }

    let mut Sm = Vec::new();

    let mut i1 = 0;
    let mut i2 = 0;
    let mut from = 0;

    while i1 < Sa.len() || i2 < Sb.len() {
        let mut val;
        let mut rmin;
        let mut rmax;

        if i1 < Sa.len() && i2 < Sb.len() {
            if Sa[i1].val < Sb[i2].val {
                val = Sa[i1].val.clone();
                from = 1;
            } else {
                val = Sb[i2].val.clone();
                from = 2;
            }
        } else if i1 < Sa.len() && i2 >= Sb.len() {
            val = Sa[i1].val.clone();
            from = 1;
        } else {
            val = Sb[i2].val.clone();
            from = 2;
        }

        if from == 1 {
            if 0 < i2 && i2 < Sb.len() {
                rmin = Sa[i1].rmin + Sb[i2 - 1].rmin;
                rmax = if Sa[i1].rmax + Sb[i2].rmax > 0 { Sa[i1].rmax + Sb[i2].rmax - 1  } else { 0 }
            } else if i2 == 0 {
                rmin = Sa[i1].rmin;
                rmax = if Sa[i1].rmax + Sb[i2].rmax > 0 { Sa[i1].rmax + Sb[i2].rmax - 1 } else { 0 }
            } else {
                rmin = Sa[i1].rmin + Sb[i2 - 1].rmin;
                rmax = Sa[i1].rmax + Sb[i2 - 1].rmax;
            }

            i1 += 1;
        } else {
            if 0 < i1 && i1 < Sa.len() {
                rmin = Sa[i1 - 1].rmin + Sb[i2].rmin;
                rmax = if Sa[i1].rmax + Sb[i2].rmax > 0 { Sa[i1].rmax + Sb[i2].rmax - 1 } else { 0 };
            } else if i1 == 0 {
                rmin = Sb[i2].rmin;
                rmax = if Sa[i1].rmax + Sb[i2].rmax > 0 { Sa[i1].rmax + Sb[i2].rmax - 1 } else { 0 };
            } else {
                rmin = Sa[i1 - 1].rmin + Sb[i2].rmin;
                rmax = Sa[i1 - 1].rmax + Sb[i2].rmax;
            }

            i2 += 1;
        }

        let rank_info = RankInfo::new(val, rmin, rmax);
        Sm.push(rank_info);
    }

    Sm
}

fn compress<T: Clone + std::fmt::Debug>(S0: &Vec<RankInfo<T>>, B: usize, epsilon: f64) -> Vec<RankInfo<T>> {
    let mut Sc = Vec::new();

    let mut S0_range = 0;
    let mut e: f64 = 0.0;

    for r in S0 {
        if S0_range < r.rmax {
            S0_range = r.rmax;
        }

        if (r.rmax - r.rmin) as f64 > e {
            e = (r.rmax - r.rmin) as f64;
        }
    }

    let epsilonN: f64 = epsilon * (S0_range as f64);
    assert!(2.0 * epsilonN >= e, "precision condition violated.");

    let mut i = 0;
    let mut j = 0;
    while i <= B && j < S0.len() {
        let r = ((i as f64) * (S0_range as f64) / (B as f64)).floor() as i64;


        while j < S0.len() {
            if S0[j].rmax >= r {
                break;
            }

            j += 1;
        }

        assert!(j < S0.len(), "unable to find the summary with precision.");
        Sc.push(S0[j].clone());
        j += 1;
        i += 1;
    }

    return Sc;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_and_compress() {
        let mut S0 = Vec::new();
        let mut S1 = Vec::new();

        S0.push(RankInfo::new(2, 1, 1));
        S0.push(RankInfo::new(4, 3, 4));
        S0.push(RankInfo::new(8, 5, 6));
        S0.push(RankInfo::new(17, 8, 8));

        S1.push(RankInfo::new(1, 1, 1));
        S1.push(RankInfo::new(7, 3, 3));
        S1.push(RankInfo::new(12, 5, 6));
        S1.push(RankInfo::new(15, 8, 8));

        let merged = merge(&S0, &S1);

        assert_eq!(merged.len(), 8);
        let merged_vals: Vec<i32> = merged.iter().map(|x| x.val ).collect();
        let merged_rmins: Vec<usize> = merged.iter().map(|x| x.rmin ).collect();
        let merged_rmaxs: Vec<usize> = merged.iter().map(|x| x.rmax ).collect();
        assert_eq!(merged_vals, vec![1,2,4,7,8,12,15,17]);
        assert_eq!(merged_rmins, vec![1,2,4,6,8,10,13,16]);
        assert_eq!(merged_rmaxs, vec![1,3,6,8,11,13,15,16]);

        let epsilon: f64 = 0.2;
        let compressed = compress(&merged, 4, epsilon);
        let compressed_vals: Vec<i32> = compressed.iter().map(|x| x.val ).collect();
        assert_eq!(compressed_vals, vec![1,4,7,12,17]);
    }

    #[test]
    fn test_query() {
        let N = 100;
        let epsilon: f64 = 0.1;
        let L = ((N as f64).ln() / 2.0_f64.ln()).floor() as usize + 2;
        let mut S = EpsilonSummary::new(L, epsilon);

        let mut records = Vec::with_capacity(N);
        let mut quantile_ans: Vec<f64> = Vec::with_capacity(N);
        for i in 0..N {
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

        for i in 0..N {
            S.update(records[i]);

            let quantile_estimated = S.query(records[i]);
            dbg!(records[i]);
            dbg!(quantile_estimated);
            dbg!(quantile_ans[i]);
            assert!((quantile_ans[i] - quantile_estimated).abs() < epsilon);
        }
    }
}
