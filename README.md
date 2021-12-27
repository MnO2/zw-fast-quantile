# zw-fast-quantile

[![CI](https://github.com/MnO2/zw-fast-quantile/actions/workflows/CI.yml/badge.svg)](https://github.com/MnO2/zw-fast-quantile/actions/workflows/CI.yml)
[![Crates.io](https://img.shields.io/crates/v/zw-fast-quantile.svg)](https://crates.io/crates/zw-fast-quantile)

Zhang Wang Fast Approximate Quantiles Algorithm in Rust. [Paper](http://web.cs.ucla.edu/~weiwang/paper/SSDBM07_2.pdf)

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
zw-fast-quantile = "0.2"
```

then you are good to go. If you are using Rust 2015 you have to ``extern crate zw-fast-quantile`` to your crate root as well.

## Example

There are two implementations in this library: `FixedSizeEpsilonSummary` is with the size of the stream known before hand, and `UnboundEpsilonSummary` is for the stream with unknown size. You can adjust the error rate `epsilon` of your own to trade-off between space and accuracy.

```rust
let epsilon = 0.1;
let n = 10;
let mut s = FixedSizeEpsilonSummary::new(n, epsilon);
for i in 1..=n {
    s.update(i);
}

let ans = s.query(0.0);
let expected = 1;
assert!(expected == ans);
```

```rust
let epsilon = 0.1;
let n = 10;
let mut s = UnboundEpsilonSummary::new(epsilon);
for i in 1..=n {
    s.update(i);
}

let ans = s.query(0.0);
let expected = 1;
assert!(expected == ans);
```

## Benchmarks

We benchmark against GK01 implemented in [quantiles](https://github.com/postmates/quantiles). To test the UPDATE, we insert 5000 values in order with error rate of 0.01. GK01 is roughly 2.6x slower than ZW.

```
zw unbound quantile update
                        time:   [60.780 us 60.855 us 60.936 us]
                        change: [-1.4032% -0.9510% -0.5005%] (p = 0.00 < 0.05)
                        Change within noise threshold.
Found 8 outliers among 100 measurements (8.00%)
  2 (2.00%) high mild
  6 (6.00%) high severe
```

```
gk quantile update      time:   [156.84 us 157.02 us 157.24 us]
                        change: [-0.1907% -0.0503% +0.0969%] (p = 0.50 > 0.05)
                        No change in performance detected.
Found 11 outliers among 100 measurements (11.00%)
  6 (6.00%) high mild
  5 (5.00%) high severe

```

As for the query time, we query 11 quantiles from 0.0, 0.1 ... 1.0. GK01 is about 1.5x slower than ZW.

```
zw unbound quantile query
                        time:   [229.62 ns 230.16 ns 230.77 ns]
                        change: [+1.3422% +1.8105% +2.2504%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 11 outliers among 100 measurements (11.00%)
  3 (3.00%) high mild
  8 (8.00%) high severe
```

```
gk quantile query       time:   [350.21 ns 350.48 ns 350.76 ns]
                        change: [-0.4638% -0.3109% -0.1670%] (p = 0.00 < 0.05)
                        Change within noise threshold.
Found 8 outliers among 100 measurements (8.00%)
  1 (1.00%) low severe
  2 (2.00%) high mild
  5 (5.00%) high severe
```