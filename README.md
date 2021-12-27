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
