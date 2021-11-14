# zw-fast-quantile

Zhang Wang Fast Approximate Quantiles Algorithm in Rust. [Paper](http://web.cs.ucla.edu/~weiwang/paper/SSDBM07_2.pdf)

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
zw-fast-quantile = "0.1"
```

then you are good to go. If you are using Rust 2015 you have to ``extern crate zw-fast-quantile`` to your crate root as well.

## Example

There are two implementations in this library: `FixedSizeEpsilonSummary` is with the size of the stream known before hand, and `UnboundEpsilonSummary` is for the stream with unknown size. You can adjust the error rate `epsilon` of your own to trade-off between space and accuracy.

```rust
let epsilon = 0.2;
let n = 10;
let mut s = FixedSizeEpsilonSummary::new(n, epsilon);
for i in 1..=n {
    s.update(i);
}

for i in 1..=n {
    let ans = s.query(i);
    let expected: f64 = ((i - 1) as f64) / (n as f64);
    let error: f64 = (expected - ans).abs();
    assert!(error < epsilon);
}
```

```rust
let epsilon = 0.01;
let n = 1000;
let mut s = UnboundEpsilonSummary::new(epsilon);
for i in 1..=n {
    s.update(i);
}

for i in 1..=n {
    let ans = s.query(i);
    let expected: f64 = ((i - 1) as f64) / (n as f64);
    let error: f64 = (expected - ans).abs();
    assert!(error < epsilon);
}
```