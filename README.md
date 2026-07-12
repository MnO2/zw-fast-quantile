# zw-fast-quantile

[![CI](https://github.com/MnO2/zw-fast-quantile/actions/workflows/CI.yml/badge.svg)](https://github.com/MnO2/zw-fast-quantile/actions/workflows/CI.yml)
[![Crates.io](https://img.shields.io/crates/v/zw-fast-quantile.svg)](https://crates.io/crates/zw-fast-quantile)
[![docs.rs](https://docs.rs/zw-fast-quantile/badge.svg)](https://docs.rs/zw-fast-quantile)

A Rust implementation of the Zhang-Wang fast approximate quantile algorithm, as described in [*An Efficient Quantile Computation Technique for Approximate Query Processing* (SSDBM 2007)](http://web.cs.ucla.edu/~weiwang/paper/SSDBM07_2.pdf). This is the original algorithm that inspired the quantile techniques used in [TensorFlow Boosted Trees](https://android.googlesource.com/platform/external/tensorflow/+/6341b8975b7660244e1ca3003bfce5371f1fd167/tensorflow/core/kernels/boosted_trees/quantiles/weighted_quantiles_stream.h#30). For more context, see this [blog post](https://blog.paulme.ng/posts/2021-12-29-zw-fast-quantile%3A-approximate-quantile-for-high-speed-data-stream.html).

## Features

- Approximate quantile computation with tunable error bound (`epsilon`)
- Two variants: **fixed-size** (known stream length) and **unbounded** (unknown stream length)
- Generic over any `Clone + Ord` element type
- Sub-linear memory usage — space grows logarithmically with stream size
- Optional `serde` support for serialization/deserialization
- Query caching for efficient repeated lookups

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
zw-fast-quantile = "1.0"
```

## Usage

### FixedSizeEpsilonSummary

Use when the total number of elements is known ahead of time. This allows tighter memory allocation.
Inserting more than the declared stream size `n` panics.

```rust
use zw_fast_quantile::FixedSizeEpsilonSummary;

let epsilon = 0.01;  // 1% error tolerance
let n = 10000;       // expected stream size
let mut summary = FixedSizeEpsilonSummary::new(n, epsilon).unwrap();

for value in 1..=n {
    summary.update(value);
}

// Query the median (rank 0.5)
let median = summary.query(0.5).unwrap();

// Query any quantile from 0.0 (min) to 1.0 (max)
let p99 = summary.query(0.99).unwrap();
```

### UnboundEpsilonSummary

Use when the stream size is not known in advance. It dynamically manages internal summaries as the stream grows.

```rust
use zw_fast_quantile::UnboundEpsilonSummary;

let epsilon = 0.01;
let mut summary = UnboundEpsilonSummary::new(epsilon).unwrap();

// Feed values from any iterator or stream
for value in 1..=10000 {
    summary.update(value);
}

let median = summary.query(0.5).unwrap();
let count = summary.size();  // number of elements inserted
```

### Error Handling

Both constructors and `query()` return `Result<_, QuantileError>`:

```rust
use zw_fast_quantile::{FixedSizeEpsilonSummary, QuantileError};

// Invalid parameters
assert!(FixedSizeEpsilonSummary::<usize>::new(0, 0.1).is_err());     // n must be > 0
assert!(FixedSizeEpsilonSummary::<usize>::new(10, -0.1).is_err());   // epsilon must be positive
assert!(FixedSizeEpsilonSummary::<usize>::new(10, f64::NAN).is_err()); // epsilon must be finite

// Query errors
let s = FixedSizeEpsilonSummary::<usize>::new(10, 0.1).unwrap();
assert!(s.query(0.5).is_err());  // EmptySummary — no values inserted
```

### Choosing Epsilon

The `epsilon` parameter controls the trade-off between accuracy and memory:

| Epsilon | Max Error | Relative Memory |
|---------|-----------|-----------------|
| 0.1     | 10%       | Lowest          |
| 0.01    | 1%        | Moderate        |
| 0.001   | 0.1%      | Higher          |

A query for rank `r` on a stream of `n` elements will return an element whose true rank is within `r ± epsilon` (i.e., the result is within `epsilon * n` positions of the exact answer).

## Benchmarks

Benchmarked against GK01 from the [quantiles](https://github.com/postmates/quantiles) crate. Inserting 5000 values with `epsilon = 0.01`:

**Update (insertion) performance** — ZW is ~2.6x faster:

```
zw unbound quantile update
                        time:   [60.780 us 60.855 us 60.936 us]
```

```
gk quantile update      time:   [156.84 us 157.02 us 157.24 us]
```

**Query performance** (10 quantiles) — ZW is ~1.5x faster:

```
zw unbound quantile query
                        time:   [229.62 ns 230.16 ns 230.77 ns]
```

```
gk quantile query       time:   [350.21 ns 350.48 ns 350.76 ns]
```

Run benchmarks yourself with:

```bash
cargo bench
```

## Optional Features

- **`serde`** — Enable serialization/deserialization support via [serde](https://serde.rs/). Cache fields are skipped during serialization and lazily rebuilt on the next query.

  ```toml
  [dependencies]
  zw-fast-quantile = { version = "1.0", features = ["serde"] }
  ```

## Minimum Supported Rust Version

The MSRV is **1.65**.

## Documentation

- [API docs on docs.rs](https://docs.rs/zw-fast-quantile)
- [Algorithm internals](docs/algorithm.md) — how the data structures and operations work
- [Architecture overview](docs/architecture.md) — code organization and design decisions

## Related Projects

- [zw-fast-quantile-py](https://github.com/MnO2/zw-fast-quantile-py) — Python bindings for this library

## License

Licensed under [Apache-2.0](LICENSE).
