#[macro_use]
extern crate criterion;
extern crate quantiles;
extern crate zw_fast_quantile;

use criterion::Criterion;
use quantiles::greenwald_khanna;
use zw_fast_quantile::{FixedSizeEpsilonSummary, UnboundEpsilonSummary};

fn bench_gk_quantile_update() -> greenwald_khanna::Stream<usize> {
    let n = 5000;
    let epsilon: f64 = 0.01;

    let mut stream = greenwald_khanna::Stream::new(epsilon);
    for i in 1..n {
        stream.insert(i);
    }

    stream
}

fn bench_unbound_zw_quantile_update() -> UnboundEpsilonSummary<usize> {
    let n = 5000;
    let epsilon: f64 = 0.01;

    let mut s = UnboundEpsilonSummary::new(epsilon);
    for i in 1..=n {
        s.update(i);
    }

    s
}

fn bench_fixed_zw_quantile_update() -> FixedSizeEpsilonSummary<usize> {
    let n = 5000;
    let epsilon: f64 = 0.01;

    let mut s = FixedSizeEpsilonSummary::new(n, epsilon);
    for i in 1..=n {
        s.update(i);
    }

    s
}

fn bench_gk_quantile_query(s: &mut greenwald_khanna::Stream<usize>) {
    for rank in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] {
        s.quantile(rank);
    }
}

fn bench_unbound_zw_quantile_query(s: &mut UnboundEpsilonSummary<usize>) {
    for rank in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] {
        s.query(rank);
    }
}

fn bench_fixed_zw_quantile_query(s: &mut FixedSizeEpsilonSummary<usize>) {
    for rank in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] {
        s.query(rank);
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("zw fixed quantile update", |b| b.iter(bench_fixed_zw_quantile_update));
    c.bench_function("zw unbound quantile update", |b| {
        b.iter(bench_unbound_zw_quantile_update)
    });
    c.bench_function("gk quantile update", |b| b.iter(bench_gk_quantile_update));

    let mut fixed_s = bench_fixed_zw_quantile_update();
    c.bench_function("zw fixed quantile query", |b| {
        b.iter(|| bench_fixed_zw_quantile_query(&mut fixed_s))
    });

    let mut s = bench_unbound_zw_quantile_update();
    c.bench_function("zw unbound quantile query", |b| {
        b.iter(|| bench_unbound_zw_quantile_query(&mut s))
    });

    let mut stream = bench_gk_quantile_update();
    c.bench_function("gk quantile query", |b| b.iter(|| bench_gk_quantile_query(&mut stream)));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
