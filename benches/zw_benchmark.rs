#[macro_use]
extern crate criterion;
extern crate quantiles;
extern crate zw_fast_quantile;

use criterion::Criterion;
use quantiles::greenwald_khanna;
use zw_fast_quantile::UnboundEpsilonSummary;

fn bench_gk_quantile() {
    let n = 1000;
    let epsilon: f64 = 0.01;

    let mut stream = greenwald_khanna::Stream::new(epsilon);
    for i in 1..n {
        stream.insert(i);
    }
}

fn bench_zw_quantile() {
    let n = 1000;
    let epsilon: f64 = 0.01;

    let mut s = UnboundEpsilonSummary::new(epsilon);
    for i in 1..=n {
        s.update(i);
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("zw quantile", |b| b.iter(bench_zw_quantile));
    c.bench_function("gk quantile", |b| b.iter(bench_gk_quantile));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
