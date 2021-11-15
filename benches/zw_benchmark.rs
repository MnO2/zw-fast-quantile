#[macro_use]
extern crate criterion;
extern crate zw_fast_quantile;

use criterion::Criterion;
use zw_fast_quantile::FixedSizeEpsilonSummary;

fn bench_zw_quantile() {
    let epsilon = 0.1;
    let n = 10000;
    let mut s = FixedSizeEpsilonSummary::new(n, epsilon);
    for i in 1..=n {
        s.update(i);
    }

    for i in 1..=n {
        let _ = s.query(i);
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("zw quantile", |b| b.iter(bench_zw_quantile));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
