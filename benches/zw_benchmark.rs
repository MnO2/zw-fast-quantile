use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
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

    let mut s = UnboundEpsilonSummary::new(epsilon).unwrap();
    for i in 1..=n {
        s.update(i);
    }

    s
}

fn bench_fixed_zw_quantile_update() -> FixedSizeEpsilonSummary<usize> {
    let n = 5000;
    let epsilon: f64 = 0.01;

    let mut s = FixedSizeEpsilonSummary::new(n, epsilon).unwrap();
    for i in 1..=n {
        s.update(i);
    }

    s
}

fn bench_gk_quantile_query(s: &greenwald_khanna::Stream<usize>) {
    for rank in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] {
        black_box(s.quantile(rank));
    }
}

fn bench_unbound_zw_quantile_query(s: &UnboundEpsilonSummary<usize>) {
    for rank in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] {
        black_box(s.query(rank).unwrap());
    }
}

fn bench_fixed_zw_quantile_query(s: &FixedSizeEpsilonSummary<usize>) {
    for rank in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] {
        black_box(s.query(rank).unwrap());
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("zw fixed quantile update", |b| b.iter(bench_fixed_zw_quantile_update));
    c.bench_function("zw unbound quantile update", |b| {
        b.iter(bench_unbound_zw_quantile_update)
    });
    c.bench_function("gk quantile update", |b| b.iter(bench_gk_quantile_update));

    let fixed_s = bench_fixed_zw_quantile_update();
    fixed_s.query(0.5).unwrap();
    c.bench_function("zw fixed quantile query (warm cache)", |b| {
        b.iter(|| bench_fixed_zw_quantile_query(&fixed_s))
    });
    let fixed_cold = bench_fixed_zw_quantile_update();
    c.bench_function("zw fixed quantile query (cold cache)", |b| {
        b.iter_batched(
            || fixed_cold.clone(),
            |summary| bench_fixed_zw_quantile_query(&summary),
            BatchSize::SmallInput,
        )
    });

    let s = bench_unbound_zw_quantile_update();
    s.query(0.5).unwrap();
    c.bench_function("zw unbound quantile query (warm cache)", |b| {
        b.iter(|| bench_unbound_zw_quantile_query(&s))
    });
    let unbound_cold = bench_unbound_zw_quantile_update();
    c.bench_function("zw unbound quantile query (cold cache)", |b| {
        b.iter_batched(
            || unbound_cold.clone(),
            |summary| bench_unbound_zw_quantile_query(&summary),
            BatchSize::SmallInput,
        )
    });

    let stream = bench_gk_quantile_update();
    c.bench_function("gk quantile query", |b| b.iter(|| bench_gk_quantile_query(&stream)));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
