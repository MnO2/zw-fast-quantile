# Audit Execution Plan (2026-07-12)

Source: project audit run on 2026-07-12 against `main` @ `53ba02a` (clean worktree).
Audience: implementing agents. Each work package (WP) is self-contained, lists its
files, exact changes, and acceptance criteria. Execute in order unless noted —
WP1–WP3 are independent of each other; WP4+ depend on WP1–WP3 being merged.

## Context snapshot (verified at audit time)

- `cargo test --all-features` → 24 pass (22 unit + 2 doc).
- `cargo clippy --all-features -- -D warnings` (CI invocation) → pass;
  with `--all-targets` → **12 errors**, all in `src/lib.rs` test code.
- `cargo fmt --all -- --check` → pass.
- crates.io latest is **0.2.3** (2021-12-30). Repo `Cargo.toml` says 1.0.0.
  No git tag for 1.0.0 exists. The 1.0 modernization is entirely unshipped —
  breaking API changes are still cheap.
- Accuracy probe (external binary, rank-error vs the ε·N contract on random
  `u64` streams):
  - FixedSize at design capacity (n=10k, ε=0.01): max rank error 0.53 × ε·N — good.
  - **FixedSize overfilled 10× (insert 100k into n=10k): max rank error 85 × ε·N,
    no panic, no error, `size()` reports 100000.** Root cause: `update`'s level
    loop (`src/lib.rs:190-199`) drops the compressed block when all levels are
    occupied.
  - Unbound ε=0.1 at 2M elements (past the `x + x >= 32` cutoff): 0.15 × ε·N — good.
  - Unbound ε=0.01 at 2M: 0.12 × ε·N — good.

Global verification command for every WP (run from repo root):

```bash
cargo fmt --all -- --check \
  && cargo clippy --all-features --all-targets -- -D warnings \
  && cargo test --all-features
```

(Until WP2 lands, the clippy step fails for pre-existing reasons; WP1 may run
tests only.)

---

## WP1 — Overfill contract on `FixedSizeEpsilonSummary::update` (CRIT)

**Problem.** `src/lib.rs:190-199`: when `s[0]` fills and every level `1..level`
is already occupied, the merged block `s_c` falls out of scope and the data is
silently discarded. Nothing guards against inserting more than the declared `n`.
Measured effect: 85× the error guarantee while `size()` keeps counting.

**Decision (default — implement this unless the owner overrides):** panic on
overfill with a clear message. Keep `update(&mut self, e: T)` — a fixed-capacity
structure receiving element `n+1` is a programmer error, analogous to
out-of-bounds indexing. The alternative (make `update` return
`Result<(), QuantileError>` with a `CapacityExceeded` variant) touches every
call site, doctest, and bench; only choose it if the owner asks.

**Changes.**
1. In `FixedSizeEpsilonSummary`, store the declared capacity (`n`) at
   construction (new private field, e.g. `capacity: usize`).
2. At the top of `update`, before mutating anything:
   `assert!(self.cnt < self.capacity, "FixedSizeEpsilonSummary capacity exceeded: constructed for n={} elements", self.capacity);`
3. `UnboundEpsilonSummary::update` (`src/lib.rs:484-504`) feeds its internal
   `FixedSizeEpsilonSummary` values beyond the "increment and continue" cutoff
   (`x + x >= 32`, reachable past ~2³⁰/ε elements). The inner summary created at
   x=15 is sized for `boundaries[30]`, so the assert must NOT fire before that.
   Verify the capacity bookkeeping keeps the existing unbound tests green and
   add a comment at the cutoff branch explaining the headroom.
4. Docs: state the contract in the `FixedSizeEpsilonSummary` rustdoc and in
   README "FixedSizeEpsilonSummary" section: inserting more than `n` elements
   panics.
5. Keep the internal level-loop as-is (the assert makes the dead branch
   unreachable), but add a `debug_assert!(false, ...)` or comment at the loop
   exit documenting why it cannot be reached.

**Tests to add** (in `src/lib.rs` tests module):
- `#[should_panic]` test: `new(10, 0.1)`, insert 11 elements.
- Rank-error contract test (this is the test class that would have caught the
  bug): insert exactly `n` random `u64` values (seeded `StdRng`), sort a copy,
  and for a grid of ranks assert the returned value's true rank is within
  ε·n of the target. Do this for FixedSize (n=10_000, ε=0.01) and Unbound
  (n=100_000, ε=0.01). Compute true rank via
  `sorted.partition_point(|x| *x < v)` / `partition_point(|x| *x <= v)`
  and assert the target rank falls within `[lo - εn, hi + εn]`.

**Acceptance.** New tests pass; all existing tests pass; README and rustdoc
mention the panic contract.

---

## WP2 — Fix clippy debt and close the CI gate (STRUCT)

**Problem.** CI (`.github/workflows/CI.yml:41`) runs clippy without
`--all-targets`, so 12 existing errors in test code are invisible:
4× `needless_range_loop` (the `for i in 0..n { s.update(records[i]) }` loops in
the four distribution tests), 8× needless borrow.

**Changes.**
1. Run `cargo clippy --all-features --all-targets --fix --allow-dirty` then
   hand-review; or fix manually (iterate `for &x in records.iter().take(n)`).
2. `.github/workflows/CI.yml` clippy job: change to
   `cargo clippy --all-features --all-targets -- -D warnings`.

**Acceptance.** Global verification command passes end-to-end.

---

## WP3 — Test hygiene (INCR)

All in `src/lib.rs` tests module. Coordinates are pre-WP1 line numbers; re-locate
by name after WP1/WP2 land.

1. Deduplicate the four 1M-element distribution tests
   (`test_{normal,pareto}_distribution_generated_seq_on_{fixed,unbound}_summary`,
   lines 688-861): extract one helper parameterized by distribution, ε, and a
   generic `update`/`query` closure or a small trait over the two summary types.
   Keep the four `#[test]` entry points (names are part of CI history).
2. Rename `test_query_unbound_summary_with_insufficient_values`
   (line 642) — it exercises `FixedSizeEpsilonSummary`, not Unbound. New name:
   `test_query_fixed_summary_with_insufficient_values`.
3. Optional: drop the `stats_alloc` global-allocator instrumentation and the
   `println!("{}k", ...)` memory printouts from the distribution tests (they are
   benchmark residue, slow every test run, and assert nothing). If memory
   tracking is worth keeping, move it to a criterion bench or an
   `#[ignore]`d test. Removing it also allows dropping the `stats_alloc`
   dev-dependency.
4. Add `#[derive(Clone)]` to `UnboundEpsilonSummary` (line 442-443) for parity
   with `FixedSizeEpsilonSummary`, plus a one-line test that clones and queries.

**Acceptance.** Global verification passes; test count does not decrease
(dedup keeps entry points; WP1 added 3+).

---

## WP4 — Hot-path performance (INCR) — after WP1 (needs the rank-error test as a safety net)

All changes in `src/lib.rs`. Justify each with a before/after `cargo bench`
run on `zw fixed quantile update` / `zw unbound quantile update`; keep only wins.

1. Line 180 (`update`), line 213 (query cache fill), line 238 (`calc_s_m`):
   `sort()` → `sort_unstable()`. Ties are between equal values; stability is
   irrelevant here.
2. Lines 187-189: replace
   `compress(self.s[0].clone(), ...); self.s[0].clear();` with
   `let block = std::mem::replace(&mut self.s[0], Vec::with_capacity(self.b));
   let mut s_c = compress(block, compressed_size, self.epsilon);`
   — avoids copying the full block every rotation.
3. Line 195: `merge(s_c, &self.s[k])` clones every element of `s[k]`, which is
   cleared on the next line. Change `merge` to take the second argument by
   value (`std::mem::take(&mut self.s[k])`) or add an owning variant; removes a
   full clone per occupied level per rotation (matters most when `T` is not
   `Copy`). Callers: `update` (195), `query` (220, 518), `calc_s_m` (246). The
   query/calc_s_m call sites must NOT consume `self.s[i]` — keep a borrowing
   path for them (e.g. `merge_owned` used only by `update`).

**Acceptance.** Global verification passes; rank-error tests from WP1 pass;
bench summary (before/after table) included in the PR description.

---

## WP5 — Benchmark quality and CI cost (INCR)

1. `benches/zw_benchmark.rs`: wrap query results in `criterion::black_box`
   (`black_box(s.query(rank).unwrap())`, same for GK) so the compiler cannot
   elide the work.
2. Add a cache-cold ZW query variant: call `update` once between queries or
   construct per-iteration, so the `RefCell` cache is not the only thing
   measured. Label the existing one `zw ... query (warm cache)`.
3. Optional: extend to a small grid (n ∈ {10⁴, 10⁵}, ε ∈ {0.1, 0.01}) using
   `BenchmarkGroup` — keeps runtime sane while giving scaling data.
4. README "Benchmarks" section: regenerate all numbers with `cargo bench` on
   current code (existing numbers predate the modernization), and add one
   sentence noting the ZW query figure benefits from result caching.
5. `.github/workflows/bench.yml`: benchmarks currently run on every push and
   PR and the output is discarded. Either
   (a) wire the output into `benchmark-action/github-action-benchmark` with
   `output-format bencher`, storing history on `gh-pages` and alerting on
   regression, or (b) restrict the trigger to `workflow_dispatch` +
   `push: branches: [main]`. Option (b) is the low-effort default.
6. `.github/workflows/CI.yml` codecov job: replace
   `cargo install cargo-tarpaulin` (compiles from source every run) with
   `uses: taiki-e/install-action@cargo-tarpaulin`.

**Acceptance.** `cargo bench` runs clean locally; CI workflows lint (actionlint
if available) and pass on a test push; README numbers match a run recorded in
the PR.

---

## WP6 — Investigate `boundaries[x + x]` sizing (STRUCT, investigation)

`src/lib.rs:496`: the sketch rotated in at boundary `x` only absorbs
`boundaries[x+1] - boundaries[x]` elements before the next rotation, but is
constructed with `n = boundaries[x + x] - cnt - 1` — exponentially oversized
(≈2^(x-1)× at boundary x). Accuracy is unaffected (probe: 0.12–0.15 × ε·N at
2M elements); the cost is memory (larger block size and level count in the
active sketch).

Task: check the Zhang-Wang paper (link in README) for the intended sizing.
- If `x + 1` is correct: change it, re-run the WP1 rank-error tests plus a
  2M-element ε=0.1 probe, and note the memory reduction (use the
  `stats_alloc` harness from WP3.3 if retained).
- If `x + x` is intentional: add a comment at line 496 citing the paper section
  so the next reader doesn't re-open this.
Do not change behavior without the rank-error tests from WP1 in place.

---

## WP7 — Ship 1.0.0 (release follow-through; needs owner approval to publish)

Blocked on WP1–WP2 minimum (WP1 is API-affecting; do not publish before it).

1. Add `CHANGELOG.md`: summarize 0.2.3 → 1.0.0 (Result-returning constructors
   and query, `QuantileError`, `query(&self)` via cache, serde feature, MSRV
   1.65, overfill panic contract from WP1, plus fixes listed in
   `docs/plans/2026-04-04-modernization-design-final.md` §1).
2. Prune `docs/plans/`: keep `2026-04-04-modernization-design-final.md` and
   this file; delete the two superseded drafts and two review snapshots (they
   are recorded in git history).
3. Preflight: clean worktree, `cargo package --list` review,
   `cargo publish --dry-run --all-features`.
4. **Stop and get explicit owner approval**, then: `git tag v1.0.0`, push tag,
   `cargo publish`. Verify: crates.io shows 1.0.0, docs.rs builds.
5. Also update the sibling project pinned in README
   (`zw-fast-quantile-py`) — out of scope here, but note it in the release PR.

**Acceptance.** Tag pushed, crate published, README badges resolve to 1.0.0.
Every step before `git tag`/`cargo publish` may run autonomously; the
tag/publish step requires the owner's explicit go-ahead in-thread.

---

## Ordering summary

```
WP1 (CRIT, API)  ─┐
WP2 (clippy/CI)  ─┼─ independent, do first, any order
WP3 (tests)      ─┘
WP4 (perf)       — after WP1 (safety net) ; WP2 (clean clippy base)
WP5 (bench/CI)   — anytime after WP2; README numbers after WP4
WP6 (x+x)        — after WP1
WP7 (release)    — last, gated on owner approval
```
