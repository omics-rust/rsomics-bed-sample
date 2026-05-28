use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use rsomics_bed_sample::sample;
use std::io::sink;

fn make_fixture(n: usize) -> Vec<u8> {
    use std::fmt::Write as _;
    let mut s = String::with_capacity(n * 28);
    let chroms = ["chr1", "chr2", "chr3", "chr4", "chr5"];
    for i in 0..n {
        let chrom = chroms[i % chroms.len()];
        let start = (i * 200) as u64;
        let end = start + 100;
        let _ = writeln!(s, "{chrom}\t{start}\t{end}");
    }
    s.into_bytes()
}

fn bench_sample(c: &mut Criterion) {
    let fixture = make_fixture(100_000);
    let mut group = c.benchmark_group("sample");
    group.throughput(Throughput::Bytes(fixture.len() as u64));
    group.bench_function("100k_sample_1000", |b| {
        b.iter(|| {
            sample(fixture.as_slice(), sink(), 1_000, Some(42)).unwrap();
        });
    });
    group.finish();
}

criterion_group!(benches, bench_sample);
criterion_main!(benches);
