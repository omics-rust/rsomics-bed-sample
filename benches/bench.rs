use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_bed_sample(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bed-sample");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bed = manifest
        .parent()
        .unwrap()
        .join("rsomics-bed-merge/tests/golden/sorted.bed");
    c.bench_function("rsomics-bed-sample golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args(["-i", bed.to_str().unwrap(), "-n", "5"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_bed_sample);
criterion_main!(benches);
