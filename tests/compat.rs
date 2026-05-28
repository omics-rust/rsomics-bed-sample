//! Correctness tests for rsomics-bed-sample.
//!
//! bedtools sample uses a different RNG so we cannot byte-match its output.
//! Instead we verify count, order preservation, and subset-of-input properties.

use rsomics_bed_sample::sample;
use std::io::Cursor;

fn run(input: &str, n: usize, seed: Option<u64>) -> Vec<String> {
    let mut out = Vec::new();
    sample(Cursor::new(input), &mut out, n, seed).unwrap();
    String::from_utf8(out)
        .unwrap()
        .lines()
        .map(|l| l.to_owned())
        .collect()
}

#[test]
fn output_is_subset_of_input() {
    let inp = "chr1\t1\t10\nchr1\t20\t30\nchr2\t1\t5\nchr3\t100\t200\nchr4\t0\t50\n";
    let records: Vec<&str> = inp.lines().collect();
    let got = run(inp, 3, Some(1));
    for line in &got {
        assert!(
            records.contains(&line.as_str()),
            "sampled record not in input: {line}"
        );
    }
}

#[test]
fn output_count_is_exact() {
    let inp = "chr1\t1\t10\nchr1\t20\t30\nchr2\t1\t5\nchr3\t100\t200\nchr4\t0\t50\n";
    for n in [1, 3, 5] {
        let got = run(inp, n, Some(42));
        assert_eq!(got.len(), n, "expected {n} records, got {}", got.len());
    }
}

#[test]
fn output_no_duplicates() {
    let inp = "chr1\t1\t10\nchr1\t20\t30\nchr2\t1\t5\nchr3\t100\t200\nchr4\t0\t50\n";
    let got = run(inp, 4, Some(7));
    let mut deduped = got.clone();
    deduped.sort_unstable();
    deduped.dedup();
    assert_eq!(deduped.len(), got.len(), "duplicate records in output");
}

#[test]
fn output_is_in_input_order() {
    let inp = "chr1\t1\t10\nchr1\t20\t30\nchr2\t1\t5\nchr3\t100\t200\nchr4\t0\t50\n";
    let records: Vec<&str> = inp.lines().collect();
    let got = run(inp, 3, Some(5));
    let indices: Vec<usize> = got
        .iter()
        .map(|l| records.iter().position(|&r| r == l.as_str()).unwrap())
        .collect();
    assert!(
        indices.windows(2).all(|w| w[0] < w[1]),
        "output not in input order"
    );
}
