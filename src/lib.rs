//! Sample random BED records — `bedtools sample` equivalent.
//!
//! Implements reservoir sampling (Vitter's Algorithm R) to draw an exact
//! sample of `n` records from a BED stream without reading the entire file
//! into memory first.  The implementation stores exactly `n` records in
//! memory at any time; only the reservoir occupies RAM — not the whole file.
//!
//! ## Determinism
//!
//! Provide `--seed` to get reproducible output.  Without a seed the output
//! is non-deterministic (but the count is always exactly `n`, or the full
//! file if there are fewer than `n` records).
//!
//! ## Output order
//!
//! The sampled records are emitted in their original input order (by line
//! number), not in sampling order.  This matches bedtools sample behaviour.

use std::io::{BufRead, BufReader, Read, Write};

use rsomics_common::{Result, RsomicsError};

/// Sample up to `n` records from a BED stream `r`, writing to `w`.
///
/// Uses reservoir sampling (Vitter's Algorithm R).  Records are emitted in
/// their original input order.  Header lines (`#`, `track`, `browser`) and
/// blank lines are skipped before counting.
///
/// `seed`: `None` → non-deterministic; `Some(s)` → deterministic LCG-seeded.
pub fn sample<R: Read, W: Write>(r: R, w: W, n: usize, seed: Option<u64>) -> Result<()> {
    if n == 0 {
        return Ok(());
    }

    let mut rdr = BufReader::new(r);
    let mut line = String::new();
    let mut rng = LcgRng::new(seed.unwrap_or_else(|| {
        // Non-deterministic seed from the system.
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(12345, |d| d.subsec_nanos().into())
    }));

    // Reservoir: (original_index, line_bytes).
    let mut reservoir: Vec<(usize, String)> = Vec::with_capacity(n);
    let mut count: usize = 0;

    loop {
        line.clear();
        let bytes = rdr.read_line(&mut line).map_err(RsomicsError::Io)?;
        if bytes == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(['\n', '\r']);
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track")
            || trimmed.starts_with("browser")
        {
            continue;
        }

        if reservoir.len() < n {
            reservoir.push((count, line.trim_end_matches(['\n', '\r']).to_owned()));
        } else {
            // Replace random element with probability n / (count + 1).
            let j = rng.next_usize(count + 1);
            if j < n {
                reservoir[j] = (count, line.trim_end_matches(['\n', '\r']).to_owned());
            }
        }
        count += 1;
    }

    // Sort by original input order before output.
    reservoir.sort_unstable_by_key(|(idx, _)| *idx);

    let mut bw = std::io::BufWriter::new(w);
    for (_, rec) in &reservoir {
        bw.write_all(rec.as_bytes()).map_err(RsomicsError::Io)?;
        bw.write_all(b"\n").map_err(RsomicsError::Io)?;
    }
    bw.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

/// Minimal LCG (linear congruential generator) for reproducible reservoir sampling.
/// Constants from Knuth's MMIX / Numerical Recipes.
struct LcgRng {
    state: u64,
}

impl LcgRng {
    fn new(seed: u64) -> Self {
        Self {
            state: seed.wrapping_add(1),
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    /// Returns a uniform random value in `0..max`.
    ///
    /// `max` must be ≤ `u64::MAX`; on platforms where `usize` is 32-bit the
    /// result is silently truncated, but `max` (a reservoir index or record
    /// count) is always < `usize::MAX` in practice.
    #[allow(clippy::cast_possible_truncation)]
    fn next_usize(&mut self, max: usize) -> usize {
        (self.next_u64() % (max as u64)) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn run(input: &str, n: usize, seed: Option<u64>) -> Vec<String> {
        let mut out = Vec::new();
        sample(Cursor::new(input), &mut out, n, seed).unwrap();
        String::from_utf8(out)
            .unwrap()
            .lines()
            .map(str::to_owned)
            .collect()
    }

    #[test]
    fn sample_all_when_n_exceeds_count() {
        let inp = "chr1\t1\t10\nchr2\t5\t15\nchr3\t20\t30\n";
        let got = run(inp, 10, Some(0));
        assert_eq!(got.len(), 3);
    }

    #[test]
    fn sample_exact_count() {
        let inp = "chr1\t1\t10\nchr1\t20\t30\nchr2\t1\t5\nchr3\t100\t200\nchr4\t0\t50\n";
        let got = run(inp, 3, Some(42));
        assert_eq!(got.len(), 3);
    }

    #[test]
    fn sample_preserves_input_order() {
        // Sampled records must appear in the order they appeared in the input.
        let inp = "chr1\t1\t10\nchr1\t20\t30\nchr2\t1\t5\nchr3\t100\t200\nchr4\t0\t50\n";
        let records: Vec<&str> = inp.lines().collect();
        let got = run(inp, 3, Some(99));
        for i in 1..got.len() {
            let a = records.iter().position(|&r| r == got[i - 1]).unwrap();
            let b = records.iter().position(|&r| r == got[i]).unwrap();
            assert!(a < b, "output not in input order: {a} >= {b}");
        }
    }

    #[test]
    fn zero_n_produces_empty_output() {
        let inp = "chr1\t1\t10\n";
        let got = run(inp, 0, None);
        assert!(got.is_empty());
    }

    #[test]
    fn skips_headers_and_blanks() {
        let inp = "# header\nchr1\t1\t10\n\nchr2\t5\t15\n";
        let got = run(inp, 10, Some(1));
        assert_eq!(got.len(), 2);
        assert!(!got.iter().any(|l| l.starts_with('#')));
    }

    #[test]
    fn deterministic_with_seed() {
        let inp = "chr1\t1\t10\nchr1\t20\t30\nchr2\t1\t5\nchr3\t100\t200\nchr4\t0\t50\n";
        let a = run(inp, 3, Some(777));
        let b = run(inp, 3, Some(777));
        assert_eq!(a, b, "same seed must produce same output");
    }

    #[test]
    fn different_seeds_may_differ() {
        let inp =
            "chr1\t1\t10\nchr1\t20\t30\nchr2\t1\t5\nchr3\t100\t200\nchr4\t0\t50\nchrX\t0\t100\n";
        let a = run(inp, 3, Some(1));
        let b = run(inp, 3, Some(999));
        // This could theoretically be equal by chance, but is astronomically unlikely
        // for 2 different seeds on 6 items choosing 3.
        assert_ne!(
            a, b,
            "different seeds should (almost certainly) produce different output"
        );
    }
}
