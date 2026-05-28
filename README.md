# rsomics-bed-sample

Sample random BED records using reservoir sampling — `bedtools sample` equivalent.

## Usage

```
rsomics-bed-sample [OPTIONS] -n <N>

Options:
  -i, --input <FILE>   BED input file (default: stdin)
  -o, --out <FILE>     Output file (default: stdout)
  -n, --num <N>        Number of records to sample (required)
      --seed <INT>     Random seed for reproducible output
  -h, --help           Print help
  -V, --version        Print version
```

## Description

Draws exactly `n` records from the BED stream using Vitter's Algorithm R
(reservoir sampling).  The sample is uniform and unbiased.  Output records
are emitted in their original input order (not sampling order).

Header lines (`#`, `track`, `browser`) and blank lines are skipped before
counting.  If the file has fewer than `n` data records, all records are output.

## Example

```
$ rsomics-bed-sample -n 3 -i intervals.bed
chr1	1	10
chr2	5	15
chr4	0	50
```

## Install

```sh
cargo install rsomics-bed-sample
```

## Origin

Rust reimplementation of `bedtools sample`.  Algorithm: Vitter, J.S. (1985).
"Random Sampling with a Reservoir." *ACM Transactions on Mathematical Software*,
11(1), 37–57.  Informed by the
[bedtools2 source](https://github.com/arq5x/bedtools2) (MIT License).

License: MIT OR Apache-2.0.
Upstream credit: bedtools2 <https://github.com/arq5x/bedtools2> (MIT License).
