use std::fs::File;
use std::io;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bed_sample::sample;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-sample", disable_help_flag = true)]
pub struct Cli {
    /// Input BED file (default: stdin).
    #[arg(short = 'i', long = "input")]
    pub input: Option<PathBuf>,

    /// Output file (default: stdout).
    #[arg(short = 'o', long = "out")]
    pub output: Option<PathBuf>,

    /// Number of records to sample.
    #[arg(short = 'n', long = "num", required = true)]
    pub num: usize,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let mut stdout_lock;
        let mut file_out;
        let out: &mut dyn io::Write = if let Some(ref p) = self.output {
            file_out = File::create(p).map_err(RsomicsError::Io)?;
            &mut file_out
        } else {
            stdout_lock = io::stdout().lock();
            &mut stdout_lock
        };
        let seed = self.common.seed;
        if let Some(ref p) = self.input {
            let f = File::open(p).map_err(RsomicsError::Io)?;
            sample(f, out, self.num, seed)
        } else {
            let stdin = io::stdin();
            sample(stdin.lock(), out, self.num, seed)
        }
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Sample random BED records using reservoir sampling (bedtools sample equivalent).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["[OPTIONS] -n <N>"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('i'),
                long: "input",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("stdin"),
                description: "Input BED file",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "out",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("stdout"),
                description: "Output file",
                why_default: None,
            },
            FlagSpec {
                short: Some('n'),
                long: "num",
                aliases: &[],
                value: Some("<int>"),
                type_hint: Some("usize"),
                required: true,
                default: None,
                description: "Number of records to sample",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Sample 100 records from a BED file",
            command: "rsomics-bed-sample -n 100 -i intervals.bed",
        },
        Example {
            description: "Reproducible sample with seed",
            command: "rsomics-bed-sample -n 100 --seed 42 -i intervals.bed",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
