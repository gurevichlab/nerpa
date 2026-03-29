use clap::Parser;
use std::path::PathBuf;

/// NRP Matching Tool
#[derive(Parser, Debug)]
#[command(name = "nrp-matching", author, version, about)]
pub struct Args {
    /// Path to the HMM/BGC JSON file
    #[arg(short = 'H', long = "hmms_json", value_name = "FILE")]
    pub hmms_json: PathBuf,

    /// Path to the NRP linearizations JSON file
    #[arg(short = 'N', long = "nrps_json", value_name = "FILE")]
    pub nrps_json: PathBuf,

    /// Maximum number of matches to keep per BGC (0 = keep all)
    #[arg(long, value_name = "N")]
    pub max_num_matches_per_bgc: usize,

    /// Maximum number of matches to keep per NRP (0 = keep all)
    #[arg(long, value_name = "N")]
    pub max_num_matches_per_nrp: usize,

    /// Minimum number of matches to keep per BGC (overwrites max)
    #[arg(long, value_name = "N")]
    pub min_num_matches_per_bgc: usize,

    /// Minimum number of matches to keep per NRP (overwrites max)
    #[arg(long, value_name = "N")]
    pub min_num_matches_per_nrp: usize,

    /// Maximum number of matches to keep in total (0 = keep all)
    #[arg(long, value_name = "N")]
    pub max_num_matches: usize,

    /// Number of threads to use (e.g., 4)
    #[arg(short = 't', long = "threads", value_name = "N")]
    pub threads: usize,

    /// Output JSON file for the resulting matches
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    pub output: PathBuf,
}
