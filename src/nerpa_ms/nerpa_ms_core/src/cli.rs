use std::path::PathBuf;

use clap::Parser;

#[derive(Debug, Parser)]
#[command(name = "hmm_dag_top_paths_per_weight")]
pub struct Cli {
    /// Path to a JSON file with a list of (HMM, Parsed_rBAN_Record, Linearizatoin) triplets
    #[arg(long)]
    pub input: PathBuf,

    /// Maximum number of edits allowed in a structure
    #[arg(long)]
    pub max_edits: usize,

    /// The number of variants generated for each number of edits
    #[arg(long)]
    pub num_variants_per_num_edits: usize,

    /// Output JSON path
    #[arg(long)]
    pub out: PathBuf,
}
