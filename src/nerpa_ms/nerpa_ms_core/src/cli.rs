use std::path::PathBuf;

use clap::Parser;

#[derive(Debug, Parser)]
#[command(name = "hmm_dag_top_paths_per_weight")]
pub struct Cli {
    /// Path to Nerpa results
    #[arg(
        long,
        conflicts_with = "input_json",
        required_unless_present = "input_json"
    )]
    pub nerpa_results: Option<PathBuf>,

    /// Alternative input source: a path to a JSON file with a list of (HMM, Parsed_rBAN_Record, Linearizatoin) triplets
    #[arg(
        long,
        conflicts_with = "nerpa_results",
        required_unless_present = "nerpa_results"
    )]
    pub input_json: Option<PathBuf>,

    /// Maximum number of Nerpa matches to check
    #[arg(long,
    	  requires = "nerpa_results")]
    pub max_nerpa_matches: Option<usize>,

    /// Consider only matches with this score or greater
    #[arg(long,
	  requires = "nerpa_results")]
    pub min_nerpa_score: Option<f64>,

    /// Maximum number of edits allowed in a structure
    #[arg(long)]
    pub max_edits: usize,

    /// The number of variants generated for each number of edits
    #[arg(long)]
    pub num_variants_per_num_edits: usize,

    /// Output JSON path
    #[arg(long)]
    pub out: PathBuf,

    /// Path to a JSON file with the monomers database
    #[arg(long)]
    pub monomers_db_json: PathBuf,

    /// Path to root Nerpa dir (used to draw molecules via calling python scripts from that repo)
    #[arg(long)]
    pub nerpa_root: Option<PathBuf>,
}
