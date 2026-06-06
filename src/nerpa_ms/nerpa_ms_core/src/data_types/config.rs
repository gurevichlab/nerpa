use std::path::PathBuf;
use crate::cli::Cli;

pub struct NerpaMS_Config {
    /// Maximum number of Nerpa matches to check
    pub max_nerpa_matches: Option<usize>,

    /// Consider only matches with this score or greater
    pub min_nerpa_score: Option<f64>,

    /// Maximum number of edits allowed in a structure
    pub max_edits: usize,

    /// The number of variants generated for each number of edits
    pub num_variants_per_num_edits: usize,
}

impl NerpaMS_Config {
    pub fn from_cli(cli: &Cli) -> Self {
	Self {
	    max_nerpa_matches: cli.max_nerpa_matches,
	    min_nerpa_score: cli.min_nerpa_score,
	    max_edits: cli.max_edits,
	    num_variants_per_num_edits: cli.num_variants_per_num_edits,
	}
    }
}

pub struct DebugConfig {
    /// Write debug log to stdout
    pub debug_stdout: bool,
    
    /// Whether to write alignments to disk
    pub write_alignments: bool,

    /// Whether to draw HMMs with optimal paths
    pub draw_hmm_opt_paths: bool,

    /// Whether to draw DAGs with optimal paths
    pub draw_hmm_dag_opt_paths: bool,

    /// Whether to draw output variants
    pub draw_output_variants: bool,

    /// Directory to write debug outputs to
    pub out: PathBuf,

    /// Path to Nerpa root dir
    pub nerpa_root: PathBuf,
}

impl DebugConfig {
    pub fn from_cli(cli: &Cli) -> Self {
	Self {
	    debug_stdout: cli.debug_stdout,
	    write_alignments: cli.write_alignments,
	    draw_hmm_opt_paths: cli.draw_hmm_opt_paths,
	    draw_hmm_dag_opt_paths: cli.draw_hmm_dag_opt_paths,
	    draw_output_variants: cli.draw_output_variants,
	    out: cli.out.join("debug"),
	    nerpa_root: cli.nerpa_root.clone(),
	}
    }
}

