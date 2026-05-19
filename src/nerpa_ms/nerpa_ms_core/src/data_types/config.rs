use std::path::PathBuf;

pub struct Nerpa_MS_Config {
    /// Path to Nerpa root dir
    pub nerpa_root: PathBuf,

    /// Maximum number of Nerpa matches to check
    pub max_nerpa_matches: Option<usize>,

    /// Consider only matches with this score or greater
    pub min_nerpa_score: Option<f64>,

    /// Maximum number of edits allowed in a structure
    pub max_edits: usize,

    /// The number of variants generated for each number of edits
    pub num_variants_per_num_edits: usize,
}

impl Nerpa_MS_Config {
    pub fn from_cli(cli: &Cli) -> Self {
	Self {
	    nerpa_root: cli.nerpa_root.clone(),
	    max_nerpa_matches: cli.max_nerpa_matches,
	    min_nerpa_score: cli.min_nerpa_score,
	    max_edits: cli.max_edits,
	    num_variants_per_num_edits: cli.num_variants_per_num_edits,
	}
    }
}
