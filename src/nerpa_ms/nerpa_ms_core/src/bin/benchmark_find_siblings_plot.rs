use nerpa_ms_core::bench::plot_find_siblings::plot;
use nerpa_ms_core::bench::find_target::TargetSearchResults;
use std::collections::HashMap;
use anyhow::Result;
use std::path::{Path, PathBuf};
use clap::Parser;

#[derive(Debug, Parser)]
#[command(name = "benchmark_find_sibglings_plot")]
pub struct Cli {
    /// Path to a JSON with a list of TargetSearchResult objects
    #[arg(long)]
    pub target_search_results: PathBuf,

    /// Path to the output SVG with the benchmarking results
    #[arg(long)]
    pub out: PathBuf,

    /// Path to nerpa root dir
    #[arg(long)]
    pub nerpa_root: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let target_search_results_json = std::fs::read_to_string(&cli.target_search_results)
	.expect(&format!("Failed to read target search results JSON: {}",
			 cli.target_search_results.display()));
    let target_search_results: Vec<TargetSearchResults> = serde_json::from_str(&target_search_results_json)
	.expect(&format!("Failed to parse target search results JSON: {}",
			 cli.target_search_results.display()));
    let nerpa_ms_root = (
	&cli.nerpa_root
	    .join("src")
	    .join("nerpa_ms")
	    .join("nerpa_ms_core")
    );
	    
    plot(
	&target_search_results,
	nerpa_ms_root,
	&cli.out
    );
    Ok(())
}
