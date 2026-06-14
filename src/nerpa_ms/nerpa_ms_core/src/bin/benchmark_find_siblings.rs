#![allow(non_snake_case, non_camel_case_types, non_upper_case_globals)]

use nerpa_ms_core::bench::find_target::TargetSearchResults;
use nerpa_ms_core::bench::bench_find_siblings::find_siblings;
use nerpa_ms_core::data_types::monomers_db::load_monomers_db;
use nerpa_ms_core::data_types::config::{DebugConfig, NerpaMS_Config};
use nerpa_ms_core::data_types::monomers_db::MonomersDB;
use nerpa_ms_core::io::input::{InputItem, get_input_from_nerpa_results};
use nerpa_ms_core::bench::plot_find_siblings::plot;

use std::collections::HashMap;
use anyhow::Result;
use std::path::{Path, PathBuf};
use clap::Parser;

#[derive(Debug, Parser)]
#[command(name = "benchmark_find_sibglings")]
pub struct Cli {
    /// Path to Nerpa results on MIBiG: each BGC against its respective NRPs in its own folder
    #[arg(long)]
    pub nerpa_results: PathBuf,

    /// Maximum number of edits allowed in a structure
    #[arg(long)]
    pub max_edits: usize,

    /// The number of variants generated for each number of edits
    #[arg(long)]
    pub num_variants_per_num_edits: usize,

    /// Path to the output tsv with the benchmarking results
    #[arg(long)]
    pub out: PathBuf,

    /// Path to nerpa root dir
    #[arg(long)]
    pub nerpa_root: PathBuf,

    /// Path to a JSON file with the monomers database
    #[arg(long)]
    pub monomers_db_json: PathBuf,
}


fn get_compound_distances(path: &PathBuf) -> Result<HashMap<(String, String), u32>> {
    let text = std::fs::read_to_string(path)?;

    let mut lines = text.lines().map(str::trim).filter(|l| !l.is_empty());

    let header_line = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("Empty compound distances file"))?;
    let header_ids: Vec<String> = header_line
        .split_whitespace()
        .map(|s| s.to_string())
        .collect();

    let mut map: HashMap<(String, String), u32> = HashMap::new();

    for line in lines {
        let mut it = line.split_whitespace();
        let row_id = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("Malformed line in compound distances: '{line}'"))?
            .to_string();

        let values: Vec<&str> = it.collect();
        if values.len() != header_ids.len() {
            return Err(anyhow::anyhow!(
                "Wrong number of distances for '{row_id}': expected {}, got {}",
                header_ids.len(),
                values.len()
            ));
        }

        for (col_id, dist_str) in header_ids.iter().zip(values.into_iter()) {
            let d: i32 = dist_str.parse().map_err(|e| {
                anyhow::anyhow!(
                    "Failed to parse distance '{dist_str}' for ({row_id}, {col_id}): {e}"
                )
            })?;
	    if d >= 0 {  // negative values indicate "infinite" distance
		map.insert((row_id.clone(), col_id.clone()), d as u32);
		map.insert((col_id.clone(), row_id.clone()), d as u32);
            }
	}
    }

    Ok(map)
}

fn process_single_output_dir(
    dir: &PathBuf,
    monomers_db: &MonomersDB,
    nerpa_ms_cfg: &NerpaMS_Config,
    debug_cfg: &DebugConfig,
) -> Vec<TargetSearchResults> {
    let input_items: Vec<InputItem> = {
	get_input_from_nerpa_results(
	    dir, None, None
	)
	    .expect(&format!("Failed to load Nerpa results from {}",
			    dir.display()))
    };

    find_siblings(
	&input_items,
	monomers_db,
	nerpa_ms_cfg,
	debug_cfg,
    )
}

fn write_results(results: &[TargetSearchResults], out: &Path) -> Result<()> {
    std::fs::create_dir_all(out.parent().unwrap())?; // Create parent directories
    let json = serde_json::to_string_pretty(results)?;
    std::fs::write(out, json)?;
    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let nerpa_ms_cfg = NerpaMS_Config {
	max_nerpa_matches: None,
	min_nerpa_score: None,
	max_edits: cli.max_edits,
	num_variants_per_num_edits: cli.num_variants_per_num_edits,
    };

    let debug_cfg = DebugConfig {
	write_alignments: false,
	debug_stdout: false,
	draw_hmm_opt_paths: false,
	draw_hmm_dag_opt_paths: false,
	draw_output_variants: false,
	out: cli.out.join("debug"),
	nerpa_root: cli.nerpa_root.clone(),
    };

    let nerpa_output_dirs: Vec<PathBuf> = {
	let root = cli.nerpa_results.clone();

	if root.join("nerpa.log").is_file() {
            vec![root]
	} else {
            std::fs::read_dir(&root)?
		.filter_map(|e| e.ok())
		.map(|e| e.path())
		.filter(|p| p.is_dir())
		.filter(|p| p.join("nerpa.log").is_file())
		.collect()
	}
    };

    let monomers_db = load_monomers_db(&cli.monomers_db_json)?;

    let all_results: Vec<TargetSearchResults> = {
	nerpa_output_dirs
	    .iter()
	    .flat_map(|dir|
		      process_single_output_dir(
			  dir,
			  &monomers_db,
			  &nerpa_ms_cfg,
			  &debug_cfg,
		      ))
	    .collect()
    };
	
    plot(&all_results, &debug_cfg.nerpa_ms_dir())?;

    Ok(())
}
				     
