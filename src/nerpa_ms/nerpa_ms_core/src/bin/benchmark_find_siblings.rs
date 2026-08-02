#![allow(non_snake_case, non_camel_case_types, non_upper_case_globals)]

use nerpa_ms_core::algo::algo_main::process_input_item;
use nerpa_ms_core::bench::find_target::TargetSearchResults;
use nerpa_ms_core::bench::bench_find_siblings::find_siblings;
use nerpa_ms_core::data_types::monomers_db::load_monomers_db;
use nerpa_ms_core::data_types::config::{DebugConfig, NerpaMS_Config};
use nerpa_ms_core::data_types::monomers_db::MonomersDB;
use nerpa_ms_core::io::input::{InputItem, get_input_from_nerpa_results};
use nerpa_ms_core::io::output::OutputItem;
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

    /// Precomputed Nerpa-MS results for the given Nerpa rsults. If not provided, Nerpa-MS will be run on the fly for each input item
    #[arg(long)]
    pub nerpa_ms_results: Option<PathBuf>,

    /// Maximum number of edits allowed in a structure
    #[arg(long, required_unless_present = "nerpa_ms_results")]
    pub max_edits: Option<usize>,

    /// The number of variants generated for each number of edits
    #[arg(long, required_unless_present = "nerpa_ms_results")]
    pub num_variants_per_num_edits: Option<usize>,

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

fn load_nerpa_matches_by_bgc_nrp(nerpa_results_dir: &PathBuf) -> HashMap<String, HashMap<String, InputItem>> {
    let mut nerpa_matches_by_bgc_nrp: HashMap<String, HashMap<String, InputItem>> = HashMap::new();

    let results_subdirs = std::fs::read_dir(nerpa_results_dir)
	.expect(&format!("Failed to read Nerpa results dir: {}", nerpa_results_dir.display()));

    for maybe_results_subdir in results_subdirs {
	let subdir_entry = match maybe_results_subdir {
	    Ok(s) => s,
	    Err(_) => continue,
	};

	if !subdir_entry.path().is_dir() {
	    continue;
	}
	let subdir = subdir_entry.path();

	let bgc_name = {
	    subdir
		.file_name()
		.unwrap()
		.to_string_lossy()
		.to_string()
	};
	let maybe_input_items: Result<Vec<InputItem>> =
	    get_input_from_nerpa_results(
		&subdir, None, None
	    );

	match maybe_input_items {
	    Ok(input_items) => {
		let input_items_by_nrp: HashMap<String, InputItem> = {
		    input_items
			.into_iter()
			.map(|item| (item.rban_record.compound_id.clone(), item))
			.collect()
		};

		nerpa_matches_by_bgc_nrp.insert(bgc_name, input_items_by_nrp);
	    },
	    Err(e) => {
		println!("Failed to load Nerpa results from {}: {}",
			 subdir.display(), e);
		continue;
	    }
	};

    }

	nerpa_matches_by_bgc_nrp
}

fn write_results(results: &[TargetSearchResults], out: &Path) -> Result<()> {
    std::fs::create_dir_all(out.parent().unwrap())?; // Create parent directories
    let json = serde_json::to_string_pretty(results)?;
    std::fs::write(out, json)?;
    Ok(())
}

fn load_nerpa_ms_results_by_bgc_nrp(path: &PathBuf) -> HashMap<String, HashMap<String, OutputItem>> {
    let mut nerpa_ms_results_by_bgc_nrp: HashMap<String, HashMap<String, OutputItem>> = HashMap::new();

    let results_subdirs = std::fs::read_dir(path)
        .expect(&format!("Failed to read Nerpa-MS results dir: {}", path.display()));

    for maybe_results_subdir in results_subdirs {
	let subdir_entry = match maybe_results_subdir {
            Ok(s) => s,
            Err(_) => continue,
	};

	if !subdir_entry.path().is_dir() {
            continue;
	}
	let subdir = subdir_entry.path();

        let bgc_name = {
	    subdir
		.file_name()
		.unwrap()
		.to_string_lossy()
		.to_string()
	};
        let results_json_path = subdir.join("new_variants.json");

        let json_str = std::fs::read_to_string(&results_json_path)
	    .expect(&format!(
		"Failed to read Nerpa-MS results file: {}",
		results_json_path.display()
            ));

        let output_items: Vec<OutputItem> = {
	    serde_json::from_str(&json_str)
		.expect(&format!(
		    "Failed to parse Nerpa-MS results JSON: {}",
		    results_json_path.display()
		))
	};

	let output_items_by_nrp: HashMap<String, OutputItem> = {
	    output_items
		.into_iter()
		.map(|item| (item.compound_id.clone(), item))
		.collect()
	};

        nerpa_ms_results_by_bgc_nrp.insert(bgc_name, output_items_by_nrp);
    }

    nerpa_ms_results_by_bgc_nrp
}


fn compute_nerpa_ms_results_by_bgc_nrp(
    nerpa_matches_by_bgc_nrp: &HashMap<String, HashMap<String, InputItem>>,
    monomers_db: &MonomersDB,
    max_edits: usize,
    num_variants_per_num_edits: usize,
    nerpa_root: &PathBuf,
) -> HashMap<String, HashMap<String, OutputItem>> {
    let nerpa_ms_cfg = NerpaMS_Config {
	max_nerpa_matches: None,
	min_nerpa_score: None,
	max_edits,
	num_variants_per_num_edits,
    };

    let debug_cfg = DebugConfig {
	write_alignments: false,
	debug_stdout: false,
	draw_hmm_opt_paths: false,
	draw_hmm_dag_opt_paths: false,
	draw_output_variants: false,
	out: PathBuf::from("debug"),
	nerpa_root: nerpa_root.clone(),
    };

    let mut nerpa_ms_variants_by_bgc_nrp: HashMap<String, HashMap<String, OutputItem>> = HashMap::new();
    for (bgc_name, nerpa_matches_by_nrp) in nerpa_matches_by_bgc_nrp {
	let nerpa_ms_variants_by_nrp: HashMap<String, OutputItem> = {
	    nerpa_matches_by_nrp
		.iter()
		.map(|(nrp_id, item)| {
		    let nerpa_ms_variants = process_input_item(
			item,
			&monomers_db,
			&nerpa_ms_cfg,
			&debug_cfg,
		    );
		    (nrp_id.clone(), nerpa_ms_variants)
		})
		.collect()
	};
	nerpa_ms_variants_by_bgc_nrp.insert(bgc_name.clone(), nerpa_ms_variants_by_nrp);
    }

    nerpa_ms_variants_by_bgc_nrp
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let nerpa_matches_by_bgc_nrp: HashMap<String, HashMap<String, InputItem>> =
	load_nerpa_matches_by_bgc_nrp(&cli.nerpa_results);

    let monomers_db: MonomersDB = {
	load_monomers_db(&cli.monomers_db_json)
	    .expect(&format!("Failed to load monomers database from {}",
			     cli.monomers_db_json.display()))
	};

    let nerpa_ms_variants_by_bgc_nrp: HashMap<String, HashMap<String, OutputItem>> = {
	if let Some(path) = cli.nerpa_ms_results {
	    load_nerpa_ms_results_by_bgc_nrp(&path)
	} else {
	    compute_nerpa_ms_results_by_bgc_nrp(
		&nerpa_matches_by_bgc_nrp,
		&monomers_db,
		cli.max_edits
		    .expect("max_edits is required if nerpa_ms_results is not provided"),
		cli.num_variants_per_num_edits
		    .expect("num_variants_per_num_edits is required if nerpa_ms_results is not provided"),
		&cli.nerpa_root,
	    )
	}
    };

    let nerpa_only_bgcs: Vec<String> = {
	nerpa_matches_by_bgc_nrp
	    .keys()
	    .filter(|bgc_name| !nerpa_ms_variants_by_bgc_nrp.contains_key(*bgc_name))
	    .cloned()
	    .collect()
    };
    if !nerpa_only_bgcs.is_empty() {
	println!("Warning: Nerpa-MS results not found for the following BGCs: {:?}",
		 nerpa_only_bgcs);
    }

    let mut all_results: Vec<TargetSearchResults> = Vec::new();
    for (i, (bgc_name, nerpa_ms_variants_by_nrp)) in nerpa_ms_variants_by_bgc_nrp.iter().enumerate() {
	println!(
	    "{}/{}. Processing BGC {}...",
	    i, nerpa_ms_variants_by_nrp.len(),
	    bgc_name
	);
	let nerpa_matches_by_nrp: &HashMap<String, InputItem> = {
	    nerpa_matches_by_bgc_nrp
		.get(bgc_name)
		.expect(&format!("No Nerpa results found for BGC {}",
				 bgc_name))
	};

	let target_search_results: Vec<TargetSearchResults> = find_siblings(
	    nerpa_matches_by_nrp,
	    nerpa_ms_variants_by_nrp,
	    &monomers_db,
	    &cli.nerpa_root,
	);

	all_results.extend(target_search_results);
    }

    write_results(&all_results, &cli.out.join("benchmarking_results.json"))?;

    Ok(())
}
				     
