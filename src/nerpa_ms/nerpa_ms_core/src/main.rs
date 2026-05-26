mod algo;
mod cli;
mod data_types;
mod io;

use std::path::PathBuf;

use algo::algo_main::process_input_item;
use clap::Parser;
use anyhow::{Context, Result};
use data_types::parsed_rban_record::Parsed_rBAN_Record;
use data_types::monomers_db::load_monomers_db;
use io::draw_molecules::draw_output_variants;
use io::output::{OutputItem, write_output};
use data_types::config::{DebugConfig, NerpaMS_Config};

		      

fn main() -> Result<()> {
    println!("Nerpa-MS variant generation has started...");
    let cli = cli::Cli::parse();
    println!("cli received: {:#?}", cli);

    let nerpa_ms_config = NerpaMS_Config::from_cli(&cli);
    let debug_config = DebugConfig::from_cli(&cli);

    let input_items = io::input::get_input(&cli)
	.context("Failed to load input")?;
    println!("Loaded {} input items", input_items.len());

    println!("Loading monomers database from {}", cli.monomers_db_json.display());
    let monomers_db = load_monomers_db(&cli.monomers_db_json)?;

    println!("Generating variants...");
    let output_items: Vec<OutputItem> = {
	input_items
	    .iter()
	    .map(|item| {
		process_input_item(
		    item,
		    &monomers_db,
		    &nerpa_ms_config,
		    &debug_config
		)
	    })
	    .collect()
    };

    write_output(&output_items,
		 &cli.out.join("new_variants.json"))?;

    if debug_config.draw_output_variants {
	let original_records: Vec<&Parsed_rBAN_Record> = {
	    input_items
		.iter()
		.map(|item| &item.rban_record)
		.collect()
	};
	let new_variants_figures_dir: PathBuf = {
	    debug_config.out
		.join("figures")
		.join("new_variants")
	};

	println!("Drawing new variants...");
	draw_output_variants(&output_items,
			     &original_records,
			     &cli.nerpa_root,
			     &new_variants_figures_dir)?;
    }

    Ok(())
}
