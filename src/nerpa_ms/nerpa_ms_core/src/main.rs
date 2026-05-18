mod algo;
mod cli;
mod data_types;
mod io;

use std::path::{Path, PathBuf};

use algo::algo_main::process_input_item;
use algo::gen_new_variants::NewVariantWithOptPaths;
use clap::Parser;
use anyhow::{Context, Result};
use data_types::dag::DAG;
use data_types::parsed_rban_record::Parsed_rBAN_Record;
use data_types::{monomer_graph::MonomerGraph, monomers_db::load_monomers_db};
use io::draw_molecules::draw_output_variants;
use io::input::InputItem;
use io::output::{OutputItem, write_output};
use io::draw_hmm_dag::{draw_hmm_dag_opt_paths};
use nerpa_ms_core::data_types::monomers_db::MonomersDB;
use crate::io::draw_nerpa_hmm::draw_nerpa_hmm_with_linearization;

use crate::data_types::alignment::Alignment;

		      

fn main() -> Result<()> {
    println!("Nerpa-MS variant generation has started...");
    let cli = cli::Cli::parse();
    println!("cli received: {:#?}", cli);
    let input_items = io::input::get_input(&cli)
	.context("Failed to load input")?;
    println!("Loaded {} input items", input_items.len());

    println!("Loading monomers database from {}", cli.monomers_db_json.display());
    let monomers_db = load_monomers_db(&cli.monomers_db_json);

    println!("Generating variants...");
    let output_items: Vec<OutputItem> = {
	input_items
	    .iter()
	    .map(|item| { process_input_item(item, &monomers_db, &cli) })
	    .collect()
    };

    write_output(&output_items,
		 &cli.out.join("new_variants.json"))?;

    if cli.draw_output_variants {
	let original_records: Vec<&Parsed_rBAN_Record> = {
	    input_items
		.iter()
		.map(|item| &item.rban_record)
		.collect()
	};
	let new_variants_figures_dir: PathBuf = {
	    cli.out
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
