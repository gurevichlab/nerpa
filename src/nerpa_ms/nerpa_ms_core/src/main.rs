mod algo;
mod cli;
mod data_types;
mod io;

use clap::Parser;
use anyhow::Result;
use data_types::monomers_db::load_monomers_db;
use io::output::{OutputItem, write_output};
use algo::algo_main::generate_new_variants_per_weight;

fn main() -> Result<()> {
    println!("Nerpa-MS variant generation has started...");
    let cli = cli::Cli::parse();
    println!("cli received: {:#?}", cli);
    let input_items = io::input::get_input(&cli)
	.map_err(|e| anyhow::anyhow!("Failed to load input:\n{e}"))?;
    println!("Loaded {} input items", input_items.len());

    println!("Loading monomers database from {}", cli.monomers_db_json.display());
    let monomers_db = load_monomers_db(&cli.monomers_db_json);

    println!("Generating variants...");
    let mut output_items: Vec<OutputItem> = Vec::with_capacity(input_items.len());
    for item in input_items {
	let variants_by_weight = generate_new_variants_per_weight(
	    &item.hmm,
	    &item.rban_record,
	    &item.linearization,
	    &monomers_db,
	    cli.max_edits,
	    cli.num_variants_per_num_edits
	);
	output_items.push(OutputItem {
	    bgc_variant_id: item.hmm.bgc_variant_id.clone(),
	    compound_id: item.rban_record.compound_id.clone(),
	    variants_by_weight,
	});
    }
	
    write_output(&output_items, &cli.out);
    Ok(())
}
