mod algo;
mod cli;
mod data_types;
mod io;

use clap::Parser;
use anyhow::Result;

fn main() -> Result<()> {
    let cli = cli::Cli::parse();
    println!("cli received: {:#?}", cli);
    let input_items = io::input::parse_input(&cli.input)
	.unwrap_or_else(|err| panic!("Failed to load (HMM, Parsed_rBAN_Record, Linearization) triplets from {}: {}",
				     cli.input.display(), err));
    println!("loaded {} input items", input_items.len());
    let fst_item = &input_items[0];
    println!("first item:");
    println!("HMM: {:?}", fst_item.hmm.bgc_variant_id);
    println!("Parsed_rBAN_Record: {:?}", fst_item.rban_record.compound_id);
    println!("Linearization: {:?}", fst_item.linearization);
    Ok(())
}
