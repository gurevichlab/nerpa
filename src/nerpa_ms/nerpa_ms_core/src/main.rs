mod algo;
mod cli;
mod data_types;
mod io;

use clap::Parser;
use anyhow::Result;
use data_types::monomers_db::load_monomers_db;
use io::draw_dag::Draw_DAG_Config;

fn main() -> Result<()> {
    let cli = cli::Cli::parse();
    println!("cli received: {:#?}", cli);
    let input_items = io::input::parse_input(&cli.input)
	.unwrap_or_else(|err| panic!("Failed to load (HMM, Parsed_rBAN_Record, Linearization) triplets from {}: {}",
				     cli.input.display(), err));
    println!("loaded {} input items", input_items.len());

    let rban_record = &input_items[0].rban_record;
    let linearization = &input_items[0].linearization;
    let monomers_db = load_monomers_db(&cli.monomers_db_json);
    let dag = algo::graph_to_dag::create_dag(rban_record,
					     linearization,
					     &monomers_db,);
    let out_path = cli.out.join("dag.svg");
    dag.draw_svg(&out_path,
		 &Draw_DAG_Config{node_indexes: false},
		 None)?;
	
    Ok(())
}
