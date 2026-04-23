mod algo;
mod cli;
mod data_types;
mod io;

use clap::Parser;
use anyhow::Result;
use data_types::{monomer_graph::MonomerGraph, monomers_db::load_monomers_db};
use io::output::{OutputItem, write_output};
use algo::{algo_main::{generate_new_variants_with_opt_paths}, graph_to_dag::create_dag};

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
	let dag = create_dag(&item.rban_record,
			     &item.linearization,
			     &monomers_db);
	let monomer_graph = MonomerGraph::from(&item.rban_record);
	let new_variants_with_opt_paths = generate_new_variants_with_opt_paths(
	    &item.hmm,
	    &monomer_graph,
	    &dag,
	    cli.max_edits,
	    cli.num_variants_per_num_edits
	);

	let bgc_id_short = &item.hmm
	    .bgc_variant_id
	    .bgc_id
	    .to_str_short();
	let figs_dir = cli.out
	    .join("figures")
	    .join("opt_paths")
	    .join(format!("{bgc_id_short}_{&item.rban_record.compound_id}"))

	for (i, new_variant_with_opt_paths) in new_variants_with_opt_paths.iter().enumerate() {
	    let res = draw_hmm_dag_opt_paths(
		&item.hmm,
		&dag,
		&new_variant_with_opt_paths.hmm_path,
		&new_variant_with_opt_paths.dag_path,
		&fig_dir.join(format!("{i}"))
	    );
	    if let Err(e) = res {
		eprintln!("Failed to draw HMM-DAG optimal paths for variant {i} of
 bgc {bgc_id_short}, compound {}: {e}", &item.rban_record.compound_id);
	    }
	}

	let new_variants = new_variants_with_opt_paths
	    .into_iter()
	    .map(|v| v.new_variant)
	    .collect::<Vec<_>>();
	output_items.push(OutputItem {
	    bgc_variant_id: item.hmm.bgc_variant_id.clone(),
	    compound_id: item.rban_record.compound_id.clone(),
	    new_variants,
	});
    }
	
    write_output(&output_items, &cli.out);
    Ok(())
}
