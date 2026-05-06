mod algo;
mod cli;
mod data_types;
mod io;

use std::path::PathBuf;

use algo::algo_main::NewVariantWithOptPaths;
use clap::Parser;
use anyhow::{Context, Result};
use data_types::dag::DAG;
use data_types::parsed_rban_record::Parsed_rBAN_Record;
use data_types::{monomer_graph::MonomerGraph, monomers_db::load_monomers_db};
use io::draw_molecules::draw_output_variants;
use io::input::InputItem;
use io::output::{OutputItem, write_output};
use io::draw_hmm_dag::{draw_hmm_dag_opt_paths, draw_hmm_opt_path_with_monomers};
use algo::{algo_main::{generate_new_variants_with_opt_paths}, graph_to_dag::create_dag};

use crate::data_types::alignment::Alignment;

fn draw_hmm_dags_optimal_paths(
    new_variants_with_opt_paths: &[NewVariantWithOptPaths],
    item: &InputItem,
    dag: &DAG<'_>,
    figs_dir: &PathBuf
) {
	let bgc_id_short = {
	    &item.hmm
	    .bgc_variant_id
	    .bgc_id
	    .to_str_short()
	};
	for (i, new_variant_with_opt_paths) in new_variants_with_opt_paths.iter().enumerate() {
	    let res = draw_hmm_dag_opt_paths(
		&item.hmm,
		&dag,
		&new_variant_with_opt_paths.hmm_path,
		&new_variant_with_opt_paths.dag_path,
		&figs_dir.join(format!("{i}"))
	    );	     
	    if let Err(e) = res {
		eprintln!("Failed to draw HMM-DAG optimal paths for variant {i} of
 bgc {bgc_id_short}, compound {}: {e}", &item.rban_record.compound_id);
	    }
	}
}

fn draw_hmms_optimal_paths(
	new_variants_with_opt_paths: &[NewVariantWithOptPaths],
	item: &InputItem,
	figs_dir: &PathBuf
) {
	let bgc_id_short = {
	    &item.hmm
	    .bgc_variant_id
	    .bgc_id
	    .to_str_short()
	};
	for (i, new_variant_with_opt_paths) in new_variants_with_opt_paths.iter().enumerate() {
	    let res = draw_hmm_opt_path_with_monomers(
		&item.hmm,
		&new_variant_with_opt_paths.hmm_path,
		&new_variant_with_opt_paths.linearization,
		&item.rban_record,
		&figs_dir.join(format!("{i}"))
	    );	     
	    if let Err(e) = res {
		eprintln!("Failed to draw HMM optimal path with monomers for variant {i} of
 bgc {bgc_id_short}, compound {}: {e}", &item.rban_record.compound_id);
	    }
	}
}

fn write_alignments(
	new_variants_with_opt_paths: &[NewVariantWithOptPaths],
	item: &InputItem,
	alignments_dir: &PathBuf
) -> anyhow::Result<()> {
    if !alignments_dir.exists() {
	std::fs::create_dir_all(alignments_dir)
            .context("Failed to create alignments directory")?;
    }

    for (i, new_variant_with_opt_paths) in new_variants_with_opt_paths.iter().enumerate() {
	let alignment = Alignment::new(
	    &new_variant_with_opt_paths.hmm_path,
	    &new_variant_with_opt_paths.linearization,
	    &item.bgc_variant,
	    &item.hmm,
	    &new_variant_with_opt_paths.new_variant.new_record,
	);
	std::fs::write(alignments_dir.join(format!("alignment_{}.txt", i)), alignment.to_tsv_string_aligned())?;
    }
    Ok(())
}

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
    let mut output_items: Vec<OutputItem> = Vec::with_capacity(input_items.len());
    for item in &input_items {
	println!("Processing BGC {} and compound {}...", item.hmm.bgc_variant_id.bgc_id.to_str_short(), item.rban_record.compound_id);
	println!("Creating monomer graph...");
	let monomer_graph = MonomerGraph::from(&item.rban_record);
	println!("Creating DAG from linearization...");
	let dag = create_dag(&monomer_graph,
			     &item.linearization,
			     &monomers_db);
	println!("Creating new variants with optimal paths...");
	let new_variants_with_opt_paths = generate_new_variants_with_opt_paths(
	    &item.hmm,
	    &monomer_graph,
	    &dag,
	    cli.max_edits,
	    cli.num_variants_per_num_edits,
	    &monomers_db
	);

	let bgc_id_short = {
	    &item.hmm
	    .bgc_variant_id
	    .bgc_id
	    .to_str_short()
	};
	let alignments_dir = {
	    cli.out
		.join("alignments")
		.join(format!("{bgc_id_short}_{}", &item.rban_record.compound_id))
	};

	write_alignments(&new_variants_with_opt_paths, item, &alignments_dir)?;

	println!("Drawing HMM with optimal path...");
	let figs_dir = {
	    cli.out
		.join("figures")
		.join("opt_paths")
		.join(format!("{bgc_id_short}_{}", &item.rban_record.compound_id))
	};
	draw_hmms_optimal_paths(&new_variants_with_opt_paths, item, &figs_dir);

	let new_variants = {
	    new_variants_with_opt_paths
	    .into_iter()
	    .map(|v| v.new_variant)
	    .collect::<Vec<_>>()
	};
	output_items.push(OutputItem {
	    bgc_variant_id: item.hmm.bgc_variant_id.clone(),
	    compound_id: item.rban_record.compound_id.clone(),
	    original_score: item.score,
	    new_variants,
	});
    }
	
    write_output(&output_items,
		 &cli.out.join("new_variants.json"))?;

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
    if let Some(nerpa_root) = &cli.nerpa_root {
	draw_output_variants(&output_items,
			     &original_records,
			     nerpa_root,
    			     &new_variants_figures_dir)?;
    } else {
	eprintln!("Nerpa root not provided, skipping drawing new variants");
    }

    Ok(())
}
