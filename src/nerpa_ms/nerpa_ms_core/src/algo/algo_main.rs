use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crate::algo::graph_to_dag::create_mod_graph;
use crate::cli;
use crate::algo::gen_new_variants::{Altered_rBAN_Record, NewVariantWithOptPaths, generate_new_variants_with_opt_paths};
use crate::data_types::config::{DebugConfig, NerpaMS_Config};
use anyhow::{Context, Result};
use crate::data_types::mod_graph::ModGraph;
use crate::data_types::parsed_rban_record::Parsed_rBAN_Record;
use crate::data_types::{monomer_graph::MonomerGraph, monomers_db::load_monomers_db};
use crate::io::input::InputItem;
use crate::io::output::OutputItem;
use crate::io::draw_hmm_dag::{draw_hmm_dag_opt_paths};
use crate::data_types::monomers_db::MonomersDB;
use crate::io::draw_nerpa_hmm::draw_nerpa_hmm_with_linearization;

use crate::data_types::alignment::Alignment;

fn draw_hmm_dags_optimal_paths(
    new_variants_with_opt_paths: &[NewVariantWithOptPaths],
    item: &InputItem,
    dag: &ModGraph<'_>,
    figs_dir: &Path,
    nerpa_root: &Path
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
		&figs_dir.join(format!("{i}.svg")),
		nerpa_root
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
    figs_dir: &Path,
    nerpa_root: &Path
) {
	let bgc_id_short = {
	    &item.hmm
	    .bgc_variant_id
	    .bgc_id
	    .to_str_short()
	};
	for (i, new_variant_with_opt_paths) in new_variants_with_opt_paths.iter().enumerate() {
	    let res = draw_nerpa_hmm_with_linearization(
		&item.hmm,
		&new_variant_with_opt_paths.new_variant.new_record,
		&new_variant_with_opt_paths.hmm_path,
		&new_variant_with_opt_paths.linearization,
		&figs_dir.join(format!("{i}.svg")),
		nerpa_root
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

pub fn process_input_item(
    item: &InputItem,
    monomers_db: &MonomersDB,
    cfg: &NerpaMS_Config,
    debug_cfg: &DebugConfig,
) -> OutputItem {
    if debug_cfg.debug_stdout {
	    println!("Processing BGC {} and compound {}...",
		         item.hmm.bgc_variant_id.bgc_id.to_str_short(),
		         item.rban_record.compound_id);
	    println!("BGC:");
	    println!("{}", item.bgc_variant);
	    println!("NRP linearization:");
	    println!("{}\n", item.rban_record.show_linearization(&item.linearization).unwrap_or_else(|_| "Failed to show linearization".to_string()));
	    println!("Creating monomer graph...");
    }

    let monomer_graph = MonomerGraph::from(&item.rban_record);

    if debug_cfg.debug_stdout {
	    println!("Creating DAG from linearization...");
    }

    let mod_graph = create_mod_graph(
	    &monomer_graph,
	    &item.linearization,
	    monomers_db
    );

    if debug_cfg.debug_stdout {
	    println!("Creating new variants with optimal paths...");
    }

    let new_variants_with_opt_paths: Vec<NewVariantWithOptPaths> = {
        generate_new_variants_with_opt_paths(
	        &item.hmm,
	        &monomer_graph,
	        &mod_graph,
	        cfg.max_edits,
	        cfg.num_variants_per_num_edits,
	        monomers_db,
	        debug_cfg.debug_stdout
        )
    };

    let bgc_id_short = {
	    &item.hmm
	        .bgc_variant_id
	        .bgc_id
	        .to_str_short()
    };

    if debug_cfg.write_alignments {
	    let alignments_dir = {
	        debug_cfg.out
		        .join("alignments")
		        .join(format!("{bgc_id_short}_{}", &item.rban_record.compound_id))
	    };
	    let r = write_alignments(&new_variants_with_opt_paths, item, &alignments_dir);
	if let Err(e) = r {
	    eprintln!("Failed to write alignments for bgc {bgc_id_short} and compound {}: {e}",
		      &item.rban_record.compound_id);
	}
    }

    if debug_cfg.draw_hmm_opt_paths {
	    let figs_dir = {
	        debug_cfg.out
		        .join("figures")
		        .join("opt_paths")
		        .join(format!("{bgc_id_short}_{}", &item.rban_record.compound_id))
	    };
	    if debug_cfg.debug_stdout {
	        println!("Drawing HMM with optimal paths...");
	    }
	    draw_hmms_optimal_paths(
	        &new_variants_with_opt_paths,
	        item,
	        &figs_dir,
	        &debug_cfg.nerpa_root
	    );
    }

    if debug_cfg.draw_hmm_dag_opt_paths {
	    if debug_cfg.debug_stdout {
	        println!("Drawing DAG with optimal paths...");
	    }
	    let figs_dir = {
	        debug_cfg.out
		        .join("figures")
		        .join("hmm_dag_opt_paths")
		        .join(format!("{bgc_id_short}_{}", &item.rban_record.compound_id))
	    };
	    draw_hmm_dags_optimal_paths(
	        &new_variants_with_opt_paths,
	        item,
	        &mod_graph,
	        &figs_dir,
	        &debug_cfg.nerpa_root
	    );
    }

    let new_variants: Vec<Altered_rBAN_Record> = {
	    new_variants_with_opt_paths
	        .into_iter()
	        .map(|v| v.new_variant)
	        .collect()
    };
    OutputItem {
	    bgc_variant_id: item.hmm.bgc_variant_id.clone(),
	    compound_id: item.rban_record.compound_id.clone(),
	    original_score: item.score,
        original_record: item.rban_record.clone(),
	    new_variants,
    }
}
