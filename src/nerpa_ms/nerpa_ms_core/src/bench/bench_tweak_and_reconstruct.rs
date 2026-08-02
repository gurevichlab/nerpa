use std::collections::HashMap;
use std::path::Path;

use crate::algo::apply_modifications::{AlteredMonomerGraph, apply_modifications, modifications_consistent};
use crate::algo::graph_to_dag::get_inserts_for_linearization;
use crate::data_types::monomer_graph::MonomerGraph;
use crate::data_types::graph_modifications::GraphModification;
use crate::bench::find_target::TargetSearchResults;
use crate::data_types::parsed_rban_record::Parsed_rBAN_Record;
use crate::data_types::common_types::MonomerIdx;
use crate::data_types::config::{DebugConfig, NerpaMS_Config};
use crate::io::input::InputItem;
use crate::io::output::OutputItem;
use crate::data_types::monomers_db::MonomersDB;
use crate::algo::algo_main::process_input_item;

use super::find_target::TargetSearchData;


pub fn get_all_modifications<'a>(
    monomer_graph: &MonomerGraph,
    linearization: &[MonomerIdx],
    monomers_db: &'a MonomersDB,
    ) -> Vec<Vec<Vec<GraphModification<'a>>>> { // type(sub/ins/del) -> place of modification -> list of possible modifications
    let max_inserts = usize::MAX; // or some smaller number if you want to limit the number of insertions per position
    let max_subs = usize::MAX; // or some smaller number if you want to limit the number of substitutions per position

    let inserts_per_position: Vec<Vec<GraphModification>> =
	get_inserts_for_linearization(
	    monomer_graph,
	    linearization,
	    monomers_db,
	    max_inserts
	);

    let mut subs_per_position: Vec<Vec<GraphModification>> = Vec::new();
    for monomer_idx in linearization.iter().cloned() {
        let subs: Vec<GraphModification> = {
            monomer_graph.possible_substitutions(monomer_idx, monomers_db, false)
		.into_iter()
                .take(max_subs)
                .map(|entry| GraphModification::Substitute {
		    monomer_idx,
		    mon_db_entry: entry,
		})
                .collect()
        };
	subs_per_position.push(subs);
    }
    
    let mut dels_per_position = Vec::new();
    for monomer_idx in linearization.iter().cloned() {
	if monomer_graph.can_remove(monomer_idx, monomers_db) {
	    let del_mod = GraphModification::Remove { monomer_idx };
	    dels_per_position.push(vec![del_mod]);
	}
    }
    
    let all_mods = vec![inserts_per_position, subs_per_position, dels_per_position];

    let mut all_mods_non_empty = vec![Vec::new(); all_mods.len()];
    for (i, mods_for_type) in all_mods.into_iter().enumerate() {
	for mods_for_position in mods_for_type.into_iter() {
	    if !mods_for_position.is_empty() {
		all_mods_non_empty[i].push(mods_for_position);
	    }
	}
	}
	
	all_mods_non_empty
}

pub fn sample_modifications_unchecked<'a>(
    all_mods: &[Vec<Vec<GraphModification<'a>>>],
    n: usize,
    rng: &mut impl rand::Rng,
) -> Vec<GraphModification<'a>> {
    let mut sampled_mods = Vec::new();
    for _ in 0..n {
	let mod_type_idx = rng.gen_range(0..all_mods.len());
	let mods_for_type = &all_mods[mod_type_idx];
	let position_idx = rng.gen_range(0..mods_for_type.len());
	let mods_for_position = &mods_for_type[position_idx];
	let mod_idx = rng.gen_range(0..mods_for_position.len());
	sampled_mods.push(mods_for_position[mod_idx].clone());
    }
    sampled_mods
}

pub fn sample_modifications<'a>(
    all_mods: &[Vec<Vec<GraphModification<'a>>>],
    n: usize,
    rng: &mut impl rand::Rng,
) -> Vec<GraphModification<'a>> {
    let max_iterations = 100;
    for _ in 0..max_iterations {
       let mods = sample_modifications_unchecked(all_mods, n, rng);
	if modifications_consistent(&mods) {
	    return mods;
	}
    }
    panic!("Failed to sample consistent modifications after {} iterations", max_iterations);
}

pub fn sample_graph_neighborhood(
    monomer_graph: &MonomerGraph,
    linearization: &[MonomerIdx],
    monomers_db: &MonomersDB,
    n_mods: usize,
    n_samples: usize,
    rng: &mut impl rand::Rng,
) -> Vec<AlteredMonomerGraph> {
	let all_mods = get_all_modifications(monomer_graph, linearization, monomers_db);
	let mut sampled_graphs = Vec::new();
	for _ in 0..n_samples {
	    let mods = sample_modifications(&all_mods, n_mods, rng);
	    let maybe_new_graph = apply_modifications(
		monomer_graph,
		&mods,
		monomers_db,
		false // debug_stdout
	    );
	    let new_graph = maybe_new_graph
		.expect("Failed to apply sampled modifications");

	    sampled_graphs.push(new_graph);
	}
	sampled_graphs
}


pub fn tweak_and_reconstruct(
    input_item: &InputItem,
    monomers_db: &MonomersDB,
    config: &NerpaMS_Config,
    debug_cfg: &DebugConfig,
    n_mods: usize,
    rng: &mut impl rand::Rng,
) -> TargetSearchResults {
    let orig_mon_graph = MonomerGraph::from(&input_item.rban_record);
    let orig_linearization = &input_item.linearization;

    let tweaked: AlteredMonomerGraph = {
	sample_graph_neighborhood(
	    &orig_mon_graph,
	    orig_linearization,
	    monomers_db,
	    n_mods,
	    1, // just one sample for now, can increase if want to test multiple tweaks
	    rng,
	)
	    .into_iter()
	    .next()
	    .expect("Failed to sample any altered graphs")
    };

    let new_input_item = InputItem {
	bgc_variant: input_item.bgc_variant.clone(),
	hmm: input_item.hmm.clone(),
	rban_record: Parsed_rBAN_Record::from(&tweaked.new_monomer_graph),
	linearization: tweaked.linearization.clone(),
	score: input_item.score, // not really meaningful for the tweaked graph, but can keep the same score for testing purposes
    };

    let nerpa_ms_out: OutputItem  = process_input_item(
	&new_input_item,
	monomers_db,
	config,
	debug_cfg
    );
	
    let template = MonomerGraph::from(&new_input_item.rban_record);
    let template_linearization = &new_input_item.linearization;
    let target = MonomerGraph::from(&input_item.rban_record);
    let target_linearization = &input_item.linearization;
    let target_search_data = TargetSearchData {
	hmm: &input_item.hmm,
	template: &template,
	template_linearization,
	target: &target,
	target_linearization,
	monomers_db,
	nerpa_root: &debug_cfg.nerpa_root
    };
	
    TargetSearchResults::new(
	&target_search_data,
	&nerpa_ms_out.new_variants,
    )
}
