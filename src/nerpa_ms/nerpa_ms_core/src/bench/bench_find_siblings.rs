use std::collections::HashMap;
use std::path::Path;

use crate::algo::gen_new_variants::Altered_rBAN_Record;
use crate::bench::find_target::{TargetSearchData, TargetSearchResults};
use crate::data_types::monomer_graph::MonomerGraph;
use crate::data_types::parsed_rban_record::Parsed_rBAN_Record;
use crate::data_types::common_types::MonomerIdx;
use crate::data_types::config::{DebugConfig, NerpaMS_Config};
use crate::io::input::InputItem;
use crate::io::output::OutputItem;
use crate::data_types::monomers_db::MonomersDB;
use crate::algo::algo_main::process_input_item;

fn lookup_distance(
    compound_distances: &HashMap<(String, String), usize>,
    template: &str,
    target: &str,
) -> Option<usize> {
    let direct_key = (template.to_string(), target.to_string());
    if let Some(&d) = compound_distances.get(&direct_key) {
	return Some(d);
    }

    let reverse_key = (target.to_string(), template.to_string());
    if let Some(&d) = compound_distances.get(&reverse_key) {
	return Some(d);
    }

    None
}

/// Find sibling compounds: for each input item as template, try to find every other item as target.
///
/// Assumes `input_items` all belong to the same BGC run.
///
/// NOTE: This function needs access to `nerpa_root` to call `find_target`.
/// I assume you have it inside `cfg` (e.g. `cfg.nerpa_root: PathBuf`).
pub fn find_siblings(
    input_items: &[InputItem],
    monomers_db: &MonomersDB,
    nerpa_ms_cfg: &NerpaMS_Config,
    debug_cfg: &DebugConfig,
) -> Vec<TargetSearchResults> {
    // If cfg doesn’t actually store this, change the function signature to accept `nerpa_root: &Path`.
    let nerpa_root: &Path = debug_cfg.nerpa_root.as_path();

    let targets_with_linearizations: Vec<(&Parsed_rBAN_Record, &[MonomerIdx])> = {
	input_items
	    .iter()
	    .map(|item| (&item.rban_record, &item.linearization[..]))
	    .collect()
    };

    let nerpa_ms_outputs: Vec<OutputItem> = {
	input_items
	    .iter()
	    .take(1) // debug
	    .map(|item|
		 process_input_item(
		     item,
		     monomers_db,
		     nerpa_ms_cfg,
		     debug_cfg,
		 ))
	    .collect()
    };


    let mut results: Vec<TargetSearchResults> = Vec::new();

    for (input_item, nerpa_ms_output) in (input_items.iter()
					  .zip(nerpa_ms_outputs.iter())) {
	println!("BGC:\n{}", input_item.bgc_variant);
	println!("NRPs:");
	for (target, target_linearization) in targets_with_linearizations.iter() {
	    println!("{}\t{}", target.compound_id, target.show_linearization(target_linearization).unwrap_or_else(|_| "Failed to show linearization".to_string()));

	}

	    
	for (target, target_linearization) in targets_with_linearizations.iter().cloned() {
	    if target.compound_id == nerpa_ms_output.compound_id {
		continue;
	    }

	    let template_mon_graph = MonomerGraph::from(&input_item.rban_record);
	    let target_mon_graph = MonomerGraph::from(target);
	    let target_search_data = TargetSearchData {
		hmm: &input_item.hmm,
		template: &template_mon_graph,
		template_linearization: &input_item.linearization,
		target: &target_mon_graph,
		target_linearization, 
		monomers_db,
		nerpa_root,
	    };
	    
	    results.push(TargetSearchResults::new(
		&target_search_data,
		&nerpa_ms_output.new_variants,
	    ));
		
	}
    }

    results
}
