use std::collections::HashMap;
use std::path::Path;

use crate::algo::gen_new_variants::Altered_rBAN_Record;
use crate::bench::find_target::{FindTargetResults, find_target};
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
    compound_distances: &HashMap<(String, String), u32>,
) -> Vec<FindTargetResults> {
    // If cfg doesn’t actually store this, change the function signature to accept `nerpa_root: &Path`.
    let nerpa_root: &Path = debug_cfg.nerpa_root.as_path();

    let targets: Vec<&Parsed_rBAN_Record> = {
	input_items
	    .iter()
	    .map(|item| &item.rban_record)
	    .collect()
    };

    let nerpa_ms_outputs: Vec<OutputItem> = {
	input_items
	    .iter()
	    .map(|item|
		 process_input_item(
		     item,
		     monomers_db,
		     nerpa_ms_cfg,
		     debug_cfg,
		 ))
	    .collect()
    };


    let mut results: Vec<FindTargetResults> = Vec::new();

    for target in &targets {
	for nerpa_ms_output in nerpa_ms_outputs.iter() {
	    if target.compound_id == nerpa_ms_output.compound_id {
		continue;
	    }

	    let distance: u32 = {
		let target_id = target.compound_id.clone();
		let template_id = nerpa_ms_output.compound_id.clone();

		*compound_distances
		    .get(&(template_id.clone(), target_id.clone()))
		    .or_else(|| compound_distances.get(&(target_id.clone(), template_id.clone())))
		    .expect(&format!("Distance not found for pair (template={}, target={})",
				     template_id, target_id))
	    };

	    let one_result = find_target(
		nerpa_ms_output,
		target,
		distance,
		nerpa_root,
	    );
	    results.push(one_result);
	}
    }

    results
}
