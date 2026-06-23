use std::collections::HashMap;
use std::path::Path;

use crate::algo::gen_new_variants::Altered_rBAN_Record;
use crate::bench::find_target::{TargetSearchData, TargetSearchResults};
use crate::data_types::bgc_variant::BGC_Variant;
use crate::data_types::monomer_graph::MonomerGraph;
use crate::data_types::parsed_rban_record::Parsed_rBAN_Record;
use crate::data_types::common_types::MonomerIdx;
use crate::data_types::config::{DebugConfig, NerpaMS_Config};
use crate::io::input::InputItem;
use crate::io::output::OutputItem;
use crate::data_types::monomers_db::MonomersDB;
use crate::algo::algo_main::process_input_item;


/// Find sibling compounds: for each input item as template, try to find every other item as target.
///
/// Assumes `input_items` all belong to the same BGC run.
///
/// NOTE: This function needs access to `nerpa_root` to call `find_target`.
pub fn find_siblings(
    nerpa_matches_by_nrp: &HashMap<String, InputItem>,
    nerpa_ms_variants_by_nrp: &HashMap<String, OutputItem>,
    monomers_db: &MonomersDB,
    nerpa_root: &Path,
) -> Vec<TargetSearchResults> {

    let targets_with_linearizations: Vec<(&Parsed_rBAN_Record, &[MonomerIdx])> = {
	nerpa_matches_by_nrp
	    .values()
	    .map(|item| (&item.rban_record, &item.linearization[..]))
	    .collect()
    };

    let bgc_variant: &BGC_Variant = {
	&nerpa_matches_by_nrp
	    .values()
	    .next()
	    .expect("No input items found for this BGC")
	    .bgc_variant
    };

    println!("BGC:\n{}", bgc_variant);

    let mut results: Vec<TargetSearchResults> = Vec::new();

    for (i, (nrp_id, input_item)) in (nerpa_matches_by_nrp.iter().enumerate()) {
	let template_lin_str = input_item.rban_record.show_linearization(&input_item.linearization)
	    .unwrap_or_else(|_| "Failed to show linearization".to_string());
	println!(
	    "{}/{}. Template:\n{}\t{}",
	    i + 1,
	    nerpa_matches_by_nrp.len(),
	    input_item.rban_record.compound_id,
	    template_lin_str
	);

	let nerpa_ms_output = nerpa_ms_variants_by_nrp
	    .get(nrp_id)
	    .unwrap_or_else(|| panic!("No Nerpa-MS output found for NRP ID {}", nrp_id));

	println!("Targets:");
	for (target, target_linearization) in targets_with_linearizations.iter() {
	    if target.compound_id == nerpa_ms_output.compound_id {
		continue;
	    }
	    let target_lin_str = target.show_linearization(target_linearization)
		.unwrap_or_else(|_| "Failed to show linearization".to_string());
	    println!(
		"{}\t{}",
		target.compound_id,
		target_lin_str
	    );

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
