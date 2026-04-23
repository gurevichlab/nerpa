use serde::Serialize;

use crate::{algo::algo_main::Altered_rBAN_Record, data_types::hmm::BGC_Variant_ID};

#[derive(Debug, Clone, Serialize)]
pub struct OutputItem {
    pub bgc_variant_id: BGC_Variant_ID,
    pub compound_id: String,
    pub new_variants: Vec<Altered_rBAN_Record>
}

use std::path::Path;

pub fn write_output(output_items: &[OutputItem], output_path: &Path) {
	// For simplicity, we will write the output as a JSON file containing a list of OutputItem
	let json_output = serde_json::to_string_pretty(output_items)
		.expect("Failed to serialize output items to JSON");
	std::fs::write(output_path, json_output)
		.expect("Failed to write output to file");
}
    
