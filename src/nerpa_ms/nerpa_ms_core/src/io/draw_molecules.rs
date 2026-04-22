use std::collections::HashMap;
use std::path::Path;
use serde::Serialize;
use crate::data_types::parsed_rban_record::Parsed_rBAN_Record;
use crate::algo::algo_main::Altered_rBAN_Record;
use crate::io::output::OutputItem;

#[derive(Debug, Clone, Serialize)]
struct ItemForDrawing {
	original: Parsed_rBAN_Record,
	new_variant: Altered_rBAN_Record,
}

pub fn draw_new_variants(
    new_variants: &[&Altered_rBAN_Record],
    original_records_by_id: HashMap<String, &Parsed_rBAN_Record>,
    nerpa_root: &Path,
    output_dir: &Path,
) {
    let items_for_drawing: Vec<ItemForDrawing> = new_variants
	.iter()
	.map(|new_variant| {
	    let original = original_records_by_id
		.get(&new_variant.new_record.compound_id)
		.expect("Original record not found for compound_id");
	    ItemForDrawing {
		original: (*original).clone(),
		new_variant: (*new_variant).clone(),
	    }
	})
	.collect();

    let json_output = serde_json::to_string(&items_for_drawing).expect("Failed to serialize items for drawing");
    let json_output_path = output_dir.join("items_for_drawing.json");
    std::fs::write(&json_output_path, json_output)
	.expect("Failed to write items for drawing to file");

    let drawing_script = nerpa_root
	.join("scripts")
	.join("draw_nerpa_ms_variants.py");
    let output = std::process::Command::new("python")
	.arg(drawing_script)
	.arg("--input_json")
	.arg(json_output_path)
	.arg("--output_dir")
	.arg(output_dir)
	.output()
	.expect("Failed to execute drawing script");
}

pub fn draw_output_variants(
	output_items: &[OutputItem],
	original_records: &[Parsed_rBAN_Record],
	nerpa_root: &Path,
	output_dir: &Path,
) {
    let new_variants: Vec<&Altered_rBAN_Record> = output_items
	.iter()
	.flat_map(|item| item.variants_by_weight.iter().flatten())
	.collect();
    let records_by_id: HashMap<String, &Parsed_rBAN_Record> = original_records
	.iter()
	.map(|record| (record.compound_id.clone(), record))
	.collect();

    draw_new_variants(&new_variants,
		      records_by_id,
		      nerpa_root,
		      output_dir);
}
