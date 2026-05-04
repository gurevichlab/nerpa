use std::collections::HashMap;
use std::path::Path;
use serde::Serialize;
use crate::data_types::common_types::LogProb;
use crate::data_types::{hmm::BGC_ID, parsed_rban_record::Parsed_rBAN_Record};
use crate::algo::algo_main::Altered_rBAN_Record;
use crate::io::output::OutputItem;

#[derive(Debug, Clone, Serialize)]
struct ItemForDrawing {
    bgc_id: BGC_ID,
    original: Parsed_rBAN_Record,
    original_score: LogProb,
    new_variant: Altered_rBAN_Record,
}

use anyhow::{anyhow, bail, Context};


pub fn draw_output_variants(
	output_items: &[OutputItem],
	original_records: &[&Parsed_rBAN_Record],
	nerpa_root: &Path,
	output_dir: &Path,
) -> anyhow::Result<()> {
    let records_by_id: HashMap<String, &Parsed_rBAN_Record> = original_records
	.iter()
	.map(|&record| (record.compound_id.clone(), record))
	.collect();

    let mut items_for_drawing: Vec<ItemForDrawing> = Vec::new();
    for output_item in output_items {
	for new_variant in &output_item.new_variants {
	    let original = records_by_id
		.get(&new_variant.new_record.compound_id)
		.with_context(|| {
		    format!(
			"Original record not found for compound_id={}",
			new_variant.new_record.compound_id
		    )
		})?;
	    items_for_drawing.push(ItemForDrawing {
		bgc_id: output_item.bgc_variant_id.bgc_id.clone(),
		original: (*original).clone(),
		original_score: output_item.original_score,
		new_variant: (*new_variant).clone(),
	    });
	}
    }

    println!("Created {} items for drawing...", items_for_drawing.len());

    let json_output = serde_json::to_string(&items_for_drawing)
        .context("Failed to serialize items for drawing")?;

    let cwd = std::env::current_dir().context("Failed to get current directory")?;
    let output_dir_abs = if output_dir.is_absolute() {
	output_dir.to_path_buf()
    } else {
	cwd.join(output_dir)
    };
    let json_output_path = output_dir_abs.join("items_for_drawing.json");
    if let Some(parent) = json_output_path.parent() {
	std::fs::create_dir_all(parent)
	    .with_context(|| format!("Failed to create parent directory for {}", json_output_path.display()))?;
    }
    std::fs::write(&json_output_path, json_output)
	.with_context(|| format!("Failed to write {}", json_output_path.display()))?;

    let drawing_script = nerpa_root
	.join("scripts")
	.join("draw_nerpa_ms_variants.py");

    let mut cmd = std::process::Command::new("python"); // keep this first
    cmd.arg(&drawing_script)
	.arg("--input_json").arg(&json_output_path)
	.arg("--output_dir").arg(output_dir_abs);

    cmd.current_dir(nerpa_root);
    cmd.env("PYTHONPATH", nerpa_root);

    let output = cmd.output()
	.with_context(|| format!("Failed to execute drawing script at {}",
				 drawing_script.display()))?;


    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);
        bail!(
            "Drawing script failed (exit={}):\n--- stdout ---\n{}\n--- stderr ---\n{}",
            output.status,
            stdout,
            stderr
        );
    }
    Ok(())
}
