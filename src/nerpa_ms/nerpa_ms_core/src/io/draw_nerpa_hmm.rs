use std::path::Path;

use crate::data_types::{common_types::MonomerIdx, hmm::{HMM, StateIdx}, parsed_rban_record::Parsed_rBAN_Record};
use anyhow::{Result, Context, bail};
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
struct ItemForDrawing {
    hmm: HMM,
    opt_path: Vec<StateIdx>,
    emission_labels: Vec<String>,
}

pub fn draw_nerpa_hmm_with_linearization(
    hmm: &HMM,
    rban_record: &Parsed_rBAN_Record,
    opt_path: &[StateIdx],
    linearization: &[MonomerIdx],
    out_file: &Path,
    nerpa_root: &Path,
) -> Result<()> {
    let emission_labels: Vec<String> = {
	linearization.iter()
	    .map(|monomer_idx| -> Result<String> {
		let monomer = rban_record
		    .monomers
		    .get(monomer_idx)
		    .with_context(|| {
			format!(
			    "Monomer index {} in linearization is absent in the rban record for {}",
			    monomer_idx, rban_record.compound_id
			)
		    })?;

		Ok(format!("{}.{}", monomer.name.0, monomer_idx.0))
	    })
	    .collect::<Result<Vec<_>>>()?
    };
		
    draw_nerpa_hmm(
	hmm,
	opt_path,
	&emission_labels,
	out_file,
	nerpa_root
    )
}
pub fn draw_nerpa_hmm(
    hmm: &HMM,
    opt_path: &[StateIdx],
    emission_labels: &[String],
    out_file: &Path,
    nerpa_root: &Path,
    ) -> Result<()> {
    let item_for_drawing = ItemForDrawing {
	hmm: hmm.clone(),
	opt_path: opt_path.to_vec(),
	emission_labels: emission_labels.to_vec(),
    };

    let json_output = serde_json::to_string(&item_for_drawing)
	.context("Failed to serialize item for drawing")?;

    let json_output_path = out_file.with_extension("json");

    if let Some(parent) = json_output_path.parent() {
	std::fs::create_dir_all(parent)
	    .with_context(|| format!("Failed to create parent directory for {}", json_output_path.display()))?;
    }
    std::fs::write(&json_output_path, json_output)
	.with_context(|| format!("Failed to write {}", json_output_path.display()))?;

    let drawing_script = nerpa_root
	.join("scripts")
	.join("draw_hmm_with_opt_path.py");

    let out_file_abs = if out_file.is_absolute() {
	out_file.to_path_buf()
    } else {
	std::env::current_dir()
	    .context("Failed to get current directory")?
	    .join(out_file)
    };
    
    let mut cmd = std::process::Command::new("python"); // keep this first
    cmd.arg(&drawing_script)
	.arg("--input").arg(&json_output_path.canonicalize()?)
	.arg("--output").arg(&out_file_abs);

    cmd.current_dir(nerpa_root);
    cmd.env("PYTHONPATH", nerpa_root);

    let output = cmd.output()
	.with_context(|| format!("Failed to execute drawing script at {}",
				 drawing_script.display()))?;


    let stderr = String::from_utf8_lossy(&output.stderr);
    let stdout = String::from_utf8_lossy(&output.stdout);
    if output.status.success() {
	println!("Drawing script completed successfully.");
	println!("--- stdout ---\n{}\n--- stderr ---\n{}", stdout, stderr);
    } else {
        bail!(
            "Drawing script failed (exit={}):\n--- stdout ---\n{}\n--- stderr ---\n{}",
            output.status,
            stdout,
            stderr
        );
    }
    Ok(())
	
}
