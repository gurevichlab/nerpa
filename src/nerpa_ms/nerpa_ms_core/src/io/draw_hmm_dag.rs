use crate::data_types::common_types::MonomerIdx;
use crate::data_types::dag::{VertexId, DAG};
use crate::data_types::hmm::{StateIdx, HMM};
use crate::data_types::parsed_rban_record::Parsed_rBAN_Record;
use crate::io::draw_dag::Draw_DAG_Config;
use crate::io::draw_nerpa_hmm::draw_nerpa_hmm;
use crate::io::join_svgs::{TitlePlacement, join_svgs_vertical};
use anyhow::Result;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

use super::join_svgs::svg_with_title;

pub fn draw_hmm_dag_opt_paths(
    hmm: &HMM,
    dag: &DAG<'_>,
    hmm_path: &[StateIdx],
    dag_path: &[VertexId],
    out_file: &Path,
) -> Result<()> {
    let out_dir = out_file
        .parent()
        .ok_or_else(|| anyhow::anyhow!("Output file {:?} has no parent directory", out_file))?;
    let hmm_svg_path = out_dir.join("hmm.svg");
    let dag_svg_path = out_dir.join("dag.svg");

    draw_nerpa_hmm(hmm,
		   Some(&hmm_svg_path),
		   Some(hmm_path),
		   false)?;

    dag.draw_svg(&dag_svg_path,
		 &Draw_DAG_Config { node_indexes: true },
		 Some(dag_path))?;
    join_svgs_vertical(&[&hmm_svg_path, &dag_svg_path], out_file)?;
    Ok(())
}

use anyhow::Context;

pub fn draw_hmm_opt_path_with_monomers(
    hmm: &HMM,
    hmm_path: &[StateIdx],
    linearization: &[MonomerIdx],
    rban_record: &Parsed_rBAN_Record,
    out_file: &Path,
) -> Result<()> {
    let emitting_state_indices_on_path: Vec<StateIdx> = {
	hmm_path.iter()
	.filter(|&&state_idx| hmm.is_emitting(state_idx))
	.cloned()
	.collect()
    };

    if emitting_state_indices_on_path.len() != linearization.len() {
	panic!("Number of emitting states on HMM path must match number of monomer names");
    }

    let monomer_names: Vec<String> = {
	linearization.iter()
	    .map(|&monomer_idx| {
		rban_record
		    .monomers
		    .get(&monomer_idx)
		    .expect("Monomer index in linearization must be valid")
		    .name
		    .0
		    .clone()
	    })
	    .collect()
    };

    let monomer_names_with_idx: Vec<String> = {
	monomer_names.iter()
	.zip(linearization.iter())
	.map(|(name, &monomer_idx)| format!("{}.{}", name, monomer_idx))
	.collect()
    };

    let labels: HashMap<StateIdx, String> = {
	let mut labels: HashMap<StateIdx, String> = HashMap::new();

	for (&state_idx, monomer_name) in (emitting_state_indices_on_path.iter()
					  .zip(monomer_names_with_idx.iter())) {
	    labels
		.entry(state_idx)
		.and_modify(|s| {
		    s.push(',');
		    s.push_str(monomer_name);
		})
    .or_insert_with(|| monomer_name.clone());
	}
	labels
    };

    let out_dir = {
	out_file
	.parent()
	.ok_or_else(|| anyhow::anyhow!("Output file {:?} has no parent directory", out_file))?
    };
    let hmm_svg_path = out_dir.join("hmm.svg");

    let hmm_with_emission_labels = {
	let mut hmm_clone = hmm.clone();
	let empty_label = "".to_string();
	for state_idx in emitting_state_indices_on_path {
	    let old_label = &hmm.state_labels[state_idx];
	    let suffix = labels.get(&state_idx).unwrap_or(&empty_label);
	    hmm_clone.state_labels[state_idx] = format!("{}*{}", old_label, suffix);
	}
	hmm_clone
    };

    draw_nerpa_hmm(&hmm_with_emission_labels,
		   Some(&hmm_svg_path),
		   Some(hmm_path),
		   true)?;

    let title = monomer_names_with_idx.join("  ->  ");
    let hmm_svg_with_title = {
	let hmm_svg_text = std::fs::read_to_string(&hmm_svg_path)
	    .with_context(|| format!("Failed to read generated HMM SVG from {:?}", hmm_svg_path))?;
	svg_with_title(
	    &hmm_svg_text,
	    &title,
	    TitlePlacement::Bottom
	)?
    };

    std::fs::write(&hmm_svg_path, hmm_svg_with_title)?;

    Ok(())
}
