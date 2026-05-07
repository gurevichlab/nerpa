use crate::algo::generic::split_when;
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
    nerpa_root: &Path,
) -> Result<()> {
    let out_dir = out_file
        .parent()
        .ok_or_else(|| anyhow::anyhow!("Output file {:?} has no parent directory", out_file))?;
    let hmm_svg_path = out_dir.join("hmm.svg");
    let dag_svg_path = out_dir.join("dag.svg");

    let emission_labels: Vec<String> = {
	dag_path.iter()
	    .filter_map(|&vertex_id| {
		if dag.labels[vertex_id].monomer_code.is_some() {
		    Some(dag.labels[vertex_id].name.clone())
		} else {
		    None
		}
	    })
	    .collect()
    };

    draw_nerpa_hmm(
	hmm,
	hmm_path,
	&emission_labels,
	&hmm_svg_path,
	nerpa_root
    )?;

    dag.draw_svg(&dag_svg_path,
		 &Draw_DAG_Config { node_indexes: true },
		 Some(dag_path))?;
    join_svgs_vertical(&[&hmm_svg_path, &dag_svg_path], out_file)?;
    Ok(())
}
