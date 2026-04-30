use crate::data_types::dag::{VertexId, DAG};
use crate::data_types::hmm::{StateIdx, HMM};
use crate::io::draw_dag::Draw_DAG_Config;
use crate::io::draw_nerpa_hmm::draw_nerpa_hmm;
use crate::io::join_svgs::join_svgs_vertical;
use anyhow::Result;
use std::path::{Path, PathBuf};

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
