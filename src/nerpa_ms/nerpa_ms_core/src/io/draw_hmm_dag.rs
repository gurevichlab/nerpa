fn draw_hmm_dag_opt_paths(hmm: &HMM,
			  dag: &DAG<'_>,
			  hmm_path: &[StateIdx],
			  dag_path: &[VertexId],
			  out_file: &Path) -> Result<()> {
    let out_dir = out_file.parent()?;
    let hmm_svg_path = out_dir.join("hmm.svg");
    let dag_svg_path = out_dir.join("dag.svg");

    draw_nerpa_hmm(hmm,
		   Some(&hmm_svg_path),
		   Some(hmm_path),
		   true)?;

    dag.draw_svg(&dag_svg_path,
    		 &Draw_DAG_Config {
		     node_indexes: true,
		 },
		 None)?;
    join_svgs_vertical(&[&hmm_svg_path, &dag_svg_path], out_file)?;
    Ok(())
}
