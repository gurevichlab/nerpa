use petgraph::algo::is_isomorphic_matching;
use petgraph::graph::UnGraph;
use petgraph::prelude::NodeIndex;
use crate::data_types::parsed_rban_record::{Parsed_rBAN_Record, NorineMonomerName};
use crate::data_types::monomer_graph::MonomerGraph;
use crate::data_types::common_types::MonomerIdx;
use std::collections::{HashMap};
use anyhow::{bail, Context, Result};
use tempfile::NamedTempFile;
use std::io::Write;
use std::path::Path;

impl Parsed_rBAN_Record {
    pub fn is_isomorphic_to(&self, other: &Parsed_rBAN_Record) -> bool {
	if self.monomers.len() != other.monomers.len() {
            return false;
	}

	if self.monomer_bonds.len() != other.monomer_bonds.len() {
            return false;
	}

	let (self_graph, _self_node_map) = self.build_pet_graph();
	let (other_graph, _other_node_map) = other.build_pet_graph();

	is_isomorphic_matching(
            &self_graph,
            &other_graph,
            |left_name, right_name| left_name == right_name, // node match: monomer name
            |_left_edge, _right_edge| true,                  // edge match: trivial
	)
    }

    fn build_pet_graph(&self) -> (UnGraph<NorineMonomerName, ()>, HashMap<MonomerIdx, NodeIndex>) {
	let mut graph: UnGraph<NorineMonomerName, ()> = UnGraph::new_undirected();
	let mut node_indices: HashMap<MonomerIdx, NodeIndex> = HashMap::new();

	for (monomer_idx, monomer_info) in &self.monomers {
            let node_index = graph.add_node(monomer_info.name.clone());
            node_indices.insert(monomer_idx.clone(), node_index);
	}

	for ((left_idx, right_idx), _edge_info) in &self.monomer_bonds {
            let left_node = *node_indices
		.get(left_idx)
		.expect("monomer_bonds references unknown left monomer idx");
            let right_node = *node_indices
		.get(right_idx)
		.expect("monomer_bonds references unknown right monomer idx");

            graph.add_edge(left_node, right_node, ());
	}

	(graph, node_indices)
    }


    pub fn get_canonical_smiles(
	records: &[&Parsed_rBAN_Record],
	nerpa_root: &Path,
    ) -> anyhow::Result<HashMap<String, String>> {
	let records_yaml: String = {
            serde_yaml::to_string(records)
		.context("Failed to serialize records to YAML")?
	};

	let mut input_yaml_file = NamedTempFile::new()
            .context("Failed to create temporary YAML file for records")?;
	input_yaml_file
            .write_all(records_yaml.as_bytes())
            .context("Failed to write records YAML to temporary file")?;
	input_yaml_file
            .flush()
            .context("Failed to flush records YAML temporary file")?;

	let output_yaml_file = NamedTempFile::new()
            .context("Failed to create temporary YAML file for id->smiles output")?;

	let records_to_smiles_script = nerpa_root
            .join("scripts")
            .join("rban_records_to_canon_smiles.py");

	// Script interface:
	//   records_to_smiles.py INPUT.yaml [-o OUTPUT.yaml] [--indent N]
	let mut cmd = std::process::Command::new("python"); // keep this first
	cmd.arg(&records_to_smiles_script)
            .arg(input_yaml_file.path())
            .arg("-o").arg(output_yaml_file.path());

	cmd.current_dir(nerpa_root);
	cmd.env("PYTHONPATH", nerpa_root);

	let output = cmd.output()
            .with_context(|| format!(
		"Failed to execute records_to_smiles script at {}",
		records_to_smiles_script.display()
            ))?;

	let stderr = String::from_utf8_lossy(&output.stderr);
	let stdout = String::from_utf8_lossy(&output.stdout);
	if !output.status.success() {
            bail!(
		"records_to_smiles.py failed (exit={}):\n--- stdout ---\n{}\n--- stderr ---\n{}",
		output.status,
		stdout,
		stderr
            );
	}

	let output_yaml_str = std::fs::read_to_string(output_yaml_file.path())
            .with_context(|| format!(
		"Failed to read records_to_smiles output YAML from {}",
		output_yaml_file.path().display()
            ))?;

	let id_to_smiles: HashMap<String, String> = serde_yaml::from_str(&output_yaml_str)
            .context("Failed to parse records_to_smiles output YAML as HashMap<String, String>")?;

	Ok(id_to_smiles)
    }
}

impl MonomerGraph {
    pub fn is_isomorphic_to(&self, other: &MonomerGraph) -> bool {
	let self_record = Parsed_rBAN_Record::from(self);
	let other_record = Parsed_rBAN_Record::from(other);
	self_record.is_isomorphic_to(&other_record)
    }
}
