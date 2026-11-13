use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{bail, Context, Result};
use serde::Deserialize;
use tempfile::NamedTempFile;

use crate::data_types::parsed_rban_record::{
    AtomId,
    Parsed_rBAN_Record,
};

#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct CanonicalSmilesData {
    pub smiles: String,

    /// Atom IDs in the traversal order used by `smiles`.
    pub atom_ids_in_smiles_order: Vec<AtomId>,
}

pub fn canonical_smiles_data(
    record: &Parsed_rBAN_Record,
) -> Result<CanonicalSmilesData> {
    let repo_root = repo_root();
    let script_path = repo_root
        .join("scripts")
        .join("for_rust_core")
        .join("rban_records_to_molecule_analysis.py");

    let records_yaml = serde_yaml::to_string(&[record])
        .context("failed to serialize rBAN record for RDKit")?;

    let mut input_file = NamedTempFile::new()
        .context("failed to create temporary RDKit input file")?;
    input_file
        .write_all(records_yaml.as_bytes())
        .context("failed to write temporary RDKit input file")?;
    input_file
        .flush()
        .context("failed to flush temporary RDKit input file")?;

    let output_file = NamedTempFile::new()
        .context("failed to create temporary RDKit output file")?;

    let output = Command::new("python")
        .arg(&script_path)
        .arg(input_file.path())
        .arg("--output")
        .arg(output_file.path())
        .current_dir(&repo_root)
        .env("PYTHONPATH", &repo_root)
        .output()
        .with_context(|| {
            format!(
                "failed to execute RDKit bridge at {}",
                script_path.display()
            )
        })?;

    if !output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout);
        let stderr = String::from_utf8_lossy(&output.stderr);

        bail!(
            "RDKit bridge failed with {}:\n--- stdout ---\n{}\n--- stderr ---\n{}",
            output.status,
            stdout.trim(),
            stderr.trim(),
        );
    }

    let output_yaml = std::fs::read_to_string(output_file.path())
        .context("failed to read RDKit bridge output")?;
    let mut analyses: HashMap<String, CanonicalSmilesData> =
        serde_yaml::from_str(&output_yaml)
            .context("failed to parse RDKit bridge output")?;

    analyses.remove(&record.compound_id).with_context(|| {
        format!(
            "RDKit bridge returned no result for compound {}",
            record.compound_id
        )
    })
}

// Temporary: Python helper scripts currently live relative to the source tree.
fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::data_types::parsed_rban_record::{
        AtomInfo,
        AtomicEdgeInfo,
        BondType,
        NRP_Metadata,
    };

    #[test]
    fn returns_canonical_smiles_and_matching_atom_order() -> Result<()> {
        let record = Parsed_rBAN_Record {
            compound_id: "ethanol".to_owned(),
            monomers: HashMap::new(),
            monomer_bonds: HashMap::new(),
            atoms: HashMap::from([
                (
                    AtomId(0),
                    AtomInfo {
                        name: "C".to_owned(),
                        hydrogens: 3,
                    },
                ),
                (
                    AtomId(1),
                    AtomInfo {
                        name: "C".to_owned(),
                        hydrogens: 2,
                    },
                ),
                (
                    AtomId(2),
                    AtomInfo {
                        name: "O".to_owned(),
                        hydrogens: 1,
                    },
                ),
            ]),
            atomic_bonds: HashMap::from([
                (
                    (AtomId(0), AtomId(1)),
                    AtomicEdgeInfo {
                        arity: "1.0".to_owned(),
                        bond_type: BondType(None),
                    },
                ),
                (
                    (AtomId(1), AtomId(2)),
                    AtomicEdgeInfo {
                        arity: "1.0".to_owned(),
                        bond_type: BondType(None),
                    },
                ),
            ]),
            metadata: NRP_Metadata {
                smiles: Some("CCO".to_owned()),
                ..NRP_Metadata::default()
            },
        };

        let data = canonical_smiles_data(&record)?;

        assert_eq!(data.smiles, "CCO");
        assert_eq!(
            data.atom_ids_in_smiles_order,
            vec![AtomId(0), AtomId(1), AtomId(2)],
        );

        Ok(())
    }
}