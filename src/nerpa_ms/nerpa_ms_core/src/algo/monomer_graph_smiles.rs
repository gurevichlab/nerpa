use std::collections::HashSet;
use std::fmt::Write as _;
use std::io::Write as _;
use std::process::{Command, Stdio};

use anyhow::{bail, Context, Result};

use crate::data_types::monomer_graph::{AtomData, AtomicEdgeData, MonomerGraph};

impl MonomerGraph {
    /// Regenerate canonical SMILES from the complete atomic graph.
    ///
    /// This does not use `metadata.smiles`.
    ///
    /// Stereochemistry and formal charges are currently NOT represented.
    pub fn to_smiles(&self) -> Result<String> {
	    let cml = self.to_cml()?;
	    let output = run_open_babel(&cml)?;

	    parse_open_babel_smiles(&output)
    }

    fn to_cml(&self) -> Result<String> {
	    let atoms = self.sorted_atoms();

	    if atoms.is_empty() {
	        bail!("cannot generate SMILES for an empty monomer graph");
	    }

	    let known_atom_ids = collect_atom_ids(&atoms)?;
	    let atomic_edges = self.collect_atomic_edges();

	    validate_atomic_edges(&atomic_edges, &known_atom_ids)?;
	    build_cml(&atoms, &atomic_edges)
    }

    fn sorted_atoms(&self) -> Vec<&AtomData> {
	    let mut atoms: Vec<&AtomData> = {
	        self.monomers
		        .values()
		        .flat_map(|monomer| monomer.atoms.iter())
		        .collect()
	    };

	    atoms.sort_by_key(|atom| atom.id.0);
	    atoms
    }

    fn collect_atomic_edges(&self) -> Vec<AtomicEdgeData> {
	    let mut atomic_edges: Vec<AtomicEdgeData> = {
	        self.monomers
		        .values()
		        .flat_map(|monomer| monomer.atomic_bonds.iter().cloned())
		        .collect()
	    };

	    for monomer_bond in &self.monomer_bonds {
	        atomic_edges.extend(monomer_bond.get_atomic_edges());
	    }

	    atomic_edges
    }
}

fn collect_atom_ids(atoms: &[&AtomData]) -> Result<HashSet<u32>> {
    let mut atom_ids = HashSet::new();

    for atom in atoms {
	    if !atom_ids.insert(atom.id.0) {
	        bail!("duplicate atom ID {}", atom.id.0);
	    }
    }

    Ok(atom_ids)
}

fn validate_atomic_edges(
    atomic_edges: &[AtomicEdgeData],
    known_atom_ids: &HashSet<u32>,
) -> Result<()> {
    for atomic_edge in atomic_edges {
	    let left_atom_id = atomic_edge.atom_ids.0.0;
	    let right_atom_id = atomic_edge.atom_ids.1.0;

	    if !known_atom_ids.contains(&left_atom_id) {
	        bail!("bond refers to unknown atom ID {left_atom_id}");
	    }

	    if !known_atom_ids.contains(&right_atom_id) {
	        bail!("bond refers to unknown atom ID {right_atom_id}");
	    }
    }

    Ok(())
}

fn build_cml(
    atoms: &[&AtomData],
    atomic_edges: &[AtomicEdgeData],
) -> Result<String> {
    let mut cml = String::new();

    writeln!(cml, r#"<?xml version="1.0" encoding="UTF-8"?>"#)?;
    writeln!(cml, r#"<cml xmlns="http://www.xml-cml.org/schema">"#)?;
    writeln!(cml, "  <molecule>")?;

    write_cml_atoms(&mut cml, atoms)?;
    write_cml_bonds(&mut cml, atomic_edges)?;

    writeln!(cml, "  </molecule>")?;
    writeln!(cml, "</cml>")?;

    Ok(cml)
}

fn write_cml_atoms(cml: &mut String, atoms: &[&AtomData]) -> Result<()> {
    writeln!(cml, "    <atomArray>")?;

    for atom in atoms {
	    let element_name = escape_xml_attribute(&atom.name);

	    writeln!(
	        cml,
	        r#"      <atom id="a{}" elementType="{}" hydrogenCount="{}"/>"#,
	        atom.id.0,
	        element_name,
	        atom.hydrogens,
	    )?;
    }

    writeln!(cml, "    </atomArray>")?;
    Ok(())
}

fn write_cml_bonds(
    cml: &mut String,
    atomic_edges: &[AtomicEdgeData],
) -> Result<()> {
    writeln!(cml, "    <bondArray>")?;

    for (bond_index, atomic_edge) in atomic_edges.iter().enumerate() {
	    let left_atom_id = atomic_edge.atom_ids.0.0;
	    let right_atom_id = atomic_edge.atom_ids.1.0;
	    let cml_bond_order = cml_bond_order(&atomic_edge.arity)?;

	    writeln!(
	        cml,
	        r#"      <bond id="b{}" atomRefs2="a{} a{}" order="{}"/>"#,
	        bond_index,
	        left_atom_id,
	        right_atom_id,
	        cml_bond_order,
	    )?;
    }

    writeln!(cml, "    </bondArray>")?;
    Ok(())
}

fn cml_bond_order(arity: &str) -> Result<&'static str> {
    match arity.trim() {
	    "1" | "1.0" | "1.00" => Ok("1"),
	    "1.5" | "1.50" => Ok("A"),
	    "2" | "2.0" | "2.00" => Ok("2"),
	    "3" | "3.0" | "3.00" => Ok("3"),
	    unsupported_arity => {
	        bail!("unsupported atomic bond arity {unsupported_arity:?}")
	    }
    }
}

fn run_open_babel(cml: &str) -> Result<String> {
    let mut open_babel = {
	    Command::new("obabel")
	        .args(["-icml", "-ocan"])
	        .stdin(Stdio::piped())
	        .stdout(Stdio::piped())
	        .stderr(Stdio::piped())
	        .spawn()
	        .context(
		        "failed to start Open Babel; install it with \
		         `sudo apt install openbabel`",
	        )?
    };

    {
	    let maybe_stdin = open_babel.stdin.take();
	    let mut stdin = maybe_stdin
	        .context("failed to open Open Babel's standard input")?;

	    stdin
	        .write_all(cml.as_bytes())
	        .context("failed to send the molecular graph to Open Babel")?;
    }

    let output = {
	    open_babel
	        .wait_with_output()
	        .context("failed while waiting for Open Babel")?
    };

    if !output.status.success() {
	    let stderr = String::from_utf8_lossy(&output.stderr);
	    bail!("Open Babel failed: {}", stderr.trim());
    }

    String::from_utf8(output.stdout)
	    .context("Open Babel returned non-UTF-8 output")
}

fn parse_open_babel_smiles(output: &str) -> Result<String> {
    let maybe_smiles = {
	    output
	        .lines()
	        .filter(|line| !line.trim().is_empty())
	        .find_map(|line| line.split_whitespace().next())
    };

    maybe_smiles
	    .map(str::to_owned)
	    .context("Open Babel returned no SMILES")
}

fn escape_xml_attribute(value: &str) -> String {
    value
	    .replace('&', "&amp;")
	    .replace('<', "&lt;")
	    .replace('>', "&gt;")
	    .replace('"', "&quot;")
	    .replace('\'', "&apos;")
}
