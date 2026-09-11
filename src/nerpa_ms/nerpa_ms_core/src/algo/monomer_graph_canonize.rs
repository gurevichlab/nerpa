use std::collections::{HashMap, HashSet};
use std::fmt::Write as FmtWrite;
use std::io::Write as IoWrite;
use std::process::{Command, Stdio};

use anyhow::{bail, Context, Result};

use crate::data_types::monomer_graph::{AtomData, AtomicEdgeData, MonomerGraph};
use crate::data_types::parsed_rban_record::AtomId;

impl MonomerGraph {
    /// Renumber atoms according to their order in Open Babel's canonical
    /// SMILES representation.
    ///
    /// New atom IDs start at zero. All references to atom IDs are updated.
    /// Stereochemistry is ignored when calculating the canonical order.
    pub fn canonize(&mut self) -> Result<()> {
	    let atom_ids_in_canonical_order = self.canonical_atom_order()?;
	    let old_to_new_atom_id =
	        build_atom_id_mapping(&atom_ids_in_canonical_order)?;

	    self.apply_atom_id_mapping(&old_to_new_atom_id)
    }

    fn canonical_atom_order(&self) -> Result<Vec<AtomId>> {
	    let atoms = self.sorted_atoms_for_canonization();

	    if atoms.is_empty() {
	        bail!("cannot canonize an empty monomer graph");
	    }

	    let class_to_atom_id: Vec<AtomId> = {
	        atoms
		        .iter()
		        .map(|atom| atom.id.clone())
		        .collect()
	    };

	    let atom_id_to_class = build_atom_class_mapping(&class_to_atom_id)?;
	    let atomic_edges = self.collect_atomic_edges_for_canonization()?;

	    let cml = build_classified_cml(
	        &atoms,
	        &atomic_edges,
	        &atom_id_to_class,
	    )?;

	    let canonical_smiles = run_open_babel_for_atom_order(&cml)?;
	    let atom_classes = parse_atom_classes(&canonical_smiles)?;

	    classes_to_atom_ids(&atom_classes, &class_to_atom_id)
    }

    fn sorted_atoms_for_canonization(&self) -> Vec<&AtomData> {
	    let mut atoms: Vec<&AtomData> = {
	        self.monomers
		        .values()
		        .flat_map(|monomer| monomer.atoms.iter())
		        .collect()
	    };

	    atoms.sort_by_key(|atom| atom.id.0);
	    atoms
    }

    fn collect_atomic_edges_for_canonization(
	    &self,
    ) -> Result<Vec<AtomicEdgeData>> {
	    let mut atomic_edges: Vec<AtomicEdgeData> = {
	        self.monomers
		        .values()
		        .flat_map(|monomer| monomer.atomic_bonds.iter().cloned())
		        .collect()
	    };

	    for monomer_bond in &self.monomer_bonds {
	        for atomic_bond_template in monomer_bond.bond_templ.as_slice() {
		        let maybe_left_atom_id = {
		            monomer_bond.label_to_atom.0
			            .get(&atomic_bond_template.atoms.0)
		        };
		        let left_atom_id = maybe_left_atom_id.with_context(|| {
		            format!(
			            "bond between monomers {} and {} has no atom for label {:?}",
			            monomer_bond.monomers.0,
			            monomer_bond.monomers.1,
			            atomic_bond_template.atoms.0,
		            )
		        })?;

		        let maybe_right_atom_id = {
		            monomer_bond.label_to_atom.1
			            .get(&atomic_bond_template.atoms.1)
		        };
		        let right_atom_id = maybe_right_atom_id.with_context(|| {
		            format!(
			            "bond between monomers {} and {} has no atom for label {:?}",
			            monomer_bond.monomers.0,
			            monomer_bond.monomers.1,
			            atomic_bond_template.atoms.1,
		            )
		        })?;

		        atomic_edges.push(AtomicEdgeData {
		            atom_ids: (
			            left_atom_id.clone(),
			            right_atom_id.clone(),
		            ),
		            arity: atomic_bond_template.arity.clone(),
		            bond_type: atomic_bond_template.bond_type.clone(),
		        });
	        }
	    }

	    Ok(atomic_edges)
    }

    fn apply_atom_id_mapping(
	    &mut self,
	    old_to_new_atom_id: &HashMap<AtomId, AtomId>,
    ) -> Result<()> {
	    for monomer in self.monomers.values_mut() {
	        for atom in &mut monomer.atoms {
		        atom.id = remap_atom_id(&atom.id, old_to_new_atom_id)?;
	        }

	        monomer.atoms.sort_by_key(|atom| atom.id.0);

	        for atomic_bond in &mut monomer.atomic_bonds {
		        atomic_bond.atom_ids.0 = remap_atom_id(
		            &atomic_bond.atom_ids.0,
		            old_to_new_atom_id,
		        )?;
		        atomic_bond.atom_ids.1 = remap_atom_id(
		            &atomic_bond.atom_ids.1,
		            old_to_new_atom_id,
		        )?;
	        }
	    }

	    for monomer_bond in &mut self.monomer_bonds {
	        for label_to_atom in [
		        &mut monomer_bond.label_to_atom.0,
		        &mut monomer_bond.label_to_atom.1,
	        ] {
		        for atom_id in label_to_atom.values_mut() {
		            *atom_id =
			            remap_atom_id(atom_id, old_to_new_atom_id)?;
		        }
	        }
	    }

	    Ok(())
    }
}

fn build_atom_class_mapping(
    class_to_atom_id: &[AtomId],
) -> Result<HashMap<AtomId, u32>> {
    let mut atom_id_to_class = HashMap::new();

    for (index, atom_id) in class_to_atom_id.iter().enumerate() {
	    let atom_class: u32 = {
	        let one_based_index = index
		        .checked_add(1)
		        .context("too many atoms to assign atom classes")?;

	        one_based_index
		        .try_into()
		        .context("too many atoms to assign atom classes")?
	    };

	    if atom_id_to_class
	        .insert(atom_id.clone(), atom_class)
	        .is_some()
	    {
	        bail!("duplicate atom ID {}", atom_id.0);
	    }
    }

    Ok(atom_id_to_class)
}

fn build_atom_id_mapping(
    atom_ids_in_canonical_order: &[AtomId],
) -> Result<HashMap<AtomId, AtomId>> {
    let mut old_to_new_atom_id = HashMap::new();

    for (new_index, old_atom_id) in
	    atom_ids_in_canonical_order.iter().enumerate()
    {
	    let new_atom_id_value: u32 = new_index
	        .try_into()
	        .context("too many atoms to represent with u32 atom IDs")?;

	    let new_atom_id = AtomId(new_atom_id_value);

	    if old_to_new_atom_id
	        .insert(old_atom_id.clone(), new_atom_id)
	        .is_some()
	    {
	        bail!(
		        "canonical atom order contains atom ID {} more than once",
		        old_atom_id.0
	        );
	    }
    }

    Ok(old_to_new_atom_id)
}

fn build_classified_cml(
    atoms: &[&AtomData],
    atomic_edges: &[AtomicEdgeData],
    atom_id_to_class: &HashMap<AtomId, u32>,
) -> Result<String> {
    let mut cml = String::new();

    writeln!(cml, r#"<?xml version="1.0" encoding="UTF-8"?>"#)?;
    writeln!(
	    cml,
	    r#"<cml xmlns="http://www.xml-cml.org/schema">"#
    )?;
    writeln!(cml, "  <molecule>")?;

    write_classified_atoms(&mut cml, atoms, atom_id_to_class)?;
    write_classified_bonds(
	    &mut cml,
	    atomic_edges,
	    atom_id_to_class,
    )?;

    writeln!(cml, "  </molecule>")?;
    writeln!(cml, "</cml>")?;

    Ok(cml)
}

fn write_classified_atoms(
    cml: &mut String,
    atoms: &[&AtomData],
    atom_id_to_class: &HashMap<AtomId, u32>,
) -> Result<()> {
    writeln!(cml, "    <atomArray>")?;

    for atom in atoms {
	    let atom_class = get_atom_class(&atom.id, atom_id_to_class)?;
	    let element_name = escape_xml_attribute(&atom.name);

	    // Open Babel interprets CML IDs of the form aN_M as:
	    //     N = atom identifier
	    //     M = SMILES atom class
	    writeln!(
	        cml,
	        r#"      <atom id="a{0}_{0}" elementType="{1}" hydrogenCount="{2}"/>"#,
	        atom_class,
	        element_name,
	        atom.hydrogens,
	    )?;
    }

    writeln!(cml, "    </atomArray>")?;

    Ok(())
}

fn write_classified_bonds(
    cml: &mut String,
    atomic_edges: &[AtomicEdgeData],
    atom_id_to_class: &HashMap<AtomId, u32>,
) -> Result<()> {
    writeln!(cml, "    <bondArray>")?;

    for (bond_index, atomic_edge) in atomic_edges.iter().enumerate() {
	    let left_atom_class = get_atom_class(
	        &atomic_edge.atom_ids.0,
	        atom_id_to_class,
	    )?;
	    let right_atom_class = get_atom_class(
	        &atomic_edge.atom_ids.1,
	        atom_id_to_class,
	    )?;
	    let bond_order = cml_bond_order(&atomic_edge.arity)?;

	    writeln!(
	        cml,
	        r#"      <bond id="b{0}" atomRefs2="a{1}_{1} a{2}_{2}" order="{3}"/>"#,
	        bond_index,
	        left_atom_class,
	        right_atom_class,
	        bond_order,
	    )?;
    }

    writeln!(cml, "    </bondArray>")?;

    Ok(())
}

fn get_atom_class(
    atom_id: &AtomId,
    atom_id_to_class: &HashMap<AtomId, u32>,
) -> Result<u32> {
    atom_id_to_class
	    .get(atom_id)
	    .copied()
	    .with_context(|| format!("unknown atom ID {}", atom_id.0))
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

fn run_open_babel_for_atom_order(cml: &str) -> Result<String> {
    let mut open_babel = {
	    Command::new("obabel")
	        .args([
		        "-icml",
		        "-ocan",
		        "-xa",
		        "-xi",
		        "-xn",
	        ])
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
	        .context("failed to send CML to Open Babel")?;
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

    let stdout = String::from_utf8(output.stdout)
	    .context("Open Babel returned non-UTF-8 output")?;

    let maybe_smiles = {
	    stdout
	        .lines()
	        .map(str::trim)
	        .find(|line| !line.is_empty())
	        .and_then(|line| line.split_whitespace().next())
    };

    maybe_smiles
	    .map(str::to_owned)
	    .context("Open Babel returned no canonical SMILES")
}

/// Extract atom classes from mapped atoms such as `[C:3]` and `[nH:7]`.
///
/// Because every input atom receives an atom class, their occurrence order
/// is the atom traversal order in the generated canonical SMILES.
fn parse_atom_classes(smiles: &str) -> Result<Vec<u32>> {
    let bytes = smiles.as_bytes();
    let mut atom_classes = Vec::new();
    let mut position = 0;

    while position < bytes.len() {
	    if bytes[position] != b':' {
	        position += 1;
	        continue;
	    }

	    let digits_start = position + 1;
	    let mut digits_end = digits_start;

	    while digits_end < bytes.len()
	        && bytes[digits_end].is_ascii_digit()
	    {
	        digits_end += 1;
	    }

	    let is_atom_class = digits_end > digits_start
	        && digits_end < bytes.len()
	        && bytes[digits_end] == b']';

	    if !is_atom_class {
	        position += 1;
	        continue;
	    }

	    let class_text = &smiles[digits_start..digits_end];
	    let atom_class: u32 = class_text.parse().with_context(|| {
	        format!("invalid atom class {class_text:?} in {smiles:?}")
	    })?;

	    atom_classes.push(atom_class);
	    position = digits_end + 1;
    }

    if atom_classes.is_empty() {
	    bail!(
	        "Open Babel returned canonical SMILES without atom classes: \
	         {smiles:?}"
	    );
    }

    Ok(atom_classes)
}

fn classes_to_atom_ids(
    atom_classes: &[u32],
    class_to_atom_id: &[AtomId],
) -> Result<Vec<AtomId>> {
    if atom_classes.len() != class_to_atom_id.len() {
	    bail!(
	        "canonical SMILES contains {} mapped atoms, but the graph \
	         contains {} atoms",
	        atom_classes.len(),
	        class_to_atom_id.len(),
	    );
    }

    let mut seen_classes = HashSet::new();
    let mut atom_ids = Vec::with_capacity(atom_classes.len());

    for atom_class in atom_classes {
	    if !seen_classes.insert(*atom_class) {
	        bail!(
		        "canonical SMILES contains atom class {} more than once",
		        atom_class
	        );
	    }

	    let zero_based_index = atom_class
	        .checked_sub(1)
	        .context("Open Babel returned atom class zero")?;

	    let class_index: usize = zero_based_index
	        .try_into()
	        .context("atom class does not fit into usize")?;

	    let maybe_atom_id = class_to_atom_id.get(class_index);
	    let atom_id = maybe_atom_id.with_context(|| {
	        format!(
		        "Open Babel returned unknown atom class {}",
		        atom_class
	        )
	    })?;

	    atom_ids.push(atom_id.clone());
    }

    Ok(atom_ids)
}

fn remap_atom_id(
    old_atom_id: &AtomId,
    old_to_new_atom_id: &HashMap<AtomId, AtomId>,
) -> Result<AtomId> {
    old_to_new_atom_id
	    .get(old_atom_id)
	    .cloned()
	    .with_context(|| {
	        format!(
		        "atom ID {} is absent from the canonical order",
		        old_atom_id.0
	        )
	    })
}

fn escape_xml_attribute(value: &str) -> String {
    value
	    .replace('&', "&amp;")
	    .replace('<', "&lt;")
	    .replace('>', "&gt;")
	    .replace('"', "&quot;")
	    .replace('\'', "&apos;")
}
