use crate::data_types::{common_types::MonomerIdx, graph_modifications::GraphModification, monomers_db::MonomersDB_Entry, parsed_rban_record::Parsed_rBAN_Record};


#[derive(Debug, Clone)]
pub struct Altered_NRP_Variant {
    pub new_molecule: Parsed_rBAN_Record,
    pub old_to_new_mon_map: Vec<(Option<MonomerIdx>, Option<MonomerIdx>)>,
}

pub fn apply_modifications(rban_record: &Parsed_rBAN_Record,
			   modifications: &[GraphModification]) -> Altered_NRP_Variant {
    let mut new_molecule = rban_record.clone();
    let mut old_to_new_mon_map: Vec<(Option<MonomerIdx>, Option<MonomerIdx>)> = Vec::new();
    for modification in modifications {
	match modification {
	    GraphModification::Substitute { monomer_idx, mon_db_entry } => {
		new_molecule.sub_inplace(monomer_idx, mon_db_entry);
		old_to_new_mon_map.push((Some(*monomer_idx), Some(*monomer_idx)));
	    },
	    _ => unimplemented!("Only substitution modifications are supported for now"),
	}
    }

    Altered_NRP_Variant {
	new_molecule,
	old_to_new_mon_map,
    }
}

impl Parsed_rBAN_Record {
    fn sub_inplace(&mut self,
		   monomer_idx: &MonomerIdx,
		   mon_db_entry_ref: &MonomersDB_Entry) {
	if !self.monomers.contains_key(monomer_idx) {
	    panic!("Monomer index {:?} not found in the molecule", monomer_idx);
	}

	let mut mon_db_entry = mon_db_entry_ref.clone();

	// 1. Shift AtomIds in the mon_db_entry by the current max AtomId in the molecule + 1 (to avoid collisions)
	let max_atom_id = self.atoms.keys()
	    .map(|atom_id| atom_id.0)
	    .max()
	    .expect("Molecule must have at least one atom");
	mon_db_entry.shift_atom_ids(max_atom_id + 1);

	// 2. update self.atoms: remove the old atoms and add the new ones
	for old_atom_id in self.monomers[monomer_idx].atoms {
	    self.atoms.remove(&old_atom_id);
	}
	for new_atom_id in mon_db_entry.monomer.atoms {
	    let atom_info = mon_db_entry
		.atoms_info
		.get(&new_atom_id)
		.expect("MonomersDB entry must have AtomInfo for all its atoms");
	    self.atoms.insert(atom_id, atom_info.clone());
    }
}
