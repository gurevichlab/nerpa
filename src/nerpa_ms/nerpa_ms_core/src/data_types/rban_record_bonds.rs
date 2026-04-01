use crate::data_types::bonds::{AtomicBondTemplate, Bond, BondSide, BondAtomLabel, BindingSiteType, BondTemplate};
use crate::data_types::common_types::MonomerIdx;
use crate::data_types::parsed_rban_record::{AtomId, MonomerEdgeInfo, Parsed_rBAN_Record};
use std::collections::HashMap;

impl Parsed_rBAN_Record {
    /// Generate symbolic bond atom labels for a given set of atom IDs.
    /// In case of multiple atoms with the same element, the assignment is non-deterministic
    /// (it depends on IDs)
    fn labels_for_atoms(&self, atom_ids: &[AtomId]) -> HashMap<AtomId, BondAtomLabel> {
        let mut labels = HashMap::new();
        let mut counts = HashMap::<&String, usize>::new();

        // Not needed for correctness, but makes the output deterministic and easier to read.
        let sorted_ids = {
            let mut ids = atom_ids.to_vec();
            ids.sort();
            ids
        };

        for atom_id in sorted_ids.iter() {
            let atom_name = &self.atoms.get(atom_id).unwrap().name;
            let count = counts.entry(atom_name).or_insert(0);
            *count += 1;
            let label = BondAtomLabel::new(format!("{}{}", atom_name, count));
            labels.insert(atom_id.clone(), label);
        }
        labels
    }

    /// Convert an (undirected) monomer edge into an internal `Bond`.
    ///
    /// # Panics
    /// Panics if `edge` is not present in `self.monomer_bonds` (in either direction),
    /// or if the stored bond info does not contain atom IDs for both monomers.
    pub fn edge_to_bond(&self, edge: (MonomerIdx, MonomerIdx)) -> Bond {
        let (left_mon, right_mon) = edge;

        let mon_edge_info: &MonomerEdgeInfo = self
            .monomer_bonds
            .get(&edge)
            .or_else(|| self.monomer_bonds.get(&(right_mon, left_mon)))
            .unwrap_or_else(|| {
		panic!(
		    "edge_to_bond: monomer edge {:?} not found in Parsed_rBAN_Record.monomer_bonds (compound_id: {:?})",
		    edge, self.compound_id
		)
            });

        let atomic_bonds = &mon_edge_info.all_edges;
        let left_atoms: Vec<AtomId> = atomic_bonds
            .iter()
            .map(|e| e.monomer_to_atom.get(&left_mon).unwrap())
            .cloned()
            .collect();
        let right_atoms: Vec<AtomId> = atomic_bonds
            .iter()
            .map(|e| e.monomer_to_atom.get(&right_mon).unwrap())
            .cloned()
            .collect();
        let left_labels: HashMap<AtomId, BondAtomLabel> = self.labels_for_atoms(&left_atoms);
        let right_labels: HashMap<AtomId, BondAtomLabel> = self.labels_for_atoms(&right_atoms);

        let atomic_bond_templates_unflipped: Vec<AtomicBondTemplate> = {
	    let mut _atomic_bond_templates: Vec<AtomicBondTemplate> =
		atomic_bonds
		.iter()
		.map(|e| AtomicBondTemplate {
		    bond_type: e.bond_type.clone(),
		    arity: e.arity.clone(),
		    atoms: (
			left_labels
			    .get(e.monomer_to_atom.get(&left_mon).unwrap())
			    .unwrap()
			    .clone(),
			right_labels
			    .get(e.monomer_to_atom.get(&right_mon).unwrap())
			    .unwrap()
			    .clone(),
		    ),
		})
            .collect();
	    _atomic_bond_templates.sort();
	    _atomic_bond_templates
	};

        let atomic_bond_templates_flipped: Vec<AtomicBondTemplate> = {
            let mut _atomic_bond_templates: Vec<AtomicBondTemplate> =
		atomic_bond_templates_unflipped
                .iter()
                .map(|b| AtomicBondTemplate {
                    bond_type: b.bond_type.clone(),
                    arity: b.arity.clone(),
                    atoms: (b.atoms.1.clone(), b.atoms.0.clone()),
                })
                .collect();
	    _atomic_bond_templates.sort();
	    _atomic_bond_templates
	};

        if atomic_bond_templates_unflipped <= atomic_bond_templates_flipped {
            Bond {
                bond_templ: BondTemplate::new(atomic_bond_templates_unflipped),
                monomers: (left_mon, right_mon),
                label_to_atom: (
                    left_labels
                        .into_iter()
                        .map(|(atom_id, label)| (label, atom_id))
                        .collect(),
                    right_labels
                        .into_iter()
                        .map(|(atom_id, label)| (label, atom_id))
                        .collect(),
                ),
            }
        } else {
            Bond {
                bond_templ: BondTemplate::new(atomic_bond_templates_flipped),
                monomers: (right_mon, left_mon),
                label_to_atom: (
                    right_labels
                        .into_iter()
                        .map(|(atom_id, label)| (label, atom_id))
                        .collect(),
                    left_labels
                        .into_iter()
                        .map(|(atom_id, label)| (label, atom_id))
                        .collect(),
                ),
            }
        }
    }

    pub fn bonds_for_monomer(&self, mon_idx: MonomerIdx) -> Vec<(BindingSiteType, Bond)> {
	let bonds: Vec<Bond> = self
	    .monomer_bonds
	    .iter()
	    .filter(|(edge, _)| edge.0 == mon_idx || edge.1 == mon_idx)
	    .map(|(edge, _)| self.edge_to_bond(*edge))
	    .collect();

	let binding_sites: Vec<BindingSiteType> = bonds
	    .iter()
	    .map(|bond| {
		if bond.monomers.0 == mon_idx {
		    BindingSiteType {side: BondSide::Left,
		                     bond_templ: bond.bond_templ.clone()}
				     
		} else if bond.monomers.1 == mon_idx {
		    BindingSiteType {side: BondSide::Right,
		                     bond_templ: bond.bond_templ.clone()}
		} else {
		    unreachable!();
		}
	    })
	    .collect();
			
	binding_sites.into_iter()
	    .zip(bonds.into_iter())
	    .collect()
    }
}
