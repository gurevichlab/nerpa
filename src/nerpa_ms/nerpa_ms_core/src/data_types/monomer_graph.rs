use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use super::{bonds::{BindingSiteType, Bond, BondSide, BondsByBSType}, common_types::{MonomerCode, MonomerIdx}, monomers_db::MonomersDB_Entry, parsed_rban_record::{AtomId, BondType, Chirality, NRP_Metadata, NerpaCoreResidue, NorineMonomerName}};


#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
pub struct MonomerFeatures {
    pub name: NorineMonomerName,
    pub nerpa_core: NerpaCoreResidue,
    pub chirality: Chirality,
    pub methylated: bool,
    pub is_pks_hybrid: bool,
    pub mon_code: MonomerCode,
}
#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
pub struct AtomData {
    pub id: AtomId,
    pub name: String, // e.g. 'C', 'N', 'O', 'Cl'
    pub hydrogens: u32,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
pub struct AtomicEdgeData {
    pub atom_ids: (AtomId, AtomId),
    pub arity: String, // "1", "1.5", "2", etc. -- use string to compare fractional arities like "1.5"

    pub bond_type: BondType,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
pub struct Monomer {
    pub atoms: Vec<AtomData>,
    pub atomic_bonds: Vec<AtomicEdgeData>,
    pub features: MonomerFeatures,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MonomerGraph {
    pub monomers: HashMap<MonomerIdx, Monomer>,
    pub monomer_bonds: Vec<Bond>,

    pub compound_id: String,
    pub metadata: NRP_Metadata,
}

impl MonomerGraph {
    pub fn bonds_by_bs_type(&self, monomer_idx: MonomerIdx) -> BondsByBSType {
        let bonds_by_bs = self.monomer_bonds
            .iter()
            .filter_map(|bond| {
                let bs_side = {
                    if bond.monomers.0 == monomer_idx {
                        Some(BondSide::Left)
                    } else if bond.monomers.1 == monomer_idx {
                        Some(BondSide::Right)
                    } else {
                        None
                    }
                };

		if let Some(side) = bs_side {
                    Some(
                        (BindingSiteType {
                            bond_templ: bond.bond_templ.clone(),
                            side,
                        },
                        bond.clone())
                    )
                }
		else {
		    None
		}
            })
            .collect();

	BondsByBSType::new(bonds_by_bs)
    }

    pub fn substitute(&mut self, monomer_idx: MonomerIdx, mon_db_entry: &MonomersDB_Entry) {
	let mon_bonds_by_bs = self.bonds_by_bs_type(monomer_idx);
	if !mon_bonds_by_bs.compatible_with(&mon_db_entry.bonds_by_bs) {
	    panic!("Attempting to substitute monomer {} with a monomer that has an incompatible binding sites profile", monomer_idx);
	}

	let new_bonds_by_mon_indices = {
	    let mut new_bonds_by_mons = HashMap::new();
	    for i in 0..mon_bonds_by_bs.len() {
		let (mon_bs, mon_bond) =
		    mon_bonds_by_bs.get(i).unwrap();
		let (db_entry_bs, db_entry_bond) =
		    mon_db_entry.bonds_by_bs.get(i).unwrap();
		debug_assert_eq!(*mon_bs, *db_entry_bs);
		let side = mon_bs.side;

		let mut new_bond = mon_bond.clone();
		match side {
		    BondSide::Left => {
			new_bond.label_to_atom.0 = db_entry_bond.label_to_atom.0.clone();
		    },
		    BondSide::Right => {
			new_bond.label_to_atom.1 = db_entry_bond.label_to_atom.1.clone();
		    }
		}

		new_bonds_by_mons.insert(mon_bond.monomers, new_bond);
	    }
	    new_bonds_by_mons
	};

	for old_bond in self.monomer_bonds.iter_mut() {
	    if let Some(new_bond) = new_bonds_by_mon_indices.get(&old_bond.monomers) {
		*old_bond = new_bond.clone();
	    }
	}

	self.monomers.insert(monomer_idx, mon_db_entry.monomer.clone());

    }
}
