use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use super::{bond_consts::{AMINO_BINDING_SITE_C, AMINO_BINDING_SITE_N, AMINO_BOND}, bonds::{BindingSiteType, BindingSitesProfile, Bond, BondSide, BondsByBSType}, common_types::{MonomerCode, MonomerIdx}, monomers_db::{MonomersDB, MonomersDB_Entry}, parsed_rban_record::{AtomId, BondType, Chirality, NRP_Metadata, NerpaCoreResidue, NorineMonomerName}};


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

    pub fn get_bond(&self, mon1: MonomerIdx, mon2: MonomerIdx) -> Option<&Bond> {
	self.monomer_bonds.iter().find(|bond| {
	    (bond.monomers.0 == mon1 && bond.monomers.1 == mon2)
		|| (bond.monomers.0 == mon2 && bond.monomers.1 == mon1)
	})
    }

    pub fn degree(&self, monomer_idx: MonomerIdx) -> usize {
	self.monomer_bonds.iter()
	    .filter(|bond| {
		bond.monomers.0 == monomer_idx
		    || bond.monomers.1 == monomer_idx
	    })
	    .count()
    }

	
}
