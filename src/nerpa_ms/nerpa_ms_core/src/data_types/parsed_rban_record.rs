use std::{collections::HashMap, sync::LazyLock};
use std::hash::Hash;
use serde::{Deserialize};

use crate::data_types::common_types::{MonomerIdx, MonomerCode};
use serde::Serialize;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub struct AtomId(pub u32);

pub type AtomicEdge = (AtomId, AtomId);
pub type MonomerEdge = (MonomerIdx, MonomerIdx);

#[derive(Debug, Clone, PartialEq, Eq, Hash, Deserialize, Serialize)]
pub struct NorineMonomerName(pub String);


#[derive(Debug, Clone, PartialEq, Eq, Hash, Deserialize, Serialize)]
pub struct NerpaCoreResidue(pub String);

impl NerpaCoreResidue {
    pub fn is_unknown(&self) -> bool {
	self.0 == "_UNKNOWN"
    }
}

pub static UNKNOWN_RESIDUE: LazyLock<NerpaCoreResidue> = LazyLock::new(|| {
    NerpaCoreResidue("_UNKNOWN".to_string())
});


#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Deserialize, Serialize)]
pub enum Chirality {
    D,
    L,
    UNKNOWN,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
pub struct AtomInfo {
    pub name: String, // e.g. 'C', 'N', 'O', 'Cl'
    pub hydrogens: u32,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Deserialize, Serialize)]
pub struct BondType(pub Option<String>);

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
pub struct AtomicEdgeInfo {
    #[serde(deserialize_with = "crate::data_types::json_helpers::de_str_or_num_to_str")]
    pub arity: String, // usually 1.0, but can be fractional for aromatic bonds
    pub bond_type: BondType,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
pub struct MonomerInfo {
    pub name: NorineMonomerName,
    pub nerpa_core: NerpaCoreResidue,
    pub chirality: Chirality,
    pub methylated: bool,
    pub is_pks_hybrid: bool,
    pub mon_code: MonomerCode,
    pub atoms: Vec<AtomId>,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
pub struct MonomerEdgeInfoSingle {
    pub monomer_to_atom: HashMap<MonomerIdx, AtomId>,
    pub atomic_edge: AtomicEdgeInfo,
}

pub type MonomerEdgeInfo = Vec<MonomerEdgeInfoSingle>;

#[derive(Debug, Clone, PartialEq, Eq, Default, Deserialize, Serialize)]
pub struct NRP_Metadata {
    pub name: Option<String>,
    pub smiles: Option<String>,
    pub origin: Option<String>,
    pub inchikey: Option<String>,
    pub source: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct Parsed_rBAN_Record {
    pub compound_id: String,
    pub monomers: HashMap<MonomerIdx, MonomerInfo>,

    #[serde(
        serialize_with = "crate::data_types::json_helpers::ser_hashmap_as_vec_pairs",
	deserialize_with = "crate::data_types::json_helpers::de_vec_pairs_to_hashmap"
    )]
    pub monomer_bonds: HashMap<MonomerEdge, MonomerEdgeInfo>,

    pub atoms: HashMap<AtomId, AtomInfo>,

    #[serde(
        serialize_with = "crate::data_types::json_helpers::ser_hashmap_as_vec_pairs",
	deserialize_with = "crate::data_types::json_helpers::de_vec_pairs_to_hashmap"
    )]
    pub atomic_bonds: HashMap<AtomicEdge, AtomicEdgeInfo>,

    pub metadata: NRP_Metadata,
}

use anyhow::Result;

impl Parsed_rBAN_Record {
    pub fn get_monomers(&self, monomer_indices: &[MonomerIdx]) -> Result<Vec<MonomerInfo>> {
	let mut result = Vec::new();
	for mon_idx in monomer_indices {
	    if let Some(monomer_info) = self.monomers.get(mon_idx) {
		result.push(monomer_info.clone());
	    } else {
		return Err(anyhow::anyhow!("Monomer index {} not found in record", mon_idx));
	    }
	}
	Ok(result)
    }
    
    pub fn show_linearization(&self, linearization: &[MonomerIdx]) -> Result<String> {
	let mut result = String::new();
	for (i, mon_idx) in linearization.iter().enumerate() {
	    if i > 0 {
		result.push_str(", ");
	    }
	    if let Some(monomer_info) = self.monomers.get(mon_idx) {
		result.push_str(
		    &format!("{}({})",
			     monomer_info.name.0,
		    	     mon_idx
		    ));
	    } else {
		return Err(anyhow::anyhow!("Monomer index {} not found in record", mon_idx));
	    }
	}
	Ok(result)
    }
}

