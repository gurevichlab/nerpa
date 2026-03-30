use std::collections::HashMap;
use std::hash::Hash;
use serde::{Deserialize};

use crate::data_types::common_types::{MonomerIdx, MonomerCode};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomId(pub u32);

pub type AtomicEdge = (AtomId, AtomId);

pub type MonomerEdge = (MonomerIdx, MonomerIdx);

#[derive(Debug, Clone, PartialEq, Eq, Hash, Deserialize)]
pub struct NorineMonomerName(pub String);

#[derive(Debug, Clone, PartialEq, Eq, Hash, Deserialize)]
pub struct NerpaCoreResidue(pub String);

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Deserialize)]
pub enum Chirality {
    D,
    L,
    Unknown,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct AtomInfo {
    pub name: char, // e.g. 'C', 'N', 'O'
    pub hydrogens: u32,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Deserialize)]
pub struct BondType(pub Option<String>);

#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct AtomicEdgeInfo {
    pub arity: f64, // usually 1.0, but can be fractional for aromatic bonds
    pub bond_type: BondType,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct MonomerInfo {
    pub name: NorineMonomerName,
    pub nerpa_core: NerpaCoreResidue,
    pub chirality: Chirality,
    pub methylated: bool,
    pub is_pks_hybrid: bool,
    pub mon_code: MonomerCode,
    pub atoms: Vec<AtomId>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MonomerEdgeInfoSingle {
    pub monomer_to_atom: HashMap<MonomerIdx, AtomId>,
    
    pub arity: String, // "1", "1.5", "2", etc. -- use string to compare fractional arities like "1.5"

    pub bond_type: BondType,
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct MonomerEdgeInfo {
    pub monomer_to_atom: HashMap<MonomerIdx, AtomId>,

    #[serde(deserialize_with = "crate::data_types::json_helpers::de_str_or_num_to_str")]
    pub arity: String, // "1", "1.5", "2", etc. -- use string to compare fractional arities like "1.5"

    pub bond_type: BondType,
    pub all_edges: Vec<MonomerEdgeInfoSingle>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default, Deserialize)]
pub struct NRP_Metadata {
    pub name: Option<String>,
    pub smiles: Option<String>,
    pub origin: Option<String>,
    pub inchikey: Option<String>,
    pub source: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct Parsed_rBAN_Record {
    pub compound_id: String,
    pub monomers: HashMap<MonomerIdx, MonomerInfo>,

    // JSON encodes this as an array of [key, value] pairs, not as an object.
    #[serde(deserialize_with = "crate::data_types::json_helpers::de_vec_pairs_to_hashmap")]
    pub monomer_bonds: HashMap<MonomerEdge, MonomerEdgeInfo>,

    pub atoms: HashMap<AtomId, AtomInfo>,

    // JSON encodes this as an array of [key, value] pairs, not as an object.
    #[serde(deserialize_with = "crate::data_types::json_helpers::de_vec_pairs_to_hashmap")]
    pub atomic_bonds: HashMap<AtomicEdge, AtomicEdgeInfo>,

    pub metadata: NRP_Metadata,
}


