use crate::data_types::common_types::MonomerIdx;
use crate::data_types::parsed_rban_record::{AtomId, BondType};

/// Symbolic atom label used inside a bond template, e.g. "C1", "N1".
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BondAtomLabel(pub String);

impl BondAtomLabel {
    pub fn new<S: Into<String>>(s: S) -> Self {
        Self(s.into())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomicBondTemplate {
    pub bond_type: BondType,
    pub arity: String, // "1", "1.5", "2", etc. -- use string to compare fractional arities like "1.5"
    pub atoms: (BondAtomLabel, BondAtomLabel),
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BondTemplate(Vec<AtomicBondTemplate>);

impl BondTemplate {
    // Verify that the input vector is sorted
    // and lexicographically smaller than its flipped version
    pub fn new(v: Vec<AtomicBondTemplate>) -> Self {
        for i in 1..v.len() {
            if v[i - 1] > v[i] {
                panic!("BondTemplate must be constructed from a sorted vector");
            }
        }
        let v_flipped: Vec<AtomicBondTemplate> = {
            let mut _vf: Vec<AtomicBondTemplate> = v
                .iter()
                .map(|b| AtomicBondTemplate {
                    bond_type: b.bond_type.clone(),
                    arity: b.arity.clone(),
                    atoms: (b.atoms.1.clone(), b.atoms.0.clone()),
                })
                .collect();
            _vf.sort();
            _vf
        };

        if v_flipped < v {
            panic!("BondTemplate must be constructed from the lexicographically smaller of the vector and its flipped version");
        }
        BondTemplate(v)
    }
      pub fn as_slice(&self) -> &[AtomicBondTemplate] {
          &self.0
      }
      pub fn into_vec(self) -> Vec<AtomicBondTemplate> {
          self.0
      }
}

use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bond {
    pub bond_templ: BondTemplate,
    pub monomers: (MonomerIdx, MonomerIdx),
    pub label_to_atom: (
        HashMap<BondAtomLabel, AtomId>,
        HashMap<BondAtomLabel, AtomId>,
    ), // for each monomer, map from bond atom labels to actual atom IDs
}

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BondSide {
    Left = 0,
    Right = 1,
}

impl From<BondSide> for u8 {
    fn from(b: BondSide) -> u8 {
        b as u8
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BindingSiteType {
    pub bond_templ: BondTemplate,
    pub side: BondSide,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BindingSitesProfile(Vec<BindingSiteType>);

impl BindingSitesProfile {
	pub fn new(binding_site_types: Vec<BindingSiteType>) -> Self {
	    let mut bs_types: Vec<BindingSiteType> = binding_site_types;
	    bs_types.sort();
	    BindingSitesProfile(bs_types)
	}
}
