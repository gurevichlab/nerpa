use crate::data_types::common_types::MonomerIdx;
use crate::data_types::parsed_rban_record::{AtomId, BondType};
use serde::{Deserialize, Serialize};

/// Symbolic atom label used inside a bond template, e.g. "C1", "N1".
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct BondAtomLabel(pub String);

impl BondAtomLabel {
    pub fn new<S: Into<String>>(s: S) -> Self {
        Self(s.into())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct AtomicBondTemplate {
    pub bond_type: BondType,
    pub arity: String, // "1", "1.5", "2", etc. -- use string to compare fractional arities like "1.5"
    pub atoms: (BondAtomLabel, BondAtomLabel),
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
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

use super::monomer_graph::AtomicEdgeData;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Bond {
    pub bond_templ: BondTemplate,
    pub monomers: (MonomerIdx, MonomerIdx),
    pub label_to_atom: (
        HashMap<BondAtomLabel, AtomId>,
        HashMap<BondAtomLabel, AtomId>,
    ), // for each monomer, map from bond atom labels to actual atom IDs
}

impl Bond {
    pub fn get_atomic_edges(&self) -> Vec<AtomicEdgeData> {
	self.bond_templ.as_slice().iter()
	    .map(|at_bond_templ| {
		let atom_id_lft = self.label_to_atom.0
		    .get(&at_bond_templ.atoms.0)
		    .unwrap()
		    .clone();
		let atom_id_rgt = self.label_to_atom.1
		    .get(&at_bond_templ.atoms.1)
		    .unwrap()
		    .clone();
		AtomicEdgeData {
		    atom_ids: (atom_id_lft, atom_id_rgt),
		    arity: at_bond_templ.arity.clone(),
		    bond_type: at_bond_templ.bond_type.clone(),
		}
	    })
	    .collect()
    }
}

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Deserialize, Serialize)]
pub enum BondSide {
    Left = 0,
    Right = 1,
}

impl From<BondSide> for u8 {
    fn from(b: BondSide) -> u8 {
        b as u8
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct BindingSiteType {
    pub bond_templ: BondTemplate,
    pub side: BondSide,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub struct BindingSitesProfile(Vec<BindingSiteType>);

impl BindingSitesProfile {
    pub fn new(binding_site_types: Vec<BindingSiteType>) -> Self {
        let mut bs_types: Vec<BindingSiteType> = binding_site_types;
        bs_types.sort();
        BindingSitesProfile(bs_types)
    }

    /// Produce a stable string form suitable for using as a JSON object key.
    /// (It is JSON text representing the inner Vec.)
    pub fn to_string_key(&self) -> String {
        serde_json::to_string(&self.0).unwrap()
    }

    /// Parse a string key produced by `to_string_key()`.
    pub fn from_string_key(s: &str) -> Result<Self, serde_json::Error> {
        let v = serde_json::from_str::<Vec<BindingSiteType>>(s)?;
        Ok(BindingSitesProfile::new(v))
    }
}

// A wrapper around Vec<(BindingSiteType, Bond)> to ensure that the vector is always sorted by BindingSiteType.
pub struct BondsByBSType(Vec<(BindingSiteType, Bond)>);

impl BondsByBSType {
    pub fn new(mut v: Vec<(BindingSiteType, Bond)>) -> Self {
	v.sort_by_key(|(bs_type, _bond)| bs_type.clone());
	BondsByBSType(v)
    }
    pub fn get_profile(&self) -> BindingSitesProfile {
	BindingSitesProfile::new(self.0.iter()
				 .map(|(bs_type, _bond)| bs_type.clone())
				 .collect())
    }
    pub fn compatible_with(&self, other: &BondsByBSType) -> bool {
	self.get_profile() == other.get_profile()
    }
}

// Custom Deserialize:
// - standard form: a JSON array (the usual newtype Vec representation)
// - string form: a JSON string containing that array (useful for JSON map keys)
impl<'de> Deserialize<'de> for BindingSitesProfile {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use serde::de::{Error as _, SeqAccess, Visitor};
        use std::fmt;

        struct BindingSitesProfileVisitor;

        impl<'de> Visitor<'de> for BindingSitesProfileVisitor {
            type Value = BindingSitesProfile;

            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(
                    f,
                    "either a sequence of BindingSiteType or a string containing JSON for that sequence"
                )
            }

            fn visit_seq<A>(self, seq: A) -> Result<Self::Value, A::Error>
            where
                A: SeqAccess<'de>,
            {
                // Deserialize the sequence into Vec<BindingSiteType>
                let v = Vec::<BindingSiteType>::deserialize(
                    serde::de::value::SeqAccessDeserializer::new(seq),
                )?;
                Ok(BindingSitesProfile(v))
            }

            fn visit_str<E>(self, s: &str) -> Result<Self::Value, E>
            where
                E: serde::de::Error,
            {
                BindingSitesProfile::from_string_key(s).map_err(E::custom)
            }

            fn visit_string<E>(self, s: String) -> Result<Self::Value, E>
            where
                E: serde::de::Error,
            {
                self.visit_str(&s)
            }
        }

        deserializer.deserialize_any(BindingSitesProfileVisitor)
    }
}

impl Bond {
    pub fn shift_atom_ids(&mut self,
			  side: BondSide,
			  shift: u32) {
	let label_to_atom = match side {
	    BondSide::Left => &mut self.label_to_atom.0,
	    BondSide::Right => &mut self.label_to_atom.1,
	};

	label_to_atom
	    .iter_mut()
	    .for_each(|(_label, atom_id)| atom_id.0 += shift);
	}
}
