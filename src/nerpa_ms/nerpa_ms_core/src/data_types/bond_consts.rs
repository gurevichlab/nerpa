use crate::data_types::bonds::{AtomicBondTemplate, BondAtomLabel, BondTemplate};
use crate::data_types::parsed_rban_record::BondType;
use std::sync::LazyLock;

use super::bonds::{BindingSiteType, BindingSitesProfile, BondSide};

pub static AMINO_ATOMIC_BOND: LazyLock<AtomicBondTemplate> = LazyLock::new(|| {
    AtomicBondTemplate {
        bond_type: BondType(Some("AMINO".to_string())),
        arity: "1".to_string(),
        atoms: (BondAtomLabel::new("C1"), BondAtomLabel::new("N1")),
    }
});

pub static ESTER_ATOMIC_BOND: LazyLock<AtomicBondTemplate> = LazyLock::new(|| {
    AtomicBondTemplate {
	bond_type: BondType(Some("ESTER".to_string())),
	arity: "1".to_string(),
	atoms: (BondAtomLabel::new("C1"), BondAtomLabel::new("O1")),
    }
});

pub static ESTER_BOND: LazyLock<BondTemplate> = LazyLock::new(|| {
    BondTemplate::new(
	vec![(*ESTER_ATOMIC_BOND).clone()]
    )
});

pub static ESTER_BINDING_SITE_C: LazyLock<BindingSiteType> = LazyLock::new(|| {
    BindingSiteType {
	bond_templ: (*ESTER_BOND).clone(),
	side: BondSide::Left,
    }
});

pub static ESTER_BINDING_SITE_O: LazyLock<BindingSiteType> = LazyLock::new(|| {
    BindingSiteType {
	bond_templ: (*ESTER_BOND).clone(),
	side: BondSide::Right,
    }
});

pub static AMINO_BOND: LazyLock<BondTemplate> = LazyLock::new(|| {
    BondTemplate::new(
	vec![(*AMINO_ATOMIC_BOND).clone()]
    )
});

pub static AMINO_BINDING_SITE_C: LazyLock<BindingSiteType> = LazyLock::new(|| {
    BindingSiteType {
	bond_templ: (*AMINO_BOND).clone(),
	side: BondSide::Left,
    }
});

pub static AMINO_BINDING_SITE_N: LazyLock<BindingSiteType> = LazyLock::new(|| {
    BindingSiteType {
	bond_templ: (*AMINO_BOND).clone(),
	side: BondSide::Right,
    }
});


pub static AMINO_MIDDLE_PROFILE: LazyLock<BindingSitesProfile> = LazyLock::new(|| {
    BindingSitesProfile::new(vec![
	AMINO_BINDING_SITE_C.clone(),
	AMINO_BINDING_SITE_N.clone(),
    ])
});

pub static AMINO_N_END_PROFILE: LazyLock<BindingSitesProfile> = LazyLock::new(|| {
    BindingSitesProfile::new(vec![
	AMINO_BINDING_SITE_N.clone(),
    ])
});

pub static AMINO_C_END_PROFILE: LazyLock<BindingSitesProfile> = LazyLock::new(|| {
    BindingSitesProfile::new(vec![
	AMINO_BINDING_SITE_C.clone(),
    ])
});

pub static AMINO_N_ESTER_C_PROFILE: LazyLock<BindingSitesProfile> = LazyLock::new(|| {
	BindingSitesProfile::new(vec![
	    AMINO_BINDING_SITE_N.clone(),
	    ESTER_BINDING_SITE_C.clone(),
	])
});

pub static AMINO_CN_ESTER_O_PROFILE: LazyLock<BindingSitesProfile> = LazyLock::new(|| {
	BindingSitesProfile::new(vec![
	    AMINO_BINDING_SITE_C.clone(),
	    AMINO_BINDING_SITE_N.clone(),
	    ESTER_BINDING_SITE_O.clone(),
	])
});

