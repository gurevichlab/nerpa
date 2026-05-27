use crate::data_types::bonds::BondTemplate;
use crate::data_types::common_types::MonomerIdx;
use crate::data_types::parsed_rban_record::{MonomerEdge, MonomerInfo, Parsed_rBAN_Record};

use super::bonds::{BindingSiteType, BindingSitesProfile};
use super::monomers_db::{MonomersDB, MonomersDB_Entry};
use std::ptr;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum InsertionSite {
    Edge(MonomerIdx, MonomerIdx),
    Leaf(MonomerIdx),
}

#[derive(Debug, Clone, Copy)]
pub enum GraphModification<'a> {
    KeepAsIs {
	monomer_idx: MonomerIdx,
    },
    Substitute {
        monomer_idx: MonomerIdx,
        mon_db_entry: &'a MonomersDB_Entry,
    },
    Remove {
        monomer_idx: MonomerIdx,
    },
    Insert {
        site: InsertionSite,
        mon_db_entry: &'a MonomersDB_Entry,
    },
}

// for faster comparison check that mon_db_entry points to the same entry in memory, instead of comparing the whole monomer info
impl PartialEq for GraphModification<'_> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::KeepAsIs { monomer_idx: a }, Self::KeepAsIs { monomer_idx: b }) => a == b,
            (
                Self::Substitute {
                    monomer_idx: a_idx,
                    mon_db_entry: a_entry,
                },
                Self::Substitute {
                    monomer_idx: b_idx,
                    mon_db_entry: b_entry,
                },
            ) => a_idx == b_idx && ptr::eq(*a_entry, *b_entry),
            (Self::Remove { monomer_idx: a }, Self::Remove { monomer_idx: b }) => a == b,
            (
                Self::Insert {
                    site: a_site,
                    mon_db_entry: a_entry,
                },
                Self::Insert {
                    site: b_site,
                    mon_db_entry: b_entry,
                },
            ) => a_site == b_site && ptr::eq(*a_entry, *b_entry),
            _ => false,
        }
    }
}

impl Eq for GraphModification<'_> {}

impl GraphModification<'_> {
    pub fn to_str_short(&self) -> String {
	match self {
	    GraphModification::KeepAsIs { monomer_idx } => {
		format!("KeepAsIs({})", monomer_idx)
	    },
	    GraphModification::Substitute { monomer_idx, mon_db_entry } => {
		format!("Substitute({}, {})", monomer_idx, mon_db_entry.monomer.features.name.0)
	    },
	    GraphModification::Remove { monomer_idx } => {
		format!("Remove({}),", monomer_idx)
	    },
	    GraphModification::Insert { site, mon_db_entry } => {
		let site_str = match site {
		    InsertionSite::Edge(mon_idx1, mon_idx2) => format!("Edge({}, {})", mon_idx1, mon_idx2),
		    InsertionSite::Leaf(mon_idx) => format!("Leaf({})", mon_idx),
		};
		format!("Insert({}, {})", site_str, mon_db_entry.monomer.features.name.0)
	    },
	}
    }
}

