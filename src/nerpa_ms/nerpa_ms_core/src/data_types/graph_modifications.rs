use crate::data_types::bonds::BondTemplate;
use crate::data_types::common_types::MonomerIdx;
use crate::data_types::parsed_rban_record::{MonomerEdge, MonomerInfo, Parsed_rBAN_Record};

use super::bonds::{BindingSiteType, BindingSitesProfile};
use super::monomers_db::{MonomersDB, MonomersDB_Entry};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InsertionSite {
    Edge(MonomerIdx, MonomerIdx),
    Leaf(MonomerIdx),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

