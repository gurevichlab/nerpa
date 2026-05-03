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

