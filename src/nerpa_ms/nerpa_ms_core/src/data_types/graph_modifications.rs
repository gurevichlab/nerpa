use crate::data_types::bonds::BondTemplate;
use crate::data_types::common_types::MonomerIdx;
use crate::data_types::parsed_rban_record::{MonomerEdge, MonomerInfo, Parsed_rBAN_Record};

use super::bonds::{BindingSiteType, BindingSitesProfile};
use super::monomers_db::{MonomersDB, MonomersDB_Entry};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GraphModification<'a> {
    Insert {
        edge: MonomerEdge,
        mon_db_entry: &'a MonomersDB_Entry,
    },
    Remove {
        monomer_idx: MonomerIdx,
    },
    Substitute {
        monomer_idx: MonomerIdx,
        mon_db_entry: &'a MonomersDB_Entry,
    },
}

pub fn get_possible_subs<'a>(
    rban_record: &Parsed_rBAN_Record,
    monomer_idx: MonomerIdx,
    monomers_db: &'a MonomersDB,
) -> Vec<GraphModification<'a>> {
    let bs_profile: BindingSitesProfile = {
	let bs_with_bonds = rban_record.bonds_for_monomer(monomer_idx);
	BindingSitesProfile::new(
	    bs_with_bonds.iter()
		.map(|(bs, _bond)| bs.clone())
		.collect::<Vec<BindingSiteType>>())
    };

    if let Some(entries) = monomers_db.get(&bs_profile) {
	entries
	.iter()
	.map(|entry| GraphModification::Substitute {
		monomer_idx,
		mon_db_entry: entry,
	})
	.collect()
    }
    else {
	Vec::new()
    }
}
