use crate::data_types::bonds::{BindingSiteType, BindingSitesProfile, Bond};
use itertools::Itertools;
use crate::data_types::parsed_rban_record::{
    MonomerInfo, NerpaCoreResidue, NorineMonomerName, Parsed_rBAN_Record,
};
use serde::{Deserialize, Serialize};
use std::{cmp::Reverse, collections::HashMap, path::{Path, PathBuf}};

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
pub struct MonomersDB_Entry {
    pub monomer: MonomerInfo,
    pub bonds_by_bs: Vec<(BindingSiteType, Bond)>,
}

pub type MonomersDB = HashMap<BindingSitesProfile, Vec<MonomersDB_Entry>>;

pub fn create_monomers_db_unfiltered(rban_records: &[Parsed_rBAN_Record]) -> MonomersDB {
    let mut entries_by_bs: MonomersDB = HashMap::new();
    for record in rban_records {
        for (mon_idx, monomer) in &record.monomers {
            let bonds_by_bs = record.bonds_for_monomer(*mon_idx);
            let bs_profile = BindingSitesProfile::new(
                bonds_by_bs
                    .iter()
                    .map(|(bs, _bond)| bs.clone())
                    .collect::<Vec<BindingSiteType>>(),
            );
            let entry = MonomersDB_Entry {
                monomer: monomer.clone(),
                bonds_by_bs,
            };
            entries_by_bs.entry(bs_profile).or_default().push(entry);
        }
    }
    entries_by_bs
}

fn filter_db_entries(entries: &[MonomersDB_Entry]) -> Vec<MonomersDB_Entry> {
    // q: 1. Discard all entries with is_pks_hybrid==true
    // 2. Group entries by (monomer.nerpa_core, monomer.methylated) pair
    // 3. For each group, keep only the entry with the most common monomer.name. In case of ties, keep any one with the shortest name.
    let mut groups: HashMap<(NerpaCoreResidue, bool), Vec<&MonomersDB_Entry>> = HashMap::new();
    entries
        .iter()
        .filter(|entry| !entry.monomer.is_pks_hybrid && !entry.monomer.nerpa_core.is_unknown())
        .for_each(|entry| {
            let key = (entry.monomer.nerpa_core.clone(), entry.monomer.methylated);
            groups.entry(key).or_default().push(entry);
        });

    let mut result: Vec<MonomersDB_Entry> = Vec::new();

    // 2. For each group (by nerpa_core, methylated), pick the entry whose monomer.name is most common,
    //    breaking ties by choosing an entry with the shortest name.
    for (_key, group_entries) in groups {
        // Count name frequencies
        let name_freq = group_entries.iter().map(|e| &e.monomer.name).counts();

        let group_repr = group_entries
            .iter()
            .sorted_by_key(|e| {
                let name = &e.monomer.name;
                (Reverse(name_freq.get(name).unwrap()), name.0.len())
            })
            .next()
            .unwrap();

        result.push((*group_repr).clone());
    }

    result
}

pub fn create_monomers_db(
    rban_records: &[Parsed_rBAN_Record],
) -> HashMap<BindingSitesProfile, Vec<MonomersDB_Entry>> {
    let unfiltered_db = create_monomers_db_unfiltered(rban_records);
    unfiltered_db
        .into_iter()
        .map(|(bs_profile, entries)| (bs_profile, filter_db_entries(&entries)))
        .collect()
}

pub fn load_monomers_db(monomers_db_json: &Path) -> MonomersDB {
    let json_str = std::fs::read_to_string(monomers_db_json).unwrap_or_else(|e| {
        panic!(
            "Failed to read rBAN records JSON file {:?}: {}",
            monomers_db_json, e
        )
    });
    serde_json::from_str(&json_str).unwrap()
}
