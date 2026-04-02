#![allow(non_snake_case, non_camel_case_types, non_upper_case_globals)]

use nerpa_ms_core::data_types::parsed_rban_record::{
    MonomerInfo, NerpaCoreResidue, NorineMonomerName, Parsed_rBAN_Record,
};
use std::cmp::Reverse;
use std::path::PathBuf;

use clap::Parser;

#[derive(Debug, Parser)]
#[command(name = "create_monomers_db")]
pub struct Cli {
    /// Path to a JSON file with a list Parsed_rBAN_Record instances
    #[arg(long)]
    pub parsed_rban_records_json: PathBuf,

    /// Path to the output JSON with the monomers database
    #[arg(long)]
    pub out: PathBuf,
}

use anyhow::Result;

use nerpa_ms_core::data_types::bonds::{BindingSiteType, BindingSitesProfile, Bond};
use serde::Serialize;
use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Mon_DB_Entry {
    pub monomer: MonomerInfo,
    pub bonds_by_bs: Vec<(BindingSiteType, Bond)>,
}
pub type MonomersDB = HashMap<BindingSitesProfile, Vec<Mon_DB_Entry>>;

fn create_monomers_db_unfiltered(
    rban_records: &[Parsed_rBAN_Record],
) -> HashMap<BindingSitesProfile, Vec<Mon_DB_Entry>> {
    let mut entries_by_bs: HashMap<BindingSitesProfile, Vec<Mon_DB_Entry>> = HashMap::new();
    for record in rban_records {
        for (mon_idx, monomer) in &record.monomers {
            let bonds_by_bs = record.bonds_for_monomer(*mon_idx);
            let bs_profile = BindingSitesProfile::new(
                bonds_by_bs
                    .iter()
                    .map(|(bs, _bond)| bs.clone())
                    .collect::<Vec<BindingSiteType>>(),
            );
            let entry = Mon_DB_Entry {
                monomer: monomer.clone(),
                bonds_by_bs,
            };
            entries_by_bs.entry(bs_profile).or_default().push(entry);
        }
    }
    entries_by_bs
}

use itertools::Itertools;

fn filter_db_entries(entries: &[Mon_DB_Entry]) -> Vec<Mon_DB_Entry> {
    // q: 1. Discard all entries with is_pks_hybrid==true
    // 2. Group entries by (monomer.nerpa_core, monomer.methylated) pair
    // 3. For each group, keep only the entry with the most common monomer.name. In case of ties, keep any one with the shortest name.
    let mut groups: HashMap<(NerpaCoreResidue, bool), Vec<&Mon_DB_Entry>> = HashMap::new();
    entries
        .iter()
        .filter(|entry| {
	    !entry.monomer.is_pks_hybrid
		&& !entry.monomer.nerpa_core.is_unknown()
	})
        .for_each(|entry| {
            let key = (entry.monomer.nerpa_core.clone(), entry.monomer.methylated);
            groups.entry(key).or_default().push(entry);
        });

    let mut result: Vec<Mon_DB_Entry> = Vec::new();

    // 2. For each group (by nerpa_core, methylated), pick the entry whose monomer.name is most common,
    //    breaking ties by choosing an entry with the shortest name.
    for (_key, group_entries) in groups {
        // Count name frequencies
	let name_freq = group_entries
	    .iter()
	    .map(|e| &e.monomer.name)
	    .counts();

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

fn create_monomers_db(
    rban_records: &[Parsed_rBAN_Record],
) -> HashMap<BindingSitesProfile, Vec<Mon_DB_Entry>> {
    let unfiltered_db = create_monomers_db_unfiltered(rban_records);
    unfiltered_db
        .into_iter()
        .map(|(bs_profile, entries)| (bs_profile, filter_db_entries(&entries)))
        .collect()
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    // q: load the Parsed_rBAN_Record instances from the input JSON. Parsed_rBAN_Record implements Deserialize
    let parsed_rban_records: Vec<Parsed_rBAN_Record> = {
        let file = std::fs::File::open(&cli.parsed_rban_records_json)?;
        let filename = cli.parsed_rban_records_json.display();
        let err_msg = format!("Failed to open file {}", filename);
        serde_json::from_reader(file).expect(&err_msg)
    };

    let monomers_db = create_monomers_db(&parsed_rban_records);
    let mon_db_str_keys: HashMap<String, Vec<Mon_DB_Entry>> = monomers_db
        .into_iter()
        .map(|(bs_profile, entries)| (bs_profile.to_string_key(), entries))
        .collect();

    let out_file = std::fs::File::create(&cli.out)?;
    serde_json::to_writer_pretty(out_file, &mon_db_str_keys)?;

    Ok(())
}
