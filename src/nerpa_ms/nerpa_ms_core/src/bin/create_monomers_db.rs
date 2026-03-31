use nerpa_ms_core::data_types::parsed_rban_record::{MonomerInfo, Parsed_rBAN_Record};
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
use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Mon_DB_Entry {
    pub monomer: MonomerInfo,
    pub bonds_by_bs: Vec<(BindingSiteType, Bond)>,
}

fn collect_tentative_db_entries(
    rban_records: &[Parsed_rBAN_Record],
) -> HashMap<BindingSitesProfile, Vec<Mon_DB_Entry>> {
    let mut entries_by_bs: HashMap<BindingSitesProfile, Vec<Mon_DB_Entry>> = HashMap::new();
    for record in rban_records {
        for (mon_idx, monomer) in &record.monomers {
            let bonds_by_bs = record.bonds_for_monomer(*mon_idx);
            let bs_profile = BindingSitesProfile::new(
                bonds_by_bs
                    .iter()
                    .map(|(bs, bond)| bs.clone())
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

use std::collections::HashSet;

fn filter_db_entries(entries: &[Mon_DB_Entry]) -> Vec<Mon_DB_Entry> {
    // keep only monomers, such that name==nerpa_core an is_pks_hybird==false
    // also, keep no more than one monomer with the same (nerpa_core, methylated) pair
    let mut seen_cores_and_methylation: HashSet<(&String, bool)> = HashSet::new();
    let mut filtered_entries = Vec::new();
    for entry in entries {
        let name = &entry.monomer.name.0;
        let core = &entry.monomer.nerpa_core.0;
        let methylated = entry.monomer.methylated;

        if name == core
            && !entry.monomer.is_pks_hybrid
            && !seen_cores_and_methylation.contains(&(core, methylated))
        {
            filtered_entries.push((*entry).clone());
            seen_cores_and_methylation.insert((core, methylated));
        }
    }
    filtered_entries
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

    Ok(())
}
