#![allow(non_snake_case, non_camel_case_types, non_upper_case_globals)]

use nerpa_ms_core::data_types::parsed_rban_record::Parsed_rBAN_Record;
use nerpa_ms_core::data_types::monomers_db::{
    MonomersDB_Entry,
    create_monomers_db
};
use std::path::PathBuf;

use clap::Parser;

#[derive(Debug, Parser)]
#[command(name = "create_monomers_db")]
pub struct Cli {
    /// Path to a YAML file with a list Parsed_rBAN_Record instances
    #[arg(long)]
    pub parsed_rban_records: PathBuf,

    /// Path to the output JSON with the monomers database
    #[arg(long)]
    pub out: PathBuf,
}

use anyhow::Result;

use std::collections::HashMap;

fn main() -> Result<()> {
    let cli = Cli::parse();

    let parsed_rban_records: Vec<Parsed_rBAN_Record> = {
        let file = std::fs::File::open(&cli.parsed_rban_records)?;
        let filename = cli.parsed_rban_records.display();
        let err_msg = format!("Failed to open file {}", filename);
        serde_yaml::from_reader(file).expect(&err_msg)
    };

    let monomers_db = create_monomers_db(&parsed_rban_records);
    let mon_db_str_keys: HashMap<String, Vec<MonomersDB_Entry>> = monomers_db
        .into_iter()
        .map(|(bs_profile, entries)| (bs_profile.to_string_key(), entries))
        .collect();

    let out_file = std::fs::File::create(&cli.out)?;
    serde_json::to_writer_pretty(out_file, &mon_db_str_keys)?;

    Ok(())
}
