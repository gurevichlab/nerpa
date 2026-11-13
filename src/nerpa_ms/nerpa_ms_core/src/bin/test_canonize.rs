use std::fs::File;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, Result};
use clap::Parser;

use nerpa_ms_core::data_types::monomer_graph::MonomerGraph;
use nerpa_ms_core::data_types::parsed_rban_record::Parsed_rBAN_Record;

#[derive(Debug, Parser)]
#[command(about = "Test MonomerGraph::canonize by comparing smiles of the orignal and canonized molecules")]
struct Cli {
    /// YAML file containing a list of Parsed_rBAN_Record objects.
    #[arg(long)]
    parsed_rban_records: PathBuf,

    /// Print records whose SMILES match.
    #[arg(short, long)]
    verbose: bool,
}

#[derive(Debug, Default)]
struct TestSummary {
    passed: usize,
    failed: usize,
    skipped: usize,
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    match run(&cli) {
	    Ok(true) => ExitCode::SUCCESS,
	    Ok(false) => ExitCode::FAILURE,
	    Err(error) => {
	        eprintln!("Error: {error:#}");
	        ExitCode::FAILURE
	    }
    }
}

fn run(cli: &Cli) -> Result<bool> {
    let records = read_records(&cli.parsed_rban_records)?;
    let total_records = records.len();

    println!("Testing {total_records} records...");

    let mut summary = TestSummary::default();

    for record in records {
	    test_record(record, cli.verbose, &mut summary);
    }

    print_summary(total_records, &summary);

    Ok(summary.failed == 0)
}

fn read_records(input_path: &Path) -> Result<Vec<Parsed_rBAN_Record>> {
    let input_file = File::open(input_path)
	    .with_context(|| format!("failed to open {}", input_path.display()))?;

    serde_yaml::from_reader(input_file)
	    .with_context(|| format!("failed to parse YAML {}", input_path.display()))
}

fn test_record(
    record: Parsed_rBAN_Record,
    verbose: bool,
    summary: &mut TestSummary,
) {
    let compound_id = &record.compound_id;

    if record
        .metadata
        .smiles
        .as_deref()
        .is_some_and(smiles_has_formal_charge)
    {
        if verbose {
            println!(
                "SKIP {compound_id}: formal charge is not stored in the atomic graph"
            );
        }

        summary.skipped += 1;
        return;
    }

    let monomer_graph = MonomerGraph::from(&record);

    let original_smiles = match monomer_graph.to_smiles_from_atomic_graph() {
        Ok(smiles) => smiles,
        Err(error) => {
            println!(
                "FAIL {compound_id}: original SMILES generation failed: {error:#}"
            );
            summary.failed += 1;
            return;
        }
    };

    let mut canonized_graph = monomer_graph.clone();
    if let Err(error) = canonized_graph.canonize() {
        println!("FAIL {compound_id}: canonize failed: {error:#}");
        summary.failed += 1;
        return;
    }

    let canonized_smiles =
        match canonized_graph.to_smiles_from_atomic_graph() {
            Ok(smiles) => smiles,
            Err(error) => {
                println!(
                    "FAIL {compound_id}: canonized SMILES generation failed: \
                     {error:#}"
                );
                summary.failed += 1;
                return;
            }
        };

    if original_smiles == canonized_smiles {
	    if verbose {
	        println!("PASS {compound_id}");
	    }

	    summary.passed += 1;
	    return;
    }

    println!("FAIL {compound_id}");
    println!("  original:             {original_smiles}");
    println!("  canonized:            {canonized_smiles}");

    summary.failed += 1;
}

fn smiles_has_formal_charge(smiles: &str) -> bool {
    let mut inside_brackets = false;

    for character in smiles.chars() {
        match character {
            '[' => inside_brackets = true,
            ']' => inside_brackets = false,
            '+' | '-' if inside_brackets => return true,
            _ => {}
        }
    }

    false
}

fn print_summary(total_records: usize, summary: &TestSummary) {
    println!();
    println!("Total:   {total_records}");
    println!("Passed:  {}", summary.passed);
    println!("Failed:  {}", summary.failed);
    println!("Skipped: {}", summary.skipped);
}
