use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode, Stdio};

use anyhow::{bail, Context, Result};
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
    check_open_babel()?;

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

fn check_open_babel() -> Result<()> {
    let status = {
	    Command::new("obabel")
	        .arg("-V")
	        .stdout(Stdio::null())
	        .stderr(Stdio::null())
	        .status()
    };

    match status {
	    Ok(status) if status.success() => Ok(()),
	    _ => bail!(
	        "Open Babel is unavailable; install it with \
	         `sudo apt install openbabel`"
	    ),
    }
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
    let monomer_graph = MonomerGraph::from(&record);
    let canonized_graph = {
        let mut canonized = monomer_graph.clone();
        let _ = canonized.canonize();
        canonized
    };


    let original_smiles = match monomer_graph.to_smiles() {
	    Ok(smiles) => smiles,
	    Err(error) => {
	        println!("FAIL {compound_id}: original to_smiles failed: {error:#}");
	        summary.failed += 1;
	        return;
	    }
    };

    let canonized_smiles = match canonized_graph.to_smiles() {
	    Ok(smiles) => smiles,
	    Err(error) => {
	        println!("FAIL {compound_id}: original to_smiles failed: {error:#}");
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

fn print_summary(total_records: usize, summary: &TestSummary) {
    println!();
    println!("Total:   {total_records}");
    println!("Passed:  {}", summary.passed);
    println!("Failed:  {}", summary.failed);
    println!("Skipped: {}", summary.skipped);
}
