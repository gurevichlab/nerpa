use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode, Stdio};

use anyhow::{bail, Context, Result};
use clap::Parser;

use nerpa_ms_core::data_types::monomer_graph::MonomerGraph;
use nerpa_ms_core::data_types::parsed_rban_record::Parsed_rBAN_Record;

#[derive(Debug, Parser)]
#[command(about = "Compare regenerated SMILES with rBAN metadata SMILES")]
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
    let compound_id = record.compound_id.clone();

    let metadata_smiles = match record.metadata.smiles.as_deref() {
	    Some(smiles) if !smiles.trim().is_empty() => smiles.to_owned(),
	    _ => {
	        if verbose {
		        println!("SKIP {compound_id}: metadata contains no SMILES");
	        }

	        summary.skipped += 1;
	        return;
	    }
    };

    let monomer_graph = MonomerGraph::from(&record);

    let generated_smiles = match monomer_graph.to_smiles() {
	    Ok(smiles) => smiles,
	    Err(error) => {
	        println!("FAIL {compound_id}: to_smiles failed: {error:#}");
	        summary.failed += 1;
	        return;
	    }
    };

    let normalized_metadata_smiles = match normalize_smiles(&metadata_smiles) {
	    Ok(smiles) => smiles,
	    Err(error) => {
	        println!("FAIL {compound_id}: invalid metadata SMILES");
	        println!("  SMILES: {metadata_smiles}");
	        println!("  Error:  {error:#}");
	        summary.failed += 1;
	        return;
	    }
    };

    let normalized_generated_smiles = match normalize_smiles(&generated_smiles) {
	    Ok(smiles) => smiles,
	    Err(error) => {
	        println!("FAIL {compound_id}: generated invalid SMILES");
	        println!("  SMILES: {generated_smiles}");
	        println!("  Error:  {error:#}");
	        summary.failed += 1;
	        return;
	    }
    };

    if normalized_metadata_smiles == normalized_generated_smiles {
	    if verbose {
	        println!("PASS {compound_id}: {normalized_generated_smiles}");
	    }

	    summary.passed += 1;
	    return;
    }

    println!("FAIL {compound_id}");
    println!("  metadata:             {metadata_smiles}");
    println!("  generated:            {generated_smiles}");
    println!("  normalized metadata:  {normalized_metadata_smiles}");
    println!("  normalized generated: {normalized_generated_smiles}");

    summary.failed += 1;
}

/// Produce canonical SMILES without stereochemical or isotope markings.
///
/// Open Babel's `i` canonical-SMILES output option removes both. The test
/// therefore ignores isotopes in addition to stereochemistry.
fn normalize_smiles(smiles: &str) -> Result<String> {
    let mut open_babel = {
	    Command::new("obabel")
	        .args(["-ismi", "-ocan", "-xi", "-xn"])
	        .stdin(Stdio::piped())
	        .stdout(Stdio::piped())
	        .stderr(Stdio::piped())
	        .spawn()
	        .context("failed to start Open Babel")?
    };

    {
	    let maybe_stdin = open_babel.stdin.take();
	    let mut stdin = maybe_stdin
	        .context("failed to open Open Babel's standard input")?;

	    writeln!(stdin, "{smiles}")
	        .context("failed to send SMILES to Open Babel")?;
    }

    let output = {
	    open_babel
	        .wait_with_output()
	        .context("failed while waiting for Open Babel")?
    };

    if !output.status.success() {
	    let stderr = String::from_utf8_lossy(&output.stderr);
	    bail!("Open Babel failed: {}", stderr.trim());
    }

    let stdout = String::from_utf8(output.stdout)
	    .context("Open Babel returned non-UTF-8 output")?;

    parse_normalized_smiles(&stdout)
}

fn parse_normalized_smiles(output: &str) -> Result<String> {
    let maybe_smiles = {
	    output
	        .lines()
	        .map(str::trim)
	        .find(|line| !line.is_empty())
	        .and_then(|line| line.split_whitespace().next())
    };

    maybe_smiles
	    .map(str::to_owned)
	    .context("Open Babel returned no SMILES")
}

fn print_summary(total_records: usize, summary: &TestSummary) {
    println!();
    println!("Total:   {total_records}");
    println!("Passed:  {}", summary.passed);
    println!("Failed:  {}", summary.failed);
    println!("Skipped: {}", summary.skipped);
}
