mod algo;
mod cli;
mod data_types;
mod io;
mod logging;

use std::path::PathBuf;
use std::time::Instant;

use algo::algo_main::process_input_item;
use clap::Parser;
use anyhow::{Context, Result};
use data_types::parsed_rban_record::Parsed_rBAN_Record;
use data_types::monomers_db::load_monomers_db;
use io::draw_molecules::draw_output_variants;
use io::output::{OutputItem, write_output};
use data_types::config::{DebugConfig, NerpaMS_Config};

fn format_command_line() -> String {
    let arguments: Vec<String> = {
        std::env::args_os()
            .map(|argument| argument.to_string_lossy().into_owned())
            .map(|argument| {
                if argument.contains(char::is_whitespace) {
                    format!("{argument:?}")
                } else {
                    argument
                }
            })
            .collect()
    };

    arguments.join(" ")
}

fn main() -> Result<()> {
    let cli = cli::Cli::parse();
    let nerpa_ms_config = NerpaMS_Config::from_cli(&cli);
    let debug_config = DebugConfig::from_cli(&cli);

    let log_paths = logging::configure(
        &cli.out,
        debug_config.debug_messages
    )?;
    let started_at = Instant::now();

    log::info!("Nerpa-MS variant generation has started...");
    log::info!("Started with command: {}", format_command_line());
    log::info!("Nerpa-MS core version: {}", env!("CARGO_PKG_VERSION"));
    log::info!(
        "Logging to {} (all) and {} (warnings and errors)",
        log_paths.full_log.display(),
        log_paths.warnings_log.display()
    );
    log::debug!("Parsed CLI arguments: {cli:#?}");

    let run_result = run(
        &cli,
        &nerpa_ms_config,
        &debug_config
    );

    match &run_result {
        Ok(()) => {
            log::info!("Nerpa-MS variant generation finished successfully.");
        },
        Err(error) => {
            log::error!("Nerpa-MS variant generation failed: {error:#}");
        },
    }

    log::info!("Elapsed time: {:.2?}", started_at.elapsed());

    run_result
}

fn run(
    cli: &cli::Cli,
    nerpa_ms_config: &NerpaMS_Config,
    debug_config: &DebugConfig,
) -> Result<()> {
    let input_items = io::input::get_input(cli)
	.context("Failed to load input")?;
    log::info!("Loaded {} input items", input_items.len());

    log::info!(
        "Loading monomers database from {}",
        cli.monomers_db_json.display()
    );
    let monomers_db = load_monomers_db(&cli.monomers_db_json)?;

    log::info!("Generating variants...");
    let output_items: Vec<OutputItem> = {
	input_items
	    .iter()
	    .map(|item| {
		process_input_item(
		    item,
		    &monomers_db,
		    nerpa_ms_config,
		    debug_config
		)
	    })
	    .collect()
    };

    let output_path = cli.out.join("nerpa-ms-output.json");
    write_output(&output_items, &output_path)?;
    log::info!(
        "Generated candidate structures written to {}",
        output_path.display()
    );

    if debug_config.draw_output_variants {
	let original_records: Vec<&Parsed_rBAN_Record> = {
	    input_items
		.iter()
		.map(|item| &item.rban_record)
		.collect()
	};
	let new_variants_figures_dir: PathBuf = {
	    debug_config.out
		.join("figures")
		.join("new_variants")
	};

	log::info!("Drawing new variants...");
	draw_output_variants(&output_items,
			     &original_records,
			     &cli.nerpa_root,
			     &new_variants_figures_dir)?;
    }

    Ok(())
}
