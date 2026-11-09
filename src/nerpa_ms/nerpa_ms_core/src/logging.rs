use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use log::LevelFilter;

const LOG_FILE_NAME: &str = "nerpa-ms-core.log";
const WARNINGS_FILE_NAME: &str = "nerpa-ms-core-warnings.log";

pub struct LogPaths {
    pub full_log: PathBuf,
    pub warnings_log: PathBuf,
}

pub fn configure(output_dir: &Path, debug_messages: bool) -> Result<LogPaths> {
    std::fs::create_dir_all(output_dir).with_context(|| {
        format!(
            "Failed to create output directory {}",
            output_dir.display()
        )
    })?;

    let enabled_level = if debug_messages {
        LevelFilter::Debug
    } else {
        LevelFilter::Info
    };

    let console_output = fern::Dispatch::new()
        .level(enabled_level)
        .chain(std::io::stdout());

    let log_path = output_dir.join(LOG_FILE_NAME);
    let file_output = fern::Dispatch::new()
        .level(enabled_level)
        .chain(
            fern::log_file(&log_path)
                .with_context(|| format!("Failed to create {}", log_path.display()))?,
        );

    let warnings_path = output_dir.join(WARNINGS_FILE_NAME);
    let warnings_output = fern::Dispatch::new()
        .level(LevelFilter::Warn)
        .format(|output, message, record| {
            output.finish(format_args!(
                "{} [{}] {}",
                chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
                record.level(),
                message
            ))
        })
        .chain(
            fern::log_file(&warnings_path)
                .with_context(|| {
                    format!("Failed to create {}", warnings_path.display())
                })?,
        );

    fern::Dispatch::new()
        .level(enabled_level)
        .chain(console_output)
        .chain(file_output)
        .chain(warnings_output)
        .apply()
        .context("Failed to configure logging")?;

    Ok(LogPaths {
        full_log: log_path,
        warnings_log: warnings_path,
    })
}