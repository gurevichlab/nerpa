// src/io/input.rs

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use serde::Deserialize;
use serde_json::Value;

use crate::cli::Cli;
use crate::data_types::common_types::MonomerIdx;
use crate::data_types::hmm::{BGC_Variant_ID, HMM};
use crate::data_types::parsed_rban_record::Parsed_rBAN_Record;

pub struct InputItem {
    pub hmm: HMM,
    pub rban_record: Parsed_rBAN_Record,
    pub linearization: Vec<MonomerIdx>,
}

pub fn get_input(cli: &Cli) -> Result<Vec<InputItem>> {
    match (&cli.nerpa_results, &cli.input_json) {
        (Some(nerpa_results_path), None) => get_input_from_nerpa_results(
            nerpa_results_path,
            cli.max_nerpa_matches,
            cli.min_nerpa_score,
        ),
        (None, Some(json_path)) => parse_input_json(json_path),
        _ => bail!("exactly one of --nerpa_results or --input_json must be provided"),
    }
}

pub fn parse_input_json(json_path: &Path) -> Result<Vec<InputItem>> {
    let text = fs::read_to_string(json_path)
        .with_context(|| format!("failed to read input file: {}", json_path.display()))?;

    let root: Value = serde_json::from_str(&text)
        .with_context(|| format!("failed to parse JSON from: {}", json_path.display()))?;

    // Expect a top-level JSON array: [ {hmm, rban_record, linearization}, ... ]
    let items = match root {
        Value::Array(arr) => arr,
        _ => bail!("top-level JSON is not an array"),
    };

    let mut out = Vec::with_capacity(items.len());

    for (i, item) in items.iter().enumerate() {
        let obj = item
            .as_object()
            .ok_or_else(|| anyhow!("item[{i}] must be an object"))?;

        // Be a little tolerant with field naming.
        let hmm_value = obj
            .get("hmm")
            .ok_or_else(|| anyhow!("item[{i}] missing required field \"hmm\""))?;

        let rban_value = obj
            .get("parsed_rban_record")
            .ok_or_else(|| anyhow!("item[{i}] missing required field \"parsed_rban_record\""))?;

        let lin_value = obj
            .get("linearization")
            .ok_or_else(|| anyhow!("item[{i}] missing required field \"linearization\""))?;

        let hmm: HMM = serde_json::from_value(hmm_value.clone())
            .with_context(|| format!("failed to deserialize item[{i}].hmm"))?;

        let rban_record: Parsed_rBAN_Record = serde_json::from_value(rban_value.clone())
            .with_context(|| format!("failed to deserialize item[{i}].rban_record"))?;

        let linearization: Vec<MonomerIdx> = serde_json::from_value(lin_value.clone())
            .with_context(|| format!("failed to deserialize item[{i}].linearization"))?;

        out.push(InputItem {
            hmm,
            rban_record,
            linearization,
        });
    }

    Ok(out)
}

#[derive(Debug, Deserialize)]
struct NerpaMatchInfo {
    bgc_variant_id: BGC_Variant_ID,
    nrp_id: String,
    nrp_linearizations: Vec<Vec<MonomerIdx>>,
    _optimal_paths: Vec<Vec<MonomerIdx>>,
    score: f64,
}

pub fn get_input_from_nerpa_results(
    nerpa_results_path: &Path,
    max_nerpa_matches: Option<usize>,
    min_nerpa_score: Option<f64>,
) -> Result<Vec<InputItem>> {
    let hmms: Vec<HMM> = {
        let hmms_json = nerpa_results_path
            .join("intermediate_files")
            .join("hmms.json");
        let hmms_text = fs::read_to_string(&hmms_json).with_context(|| {
            format!(
                "failed to read hmms.json from nerpa_results at {}",
                hmms_json.display()
            )
        })?;
        serde_json::from_str(&hmms_text)
            .with_context(|| format!("failed to parse JSON from {}", hmms_json.display()))?
    };
    let hmms_by_bgc_variant_id: HashMap<BGC_Variant_ID, HMM> = hmms
        .into_iter()
        .map(|hmm| (hmm.bgc_variant_id.clone(), hmm))
        .collect::<std::collections::HashMap<_, _>>();

    let rban_records: Vec<Parsed_rBAN_Record> = {
        let rban_yaml = nerpa_results_path
            .join("preprocessed_input")
            .join("parsed_rban_records.yaml");
        let rban_text = fs::read_to_string(&rban_yaml).with_context(|| {
            format!(
                "failed to read parsed_rban_records.json from nerpa_results at {}",
                rban_yaml.display()
            )
        })?;
        serde_yaml::from_str(&rban_text)
            .with_context(|| format!("failed to parse YAML from {}", rban_yaml.display()))?
    };
    let rban_records_by_compound_id: HashMap<String, Parsed_rBAN_Record> = rban_records
        .into_iter()
        .map(|rban| (rban.compound_id.clone(), rban))
        .collect::<std::collections::HashMap<_, _>>();

    let matches_info: Vec<NerpaMatchInfo> = {
        let matches_json = nerpa_results_path
            .join("intermediate_files")
            .join("cpp_matching_output.json");
        let matches_text = fs::read_to_string(&matches_json).with_context(|| {
            format!(
                "failed to read cpp_matching_output.json from nerpa_results at {}",
                matches_json.display()
            )
        })?;
        serde_json::from_str(&matches_text)
            .with_context(|| format!("failed to parse JSON from {}", matches_json.display()))?
    };

    matches_info.into_iter()
        .take(if let Some(n) = max_nerpa_matches { n } else { usize::MAX })
	.take_while(|m| if let Some(s) = min_nerpa_score {m.score >= s } else {true})
	.map(|match_info| {
	    let hmm = hmms_by_bgc_variant_id.get(&match_info.bgc_variant_id)
		.ok_or_else(|| anyhow!("no HMM found for bgc_variant_id {} in nerpa results", match_info.bgc_variant_id))?;

	    let rban_record = rban_records_by_compound_id.get(&match_info.nrp_id)
		.ok_or_else(|| anyhow!("no Parsed_rBAN_Record found for nrp_id {} in nerpa results", match_info.nrp_id))?;

	    // For simplicity, we will just take the first linearization among the optimal paths.
	    let linearization = match_info.nrp_linearizations.into_iter().next()
		.ok_or_else(|| anyhow!("no linearizations found for match with bgc_variant_id {} and nrp_id {} in nerpa results", match_info.bgc_variant_id, match_info.nrp_id))?;

	    Ok(InputItem {
		hmm: hmm.clone(),
		rban_record: rban_record.clone(),
		linearization,
	    })
	})
		.collect()
}
