// src/io/input.rs

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use serde::Deserialize;
use serde_json::Value;

use crate::cli::Cli;
use crate::data_types::bgc_variant::BGC_Variant;
use crate::data_types::common_types::{LogOdds, MonomerIdx};
use crate::data_types::hmm::{BGC_Variant_ID, HMM};
use crate::data_types::parsed_rban_record::Parsed_rBAN_Record;

pub struct InputItem {
    pub bgc_variant: BGC_Variant, // for reporting alignments
    pub hmm: HMM,
    pub rban_record: Parsed_rBAN_Record,
    pub linearization: Vec<MonomerIdx>,
    pub score: LogOdds,
}

pub fn parse_nerpa_match_ids(nerpa_match_ids: &str) -> Result<Vec<usize>> {
    /// Parse a string of comma-separated integers and inclusive ranges.
    ///
    /// For example, `"1,3-5,8"` becomes `[1, 3, 4, 5, 8]`.
    let mut parsed_ids = Vec::new();

    for item in nerpa_match_ids.split(',') {
        let item = item.trim();

        if item.is_empty() {
            bail!("nerpa match ID list contains an empty item");
        }

        if let Some((range_start, range_end)) = item.split_once('-') {
            let range_start: usize = range_start.trim().parse().with_context(|| {
                format!("invalid start of nerpa match ID range: {item:?}")
            })?;
            let range_end: usize = range_end.trim().parse().with_context(|| {
                format!("invalid end of nerpa match ID range: {item:?}")
            })?;

            if range_start > range_end {
                bail!(
                    "invalid descending nerpa match ID range: {range_start}-{range_end}"
                );
            }

            parsed_ids.extend(range_start..=range_end);
        } else {
            let match_id: usize = item
                .parse()
                .with_context(|| format!("invalid nerpa match ID: {item:?}"))?;

            parsed_ids.push(match_id);
        }
    }

    Ok(parsed_ids)
}



pub fn get_input(cli: &Cli) -> Result<Vec<InputItem>> {
    if cli.nerpa_results.is_none() {
	    bail!("--nerpa_results is required. --input_json is currently not implemented");
    }
    let nerpa_match_ids: Option<Vec<usize>> = if let Some(ref ids_str) = cli.nerpa_match_ids {
        Some(parse_nerpa_match_ids(ids_str)?)
    } else {
        None
    };

    get_input_from_nerpa_results(
        &(cli.nerpa_results.clone().unwrap()),
        cli.max_nerpa_matches,
        cli.min_nerpa_score,
        nerpa_match_ids.as_deref()
    )
}


#[derive(Debug, Deserialize)]
struct NerpaMatchInfo {
    bgc_variant_id: BGC_Variant_ID,
    nrp_id: String,
    nrp_linearizations: Vec<Vec<MonomerIdx>>,
    // _optimal_paths: Vec<Vec<MonomerIdx>>,
    score: f64,
}

fn input_item_from_nerpa_match(
    match_info: NerpaMatchInfo,
    bgc_variants_by_id: &HashMap<BGC_Variant_ID, BGC_Variant>,
    hmms_by_bgc_variant_id: &HashMap<BGC_Variant_ID, HMM>,
    rban_records_by_compound_id: &HashMap<String, Parsed_rBAN_Record>,
) -> Result<InputItem> {
    let bgc_variant = bgc_variants_by_id
        .get(&match_info.bgc_variant_id)
        .ok_or_else(|| {
            anyhow!(
                "no BGC_Variant found for bgc_variant_id {} in nerpa results",
                match_info.bgc_variant_id
            )
        })?;

    let hmm = hmms_by_bgc_variant_id
        .get(&match_info.bgc_variant_id)
        .ok_or_else(|| {
            anyhow!(
                "no HMM found for bgc_variant_id {} in nerpa results",
                match_info.bgc_variant_id
            )
        })?;

    let rban_record = rban_records_by_compound_id
        .get(&match_info.nrp_id)
        .ok_or_else(|| {
            anyhow!(
                "no Parsed_rBAN_Record found for nrp_id {} in nerpa results",
                match_info.nrp_id
            )
        })?;

    // For simplicity, take the first linearization among the optimal paths.
    let linearization = {
        match_info
            .nrp_linearizations
            .into_iter()
            .next()
            .ok_or_else(|| {
                anyhow!(
                    "no linearizations found for match with bgc_variant_id {} and nrp_id {} in nerpa results",
                    match_info.bgc_variant_id,
                    match_info.nrp_id
                )
            })?
    };

    Ok(InputItem {
        bgc_variant: bgc_variant.clone(),
        hmm: hmm.clone(),
        rban_record: rban_record.clone(),
        linearization,
        score: match_info.score,
    })
}

pub fn get_input_from_nerpa_results(
    nerpa_results_path: &Path,
    max_nerpa_matches: Option<usize>,
    min_nerpa_score: Option<f64>,
    nerpa_match_ids: Option<&[usize]>
) -> Result<Vec<InputItem>> {
    let hmms: Vec<HMM> = {
        let hmms_json = nerpa_results_path
            .join("preprocessed_input")
            .join("hmms_nerpa_ms.json");
        let hmms_text = fs::read_to_string(&hmms_json).with_context(|| {
            format!(
                "failed to read hmms.json from nerpa_results at {}",
                hmms_json.display()
            )
        })?;
        serde_json::from_str(&hmms_text)
            .with_context(|| format!("failed to parse JSON from {}", hmms_json.display()))?
    };

    for hmm in hmms.iter() {
	    hmm.validate()?;
    }

    let hmms_by_bgc_variant_id: HashMap<BGC_Variant_ID, HMM> = {
	    hmms
            .into_iter()
            .map(|hmm| (hmm.bgc_variant_id.clone(), hmm))
            .collect::<std::collections::HashMap<_, _>>()
    };

    let bgc_variants: Vec<BGC_Variant> = {
	    let bgc_variants_yaml = {
	        nerpa_results_path
	            .join("preprocessed_input")
	            .join("BGC_variants.yaml")
	    };
	    let bgc_variants_text = {
	        fs::read_to_string(&bgc_variants_yaml)
	            .with_context(|| { 
		            format!(
		                "failed to read BGC_Variants.yaml from nerpa_results at {}",
		                bgc_variants_yaml.display()
		            )
	            })?
	    };
	    serde_yaml::from_str(&bgc_variants_text)
	        .with_context(|| format!("failed to parse YAML from {}",
				                     bgc_variants_yaml.display()))?
    };
    let bgc_variants_by_id: HashMap<BGC_Variant_ID, BGC_Variant> = {
	    bgc_variants
	        .into_iter()
	        .map(|bgc_variant| (bgc_variant.bgc_variant_id.clone(), bgc_variant))
	        .collect::<std::collections::HashMap<_, _>>()
    };

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

    let selected_matches: Box<dyn Iterator<Item = NerpaMatchInfo> + '_> =
        if let Some(match_ids) = nerpa_match_ids {
            Box::new({
                matches_info
                    .into_iter()
                    .enumerate()
                    .filter_map(|(match_id, match_info)| {
                        match_ids.contains(&match_id).then_some(match_info)
                    })
            })
        } else {
            Box::new({
                matches_info
                    .into_iter()
                    .take(max_nerpa_matches.unwrap_or(usize::MAX))
                    .take_while(move |match_info| {
                        min_nerpa_score
                            .map(|minimum_score| match_info.score >= minimum_score)
                            .unwrap_or(true)
                    })
            })
        };

    selected_matches
        .map(|match_info| {
            input_item_from_nerpa_match(
                match_info,
                &bgc_variants_by_id,
                &hmms_by_bgc_variant_id,
                &rban_records_by_compound_id,
            )
        })
        .collect()
}
