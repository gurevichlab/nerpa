// src/io/input.rs

use std::fs;
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use serde_json::Value;

use crate::data_types::hmm::HMM;
use crate::data_types::parsed_rban_record::Parsed_rBAN_Record;
use crate::data_types::common_types::{MonomerIdx};

pub struct InputItem {
    pub hmm: HMM,
    pub rban_record: Parsed_rBAN_Record,
    pub linearization: Vec<MonomerIdx>,
}

pub fn parse_input(json_path: &Path) -> Result<Vec<InputItem>> {
    let text = fs::read_to_string(json_path)
        .with_context(|| format!("failed to read input file: {}", json_path.display()))?;

    let root: Value = serde_json::from_str(&text)
        .with_context(|| format!("failed to parse JSON from: {}", json_path.display()))?;

	// Expect a top-level JSON array: [ {hmm, rban_record, linearization}, ... ]
    let items = match root {
	Value::Array(arr) => arr,
	other => bail!("top-level JSON is not an array"),
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
            .ok_or_else(|| {
                anyhow!("item[{i}] missing required field \"parsed_rban_record\"")
            })?;

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

