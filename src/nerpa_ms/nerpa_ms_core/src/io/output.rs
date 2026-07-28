use serde::{Serialize, Deserialize};

use crate::{algo::gen_new_variants::Altered_rBAN_Record, data_types::{common_types::LogOdds, hmm::BGC_Variant_ID, parsed_rban_record::Parsed_rBAN_Record}};
use crate::data_types::json_helpers::serialize_new_variants_as_id_map;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputItem {
    pub bgc_variant_id: BGC_Variant_ID,
    pub compound_id: String,
    pub original_score: LogOdds,
    pub original_record: Parsed_rBAN_Record,

    #[serde(serialize_with = "serialize_new_variants_as_id_map")]
    pub new_variants: Vec<Altered_rBAN_Record>,
}

impl OutputItem {
    pub fn id(&self) -> String {
        format!(
            "BGC_ID-{}___NRP_ID-{}",
            self.bgc_variant_id,
            self.compound_id
        )
    }
}

use std::path::Path;
use std::collections::HashMap;

pub fn write_output(
    output_items: &[OutputItem],
    output_path: &Path
) -> anyhow::Result<()> {
    let output_items_by_id: HashMap<String, OutputItem> = {
        output_items.iter()
            .map(|item| (item.id(), item.clone()))
            .collect()
    };

    let json_output = serde_json::to_string_pretty(&output_items_by_id)?;
    std::fs::create_dir_all(output_path.parent().unwrap())?;
    std::fs::write(output_path, json_output)?;
    Ok(())
}
    
