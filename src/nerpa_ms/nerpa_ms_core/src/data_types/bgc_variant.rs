use std::collections::HashMap;
use serde::{Deserialize, Serialize};

use crate::data_types::{hmm::BGC_Variant_ID, parsed_rban_record::NerpaCoreResidue};
use super::common_types::LogProb;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct BGC_Module {
    pub gene_id: String,
    pub a_domain_idx: usize,
    pub residue_score: HashMap<NerpaCoreResidue, LogProb>,
    pub modifications: Vec<String>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct BGC_Variant {
    pub bgc_variant_id: BGC_Variant_ID,
    pub modules: Vec<BGC_Module>,
}


