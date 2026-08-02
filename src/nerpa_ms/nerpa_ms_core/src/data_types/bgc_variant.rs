use crate::data_types::{hmm::BGC_Variant_ID, parsed_rban_record::NerpaCoreResidue};
use super::common_types::LogOdds;

use std::collections::HashMap;
use serde::{Deserialize, Serialize};
use std::fmt;


#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct BGC_Module {
    pub gene_id: String,
    pub a_domain_idx: usize,
    pub residue_score: HashMap<NerpaCoreResidue, LogOdds>,
    pub modifications: Vec<String>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct BGC_Variant {
    pub bgc_variant_id: BGC_Variant_ID,
    pub modules: Vec<BGC_Module>,
}

impl fmt::Display for BGC_Variant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
	let mut s = format!("{}\n", self.bgc_variant_id);
	for module in self.modules.iter() {
	    let module_id = format!("{}:{}", module.gene_id, module.a_domain_idx);
	    let top_residues = {
		let mut residues_and_scores: Vec<(&NerpaCoreResidue, &LogOdds)> = {
		    module.residue_score
			.iter()
			.collect()
		};
		residues_and_scores.sort_by(
		    |a, b| b.1.partial_cmp(a.1).unwrap_or(std::cmp::Ordering::Equal)
		);
		residues_and_scores
		    .into_iter()
		    .take(5)
		    .map(|(res, score)|
			 {
			     let res_name = if res.is_unknown() {
				 "unk".to_string()
			     }
			     else { res.0.clone() };
			     format!("{}({:.1})", res_name, score)
			 })
		    .collect::<Vec<String>>()
		    .join(",")
	    };

	    let modification_str = format!(
		"[{}]",
		module.modifications.join(",")
	    );

	    let module_str = vec![module_id, top_residues, modification_str].join("\t");

	    s.push_str(&module_str);
	    s.push_str("\n");
	}

	write!(f, "{}", s)
    }
}



