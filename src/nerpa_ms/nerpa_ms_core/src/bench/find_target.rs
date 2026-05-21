use std::collections::HashMap;
use std::path::Path;
use crate::io::output::OutputItem;
use crate::data_types::parsed_rban_record::Parsed_rBAN_Record;
use crate::algo::gen_new_variants::Altered_rBAN_Record;

pub enum MatchType {
    SmilesIdentical,
    MonomerGraphIsomorphic,
}

pub struct FindTargetResults {
    // the name of the compound used as the template for NerpaMS
    template: String,

    // the name of the compound we want to find among the candidates
    target: String,

    // edit distance between template and target (for filtering into close siblings)
    distance: u32,

    // rank of the matched variant among all new variants (if any)
    matched_variant_rank: Option<usize>, 

    // score of the matched variant (if any)
    matched_variant_score: Option<f64>, 

    match_type: Option<MatchType>,
}

impl FindTargetResults {
    pub fn to_tsv_row(&self) -> String {
	let rank_str = match self.matched_variant_rank {
	    Some(rank) => rank.to_string(),
	    None => "None".to_string(),
	};
	let score_str = match self.matched_variant_score {
	    Some(score) => score.to_string(),
	    None => "None".to_string(),
	};
	let match_type_str = match self.match_type {
	    Some(MatchType::SmilesIdentical) => "SmilesIdentical",
	    Some(MatchType::MonomerGraphIsomorphic) => "MonomerGraphIsomorphic",
	    None => "None",
	};

	format!("{}\t{}\t{}\t{}\t{}\t{}",
		self.template,
		self.target,
		self.distance,
		rank_str,
		score_str,
		match_type_str
	)
    }

    pub fn to_tsv(results: &[FindTargetResults]) -> String {
	let mut tsv = String::from("template\ttarget\tdistance\tmatched_rank\tmatched_score\texact_match\n");
	for result in results {
	    tsv.push_str(result.to_tsv_row().as_str());
	    tsv.push('\n');
	}
	tsv
    }
}

pub fn find_target(
    nerpa_ms_output: &OutputItem,
    target: &Parsed_rBAN_Record,
    distance: u32,
    nerpa_root: &Path,
) -> FindTargetResults {
    let variants = &nerpa_ms_output.new_variants;
    let template_id: String = {
	variants.iter().next()
	    .expect("No variants found in NerpaMS output")
	    .original_id.clone()
    };
    let target_id: String = target.compound_id.clone();

    let rank: Option<usize>;
    let score: Option<f64>;
    let match_type: Option<MatchType>;

    if let Some(variant) = find_target_smiles_identical(variants, target, nerpa_root) {
	rank = Some(variant.rank);
	score = Some(variant.score);
	match_type = Some(MatchType::SmilesIdentical);
    }
    else if let Some(variant) = find_target_monomer_graph_isomorphic(variants, target) {
	rank = Some(variant.rank);
	score = Some(variant.score);
	match_type = Some(MatchType::MonomerGraphIsomorphic);
    }
    else {
	rank = None;
	score = None;
	match_type = None;
    }

    FindTargetResults {
	template: template_id,
	target: target_id,
	distance,
	matched_variant_rank: rank,
	matched_variant_score: score,
	match_type,
    }
}

fn find_target_smiles_identical<'a>(
    variants: &'a [Altered_rBAN_Record],
    target: &Parsed_rBAN_Record,
    nerpa_root: &Path,
) -> Option<&'a Altered_rBAN_Record> {
    let variants_by_id: HashMap<String, &Altered_rBAN_Record> = {
	variants
	    .iter()
	    .map(|variant| (variant.new_record.compound_id.clone(), variant))
	    .collect()
    };

    let canon_smiles_to_id = {
	let rban_records: Vec<&Parsed_rBAN_Record> =
	    variants_by_id.values()
	    .map(|variant| &variant.new_record)
	    .collect();
	let canon_smiles_by_id = Parsed_rBAN_Record::get_canonical_smiles(&rban_records, nerpa_root)
	    .expect("Failed to compute canon SMILES for records");

	canon_smiles_by_id
	    .iter()
	    .map(|(id, canon_smiles)| (canon_smiles.clone(), id.clone()))
	    .collect::<HashMap<String, String>>()
    };

    let target_canon_smiles = {
	Parsed_rBAN_Record::get_canonical_smiles(&[target], nerpa_root)
	    .expect("Failed to compute canon SMILES for target record")
	    .values()
	    .next()
	    .expect("No canon SMILES found for target record")
	    .clone()
    };

    if let Some(variant_id) = canon_smiles_to_id.get(&target_canon_smiles) {
	Some(variants_by_id[variant_id])
    } else {
	None
    }
}

fn find_target_monomer_graph_isomorphic<'a>(
    variants: &'a [Altered_rBAN_Record],
    target: &Parsed_rBAN_Record,
) -> Option<&'a Altered_rBAN_Record> {
	variants
	    .iter()
	    .find(|variant| {
		variant.new_record.is_isomorphic_to(target)
	    })
}
