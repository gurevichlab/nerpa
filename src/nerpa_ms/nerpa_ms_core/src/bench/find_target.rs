pub struct MatchedVariantData {
    rank: usize,
    score: f64,
	
    // inexact means that monomer graphs match
    // but molecules do not
    exact_match: bool
}

pub struct FindTargetResults {
	// the name of the compound used as the template for NerpaMS
	template: String,

	// the name of the compound we want to find among the candidates
	target: String,

	// edit distance between template and target (for filtering into close siblings)
	distance: u32,

	matched_data: Option<MatchedVariantData>,
}

pub fn find_target(
    nerpa_ms_output: &OutputItem,
    target: &Parsed_rBAN_Record,
    distance: usize,
    nerpa_root: &Path,
) -> FindTargetResults {
    let variants = &nerpa_ms_output.new_variants;
    let matched_data = {
	find_target_smiles_identical(
	    original_record,
	    variants,
	    target,
	    nerpa_root)
	    .or_else(|| find_target_monomer_graph_isomorphic(
		original_record,
		variants,
		target,
		nerpa_root))
    };
    FindTargetResults {
	template: original_record.compound_id.clone(),
	target: target.compound_id.clone(),
	distance,
	matched_data,
    }
}

fn find_target_smiles_identical(
    original_record: &Parsed_rBAN_Record,
    variants: &[Altered_rBAN_Record],
    target: &Parsed_rBAN_Record,
    nerpa_root: &Path,
) -> Option<MatchedVariantData> {
    let randed_variants_by_id: HashMap<String, (usize, &Altered_rBAN_Record)> = {
	variants
	    .iter()
	    .enumerate()
	    .map(|(rank, variant)| (variant.new_record.compound_id.clone(), (rank, variant)))
	    .collect()
    };

    let canon_smiles_to_id = {
	let rban_records: Vec<&Parsed_rBAN_Record> =
	    ranked_variants_by_id.values()
	    .map(|(rank, variant)| &variant.new_record)
	    .collect();
	let canon_smiles_by_id = Parsed_rBAN_Record::get_canon_smiles(rban_records, nerpa_root)
	    .expect("Failed to compute canon SMILES for records");

	canon_smiles_by_id
	    .iter()
	    .map(|(id, canon_smiles)| (canon_smiles.clone(), id.clone()))
	    .collect::<HashMap<String, String>>()
    };

    let target_canon_smiles = {
	Parsed_rBAN_Records::get_canon_smiles([target], nerpa_root)
	    .expect("Target compound ID not found in records")
    };

    if let Some(variant_id) = canon_smiles_to_id.get(&target_canon_smiles) {
	let (rank, variant) = ranked_variants_by_id.get(variant_id)
	    .expect("Variant ID from canon SMILES not found in ranked variants");
	MatchedVariantData {
	    rank: *rank,
	    score: variant.score,
	    exact_match: true,
	}
    } else {
	None
    }
}

fn find_target_monomer_graph_isomorphic(
    original_record: &Parsed_rBAN_Record,
    variants: &[Altered_rBAN_Record],
    target: &Parsed_rBAN_Record,
    nerpa_root: &Path,
) -> Option<MatchedVariantData> {
    if let Some((rank, variant)) = {
	variants
	    .iter()
	    .enumerate()
	    .find(|(rank, variant)| {
		variant.new_record.is_isomorphic_to(target)
	    })
    } {
	MatchedVariantData {
	    rank: *rank,
	    score: variant.score,
	    exact_match: false,
	}
    } else {
	None
    }
}
