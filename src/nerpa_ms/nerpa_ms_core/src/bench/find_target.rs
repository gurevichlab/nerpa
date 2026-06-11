use std::collections::HashMap;
use std::path::Path;
use crate::algo::generic::f64_cmp;
use crate::algo::apply_modifications::apply_modifications;
use crate::algo::compare_linearizations::{LinearizationComparisonMode, LinearizationModification, linearization_edits};
use crate::data_types::common_types::MonomerIdx;
use crate::data_types::graph_modifications::GraphModification;
use crate::data_types::hmm::HMM;
use crate::data_types::monomer_graph::MonomerGraph;
use crate::data_types::monomers_db::MonomersDB;
use crate::data_types::parsed_rban_record::{MonomerInfo, NorineMonomerName, Parsed_rBAN_Record};
use crate::algo::gen_new_variants::Altered_rBAN_Record;
use crate::algo::rban_records_isomorphic::MonomerKey;
use serde::Serialize;
use itertools::Itertools;

fn cartesian_product_of_rows<T: Clone>(
    matrix: &[Vec<T>],
    max_combinations: usize, // optional limit on the number of combinations
) -> Vec<Vec<T>> {
    // If any row is empty, there are no combinations.
    if matrix.iter().any(|row| row.is_empty()) {
        return Vec::new();
    }

    // By convention, the product of zero rows is one empty combination.
    if matrix.is_empty() {
        return vec![vec![]];
    }

    let result = {
    matrix
        .iter()
        .map(|row| row.iter())
        .multi_cartesian_product()
        .take(max_combinations) // limit the number of combinations if needed
        .map(|combo| combo.into_iter().cloned().collect::<Vec<T>>())
        .collect::<Vec<Vec<T>>>()
    };
    if result.len() == max_combinations {
	eprintln!("Warning: reached the maximum number of combinations ({}). Consider increasing the limit if you want to explore more combinations.", max_combinations);
    }

    result
}

#[derive(Debug, Clone, Serialize)]
pub enum MatchLevel {
    LINEARIZATION,
    MONOMER_GRAPH,
    SMILES,
}

#[derive(Debug, Clone, Serialize)]
pub struct TargetSearchFeatures {
    // the name of the compound used as the template for NerpaMS
    template: String,

    // the name of the compound we want to find among the candidates
    target: String,

    // HMM score for the target linearization
    target_score: f64,


    // edit distance between template and target linearizations
    linearization_distance: u32,

    // edits to reconstruct the target involve only supported monomers
    new_monomers_supported: bool,

    // edits to reconstruct the target can be applied to the template monomer graph
    edits_applicable: bool,

    // the target monomer graph can be reconstructed from the correct linearization
    target_mon_graph_reconstructable: bool,

    // the monomers of the reconstructed monomer graph have correct masses
    target_mon_masses_reconstructable: bool,

    // the target SMILES can be reconstructed from the correct monomer graph
    target_smiles_reconstructable: bool,
}

pub struct TargetSearchData <'a>{
    pub hmm: &'a HMM,
    pub template: &'a MonomerGraph,
    pub template_linearization: &'a [MonomerIdx],
    pub target: &'a MonomerGraph,
    pub target_linearization: &'a [MonomerIdx],
    pub monomers_db: &'a MonomersDB,
    pub nerpa_root: &'a Path,
}

impl <'a>TargetSearchData <'a> {
    pub fn linearization_distance(&self) -> u32 {
	let maybe_edits = linearization_edits(
	    self.template,
	    self.template_linearization,
	    self.target,
	    self.target_linearization,
	    self.monomers_db,
	    LinearizationComparisonMode::Free,
	);
	maybe_edits.expect(
	    "In free mode editting one linearization into another \
	     should be always possible"
	)
	.len() as u32
    }

    pub fn new_monomers_supported(&self) -> bool {
	let maybe_edits = linearization_edits(
	    self.template,
	    self.template_linearization,
	    self.target,
	    self.target_linearization,
	    self.monomers_db,
	    LinearizationComparisonMode::NewMonomerNamesExist,
	);
	maybe_edits.is_some()
    }

    pub fn get_applicable_edits(&self) -> Option<Vec<LinearizationModification>> {
	linearization_edits(
	    self.template,
	    self.template_linearization,
	    self.target,
	    self.target_linearization,
	    self.monomers_db,
	    LinearizationComparisonMode::EditsCanBePerformed,
	)
    }

    pub fn reconstructed_graphs(
	&self,
	lin_mods: &[LinearizationModification]
    ) -> Vec<MonomerGraph> {
	// each lin_mod corresponds to a set of possible graph modifications
	let graph_mods: Vec<Vec<GraphModification>> = {
	    lin_mods
		.iter()
		.map(|lin_mod| {
		    LinearizationModification::to_graph_modifications(
			lin_mod,
			self.template,
			self.template_linearization,
			self.monomers_db
		    )
		})
		.collect()
	};

	// combine all possible pre-image graph modifications for each linearization modification
	cartesian_product_of_rows(&graph_mods, 100) // choose a graph modification for each linearization modification
	    .iter()
	    .filter_map(|graph_mods_vec|
		 apply_modifications(
		     self.template,
		     graph_mods_vec,
		     self.monomers_db,
		     false
		 ))
	    .map(|altered_graph| altered_graph.new_monomer_graph)
	    .collect()
    }

    pub fn reconstructed_smiles(
	&self,
	reconstructed_graph: &MonomerGraph,
    ) -> Option<String> {
	let reconstructed_record = Parsed_rBAN_Record::from(reconstructed_graph);
	let maybe_smiles = Parsed_rBAN_Record::get_canonical_smiles(
	    &[&reconstructed_record],
	    self.nerpa_root
	).map_err(|e| eprintln!("Error getting canonical SMILES: {}", e)); 

	maybe_smiles.ok()?
	    .values()
	    .next()
	    .cloned()
    }

    pub fn target_smiles(&self) -> Option<String> {
	let target_record = Parsed_rBAN_Record::from(self.target);
	let maybe_smiles = Parsed_rBAN_Record::get_canonical_smiles(
	    &[&target_record],
	    self.nerpa_root
	).map_err(|e| eprintln!("Error getting canonical SMILES: {}", e)); 

	maybe_smiles.ok()?
	    .values()
	    .next()
	    .cloned()
    }

}

impl TargetSearchFeatures {
    pub fn new(
	data: &TargetSearchData,
	debug_output: bool
    ) -> Self {
	if debug_output {
	    println!("Template: {}", data.template.compound_id);
	    println!("Target: {}", data.target.compound_id);
	}

	if debug_output {
	    println!("Computing target score...");
	}
	let target_score = data.hmm.score(
	    &data.target,
	    data.target_linearization
	);

	if debug_output {
	    println!("Computing linearization distance...");
	}
	let linearization_distance: u32 = data.linearization_distance();

	if debug_output {
	    println!("Checking if new monomers are supported...");
	}
	let new_monomers_supported: bool = data.new_monomers_supported();

	if debug_output {
	    println!("Checking if edits are applicable...");
	}
	let maybe_edits = data.get_applicable_edits();
	let edits_applicable: bool = maybe_edits.is_some();

	if debug_output {
	    println!("Reconstructing monomer graphs from edits...");
	}
	let reconstructed_graphs: Vec<MonomerGraph> = {
	    if let Some(edits) = maybe_edits {
		data.reconstructed_graphs(&edits)
	    }
	    else { vec![] }
	};

	if debug_output {
	    println!("Checking if any reconstructed graph is isomorphic to the target...");
	}
	let maybe_correct_reconstructed_graph = reconstructed_graphs
	    .iter()
	    .find(|reconstructed_graph|
		  reconstructed_graph.is_isomorphic_to(data.target, MonomerKey::Name));
		

	let target_mon_graph_reconstructable: bool = 
	    maybe_correct_reconstructed_graph.is_some();

	if debug_output {
	    println!("Checking if the monomers of the reconstructed graph have correct masses...");
	}
	let target_mon_masses_reconstructable: bool = 
	    maybe_correct_reconstructed_graph 
            .is_some_and(|g| g.is_isomorphic_to(data.target, MonomerKey::Name_Mass));

	if debug_output {
	    println!("Checking if the target SMILES can be reconstructed from the correct monomer graph...");
	}
	let maybe_reconstructed_smiles =
	    maybe_correct_reconstructed_graph
	    .and_then(|g| data.reconstructed_smiles(g));


	let maybe_target_smiles = data.target_smiles();
	let target_smiles_reconstructable: bool =
	    if let (Some(reconstructed_smiles), Some(target_smiles)) = (&maybe_reconstructed_smiles, &maybe_target_smiles) {
		reconstructed_smiles == target_smiles
	    } else {
		false
	    };


	TargetSearchFeatures {
	    template: data.template.compound_id.clone(),
	    target: data.target.compound_id.clone(),
	    target_score,
	    linearization_distance,
	    new_monomers_supported,
	    edits_applicable,
	    target_mon_graph_reconstructable,
	    target_mon_masses_reconstructable,
	    target_smiles_reconstructable,
	}
    }
}
	    

#[derive(Debug, Clone, Serialize)]
pub struct TargetSearchMatch {
    matched_variant_rank: usize, 
    match_level: MatchLevel,
}


impl TargetSearchMatch {
    pub fn find_target(
	target_record: &Parsed_rBAN_Record,
	target_linearization: &[MonomerIdx],
	generated_variants: &[Altered_rBAN_Record],
	nerpa_root: &Path,
    ) -> Option<TargetSearchMatch> {
	let mut matched_variant: Option<&Altered_rBAN_Record> = None;
	let mut match_level: Option<MatchLevel> = None;

	// trying increasing levels of match strictness
	// if a match at a level is not found, skip trying even stricter levels
	if let Some(v1) = Self::find_target_linearization_identical(
	    generated_variants,
	    target_record,
	    target_linearization
	) {
	    matched_variant = Some(v1);
	    match_level = Some(MatchLevel::LINEARIZATION);

	    if let Some(v2) = Self::find_target_monomer_graph_isomorphic(
		generated_variants,
		target_record
	    ) {
		matched_variant = Some(v2);
		match_level = Some(MatchLevel::MONOMER_GRAPH);

		if let Some(v3) = Self::find_target_smiles_identical(
		    generated_variants,
		    target_record,
		    nerpa_root
		) {
		    matched_variant = Some(v3);
		    match_level = Some(MatchLevel::SMILES);
		}
	    }
	}

	matched_variant
	    .map(|variant|
		 TargetSearchMatch {
		     matched_variant_rank: variant.rank,
		     match_level: match_level.expect("If a match is found, match level should be set"),
		 })
    }

    pub fn find_target_smiles_identical<'a>(
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
		variant.new_record.is_isomorphic_to(
		    target, MonomerKey::Name
		)
	    })
    }

    fn find_target_linearization_identical<'a>(
	variants: &'a [Altered_rBAN_Record],
	target: &Parsed_rBAN_Record,
	target_linearization: &[MonomerIdx],
    ) -> Option<&'a Altered_rBAN_Record> {
	let target_mons: Vec<NorineMonomerName> = {
	    target.get_monomers(target_linearization)
		.expect("Invalid target linearization")
		.iter()
		.map(|mon_info| mon_info.name.clone())
		.collect()
	};

	variants
	    .iter()
	    .find(|variant| {
		let variant_names: Vec<NorineMonomerName> = {
		    variant.linearization_monomers()
			.iter()
			.map(|mon_info| mon_info.name.clone())
			.collect()
		};
		variant_names == target_mons
	    })
    }

}


#[derive(Debug, Clone, Serialize)]
pub struct TargetSearchResults {
    features: TargetSearchFeatures,
    found_match: Option<TargetSearchMatch>,

    // minimum HMM score among the generated variants
    variants_min_score: f64,

    // maximum HMM score among the generated variants
    variants_max_score: f64,
}

impl TargetSearchResults {
    pub fn new(
	search_data: &TargetSearchData,
	generated_variants: &[Altered_rBAN_Record],
    ) -> TargetSearchResults {
	// println!("Generated {} variants", generated_variants.len());
	// for variant in generated_variants {
	//     let linearization_names = {
	// 	variant.linearization_monomers()
	// 	    .iter()
	// 	    .map(|mon_info| mon_info.name.clone())
	// 	    .collect::<Vec<NorineMonomerName>>()
	//     };
	//     if linearization_names.last().unwrap().0 != "Arg" {
	// 	println!("{}", linearization_names.iter().map(|name| name.0.clone()).join(", "));
	//     }	    
	// }

	let features = TargetSearchFeatures::new(search_data, false);
	let target_record = Parsed_rBAN_Record::from(search_data.target);
	let found_match = TargetSearchMatch::find_target(
	    &target_record,
	    search_data.target_linearization,
	    generated_variants,
	    search_data.nerpa_root
	);

	let variant_scores: Vec<f64> = generated_variants
	    .iter()
	    .map(|variant| variant.score)
	    .collect();

	let variants_min_score: f64 = {
	    variant_scores
		.iter()
		.min_by(f64_cmp)
		.cloned()
		.unwrap_or(f64::NAN)
	};
	let variants_max_score: f64 = {
	    variant_scores
		.iter()
		.max_by(f64_cmp)
		.cloned()
		.unwrap_or(f64::NAN)
	};
	
	TargetSearchResults {
	    features,
	    found_match,
	    variants_min_score,
	    variants_max_score,
	}
    }
}
