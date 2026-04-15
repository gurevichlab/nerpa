use crate::{algo::dp_backtrack::backtrack_solutions, data_types::{common_types::{LogProb, MonomerIdx}, hmm::HMM, monomers_db::MonomersDB, parsed_rban_record::Parsed_rBAN_Record}};

use crate::algo::graph_to_dag::create_dag;

use crate::algo::dp::compute_dp_table;

pub struct NewVariantWithScore {
	pub score: LogProb,
	pub new_molecule: Parsed_rBAN_Record,
	pub old_to_new_mon_map: Vec<(Option<MonomerIdx>, Option<MonomerIdx>)>,
}


pub fn generate_new_variants_per_weight(
    hmm: &HMM,
    rban_record: &Parsed_rBAN_Record,
    linearization: &Vec<MonomerIdx>,
    monomers_db: &MonomersDB,
    max_weight: usize,
    max_variants_per_weight: usize,
) -> Vec<Vec<NewVariantWithScore>> {
    let dag = create_dag(rban_record, linearization, monomers_db);
    let dp_table = compute_dp_table(hmm, &dag, max_weight);
    let mut variants_per_weight: Vec<Vec<Altered_NRP_Variant>> = vec![Vec::new(); max_weight + 1];

    for weight in 0..=max_weight {
	let max_solutions = if weight > 0 {max_variants_per_weight} else {1}; // for weight 0, we only want the original molecule, so we take 1 solution
	let solutions = backtrack_solutions(weight, &dp_table, &dag);

	for sol in solutions.into_iter().take(max_variants_per_weight) {
	    let mods = sol
		.dag_edges.iter()
		.filter_map(|e| e.modification)
		.collect::<Vec<_>>();
	    let new_variant = apply_modifications(rban_record, mods);
	    let variant = NewVariantWithScore {
		score: sol.dlp.to_logprob(),
		new_variant.new_molecule,
		new_variant.old_to_new_mon_map,
	    };
	    variants_per_weight[weight].push(variant);
	}
    }

    variants_per_weight
}

