use crate::{algo::{apply_modifications::apply_modifications, dp_backtrack::backtrack_solutions}, data_types::{common_types::{LogOdds, MonomerIdx}, dag::{DAG, VertexId}, graph_modifications::GraphModification, hmm::{HMM, StateIdx}, monomer_graph::MonomerGraph, monomers_db::{MonomerOrigin, MonomersDB}, parsed_rban_record::Parsed_rBAN_Record}};

use crate::algo::graph_to_dag::create_dag;
use serde::Serialize;

use crate::algo::dp::compute_dp_table;

use super::dp_backtrack::Solution;

#[derive(Debug, Clone, Serialize)]
pub struct Altered_rBAN_Record {
    pub score: LogOdds,
    pub new_record: Parsed_rBAN_Record,
    pub old_to_new_mon_map: Vec<(Option<MonomerIdx>, Option<MonomerIdx>)>,
    pub monomer_origins: Vec<MonomerOrigin>,
}

pub struct NewVariantWithOptPaths {
    pub new_variant: Altered_rBAN_Record,
    pub linearization: Vec<MonomerIdx>,
    pub hmm_path: Vec<StateIdx>,
    pub dag_path: Vec<VertexId>,
}


pub fn generate_new_variants_with_opt_paths<'mon_db>(
    hmm: &HMM,
    monomer_graph: &MonomerGraph,
    dag: &DAG<'mon_db>,
    max_weight: usize,
    max_variants_per_weight: usize,
    monomers_db: &'mon_db MonomersDB,
) -> Vec<NewVariantWithOptPaths> {
    println!("Computing DP table...");
    let dp_table = compute_dp_table(hmm, &dag, max_weight);
    println!("DP table computed. Retrieving new variants with optimal paths...");
    let mut new_variants_with_opt_paths: Vec<NewVariantWithOptPaths> = Vec::new();

    for weight in 0..=max_weight {
	let max_solutions = if weight > 0 {max_variants_per_weight} else { 1 }; // for weight 0, we only want the original molecule, so we take 1 solution
	let solutions = backtrack_solutions(weight, &dp_table, &dag);

	let mut variants_collected = 0;
	for sol in solutions.into_iter() {
	    let mods = {
		sol
		    .dag_edges.iter()
		    .filter_map(|e| e.modification)
		    .collect::<Vec<_>>()
	    };

	    if let Some(new_variant) = apply_modifications(monomer_graph, &mods, monomers_db) {
		let monomer_origins: Vec<MonomerOrigin> = {
		    let mon_db_entries = mods
			.iter()
			.filter_map(|m| {
			    match m {
			        GraphModification::Insert { site: _, mon_db_entry } => Some(mon_db_entry),
			        GraphModification::Substitute { monomer_idx: _, mon_db_entry } => Some(mon_db_entry),
			        GraphModification::Remove { monomer_idx: _ } => None,
			        GraphModification::KeepAsIs { monomer_idx: _ } => None,
			    }
			})
			.collect::<Vec<_>>();

		    mon_db_entries.iter()
			.map(|entry| entry.monomer_origin.clone())
			.collect()
		};
			
		let variant = Altered_rBAN_Record {
		    score: sol.dlo.to_logodds(),
		    new_record: Parsed_rBAN_Record::from(&new_variant.new_monomer_graph),
		    old_to_new_mon_map: new_variant.old_to_new_mon_map.clone(),
		    monomer_origins,
		};

		new_variants_with_opt_paths.push(NewVariantWithOptPaths {
		    new_variant: variant,
		    linearization: new_variant.linearization.clone(),
		    hmm_path: sol.states.clone(),
		    dag_path: {
			let mut path = vec![dag.start];
			for edge in &sol.dag_edges {
			    path.push(edge.to);
			}
			path
		    }
		});
		variants_collected += 1;
	    }

	    if variants_collected >= max_solutions {
		break;
	    }
	}

    }

    new_variants_with_opt_paths
}

