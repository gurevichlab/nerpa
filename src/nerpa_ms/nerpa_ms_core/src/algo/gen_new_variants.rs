use crate::{algo::{apply_modifications::{AlteredMonomerGraph, apply_modifications}, dp_backtrack::backtrack_solutions}, data_types::{common_types::{LogOdds, MonomerIdx}, config::DebugConfig, mod_graph::{ModGraph, VertexId}, graph_modifications::GraphModification, hmm::{HMM, StateIdx}, monomer_graph::MonomerGraph, monomers_db::{MonomerOrigin, MonomersDB}, parsed_rban_record::{MonomerInfo, Parsed_rBAN_Record}}};

use crate::algo::graph_to_dag::create_mod_graph;
use serde::Serialize;

use crate::algo::dp::compute_dp_table;

use super::dp_backtrack::Solution;
use itertools::Itertools;

#[derive(Debug, Clone, Serialize)]
pub struct Altered_rBAN_Record {
    pub score: LogOdds,
    pub rank: usize,
    pub original_id: String,
    pub new_record: Parsed_rBAN_Record,
    pub linearization: Vec<MonomerIdx>,
    pub old_to_new_mon_map: Vec<(Option<MonomerIdx>, Option<MonomerIdx>)>,
    pub monomer_origins: Vec<MonomerOrigin>,
}

impl Altered_rBAN_Record {
    pub fn linearization_monomers(&self) -> Vec<MonomerInfo> {
	self.new_record
	    .get_monomers(&self.linearization)
	    .expect("Linearization contains invalid monomer indices")
    }
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
    dag: &ModGraph<'mon_db>,
    max_weight: usize,
    max_variants_per_weight: usize,
    monomers_db: &'mon_db MonomersDB,
    debug_stdout: bool,
) -> Vec<NewVariantWithOptPaths> {
    if debug_stdout {
	println!("Computing DP table...");
    }

    let dp_table = compute_dp_table(hmm, &dag, max_weight);
    if debug_stdout {
	println!("DP table computed. Retrieving new variants with optimal paths...");
    }

    let mut new_variants_with_opt_paths: Vec<NewVariantWithOptPaths> = Vec::new();

    for weight in 0..=max_weight {
	let max_solutions = if weight > 0 {max_variants_per_weight} else { 1 }; // for weight 0, we only want the original molecule, so we take 1 solution
	let solutions_with_mods = {
	    backtrack_solutions(weight, &dp_table, dag)
		.map(|sol| {
		    let mods = sol
			.dag_edges
			.iter()
			.filter_map(|e| e.modification.clone())
			.collect::<Vec<_>>();
		    (sol, mods)
		})
	};
	let solutions_with_mods_unique = {
	    solutions_with_mods
		.take(max_solutions * 2)  // to avoid the case where there're no more unique solutions, and solutions_with_mods just keeps yielding the same (differing only by eps cycles in HMM)
		.unique_by(|(_sol, mods)| {
		    mods.iter()
			.map(|m| m.to_str_short())
			.sorted()
			.join(";")
		})
	};
		
	for (solutions_fetched, (sol, mods)) in solutions_with_mods_unique.enumerate() {
	    let new_variant: AlteredMonomerGraph;
	    if let Some(variant) = apply_modifications(
		monomer_graph,
		&mods,
		monomers_db,
		debug_stdout,
	    ) {
		new_variant = variant;
	    }
	    else { continue; } // skip this solution if modifications are inconsistent

	    if debug_stdout {
		println!("Getting monomer origins for new variant...");
	    }
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
	    let rank = new_variants_with_opt_paths.len();
			
	    if debug_stdout {
		println!("Constructing Altered_rBAN_Record");
	    }
	    let variant = Altered_rBAN_Record {
		score: sol.dlo.to_logodds(),
		rank: rank,
		original_id: monomer_graph.compound_id.clone(),
		new_record: Parsed_rBAN_Record::from(&new_variant.new_monomer_graph),
		linearization: new_variant.linearization.clone(),
		old_to_new_mon_map: new_variant.old_to_new_mon_map.clone(),
		monomer_origins,
	    };
	    if debug_stdout {
		println!("Constructed Altered_rBAN_Record for rank {rank} with score {} and linearization {:?}.", variant.score, variant.linearization);
	    }

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

	    if rank + 1 >= max_solutions {
		break;
	    }
	}

    }

    new_variants_with_opt_paths
}

