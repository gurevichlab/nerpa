use crate::{algo::{apply_modifications::{AlteredMonomerGraph, apply_modifications}, dp_backtrack::backtrack_solutions}, data_types::{common_types::{LogOdds, MonomerIdx}, config::DebugConfig, graph_modifications::GraphModification, hmm::{BGC_Variant_ID, HMM, StateIdx}, mod_graph::{ModGraph, VertexId}, monomer_graph::MonomerGraph, monomers_db::{MonomerOrigin, MonomersDB}, parsed_rban_record::{MonomerInfo, Parsed_rBAN_Record}}};

use crate::algo::graph_to_dag::create_mod_graph;
use serde::{Serialize, Deserialize};

use crate::algo::dp::compute_dp_table;

use super::dp_backtrack::Solution;
use itertools::Itertools;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Altered_rBAN_Record {
    pub score: LogOdds,
    pub rank: usize,
    pub original_id: String,
    pub bgc_variant_id: BGC_Variant_ID,
    pub new_record: Parsed_rBAN_Record,
    pub linearization: Vec<MonomerIdx>,

    // (Some(old_idx), None) -- removal
    // (None, Some(new_idx)) -- insertion
    // (Some(old_idx), Some(new_idx)) -- substitution (here always old_idx == new_idx)
    // old_idx absent -- keep as is
    pub old_to_new_mon_map: Vec<(Option<MonomerIdx>, Option<MonomerIdx>)>,
    pub monomer_origins: Vec<MonomerOrigin>,  // for debugging record the database monomers used to create the new variant
}

impl Altered_rBAN_Record {
    pub fn linearization_monomers(&self) -> Vec<MonomerInfo> {
	    self.new_record
	        .get_monomers(&self.linearization)
	        .expect("Linearization contains invalid monomer indices")
    }

    pub fn num_modifications(&self) -> usize {
        self.old_to_new_mon_map.len()
    }

    pub fn id(&self) -> String {
        format!(
            "BGC_ID-{}___NRP_ID-{}___NUM_MODS-{}___RANK-{}",
            self.bgc_variant_id,
            self.original_id,
            self.num_modifications(),
            self.rank
        )
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
		
        // rank is different from solutions_fetched because we skip some solutions if they are inconsistent (e.g. remove + insert on the same monomer)
        let mut rank = 0;
	    for (solutions_fetched, (sol, mods)) in solutions_with_mods_unique.enumerate() {

            let new_variant: AlteredMonomerGraph;
	        let apply_res = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
		        apply_modifications(
		            monomer_graph,
		            &mods,
		            monomers_db,
		            debug_stdout,
		        )
	        }));

	        let maybe_variant = match apply_res {
		        Ok(v) => v,
		        Err(_panic) => {
		            eprintln!(
			            "apply_modifications panicked for compound_id={} mods={:?}",
			            monomer_graph.compound_id,
			            mods
		            );
		            panic!("apply_modifications panicked");
		        }
	        };

	        if let Some(variant) = maybe_variant {
		        new_variant = variant;
	        } else {
		        continue; // skip this solution if modifications are inconsistent
	        }

            debug_assert_eq!(
                new_variant.old_to_new_mon_map.len(),
                weight,
                "Number of modifications in old_to_new_mon_map does not match weight"
            );

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

            rank += 1;
			
	        if debug_stdout {
		        println!("Constructing Altered_rBAN_Record");
	        }

	        let variant = {
                // not ideal: new_record id should be updated inside Altered_rBAN_Record::new
                let mut variant = Altered_rBAN_Record {
                    bgc_variant_id: hmm.bgc_variant_id.clone(),
		            score: sol.dlo.to_logodds(),
		            rank: rank,
		            original_id: monomer_graph.compound_id.clone(),
		            new_record: Parsed_rBAN_Record::from(&new_variant.new_monomer_graph),
		            linearization: new_variant.linearization.clone(),
		            old_to_new_mon_map: new_variant.old_to_new_mon_map.clone(),
		            monomer_origins,
	            };
                variant.new_record.compound_id = variant.id();
                variant
            };

	        if debug_stdout {
		        println!(
                    "Constructed Altered_rBAN_Record for rank {rank} with score {} and linearization {:?}.",
                    variant.score,
                    variant.linearization
                );
	        }

	        new_variants_with_opt_paths.push(
                NewVariantWithOptPaths {
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

	        if rank >= max_solutions {
		        break;
	        }
	    }
    }

    new_variants_with_opt_paths
}

