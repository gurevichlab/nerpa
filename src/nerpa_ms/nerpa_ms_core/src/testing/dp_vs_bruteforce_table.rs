use std::collections::{HashMap, HashSet};

use crate::{
    algo::solve_brute_force::{compute_dp_table_brute_force, PathsToCoords},
    data_types::{
        common_types::LogProb,
        dag::DAG,
        discrete_log_prob::{DiscreteLogProb, DiscreteLogProbSet},
        dp_table::{DP_Coords, DP_Table},
        hmm::HMM,
    },
};

pub struct DP_Table_Mismatch {
    pub coords: DP_Coords,
    pub dp_only_no_near: Vec<LogProb>,
    pub brute_only_no_near: Vec<LogProb>,
    pub message: String,
}

fn sets_almost_equal(
    set1: &HashSet<DiscreteLogProb>,
    set2: &HashSet<DiscreteLogProb>,
    tol: usize,
) -> bool {
    if set1.len() != set2.len() {
        return false;
    }
    for lp1 in set1 {
        if !set2.iter().any(|lp2| lp1.0.abs_diff(lp2.0) <= tol) {
            return false;
        }
    }
    true
}

fn find_frontier_mismatch<'a>(
    dp_results: &DP_Table<'a>,
    brute_results: &HashMap<DP_Coords, HashSet<DiscreteLogProb>>,
    tol_dlp: usize,
    matches: &Vec<bool>,
) -> DP_Table_Mismatch {
    // Find a "frontier" mismatch: mismatching cell whose all parents match.
    let empty_set = HashSet::new(); // for when I need a reference to an empty set
    for idx in 0..matches.len() {
        if matches[idx] {
            continue;
        }

        let c = dp_results.idx_to_coordinates(idx);
        let parents = dp_results.get_parents(&c);

        let all_parents_match = parents.iter().all(|p| matches[dp_results.idx(p)]);

        if !all_parents_match {
            continue;
        }

        let brute_set = brute_results.get(&c).unwrap_or(&empty_set);
        let dp_set: HashSet<DiscreteLogProb> = dp_results.get(&c).iter_desc().collect();

        let dp_only_no_near: Vec<LogProb> = dp_set
            .iter()
            .filter(|lp1| !brute_set.iter().any(|lp2| lp1.0.abs_diff(lp2.0) <= tol_dlp))
            .map(|lp| lp.to_logprob())
            .collect();
        let brute_only_no_near: Vec<LogProb> = brute_set
            .iter()
            .filter(|lp2| !dp_set.iter().any(|lp1| lp2.0.abs_diff(lp1.0) <= tol_dlp))
            .map(|lp| lp.to_logprob())
            .collect();

	let message = format!(
                "Mismatch at coords {:?}\n\nDP-only {:?}\nBrute-only {:?}",
                c, dp_only_no_near, brute_only_no_near
            );
        return DP_Table_Mismatch {
            coords: c,
            dp_only_no_near,
            brute_only_no_near,
            message,
        };
    }
    panic!("No frontier mismatch found, so likely parent links are missing in DP_Table.");
}

pub fn find_dp_table_mismatch<'a>(
    hmm: &HMM,
    dag: &DAG<'a>,
    max_weight: usize,
    dp_results: &DP_Table<'a>,
    brute_results: &HashMap<DP_Coords, Vec<PathsToCoords<'a>>>,
    tol: f64,
) -> Option<DP_Table_Mismatch> {
    let empty_set: HashSet<DiscreteLogProb> = HashSet::new();
    let tol_dlp = DiscreteLogProb::logprob_to_centered_discrete_lp(tol) as usize;
    let brute = compute_dp_table_brute_force(brute_results);

    let n_vertices = dag.labels.len();
    let n_states = hmm.num_states();

    // Precompute which cells match (so parent checks are cheap).
    let total = n_vertices * (max_weight + 1) * n_states;
    let mut matches = vec![true; total];

    for v in 0..n_vertices {
        for w in 0..=max_weight {
            for s in 0..n_states {
                let c = DP_Coords {
                    vertex: v,
                    weight: w,
                    state: s,
                };
                let brute_set = brute.get(&c).unwrap_or(&empty_set);
                let dp_set: HashSet<DiscreteLogProb> = dp_results.get(&c).iter_desc().collect();
                let ok = sets_almost_equal(&dp_set, &brute_set, tol_dlp);
                matches[dp_results.idx(&c)] = ok;
            }
        }
    }

    if matches.iter().all(|&x| x) {
        return None;
    } else {
        Some(find_frontier_mismatch(
            dp_results, &brute, tol_dlp, &matches,
        ))
    }
}
