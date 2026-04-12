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

fn display_lp_with_paths(lp: LogProb, ptc: &PathsToCoords) -> String {
    let mut dag_vertices = vec![0usize]; // start vertex
    for &edge in &ptc.dag_path {
        dag_vertices.push(edge.to);
    }
    format!(
        "LogProb: {:.1}\n\t DAG path: {:?}\n\t HMM path: {:?}",
        lp, dag_vertices, ptc.hmm_path
    )
}

fn build_table_mismatch_error_message(
    coords: &DP_Coords,
    missing_in_brute: &Vec<LogProb>,
    missing_in_dp: &Vec<(LogProb, &PathsToCoords)>,
    parents: &Vec<DP_Coords>,
    tol: f64,
) -> String {
    let mut message = format!(
        "DP table mismatch at coords {:?} (tolerance: {})\n",
        coords, tol
    );
    if !missing_in_brute.is_empty() {
        message.push_str(&format!(
            "LogProbs in DP but not in brute force: {:?}\n",
            missing_in_brute
        ));
    }
    if !missing_in_dp.is_empty() {
        message.push_str(&format!(
            "LogProbs in brute force but not in DP:\n{}\n",
            missing_in_dp
                .iter()
                .map(|(lp, ptc)| display_lp_with_paths(*lp, ptc))
                .collect::<Vec<String>>()
                .join("\n")
        ));
    }

    message.push_str(&format!(
		"\nDP parents:\n{:?}\n",
		parents
	));

    message
}

pub fn check_dp_table_coords<'a>(
    coords: &DP_Coords,
    dp_table: &DP_Table<'a>,
    brute_results: &HashMap<DP_Coords, Vec<PathsToCoords<'a>>>,
    tol: f64,
) -> Option<DP_Table_Mismatch> {
    let empty_vec: Vec<PathsToCoords> = Vec::new(); // for unwrap_or

    let dp_lps: Vec<LogProb> = dp_table
        .get(coords)
        .iter_desc()
        .map(|dlp| dlp.to_logprob())
        .collect();

    let brute_lps_with_paths: Vec<(LogProb, &PathsToCoords)> = brute_results
        .get(&coords)
        .unwrap_or(&empty_vec)
        .iter()
        .map(|ptc| (ptc.lp, ptc))
        .collect();

    let missing_in_brute = dp_lps
        .iter()
        .filter(|&dp_lp| {
            !brute_lps_with_paths
                .iter()
                .any(|(br_lp, _)| (dp_lp - br_lp).abs() < tol)
        })
        .cloned()
        .collect::<Vec<_>>();

    let missing_in_dp = brute_lps_with_paths
        .iter()
        .filter(|(br_lp, _)| !dp_lps.iter().any(|dp_lp| (dp_lp - br_lp).abs() < tol))
        .cloned()
        .collect::<Vec<_>>();

    if missing_in_brute.is_empty() && missing_in_dp.is_empty() {
        None
    } else {
        let message =
            build_table_mismatch_error_message(&coords, &missing_in_brute, &missing_in_dp, &dp_table.get_parents(coords), tol);
        Some(DP_Table_Mismatch {
            coords: coords.clone(),
            message,
        })
    }
}

pub fn find_dp_table_mismatch<'a>(
    hmm: &HMM,
    dag: &DAG<'a>,
    dp_table: &DP_Table<'a>,
    brute_results: &HashMap<DP_Coords, Vec<PathsToCoords<'a>>>,
    tol: f64,
) -> Option<DP_Table_Mismatch> {
    let empty_set: HashSet<DiscreteLogProb> = HashSet::new();
    let tol_dlp = DiscreteLogProb::logprob_to_centered_discrete_lp(tol) as usize;

    let (max_weight, n_vertices, n_states) = (
        dp_table.max_weight(),
        dp_table.n_vertices(),
        dp_table.n_states(),
    );

    for w in 0..=max_weight {
        for v in 0..n_vertices {
            for s in 0..n_states {
                let coords = DP_Coords {
                    weight: w,
                    vertex: v,
                    state: s,
                };
                let dp_parents = dp_table.get_parents(&coords);
                if dp_parents.is_empty() {
                    // DP state is unreachable, skip
                    continue;
                }
		if let Some(mismatch) = check_dp_table_coords(&coords,
							      dp_table,
							      brute_results,
							      tol) {
		    return Some(mismatch)
		}
            }
        }
    }
    None

}
