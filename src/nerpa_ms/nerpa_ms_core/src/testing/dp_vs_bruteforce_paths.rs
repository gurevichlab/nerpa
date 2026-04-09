use std::collections::HashMap;

use crate::{
    algo::solve_brute_force::{PathsToCoords, retrieve_paths_brute_force},
    data_types::{
        common_types::LogProb,
        dag::{DAG, Edge, VertexId},
        dp_table::{DP_Coords, DP_Table},
        hmm::{HMM, StateIdx},
    },
};

// adjust this import path if your backtracking lives elsewhere
use crate::algo::dp_backtrack::backtrack_solutions;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct PathKey {
    dag_vertices: Vec<VertexId>,
    hmm_states: Vec<StateIdx>,
}

#[derive(Debug, Clone)]
pub struct PathsMismatch {
    pub weight: usize,
    pub message: String,
}

fn dag_vertices_from_edges<'a>(dag: &DAG<'a>, edges: &[Edge<'a>]) -> Vec<usize> {
    let mut vertices = Vec::with_capacity(edges.len() + 1);
    vertices.push(dag.start);
    for e in edges {
        vertices.push(e.to);
    }
    vertices
}

/// Compare *all* paths returned by DP backtracking vs brute force for a given `weight`.
///
/// Path identity = (DAG vertex sequence, HMM state sequence).
/// Logprobs must match within `tol`.
pub fn find_paths_mismatch<'a>(
    hmm: &HMM,
    dag: &DAG<'a>,
    dp: &'a DP_Table<'a>,
    brute_paths_to_coords: &HashMap<DP_Coords, Vec<PathsToCoords<'a>>>,
    weight: usize,
    tol: f64,
) -> Option<PathsMismatch> {
    // ---- collect DP backtracked paths
    let mut dp_pairs: Vec<(PathKey, LogProb)> = Vec::new();
    for sol in backtrack_solutions(weight, dp, dag) {
        let key = PathKey {
            dag_vertices: dag_vertices_from_edges(dag, &sol.dag_edges),
            hmm_states: sol.states.clone(),
        };
        dp_pairs.push((key, sol.dlp.to_logprob()));
    }
    dp_pairs.sort_by(|a, b| a.0.cmp(&b.0));

    // ---- collect brute force paths
    let brute_paths = retrieve_paths_brute_force(hmm, dag, brute_paths_to_coords, weight);

    let mut brute_pairs: Vec<(PathKey, LogProb)> = brute_paths
        .into_iter()
        .map(|p| {
            let key = PathKey {
                dag_vertices: dag_vertices_from_edges(dag, &p.dag_path),
                hmm_states: p.hmm_path,
            };
            (key, p.lp)
        })
        .collect();

    brute_pairs.sort_by(|a, b| a.0.cmp(&b.0));

    // ---- zip-longest compare
    let mut i = 0usize;
    while i < dp_pairs.len() || i < brute_pairs.len() {
        match (dp_pairs.get(i), brute_pairs.get(i)) {
            (Some((dp_key, dp_lp)), Some((br_key, br_lp))) => {
                if dp_key != br_key {
                    return Some(PathsMismatch {
                        weight,
                        message: format!(
                            "Path key mismatch at weight {weight} index {i}.\nDP:   {:?}\nBrut: {:?}",
                            dp_key, br_key
                        ),
                    });
                }

                let diff = (*dp_lp - *br_lp).abs();
                if diff > tol {
                    return Some(PathsMismatch {
                        weight,
                        message: format!(
                            "LogProb mismatch at weight {weight} index {i} for path {:?}\nDP:   {}\nBrut: {}\n|diff|={} > tol={}",
                            dp_key, dp_lp, br_lp, diff, tol
                        ),
                    });
                }
            }

            (Some((dp_key, dp_lp)), None) => {
                return Some(PathsMismatch {
                    weight,
                    message: format!(
                        "Missing brute_force path at weight {weight} index {i}.\nDP has extra: {:?} with lp={}",
                        dp_key, dp_lp
                    ),
                });
            }

            (None, Some((br_key, br_lp))) => {
                return Some(PathsMismatch {
                    weight,
                    message: format!(
                        "Missing DP(backtrack) path at weight {weight} index {i}.\nbrute_force has extra: {:?} with lp={}",
                        br_key, br_lp
                    ),
                });
            }

            (None, None) => unreachable!("loop condition prevents this"),
        }

        i += 1;
    }

    None
}
