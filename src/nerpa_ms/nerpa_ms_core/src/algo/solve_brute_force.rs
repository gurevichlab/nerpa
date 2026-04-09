use std::collections::{HashMap, HashSet};

use crate::data_types::common_types::{LogProb, MonomerCode};
use crate::data_types::dag::{Edge, DAG};
use crate::data_types::discrete_log_prob::{DiscreteLogProb, DiscreteLogProbSet};
use crate::data_types::dp_table::DP_Coords;
use crate::data_types::hmm::{StateIdx, HMM};

pub fn all_dag_paths_until<'a>(v: usize, w: usize, dag: &DAG<'a>) -> Vec<Vec<Edge<'a>>> {
    let mut paths = Vec::new();
    let mut stack = vec![(dag.start, 0, Vec::new())]; // (vertex, weight, path)
    while let Some((current_vertex, current_weight, current_path)) = stack.pop() {
        if (current_vertex, current_weight) == (v, w) {
            paths.push(current_path);
            continue;
        }
        for edge in &dag.out_edges[current_vertex] {
            let new_weight = current_weight + edge.weight as usize;
            // For "paths until (v,w)" we only need to explore up to w, not max_weight.
            if new_weight <= w {
                let mut new_path = current_path.clone();
                new_path.push((*edge).clone());
                stack.push((edge.to, new_weight, new_path));
            }
        }
    }
    paths
}

pub fn all_hmm_paths_until(
    s: usize,
    num_emissions: usize,
    hmm: &HMM,
) -> Vec<Vec<(usize, LogProb)>> {
    let mut paths = Vec::new();
    let mut stack = vec![(0, 0, Vec::new())]; // (state, path)
    while let Some((current_state, current_num_emissions, current_path)) = stack.pop() {
        if (current_state, current_num_emissions) == (s, num_emissions) {
            paths.push(current_path);
            continue;
        }

        for &(to, edge_lp) in &hmm.transitions[current_state] {
            let new_emission = if hmm.emissions[current_state].is_empty() {
                0
            } else {
                1
            };
            let new_num_emissions = current_num_emissions + new_emission;

            if new_num_emissions <= num_emissions {
                let mut new_path = current_path.clone();
                new_path.push((to, edge_lp));
                stack.push((to, new_num_emissions, new_path));
            }
        }
    }
    paths
}

#[derive(Debug, Clone)]
pub struct PathsToCoords<'a> {
    pub dag_path: Vec<Edge<'a>>,
    pub hmm_path: Vec<StateIdx>,
    pub lp: LogProb,
}

pub fn dp_brute_force<'a>(
    hmm: &HMM,
    dag: &DAG<'a>,
    max_weight: usize,
) -> HashMap<DP_Coords, Vec<PathsToCoords<'a>>> {
    let mut results = HashMap::new();
    for v in 0..dag.num_nodes() {
        for w in 0..=max_weight {
            let dag_paths = all_dag_paths_until(v, w, dag);
            for dag_path in &dag_paths {
                // DP semantics: emissions are labels of vertices along the DAG path
                // EXCLUDING the final vertex `v`.
                let emissions: Vec<MonomerCode> = dag_path
                    .iter()
                    .take(dag_path.len().saturating_sub(1))
                    .filter_map(|edge| dag.labels[edge.to].monomer_code)
                    .collect();

                for s in 0..hmm.num_states() {
                    let hmm_paths = all_hmm_paths_until(s, emissions.len(), hmm);

                    for hmm_path in &hmm_paths {
                        // Mirror the same convention: don't include emission of the final state `s`.
                        let emitting_states: Vec<usize> = hmm_path
                            .iter()
                            .take(hmm_path.len().saturating_sub(1))
                            .filter(|&&(state, _)| !hmm.emissions[state].is_empty())
                            .map(|&(state, _)| state)
                            .collect();
                        let total_transition_lp: LogProb = hmm_path.iter().map(|&(_, lp)| lp).sum();
                        let total_emission_lp: LogProb = emitting_states
                            .iter()
                            .zip(emissions.iter())
                            .map(|(&state, &mon_code)| hmm.emissions[state][mon_code.as_usize()])
                            .sum();

                        let hmm_states: Vec<StateIdx> = std::iter::once(0)  // prepend start state
			    .chain(hmm_path.iter()
				   .map(|&(state, _)| state)).collect();
                        results
                            .entry(DP_Coords {
                                vertex: v,
                                weight: w,
                                state: s,
                            })
                            .or_insert_with(Vec::new)
                            .push(PathsToCoords {
                                dag_path: dag_path.clone(),
                                hmm_path: hmm_states,
                                lp: total_transition_lp + total_emission_lp,
                            });
                    }
                }
            }
        }
    }

    results
}

pub fn compute_dp_table_brute_force(
    paths_to_coords: &HashMap<DP_Coords, Vec<PathsToCoords>>,
) -> HashMap<DP_Coords, HashSet<DiscreteLogProb>> {
    let mut dp_table = HashMap::new();
    for (coords, paths) in paths_to_coords {
        let dlps: HashSet<DiscreteLogProb> = paths
            .iter()
            .map(|p| DiscreteLogProb::from_logprob(p.lp))
            .collect();
        dp_table.insert(coords.clone(), dlps);
    }
    dp_table
}

pub fn retrieve_paths_brute_force<'a>(
    hmm: &HMM,
    dag: &DAG<'a>,
    paths_to_coords: &HashMap<DP_Coords, Vec<PathsToCoords<'a>>>,
    target_weight: usize,
) -> Vec<PathsToCoords<'a>> {
    let target_coords = DP_Coords {
	vertex: dag.finish,
	weight: target_weight,
	state: hmm.num_states() - 1, // assuming last state is finish
    };
    paths_to_coords
	.get(&target_coords)
	.unwrap_or(&Vec::new())
	.iter()
	.cloned()
	.collect()
}
