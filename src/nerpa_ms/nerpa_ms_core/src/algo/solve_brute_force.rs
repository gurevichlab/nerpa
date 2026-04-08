use std::collections::{HashMap, HashSet};

use crate::data_types::common_types::{LogProb, MonomerCode};
use crate::data_types::dag::{Edge, DAG};
use crate::data_types::discrete_log_prob::{DiscreteLogProb, DiscreteLogProbSet};
use crate::data_types::hmm::HMM;

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
	    let new_emission = if hmm.emissions[current_state].is_empty() {0} else {1};
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

pub fn compute_dp_brute_force(
    hmm: &HMM,
    dag: &DAG<'_>,
    max_weight: usize,
) -> HashMap<(usize, usize, usize), HashSet<DiscreteLogProb>> {
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
			let total_emission_lp: LogProb = emitting_states.iter()
			    .zip(emissions.iter())
			    .map(|(&state, &mon_code)| hmm.emissions[state][mon_code.as_usize()])
			    .sum();

			results.entry((v, w, s))
			    .or_insert_with(HashSet::new)
			    .insert(DiscreteLogProb::from_logprob(total_transition_lp + total_emission_lp));
		    }
		}
	    }
	}
    }

    results
}
