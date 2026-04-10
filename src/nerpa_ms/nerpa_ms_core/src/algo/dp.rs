use crate::data_types::common_types::LogProb;
use crate::data_types::dag::{Edge, DAG};
use crate::data_types::discrete_log_prob::{DiscreteLogProb, DiscreteLogProbSet};
use crate::data_types::dp_table::{DP_Coords, DP_Table};
use crate::data_types::hmm::{StateIdx, HMM};

/// For a fixed (v, w), propagate mass through HMM epsilon transitions (non-emitting states),
/// staying at the same DAG vertex and weight.
/// Two passes handle occasional back-edges (between each two back edges, there must be at least one emitting state, so 2 passes suffice).
fn relax_non_emitting_states(dp: &mut DP_Table, hmm: &HMM, v: usize, w: usize) {
    for _pass in 0..2 {
        for s in (0..hmm.num_states()).filter(|&s| hmm.emissions[s].is_empty()) {
            let cur_coords = DP_Coords {
                vertex: v,
                weight: w,
                state: s,
            };
            if dp.get(&cur_coords).is_empty() {
                // dp state is unreachable, skip
                continue;
            }

            for &(to, edge_lp) in &hmm.transitions[s] {
                let new_coords = DP_Coords {
                    vertex: v,
                    weight: w,
                    state: to,
                };
                dp.update(&new_coords, &cur_coords, Some(edge_lp), None);
            }
        }
    }
}

/// Unlabeled DAG vertex (label=None): advance along DAG edges without consuming an HMM emission.
/// HMM state stays the same; only (vertex, weight) changes.
fn advance_dag_unlabeled<'a>(
    dp: &mut DP_Table<'a>,
    hmm: &HMM,
    dag: &DAG<'a>,
    v: usize,
    w: usize,
    max_weight: usize,
) {
    debug_assert!(
        dag.labels[v].monomer_code.is_none(),
        "advance_dag_unlabeled should be called only on DAG vertices without a monomer code"
    );
    for s in 0..hmm.num_states() {
        let cur_coords = DP_Coords {
            vertex: v,
            weight: w,
            state: s,
        };
        if dp.get(&cur_coords).is_empty() {
            // dp state is unreachable, skip
            continue;
        }

        for dag_edge in &dag.out_edges[v] {
            let new_weight = w + dag_edge.weight as usize;
            if new_weight <= max_weight {
                let new_coords = DP_Coords {
                    vertex: dag_edge.to,
                    weight: new_weight,
                    state: s,
                };
                dp.update(&new_coords, &cur_coords, None, Some(dag_edge.clone()));
            }
        }
    }
}

/// Labeled DAG vertex (label=Some): consume label[v] by emitting from the current emitting HMM state
/// when taking an HMM transition, then advance along a DAG edge.
fn advance_dag_labeled<'a>(
    dp: &mut DP_Table<'a>,
    hmm: &HMM,
    dag: &DAG<'a>,
    v: usize,
    w: usize,
    max_weight: usize,
) {
    let mon_code = dag.labels[v]
        .monomer_code
        .expect("advance_dag_labeled should be called on vertex with a monomer code");
    for s in (0..hmm.num_states()).filter(|&s| !hmm.emissions[s].is_empty()) {
        let cur_coords = DP_Coords {
            vertex: v,
            weight: w,
            state: s,
        };
        if dp.get(&cur_coords).is_empty() {
            // dp state is unreachable, skip
            continue;
        }

        let emission_lp = hmm.emissions[s][mon_code.as_usize()];

        for &(new_state, edge_lp) in &hmm.transitions[s] {
            for dag_edge in &dag.out_edges[v] {
                let new_weight = w + dag_edge.weight as usize;
                if new_weight <= max_weight {
                    // shifted (v, w, s) is recomputed here but that's
                    // fine as most of the time there's just one dag edge
                    let new_coords = DP_Coords {
                        vertex: dag_edge.to,
                        weight: new_weight,
                        state: new_state,
                    };
                    dp.update(
                        &new_coords,
                        &cur_coords,
                        Some(edge_lp + emission_lp),
                        Some(dag_edge.clone()),
                    );
                }
            }
        }
    }
}

/// dp[v][w][s] stores reachable discrete log-probabilities at DAG vertex v with total deviation
/// weight w and HMM state s, before consuming label[v] (if any).
pub fn compute_dp_table<'a>(hmm: &HMM, dag: &DAG<'a>, max_weight: usize) -> DP_Table<'a> {
    let n_vertices = dag.num_nodes();
    let n_states = hmm.num_states();

    // Base: DAG START, 0 edits, HMM START, log(1)=0.
    // is set in DP_Table::new()
    let mut dp = DP_Table::new(n_vertices, max_weight, n_states);

    for v in 0..n_vertices {
        for w in 0..=max_weight {
            if dag.labels[v].monomer_code.is_some() || v == dag.finish {
                relax_non_emitting_states(&mut dp, hmm, v, w);
            }

            match dag.labels[v].monomer_code {
                None => {
                    advance_dag_unlabeled(&mut dp, hmm, dag, v, w, max_weight);
                }
                Some(_) => {
                    advance_dag_labeled(&mut dp, hmm, dag, v, w, max_weight);
                }
            }
        }
    }

    dp
}
