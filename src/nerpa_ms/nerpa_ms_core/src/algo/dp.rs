use crate::data_types::common_types::LogProb;
use crate::data_types::dag::DAG;
use crate::data_types::discrete_log_prob::DiscreteLogProbSet;
use crate::data_types::dp_table::DP_Table;
use crate::data_types::hmm::HMM;

/// For a fixed (v, w), propagate mass through HMM epsilon transitions (non-emitting states),
/// staying at the same DAG vertex and weight.
/// Two passes handle occasional back-edges (between each two back edges, there must be at least one emitting state, so 2 passes suffice).
fn relax_non_emitting_states(dp: &mut DP_Table, hmm: &HMM, v: usize, w: usize) {
    for _pass in 0..2 {
        for s in (0..hmm.num_states()).filter(|&s| hmm.emissions[s].is_empty()) {
            if dp.get(v, w, s).is_empty() {
                // dp state is unreachable, skip
                continue;
            }

            for &(to, edge_lp) in &hmm.transitions[s] {
                let shifted_dlp_set = dp.get(v, w, s).add_to_all(edge_lp);
                dp.get_mut(v, w, to).union_inplace(&shifted_dlp_set);
            }
        }
    }
}

/// Unlabeled DAG vertex (label=None): advance along DAG edges without consuming an HMM emission.
/// HMM state stays the same; only (vertex, weight) changes.
fn advance_dag_unlabeled(
    dp: &mut DP_Table,
    hmm: &HMM,
    dag: &DAG<'_>,
    v: usize,
    w: usize,
    max_weight: usize,
) {
    assert!(
        dag.labels[v].monomer_code.is_none(),
        "advance_dag_unlabeled should be called only on DAG vertices without a monomer code"
    );
    for s in 0..hmm.num_states() {
        if dp.get(v, w, s).is_empty() {
            // dp state is unreachable, skip
            continue;
        }

        for dag_edge in &dag.out_edges[v] {
            let new_weight = w + dag_edge.weight as usize;
            if new_weight <= max_weight {
                dp.union_cell_into_cell(dp.idx(v, w, s), dp.idx(dag_edge.to, new_weight, s));
            }
        }
    }
}

/// Labeled DAG vertex (label=Some): consume label[v] by emitting from the current emitting HMM state
/// when taking an HMM transition, then advance along a DAG edge.
fn advance_dag_labeled(
    dp: &mut DP_Table,
    hmm: &HMM,
    dag: &DAG<'_>,
    v: usize,
    w: usize,
    max_weight: usize,
) {
    let mon_code = dag.labels[v]
        .monomer_code
        .expect("advance_dag_labeled should be called on vertex with a monomer code");
    for s in (0..hmm.num_states()).filter(|&s| !hmm.emissions[s].is_empty()) {
        if dp.get(v, w, s).is_empty() {
            // dp state is unreachable, skip
            continue;
        }

        let emission_lp = hmm.emissions[s][mon_code.as_usize()];

        for &(new_state, edge_lp) in &hmm.transitions[s] {
            let shifted_dlp_set = dp.get(v, w, s).add_to_all(edge_lp + emission_lp);

            for dag_edge in &dag.out_edges[v] {
                let new_weight = w + dag_edge.weight as usize;
                if new_weight <= max_weight {
                    dp.get_mut(dag_edge.to, new_weight, new_state)
                        .union_inplace(&shifted_dlp_set);
                }
            }
        }
    }
}

/// dp[v][w][s] stores reachable discrete log-probabilities at DAG vertex v with total deviation
/// weight w and HMM state s, before consuming label[v] (if any).
pub fn compute_dp_table(hmm: &HMM, dag: &DAG<'_>, max_weight: usize) -> DP_Table {
    let n_vertices = dag.labels.len();
    let n_states = hmm.transitions.len();

    let mut dp = DP_Table::new(n_vertices, max_weight, n_states);

    // Base: DAG START, 0 edits, HMM START, log(1)=0.
    *dp.get_mut(dag.start, 0, 0) = DiscreteLogProbSet::from_logprob_vec(vec![0.0]);

    for v in 0..n_vertices {
        for w in 0..=max_weight {
            relax_non_emitting_states(&mut dp, hmm, v, w);

            match dag.labels[v].monomer_code {
                None => advance_dag_unlabeled(&mut dp, hmm, dag, v, w, max_weight),
                Some(_) => advance_dag_labeled(&mut dp, hmm, dag, v, w, max_weight),
            }
        }
    }

    dp
}
