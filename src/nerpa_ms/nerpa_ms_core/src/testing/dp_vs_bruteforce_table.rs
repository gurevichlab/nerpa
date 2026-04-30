use std::collections::{HashMap, HashSet};

use crate::{
    algo::solve_brute_force::{PathsToCoords, compute_dp_table_brute_force},
    data_types::{
        common_types::LogProb,
        dag::DAG,
        discrete_log_prob::{DiscreteLogProb, DiscreteLogProbSet, MAX_LOG_PROB, MIN_LOG_PROB},
        dp_table::{DP_Coords, DP_Table},
        hmm::HMM,
    },
};

pub struct DP_Table_Mismatch {
    pub coords: DP_Coords,
    pub message: String,
}

fn is_dag_emitting(dag: &DAG<'_>, v: usize) -> bool {
    dag.labels[v].monomer_code.is_some()
}

fn is_hmm_emitting(hmm: &HMM, s: usize) -> bool {
    !hmm.emissions[s].is_empty()
}

/// Build aligned "columns". Each column may have a DAG vertex, an HMM state, or both.
/// We align emitting steps with emitting steps; epsilon steps create a `None` on the other side.
fn align_paths(
    hmm: &HMM,
    dag: &DAG<'_>,
    dag_vertices: &[usize],
    hmm_states: &[usize],
) -> (Vec<Option<usize>>, Vec<Option<usize>>) {
    let mut hmm_cols: Vec<Option<usize>> = Vec::new();
    let mut dag_cols: Vec<Option<usize>> = Vec::new();

    // Turn slices into iterators we can "peek" at without consuming.
    let mut hmm_it = hmm_states.iter().copied().peekable();
    let mut dag_it = dag_vertices.iter().copied().peekable();

    while hmm_it.peek().is_some() && dag_it.peek().is_some() {
        let hmm_state = *hmm_it.peek().unwrap();
        let dag_vertex = *dag_it.peek().unwrap();

        let hmm_emitting = is_hmm_emitting(hmm, hmm_state);
        let dag_emitting = is_dag_emitting(dag, dag_vertex);

        match (hmm_emitting, dag_emitting) {
            // DAG is non-emitting: consume DAG only, put a gap in HMM.
            (_, false) => {
                hmm_cols.push(None);
                dag_cols.push(Some(dag_it.next().unwrap()));
            }
            // HMM is non-emitting, DAG is emitting: consume HMM only, put a gap in DAG.
            (false, true) => {
                hmm_cols.push(Some(hmm_it.next().unwrap()));
                dag_cols.push(None);
            }
            // Both emitting: consume both.
            (true, true) => {
                hmm_cols.push(Some(hmm_it.next().unwrap()));
                dag_cols.push(Some(dag_it.next().unwrap()));
            }
        }
    }

    // Drain the remainder (same logic as your trailing while-loops, just iterator style).
    while let Some(hmm_state) = hmm_it.next() {
        hmm_cols.push(Some(hmm_state));
        dag_cols.push(None);
    }
    while let Some(dag_vertex) = dag_it.next() {
        hmm_cols.push(None);
        dag_cols.push(Some(dag_vertex));
    }

    (hmm_cols, dag_cols)
}

/// Print aligned columns "compactly":
/// - we only print actual numbers (Some)
/// - gaps (None) become extra spaces *before the next printed number* (no extra commas)
fn format_aligned_compact(cols: &[Option<usize>]) -> String {
    let mut out = String::new();
    out.push('[');

    for col in cols {
        match col {
            Some(v) => {
                out.push_str(&format!("{v}, "));
            }
            None => {
		out.push_str("   "); // 3 spaces for the gap (adjust as needed)
            }
        }
    }

    out.push(']');
    out
}

fn display_lp_with_paths(lp: LogProb, ptc: &PathsToCoords, hmm: &HMM, dag: &DAG<'_>) -> String {
    let mut dag_vertices = vec![0usize]; // start vertex
    for &edge in &ptc.dag_path {
        dag_vertices.push(edge.to);
    }
    let hmm_states: Vec<usize> = ptc.hmm_path.clone();

    let (hmm_cols, dag_cols) = align_paths(hmm, dag, &dag_vertices, &hmm_states);

    let hmm_str = format_aligned_compact(&hmm_cols);
    let dag_str = format_aligned_compact(&dag_cols);

    format!(
        "LogProb: {:.1}\n\t HMM path: {}\n\t DAG path: {}",
        lp, hmm_str, dag_str
    )
}

fn build_table_mismatch_error_message(
    coords: &DP_Coords,
    hmm: &HMM,
    dag: &DAG<'_>,
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
                .map(|(lp, ptc)| display_lp_with_paths(*lp, ptc, hmm, dag))
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
    hmm: &HMM,
    dag: &DAG<'a>,
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
        .filter_map(|ptc| if ptc.lp > MIN_LOG_PROB { Some((ptc.lp, ptc))} else { None })
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
            build_table_mismatch_error_message(&coords,
					       hmm, dag,
					       &missing_in_brute, &missing_in_dp, &dp_table.get_parents(coords), tol);
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
		if !((v == dag.finish || dag.labels[v].monomer_code.is_some())
		     && !hmm.emissions[s].is_empty()) {
		    // check only "emitting" DP states (those that consume a DAG vertex and an HMM state)
		    // intermittent states may be unreachable in DP but reachable in brute force, which would cause false positives
		    continue;
		}

		let dp_parents = dp_table.get_parents(&coords);
                if dp_parents.is_empty() {
                    // DP state is unreachable, skip
                    continue;
                }
		if let Some(mismatch) = check_dp_table_coords(&coords,
							      hmm,
							      dag,
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

