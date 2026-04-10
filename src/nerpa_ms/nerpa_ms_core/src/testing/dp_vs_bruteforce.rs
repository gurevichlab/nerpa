use std::{fs, path::Path};

use serde::Deserialize;

use crate::{
    algo::{
        dp::compute_dp_table, solve_brute_force::dp_brute_force
    },
    data_types::{dag::DAG, hmm::HMM}, io::{draw_dag::Draw_DAG_Config, draw_hmm::Draw_HMM_Config, join_svgs::join_svgs_vertical},
};

use crate::testing::dp_vs_bruteforce_table::find_dp_table_mismatch;
use crate::testing::dp_vs_bruteforce_paths::find_paths_mismatch;

fn draw_hmm_and_dag(hmm: &HMM, dag: &DAG<'_>, out_dir: &Path) {
    let hmm_svg_path = out_dir.join("debug_hmm.svg");
    let dag_svg_path = out_dir.join("debug_dag.svg");
    let joined_svg_path = out_dir.join("debug_joined.svg");
    let _ = hmm.draw_svg(&hmm_svg_path,
			 &Draw_HMM_Config {
			     state_indexes: true,
			 },
			 None);
    let _ = dag.draw_svg(&dag_svg_path,
    			 &Draw_DAG_Config {
			     node_indexes: true,
			 },
			 None);
    let _ = join_svgs_vertical(&[&hmm_svg_path, &dag_svg_path], &joined_svg_path);
}

pub fn dp_vs_bruteforce<'a>(hmm: &HMM, dag: &DAG<'a>, max_weight: usize, tol: f64) -> Option<String> {
    // Brute force once (paths per DP cell).
    let brute_paths_to_coords = dp_brute_force(hmm, dag, max_weight);

    // Real DP once.
    let dp = compute_dp_table(hmm, dag, max_weight);

    for w in 0..=max_weight {
        // 1) Compare DP tables (only weights 0..=w)
        if let Some(m) = find_dp_table_mismatch(
            hmm,
            dag,
            w,
            &dp,
            &brute_paths_to_coords,
            tol,
        ) {
            return Some(m.message);
        }

        // 2) Compare actual reconstructed paths at exactly weight w
        if let Some(m) = find_paths_mismatch(
            hmm,
            dag,
            &dp,
            &brute_paths_to_coords,
            w,
            tol,
        ) {
            return Some(m.message);
        }
    }

    None
}

#[derive(Debug, Clone, Deserialize)]
struct DPvsBruteforceTestCase {
    pub hmm: HMM,
    pub dag: DAG<'static>,
    pub max_weight: usize,
    pub description: String,
}

/// Loads tests from JSON and runs `dp_vs_bruteforce` on each.
/// On first failure:
/// - prints error message
/// - writes an SVG (HMM + DAG) into `out_dir`
/// - panics (so tests stop immediately)
pub fn test_dp_vs_bruteforce(tests_json: &Path, out_dir: &Path) {
    fs::create_dir_all(out_dir).expect("failed to create out_dir");

    let json_str = fs::read_to_string(tests_json)
        .unwrap_or_else(|e| panic!("failed to read tests JSON at {tests_json:?}: {e}"));

    let tests: Vec<DPvsBruteforceTestCase> = serde_json::from_str(&json_str)
        .unwrap_or_else(|e| panic!("failed to parse tests JSON at {tests_json:?}: {e}"));

    for (idx, tc) in tests.iter().enumerate() {
        if let Some(msg) = dp_vs_bruteforce(&tc.hmm, &tc.dag, tc.max_weight, 1e-3) {
            eprintln!(
                "dp_vs_bruteforce FAILED for test #{idx}: {}\n{}",
                tc.description, msg
            );

            // dump SVG next to help debugging
            let _ = draw_hmm_and_dag(&tc.hmm, &tc.dag, out_dir);
            panic!("dp_vs_bruteforce failed; wrote SVG to {out_dir:?}");
        }
    }
}

