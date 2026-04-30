use std::{fs, path::Path, io::Write as _};
use tabwriter::TabWriter;

use serde::Deserialize;

use crate::{
    algo::{
        dp::compute_dp_table, solve_brute_force::dp_brute_force
    },
    data_types::{dag::DAG, hmm::HMM}, io::{draw_dag::Draw_DAG_Config, draw_hmm::Draw_HMM_Config, join_svgs::join_svgs_vertical},
};

use crate::testing::dp_vs_bruteforce_table::find_dp_table_mismatch;
use crate::testing::dp_vs_bruteforce_paths::find_paths_mismatch;
use anyhow::{Error, Result};

fn draw_hmm_and_dag(hmm: &HMM, dag: &DAG<'_>, out_dir: &Path) -> Result<()> {
    let hmm_svg_path = out_dir.join("debug_hmm.svg");
    let dag_svg_path = out_dir.join("debug_dag.svg");
    let joined_svg_path = out_dir.join("debug_joined.svg");

    hmm.draw_svg(&hmm_svg_path,
		 &Draw_HMM_Config {
		     state_indexes: true,
		 },
		 None)?;

    dag.draw_svg(&dag_svg_path,
    		 &Draw_DAG_Config {
		     node_indexes: true,
		 },
		 None)?;
    join_svgs_vertical(&[&hmm_svg_path, &dag_svg_path], &joined_svg_path)?;
    Ok(())
}

pub fn dp_vs_bruteforce<'a>(hmm: &HMM, dag: &DAG<'a>, max_weight: usize, tol: f64) -> Option<String> {
    // Brute force once (paths per DP cell).
    let brute_paths_to_coords = dp_brute_force(hmm, dag, max_weight);

    // Real DP once.
    let dp = compute_dp_table(hmm, dag, max_weight);

    // 1) Compare DP tables
    if let Some(m) = find_dp_table_mismatch(
	hmm,
	dag,
	&dp,
	&brute_paths_to_coords,
	tol,
    ) {
	return Some(m.message);
    }

    // 2) Compare actual reconstructed paths
    for w in 0..=max_weight {
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

pub fn write_hmm_emissions(hmm: &HMM, out_path: &Path) {
    // Number of monomer-code columns = max row length (ignoring non-emitting empty rows)
    let emissions_alphabet = hmm.emissions
	.iter()
	.map(|row| row.len())
	.max()
	.unwrap_or(0);

    let mut buf: Vec<u8> = Vec::new();
    {
        let mut tw = TabWriter::new(&mut buf).padding(2);

        // Header: State, then codes 0..m-1
        write!(tw, "State").unwrap();
        for code in 0..emissions_alphabet {
            write!(tw, "\t{code}").unwrap();
        }
        writeln!(tw).unwrap();

        // Rows
        for s in 0..hmm.num_states() {
            let row = &hmm.emissions[s];

            write!(tw, "S{s}").unwrap();

            if row.is_empty() {
                // Non-emitting state: print blank cells so columns align
                for _ in 0..emissions_alphabet {
                    write!(tw, "\t").unwrap();
                }
            } else {
                for lp in row.iter() {
                    write!(tw, "\t{lp:.1}").unwrap();
                }
            }

            writeln!(tw).unwrap();
        }

        tw.flush().unwrap();
    }

    fs::write(out_path, buf).expect("failed to write HMM emissions table");
}


/// Loads tests from JSON and runs `dp_vs_bruteforce` on each.
/// On first failure:
/// - prints error message
/// - writes an SVG (HMM + DAG) into `out_dir`
/// - panics (so tests stop immediately)
pub fn test_dp_vs_bruteforce(tests_json: &Path, out_dir: &Path) {
    let tol = 0.05;
    fs::create_dir_all(out_dir).expect("failed to create out_dir");

    let json_str = fs::read_to_string(tests_json)
        .unwrap_or_else(|e| panic!("failed to read tests JSON at {tests_json:?}: {e}"));

    let tests: Vec<DPvsBruteforceTestCase> = serde_json::from_str(&json_str)
        .unwrap_or_else(|e| panic!("failed to parse tests JSON at {tests_json:?}: {e}"));

    for (idx, tc) in tests.iter().enumerate() {
        // dump SVG next to help debugging
	println!("\n======== Test {}: {}\n", idx, &tc.description);
        let res = draw_hmm_and_dag(&tc.hmm, &tc.dag, out_dir);
	if let Err(e) = res {
	    eprintln!("failed to draw HMM and DAG for test #{idx}: {e}");
	}
	write_hmm_emissions(&tc.hmm, &out_dir.join("hmm_emissions.tsv"));

        if let Some(msg) = dp_vs_bruteforce(&tc.hmm, &tc.dag, tc.max_weight, tol) {
            eprintln!(
                "dp_vs_bruteforce FAILED for test #{idx}: {}\n{}",
                tc.description, msg
            );

            panic!("dp_vs_bruteforce failed; wrote SVG to {out_dir:?}");
        }
    }
}

