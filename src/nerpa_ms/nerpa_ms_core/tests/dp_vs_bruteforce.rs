use std::collections::{HashMap, HashSet};
use std::path::Path;

use nerpa_ms_core::data_types::common_types::{LogProb, MonomerCode};
use nerpa_ms_core::data_types::dag::{DAG, Edge, VertexLabel};
use nerpa_ms_core::data_types::discrete_log_prob::{DiscreteLogProb, DiscreteLogProbSet};
use nerpa_ms_core::data_types::dp_table::DP_Coords;
use nerpa_ms_core::data_types::hmm::{BGC_ID, BGC_Variant_ID, HMM};

use nerpa_ms_core::algo::dp::compute_dp_table;
use nerpa_ms_core::algo::solve_brute_force::compute_dp_brute_force;

fn dp_cell_to_hashset(cell: &DiscreteLogProbSet) -> HashSet<DiscreteLogProb> {
    cell.iter_desc()
        .collect()
}


fn make_tiny_dag<'a>() -> DAG<'a> {
    // 0 = start, 1 = emits A, 2 = finish
    let labels = vec![
        VertexLabel { monomer_code: None, name: "START".into() },
        VertexLabel { monomer_code: Some(MonomerCode(0)), name: "A".into() },
        VertexLabel { monomer_code: None, name: "FINISH".into() },
    ];

    let mut out_edges: Vec<Vec<Edge<'a>>> = vec![vec![], vec![], vec![]];
    out_edges[0].push(Edge { to: 1, weight: 0, modification: None });
    out_edges[1].push(Edge { to: 2, weight: 0, modification: None });

    DAG {
        nrp_variant_id: "test_dag".into(),
        labels,
        out_edges,
        start: 0,
        finish: 2,
    }
}

fn make_tiny_dag_with_deviation<'a>() -> DAG<'a> {
    // 0 = start, 1 = emits A, 2 = emits B, 3 = finish
    // Two ways to reach B:
    // 0->1 (w0), 1->2 (w0)
    // 0->2 (w1 deviation)
    let labels = vec![
        VertexLabel { monomer_code: None, name: "START".into() },
        VertexLabel { monomer_code: Some(MonomerCode(0)), name: "A".into() },
        VertexLabel { monomer_code: Some(MonomerCode(1)), name: "B".into() },
        VertexLabel { monomer_code: None, name: "FINISH".into() },
    ];

    let mut out_edges: Vec<Vec<Edge<'a>>> = vec![vec![], vec![], vec![], vec![]];
    out_edges[0].push(Edge { to: 1, weight: 0, modification: None });
    out_edges[1].push(Edge { to: 2, weight: 0, modification: None });
    out_edges[0].push(Edge { to: 2, weight: 1, modification: None });
    out_edges[2].push(Edge { to: 3, weight: 0, modification: None });

    DAG {
        nrp_variant_id: "test_dag_dev".into(),
        labels,
        out_edges,
        start: 0,
        finish: 3,
    }
}

fn make_tiny_hmm() -> HMM {
    // States: 0 START (non-emit), 1 emit, 2 FINISH (non-emit)
    // Alphabet size = 2 (A,B)
    let bgc_variant_id = BGC_Variant_ID{
	bgc_id: BGC_ID{
	    antiSMASH_file: "test".into(),
	    contig_idx: 0,
	    bgc_idx: 0,
	},
	variant_idx: 0,
    };
    HMM {
        bgc_variant_id,
        transitions: vec![
            vec![(1, -0.1), (2, -1.0)], // 0 -> 1 or 2
            vec![(2, -0.2)],            // 1 -> 2
            vec![],                     // 2
        ],
        emissions: vec![
            vec![],                     // 0 non-emit
            vec![-0.3, -1.2],            // 1 emits A/B
            vec![],                     // 2 non-emit
        ],
    }
}

fn make_hmm_with_silent_middle() -> HMM {
    // 0 START non-emit
    // 1 silent
    // 2 emit
    // 3 FINISH silent
    let bgc_variant_id = BGC_Variant_ID{
	bgc_id: BGC_ID{
	    antiSMASH_file: "test_hmm_silent".into(),
	    contig_idx: 0,
	    bgc_idx: 0,
	},
	variant_idx: 0,
    };
    HMM {
        bgc_variant_id,
        transitions: vec![
            vec![(1, -0.1)],            // 0->1
            vec![(2, -0.2), (3, -2.0)],  // 1->2 or 1->3
            vec![(3, -0.3)],            // 2->3
            vec![],
        ],
        emissions: vec![
            vec![],                      // 0
            vec![],                      // 1 silent
            vec![-0.4, -0.6],             // 2 emits A/B
            vec![],                      // 3
        ],
    }
}

fn sets_almost_equal(set1: &HashSet<DiscreteLogProb>, set2: &HashSet<DiscreteLogProb>, tol: usize) -> bool {
    if set1.len() != set2.len() {
	return false;
    }
    for lp1 in set1 {
	if !set2.iter().any(|lp2| lp1.0.abs_diff(lp2.0) <= tol){
	    return false;
	}
    }
    true
}

use nerpa_ms_core::io::draw_hmm::Draw_HMM_Config;
use nerpa_ms_core::io::draw_dag::Draw_DAG_Config;
use nerpa_ms_core::io::join_svgs::join_svgs_vertical;

fn draw_hmm_and_dag(hmm: &HMM, dag: &DAG<'_>) {
    let hmm_svg_path = Path::new("debug_hmm.svg");
    let dag_svg_path = Path::new("debug_dag.svg");
    let joined_svg_path = Path::new("debug_joined.svg");
    let _ = hmm.draw_svg(Path::new("debug_hmm.svg"),
			 &Draw_HMM_Config {
			     state_indexes: true,
			 },
			 None);
    let _ = dag.draw_svg(Path::new("debug_dag.svg"),
    			 &Draw_DAG_Config {
			     node_indexes: true,
			 },
			 None);
    let _ = join_svgs_vertical(&[hmm_svg_path, dag_svg_path], joined_svg_path);
}

fn assert_dp_matches_bruteforce(hmm: &HMM, dag: &DAG<'_>, max_weight: usize, tol: f64) {
    let tol_dlp = DiscreteLogProb::logprob_to_centered_discrete_lp(tol) as usize;
    let brute = compute_dp_brute_force(hmm, dag, max_weight);
    let dp = compute_dp_table(hmm, dag, max_weight);

    let n_vertices = dag.labels.len();
    let n_states = hmm.num_states();

    let get_cell_sets = |c: &DP_Coords| {
        let brute_set = brute
            .get(c)
            .cloned()
            .unwrap_or_default();
        let dp_set = dp_cell_to_hashset(dp.get(c));
        (dp_set, brute_set)
    };

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
                let (dp_set, brute_set) = get_cell_sets(&c);
                let ok = sets_almost_equal(&dp_set, &brute_set, tol_dlp);
                matches[dp.idx(&c)] = ok;
            }
        }
    }

    if matches.iter().all(|&x| x) {
        return;
    }

    // Find a "frontier" mismatch: mismatching cell whose all parents match.
    for idx in 0..total {
        if matches[idx] {
            continue;
        }

        let c = dp.idx_to_coordinates(idx);
        let parents = dp.get_parents(&c);

        let all_parents_match = parents
            .iter()
            .all(|p| matches[dp.idx(p)]);

        if !all_parents_match {
            continue;
        }

        let (dp_set, brute_set) = get_cell_sets(&c);
	draw_hmm_and_dag(hmm, dag);
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

        panic!(
            "Frontier mismatch at cell (v={}, w={}, s={})\n\
             In DP but no near match in Brute (tol={}): {:?}\n\
             In Brute but no near match in DP (tol={}): {:?}",
            c.vertex,
            c.weight,
            c.state,
            tol,
            dp_only_no_near,
            tol,
            brute_only_no_near,
        );
    }

    // Fallback: no frontier mismatch found (e.g. parent links missing); report first mismatch.
    for idx in 0..total {
        if !matches[idx] {
	    let c = dp.idx_to_coordinates(idx);
	    panic!(
		"Mismatch at cell (v={}, w={}, s={})\nNo frontier mismatch found, so likely parent links are missing in DP_Table.",
		c.vertex, c.weight, c.state,
	    );
	}
    }
}

#[test]
fn dp_matches_bruteforce_on_tiny_cases() {
    // Case 1: simplest possible DAG + HMM
    let tol = 0.01; // small tolerance for floating point differences
    let hmm = make_tiny_hmm();
    let dag = make_tiny_dag();
    assert_dp_matches_bruteforce(&hmm, &dag, 0, tol);

    // Case 2: DAG has a deviation edge (weight=1), test max_weight=1
    let hmm = make_tiny_hmm();
    let dag = make_tiny_dag_with_deviation();
    assert_dp_matches_bruteforce(&hmm, &dag, 1, tol);

    // Case 3: HMM has a silent middle state
    let hmm = make_hmm_with_silent_middle();
    let dag = make_tiny_dag_with_deviation();
    assert_dp_matches_bruteforce(&hmm, &dag, 1, tol);
}

fn make_medium_dag<'a>() -> DAG<'a> {
    // 0 START
    // 1 emits A
    // 2 epsilon (None)        (internal silent vertex in DAG)
    // 3 emits B
    // 4 emits C
    // 5 FINISH
    //
    // Paths:
    // 0->1->3->4->5                      (w=0) emits A B C
    // 0->1->2->3->4->5                   (w=1) emits A B C   (extra epsilon step)
    // 0->1->3->5                         (w=1) emits A B     (skip C with deviation)
    // 0->1->2->3->5                      (w=2) emits A B
    let labels = vec![
        VertexLabel { monomer_code: None, name: "START".into() },
        VertexLabel { monomer_code: Some(MonomerCode(0)), name: "A".into() },
        VertexLabel { monomer_code: None, name: "eps".into() },
        VertexLabel { monomer_code: Some(MonomerCode(1)), name: "B".into() },
        VertexLabel { monomer_code: Some(MonomerCode(2)), name: "C".into() },
        VertexLabel { monomer_code: None, name: "FINISH".into() },
    ];

    let mut out_edges: Vec<Vec<Edge<'a>>> = vec![vec![]; labels.len()];

    out_edges[0].push(Edge { to: 1, weight: 0, modification: None });

    out_edges[1].push(Edge { to: 3, weight: 0, modification: None });
    out_edges[1].push(Edge { to: 2, weight: 1, modification: None }); // deviation to epsilon

    out_edges[2].push(Edge { to: 3, weight: 0, modification: None });

    out_edges[3].push(Edge { to: 4, weight: 0, modification: None });
    out_edges[3].push(Edge { to: 5, weight: 1, modification: None }); // deviation: skip C

    out_edges[4].push(Edge { to: 5, weight: 0, modification: None });

    DAG {
        nrp_variant_id: "medium_dag".into(),
        labels,
        out_edges,
        start: 0,
        finish: 5,
    }
}

fn make_medium_hmm_abc() -> HMM {
    // States:
    // 0 START silent
    // 1 silent
    // 2 emit
    // 3 emit
    // 4 FINISH silent
    //
    // Alphabet size = 3 (A,B,C) => MonomerCode(0..=2)
    let bgc_variant_id = BGC_Variant_ID {
        bgc_id: BGC_ID {
            antiSMASH_file: "test_medium_hmm".into(),
            contig_idx: 0,
            bgc_idx: 0,
        },
        variant_idx: 0,
    };

    HMM {
        bgc_variant_id,
        transitions: vec![
            vec![(1, -0.1), (2, -0.9)], // 0 -> 1 or jump to first emitter
            vec![(2, -0.2), (4, -2.0)], // 1 -> 2 or early finish
            vec![(3, -0.3), (4, -1.5)], // 2 -> 3 or finish after 1 emission
            vec![(4, -0.4)],            // 3 -> 4
            vec![],
        ],
        emissions: vec![
            vec![], // 0 silent
            vec![], // 1 silent
            vec![-0.1, -1.0, -2.0], // 2 emits A/B/C (prefers A)
            vec![-2.0, -0.2, -1.2], // 3 emits A/B/C (prefers B)
            vec![], // 4 silent
        ],
    }
}

fn make_hmm_with_back_edge() -> HMM {
    // Like medium HMM, but with a back-edge 3 -> 2 (u > v).
    // This is exactly the scenario your DP implementation mentions (needs 2 passes).
    let bgc_variant_id = BGC_Variant_ID {
        bgc_id: BGC_ID {
            antiSMASH_file: "test_hmm_backedge".into(),
            contig_idx: 0,
            bgc_idx: 0,
        },
        variant_idx: 0,
    };

    HMM {
        bgc_variant_id,
        transitions: vec![
            vec![(1, -0.1)],                 // 0 -> 1
            vec![(2, -0.2)],                 // 1 -> 2
            vec![(3, -0.3), (4, -2.0)],      // 2 -> 3 or finish
            vec![(2, -0.7), (4, -0.4)],      // 3 -> 2 back-edge OR finish
            vec![],
        ],
        emissions: vec![
            vec![],                          // 0 silent
            vec![],                          // 1 silent
            vec![-0.2, -0.5, -1.5],          // 2 emits A/B/C
            vec![-1.0, -0.1, -1.2],          // 3 emits A/B/C (prefers B)
            vec![],                          // 4 silent
        ],
    }
}

#[test]
fn dp_matches_bruteforce_on_medium_cases() {
    // Medium DAG + medium HMM, max_weight=2
    // (Still small enough that brute force shouldn't melt your laptop.)
    let dag = make_medium_dag();
    let hmm = make_medium_hmm_abc();
    let tol = 0.01; // small tolerance for floating point differences
    assert_dp_matches_bruteforce(&hmm, &dag, 2, tol);
}

#[test]
fn dp_matches_bruteforce_with_hmm_back_edge() {
    // Same medium DAG, but HMM includes a back-edge (3 -> 2).
    // This specifically tests the "iterate states twice" assumption.
    let dag = make_medium_dag();
    let hmm = make_hmm_with_back_edge();
    let tol = 0.01; // small tolerance for floating point differences
    assert_dp_matches_bruteforce(&hmm, &dag, 2, tol);
}

fn make_bigger_dag<'a>() -> DAG<'a> {
    // 0 START
    // 1 A
    // 2 B
    // 3 C
    // 4 D
    // 5 X (insertion-ish)
    // 6 E
    // 7 FINISH
    //
    // Edges:
    // 0->1 (0)
    // 1->2 (0)
    // 1->3 (1)  skip B (deviation)
    // 2->3 (0)
    // 2->4 (1)  skip C (deviation)
    // 3->4 (0)
    // 4->6 (0)  no insertion
    // 4->5 (1)  insert X (deviation)
    // 5->6 (0)
    // 6->7 (0)
    //
    // This gives a small family of paths with weights 0..2 and emissions length ~4..6.
    let labels = vec![
        VertexLabel { monomer_code: None, name: "START".into() },            // 0
        VertexLabel { monomer_code: Some(MonomerCode(0)), name: "A".into() }, // 1
        VertexLabel { monomer_code: Some(MonomerCode(1)), name: "B".into() }, // 2
        VertexLabel { monomer_code: Some(MonomerCode(2)), name: "C".into() }, // 3
        VertexLabel { monomer_code: Some(MonomerCode(3)), name: "D".into() }, // 4
        VertexLabel { monomer_code: Some(MonomerCode(5)), name: "X".into() }, // 5
        VertexLabel { monomer_code: Some(MonomerCode(4)), name: "E".into() }, // 6
        VertexLabel { monomer_code: None, name: "FINISH".into() },           // 7
    ];

    let mut out_edges: Vec<Vec<Edge<'a>>> = vec![vec![]; labels.len()];

    out_edges[0].push(Edge { to: 1, weight: 0, modification: None });

    out_edges[1].push(Edge { to: 2, weight: 0, modification: None });
    out_edges[1].push(Edge { to: 3, weight: 1, modification: None }); // skip B

    out_edges[2].push(Edge { to: 3, weight: 0, modification: None });
    out_edges[2].push(Edge { to: 4, weight: 1, modification: None }); // skip C

    out_edges[3].push(Edge { to: 4, weight: 0, modification: None });

    out_edges[4].push(Edge { to: 6, weight: 0, modification: None });
    out_edges[4].push(Edge { to: 5, weight: 1, modification: None }); // insert X

    out_edges[5].push(Edge { to: 6, weight: 0, modification: None });

    out_edges[6].push(Edge { to: 7, weight: 0, modification: None });

    DAG {
        nrp_variant_id: "bigger_dag".into(),
        labels,
        out_edges,
        start: 0,
        finish: 7,
    }
}

fn make_bigger_hmm_emit_on_leave() -> HMM {
    // Emit-on-leave convention:
    // If emissions[state] is non-empty, leaving that state consumes 1 emission.
    //
    // Alphabet size = 6: A,B,C,D,E,X (0..=5).
    //
    // States:
    // 0 START silent
    // 1 emit (wants A)
    // 2 emit (wants B/C)
    // 3 silent (router)
    // 4 emit (wants C/D)
    // 5 emit (wants D/X)
    // 6 emit (wants X/E)
    // 7 emit (wants E)
    // 8 FINISH silent
    //
    // Branches are limited (mostly 1-2 options), so brute force stays reasonable.
    let bgc_variant_id = BGC_Variant_ID {
        bgc_id: BGC_ID {
            antiSMASH_file: "test_bigger_hmm".into(),
            contig_idx: 0,
            bgc_idx: 0,
        },
        variant_idx: 0,
    };

    let a = 0usize;
    let b = 1usize;
    let c = 2usize;
    let d = 3usize;
    let e = 4usize;
    let x = 5usize;

    // Helper: make a length-6 emission vector with a simple preference.
    fn emis(pref: usize, second: usize) -> Vec<LogProb> {
        let mut v = vec![-3.0; 6];
        v[pref] = -0.1;
        v[second] = -0.5;
        v
    }

    HMM {
        bgc_variant_id,
        transitions: vec![
            vec![(1, -0.1), (3, -1.2)], // 0 -> 1 or skip to router (fewer emissions)
            vec![(2, -0.2), (3, -1.0)], // 1 -> 2 or router
            vec![(3, -0.3)],            // 2 -> router
            vec![(4, -0.2), (7, -1.0)], // 3 -> 4 (longer) or jump to near-end (shorter)
            vec![(5, -0.2)],            // 4 -> 5
            vec![(6, -0.2), (7, -0.8)], // 5 -> 6 (extra emission) or -> 7
            vec![(7, -0.2)],            // 6 -> 7
            vec![(8, -0.2)],            // 7 -> finish
            vec![],                     // 8
        ],
        emissions: vec![
            vec![],                // 0 silent
            emis(a, b),            // 1
            emis(b, c),            // 2
            vec![],                // 3 silent
            emis(c, d),            // 4
            emis(d, x),            // 5
            emis(x, e),            // 6
            emis(e, d),            // 7
            vec![],                // 8 silent
        ],
    }
}

#[test]
fn dp_matches_bruteforce_on_bigger_case() {
    let dag = make_bigger_dag();
    let hmm = make_bigger_hmm_emit_on_leave();
    let tol = 0.01; // small tolerance for floating point differences

    // Keep this at 2; with brute force, 3 can get annoying depending on branching.
    assert_dp_matches_bruteforce(&hmm, &dag, 2, tol);
}
