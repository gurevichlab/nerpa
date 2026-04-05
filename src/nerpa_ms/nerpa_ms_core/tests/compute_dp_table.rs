use nerpa_ms_core::data_types::common_types::LogProb;
use nerpa_ms_core::data_types::dag::{DAG, Edge, VertexLabel};
use nerpa_ms_core::data_types::discrete_log_prob::DiscreteLogProbSet;
use nerpa_ms_core::data_types::hmm::{HMM, BGC_Variant_ID, BGC_ID};
use nerpa_ms_core::algo::dp::compute_dp_table;

fn make_unlabeled_dag(edge_weight: u8) -> DAG<'static> {
    let labels = vec![
        VertexLabel { monomer_code: None, name: "START".into() },
        VertexLabel { monomer_code: None, name: "FINISH".into() },
    ];
    let out_edges = vec![
        vec![Edge { to: 1, weight: edge_weight, modification: None }],
        vec![],
    ];

    DAG {
        nrp_variant_id: "test".into(),
        labels,
        out_edges,
        start: 0,
        finish: 1,
    }
}

fn make_non_emitting_hmm() -> HMM {
    HMM {
        bgc_variant_id: BGC_Variant_ID{
	    bgc_id: BGC_ID{
		antiSMASH_file: "test".into(),
		contig_idx: 0,
		bgc_idx: 0,
	    },
	    variant_idx: 0,
	},
        transitions: vec![vec![(1, 0.0)], vec![]],
        emissions: vec![vec![], vec![]],
    }
}

#[test]
fn dp_propagates_on_unlabeled_dag() {
    let dag = make_unlabeled_dag(0);
    let hmm = make_non_emitting_hmm();

    let dp = compute_dp_table(&hmm, &dag, 0);

    assert_eq!(
        dp.get(dag.finish, 0, 1),
        &DiscreteLogProbSet::from_logprob_vec(vec![0.0])
    );
}

#[test]
fn dp_respects_weight_budget() {
    let dag = make_unlabeled_dag(1);
    let hmm = make_non_emitting_hmm();

    let dp0 = compute_dp_table(&hmm, &dag, 0);
    assert!(dp0.get(dag.finish, 0, 1).is_empty());

    let dp1 = compute_dp_table(&hmm, &dag, 1);
    assert!(!dp1.get(dag.finish, 1, 1).is_empty());
}
