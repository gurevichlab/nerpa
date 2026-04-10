use std::{fs, path::PathBuf};

use clap::Parser;
use nerpa_ms_core::data_types::hmm::{BGC_ID, BGC_Variant_ID};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use serde_json::json;

#[derive(Parser, Debug)]
struct Args {
    /// Where to write the JSON file
    #[arg(long)]
    out: PathBuf,

    /// Number of test cases to generate
    #[arg(long, default_value_t = 100)]
    n: usize,

    /// RNG seed (deterministic output)
    #[arg(long, default_value_t = 1)]
    seed: u64,

    /// Max weight to store in each test (DAG path weight budget)
    #[arg(long, default_value_t = 3)]
    max_weight: usize,

    /// Emission alphabet size (MonomerCode range will be 0..alphabet)
    #[arg(long, default_value_t = 4)]
    alphabet: usize,
}

fn random_logprobs(rng: &mut impl Rng, n: usize) -> Vec<f64> {
    // Random categorical distribution -> ln(prob)
    let mut xs: Vec<f64> = (0..n).map(|_| rng.gen_range(0.001..1.0)).collect();
    let sum: f64 = xs.iter().sum();
    for x in &mut xs {
        let lp = (*x / sum).ln();
        *x = (lp * 10.0).round() / 10.0;
    }
    xs
}

fn gen_hmm(rng: &mut impl Rng, alphabet: usize, test_idx: usize) -> serde_json::Value {
    let n_states = rng.gen_range(3..=7);

    // emitting[i] indicates whether state i emits.
    // Start (0) and finish (n-1) are non-emitting.
    let mut emitting = vec![false; n_states];
    for i in 1..(n_states - 1) {
        emitting[i] = rng.gen_bool(0.7);
    }

    let emissions: Vec<serde_json::Value> = (0..n_states)
        .map(|i| {
            if emitting[i] {
                json!(random_logprobs(rng, alphabet))
            } else {
                json!([]) // non-emitting convention
            }
        })
        .collect();

    // transitions[i] is Vec<(to, logprob)>
    // Rule to avoid epsilon-cycles: if state i is NON-emitting, only allow forward transitions (to > i).
    // Emitting states may have back edges, self loops, whatever (bounded by num_emissions).
    let mut transitions: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_states];

    for i in 0..(n_states - 1) {
        // Ensure at least one way to progress:
        let mut tos = vec![i + 1];

        let extra = rng.gen_range(0..=2);
        for _ in 0..extra {
            let to = if !emitting[i] {
                // forward only
                rng.gen_range((i + 1)..n_states)
            } else {
                // anywhere
                rng.gen_range(0..n_states)
            };
            tos.push(to);
        }

        // de-dup tos
        tos.sort_unstable();
        tos.dedup();

        let lps = random_logprobs(rng, tos.len());
        transitions[i] = tos.into_iter().zip(lps.into_iter()).collect();
    }

    // finish has no outgoing transitions
    transitions[n_states - 1] = Vec::new();

    let bgc_variant_id = BGC_Variant_ID{
	bgc_id: BGC_ID{antiSMASH_file: "test".into(), contig_idx: 0, bgc_idx: test_idx},
	variant_idx: 0,
    };
    json!({
        "bgc_variant_id": bgc_variant_id,
        "transitions": transitions.iter().map(|row| {
            row.iter().map(|(to, lp)| json!([to, lp])).collect::<Vec<_>>()
        }).collect::<Vec<_>>(),
        "emissions": emissions,
    })
}

fn gen_dag(rng: &mut impl Rng, alphabet: usize, test_idx: usize) -> serde_json::Value {
    let n_vertices = rng.gen_range(3..=9);
    let start = 0usize;
    let finish = n_vertices - 1;

    let labels: Vec<serde_json::Value> = (0..n_vertices)
        .map(|v| {
            let monomer_code = if v == start || v == finish {
                serde_json::Value::Null
            } else if rng.gen_bool(0.75) {
                json!(rng.gen_range(0..alphabet))
            } else {
                serde_json::Value::Null
            };

            json!({
                "monomer_code": monomer_code,
                "name": format!("v{v}"),
            })
        })
        .collect();

    let mut out_edges: Vec<Vec<serde_json::Value>> = vec![Vec::new(); n_vertices];

    // Backbone so finish is reachable with weight 0
    for v in 0..(n_vertices - 1) {
        out_edges[v].push(json!({
            "to": v + 1,
            "weight": 0,
            "modification": null
        }));
    }

    // Extra edges
    for from in 0..n_vertices {
        let k = rng.gen_range(0..=3);
        for _ in 0..k {
            let mut to = rng.gen_range(0..n_vertices);
            if from == finish {
                continue; // no outgoing from finish usually
            }

            // choose weight; enforce: weight=0 edges must be forward (to > from)
            let mut weight: u8 = if rng.gen_bool(0.5) { 0 } else { 1 };

            if weight == 0 {
                if from + 1 >= n_vertices {
                    weight = 1;
                } else {
                    to = rng.gen_range((from + 1)..n_vertices);
                }
            }

            out_edges[from].push(json!({
                "to": to,
                "weight": weight,
                "modification": null
            }));
        }
    }

    json!({
        "nrp_variant_id": format!("test_dag_{}", test_idx),
        "labels": labels,
        "out_edges": out_edges,
        "start": start,
        "finish": finish
    })
}

fn main() {
    let args = Args::parse();
    let mut rng = ChaCha8Rng::seed_from_u64(args.seed);

    let mut tests = Vec::with_capacity(args.n);

    for i in 0..args.n {
        let hmm = gen_hmm(&mut rng, args.alphabet, i);
        let dag = gen_dag(&mut rng, args.alphabet, i);

        tests.push(json!({
            "hmm": hmm,
            "dag": dag,
            "max_weight": args.max_weight,
            "description": format!("random_{i}_seed{}_mw{}", args.seed, args.max_weight),
        }));
    }

    let out = serde_json::to_string_pretty(&tests).expect("serialize tests");
    // q: create parent dir if it doesn't exist
    if let Some(parent) = args.out.parent() {
        if !parent.exists() {
            fs::create_dir_all(parent).unwrap_or_else(|e| panic!("create parent dir {}: {}", parent.display(), e));
        }
    }
    fs::write(&args.out, out).expect("write tests JSON");
}
