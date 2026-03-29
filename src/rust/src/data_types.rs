// data_types.rs — line-for-line Rust rewrite of the given C++ header
// Names and field order are preserved. Notes inline where Rust differs.
use ordered_float::NotNan;
pub const ZERO: NotNan<f64> = unsafe { NotNan::new_unchecked(0.0) };
pub const ONE: NotNan<f64> = unsafe { NotNan::new_unchecked(1.0) };
pub const INFINITY: NotNan<f64> = unsafe { NotNan::new_unchecked(f64::INFINITY) };
pub const NEG_INFINITY: NotNan<f64> = unsafe { NotNan::new_unchecked(f64::NEG_INFINITY) };

use crate::command_line_args::Args;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::path::Path;

// ========== Basic Type Aliases ==========
pub type MonCode = usize;
pub type rBAN_idx = usize;
pub type LogProb = NotNan<f64>;
pub type StateIdx = usize;
pub type NRP_ID = String;

// ========== Data Structures ==========

#[derive(Debug, Clone, Hash, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct BGC_ID {
    pub genome_id: String,
    pub contig_idx: usize,
    pub bgc_idx: usize,
}

#[derive(Debug, Clone, Hash, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct BGC_Variant_ID {
    pub bgc_id: BGC_ID,
    pub variant_idx: usize,
}

// For each linearization, store sequences of monomers (first) and rBAN indices (second).
// C++: struct NRP_Linearization { std::vector<MonCode> mon_codes; std::vector<rBAN_idx> rban_idxs; };
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct NRP_Linearization {
    pub mon_codes: Vec<MonCode>,
    pub rban_idxs: Vec<rBAN_idx>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct NRP_Linearizations_Info {
    pub nrp_id: NRP_ID,
    pub non_iterative: Vec<NRP_Linearization>,
    pub iterative: Vec<Vec<Vec<NRP_Linearization>>>,
    pub score_vs_avg_bgc: LogProb, // used for ranking matches
}

impl NRP_Linearizations_Info {
    pub fn load_nrps_linearizations(path: &Path) -> Vec<NRP_Linearizations_Info> {
        let data = std::fs::read_to_string(path).expect("Failed to read NRP linearizations file");
        let nrps: Vec<NRP_Linearizations_Info> =
            serde_json::from_str(&data).expect("Failed to parse NRP linearizations JSON");
        nrps
    }
}

// HMM data structure (names preserved)
#[derive(Debug, Clone, Deserialize)]
pub struct HMM {
    pub bgc_variant_id: BGC_Variant_ID,
    // transitions[i]: vector of (next_state, log_prob)
    pub transitions: Vec<Vec<(StateIdx, LogProb)>>,
    // emissions[i][j]: log probability of emitting monomer j from state i
    pub emissions: Vec<Vec<LogProb>>,
}

impl HMM {
    pub fn num_states(&self) -> usize {
        self.transitions.len()
    }
    pub fn check_shapes(&self) {
        assert_eq!(
            self.transitions.len(),
            self.emissions.len(),
            "states mismatch"
        );
        let mut a = 0usize;
        for e in &self.emissions {
            if !e.is_empty() {
                a = e.len();
                break;
            }
        }
        for e in &self.emissions {
            if !e.is_empty() {
                assert_eq!(e.len(), a);
            }
        }
    }
    pub fn load_hmms(path: &Path) -> Vec<HMM> {
        let data = std::fs::read_to_string(path).expect("Failed to read HMMs file");
        let hmms: Vec<HMM> = serde_json::from_str(&data).expect("Failed to parse HMMs JSON");
        for hmm in &hmms {
            hmm.check_shapes();
        }
        hmms
    }
}

pub fn max_if_0(value: usize) -> usize {
    if value > 0 {
        value
    } else {
        usize::MAX
    }
}

#[derive(Debug, Clone, Copy, Deserialize)]
pub struct MatchingConfig {
    pub max_num_matches_per_bgc: usize,
    pub max_num_matches_per_nrp: usize,
    pub min_num_matches_per_bgc: usize,
    pub min_num_matches_per_nrp: usize,
    pub max_num_matches: usize,
}
impl MatchingConfig {
    pub fn load_from_args(args: &Args) -> MatchingConfig {
        MatchingConfig {
            max_num_matches_per_bgc: max_if_0(args.max_num_matches_per_bgc),
            max_num_matches_per_nrp: max_if_0(args.max_num_matches_per_nrp),
            min_num_matches_per_bgc: args.min_num_matches_per_bgc,
            min_num_matches_per_nrp: args.min_num_matches_per_nrp,
            max_num_matches: max_if_0(args.max_num_matches),
        }
    }
}

// ========== MatchInfo Structures ==========
#[derive(Debug, Clone, Serialize)]
pub struct MatchInfo<'a> {
    pub raw_score: LogProb,
    pub score: LogProb,
    pub bgc_variant_id: &'a BGC_Variant_ID,
    pub nrp_id: &'a NRP_ID,
    pub linearizations: Vec<&'a NRP_Linearization>,
    pub optimal_paths: Vec<Vec<StateIdx>>, // one path per linearization
}
impl MatchInfo<'_> {
    pub fn dump_json_list(matches: &[MatchInfo], path: &Path) {
        let data =
            serde_json::to_string_pretty(matches).expect("Failed to serialize matches to JSON");
        std::fs::write(path, data).expect("Failed to write matches JSON to file");
    }
}

#[derive(Debug, Clone)]
pub struct MatchInfoLight<'a> {
    pub raw_score: LogProb,
    pub score: LogProb, // used for ranking matches
    pub bgc_variant_id: &'a BGC_Variant_ID,
    pub nrp_id: &'a NRP_ID,
    pub linearizations: Vec<&'a NRP_Linearization>,
}
