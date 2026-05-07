use serde::{Deserialize, Serialize};

use crate::data_types::common_types::{LogProb};

pub type StateIdx = usize;

use anyhow::{anyhow, bail};
use derive_more::Display;

#[derive(Debug, Clone, Display, Hash, Default, PartialEq, Eq, Deserialize, Serialize)]
#[display("{antiSMASH_file}:{contig_idx}:{bgc_idx}")]
pub struct BGC_ID {
    pub antiSMASH_file: String,
    pub contig_idx: usize,
    pub bgc_idx: usize,
}

impl BGC_ID {
    pub fn to_str_short(&self) -> String {
	let path = std::path::Path::new(&self.antiSMASH_file);
	let genome_id = path.file_stem()
	    .and_then(|s| s.to_str())
	    .unwrap_or(&self.antiSMASH_file)
	    .to_string();
	format!("{}:{}:{}",
		genome_id,
		self.contig_idx,
		self.bgc_idx)
    }
}

#[derive(Debug, Clone, Display, Hash, Default, PartialEq, Eq, Deserialize, Serialize)]
#[display("{bgc_id}:{variant_idx}")]
pub struct BGC_Variant_ID {
    pub bgc_id: BGC_ID,
    pub variant_idx: usize,
}


// HMM data structure (names preserved)
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct HMM {
    pub bgc_variant_id: BGC_Variant_ID,

    // transitions[i]: vector of (next_state, log_prob)
    #[serde(deserialize_with = "crate::data_types::json_helpers::de_transitions_null_lp_as_neg_inf",
    	    serialize_with = "crate::data_types::json_helpers::ser_transitions_neg_inf_lp_as_null")]
    pub transitions: Vec<Vec<(StateIdx, LogProb)>>,
    // emissions[i][j]: log probability of emitting monomer j from state i
    #[serde(deserialize_with = "crate::data_types::json_helpers::de_vec_vec_logprob_null_as_neg_inf",
    	    serialize_with = "crate::data_types::json_helpers::ser_vec_vec_logprob_neg_inf_as_null")]
    pub emissions: Vec<Vec<LogProb>>,

    // for debugging and interpretability.
    pub state_types: Vec<String>,

    // module_match_states[i] = state index of the emitting state corresponding to the i-th module in the BGC.
    pub module_match_states: Vec<StateIdx>,

    // not used. Needed for compatibility with python code
    pub module_start_states: Vec<StateIdx>,

    pub module_names: Vec<String>,
}

// Basic helpers
impl HMM {
    pub fn num_states(&self) -> usize {
        self.transitions.len()
    }

    /// Infer number of monomers M from the first emitting state.
    /// Returns None if there are no emitting states.
    pub fn num_monomers(&self) -> Option<usize> {
        self.emissions
            .iter()
            .find(|row| !row.is_empty())
            .map(|row| row.len())
    }

    pub fn is_emitting(&self, state: StateIdx) -> bool {
        !self.emissions[state].is_empty()
    }

    pub fn start_state_idx(&self) -> StateIdx {
	0
    }

    pub fn finish_state_idx(&self) -> StateIdx {
	self.num_states() - 1
    }
}

// validate HMM from JSON
impl HMM {
    pub fn validate(&self) -> Result<(), anyhow::Error> {
        let s = self.transitions.len();
        if s == 0 {
            bail!("HMM has zero states");
        }
        if self.emissions.len() != s {
            bail!(
                "HMM transitions/emissions length mismatch: transitions={}, emissions={}",
                s,
                self.emissions.len()
            );
        }

        // START and FINISH by convention
        let start = 0usize;
        let finish = s - 1;

        if !self.emissions[start].is_empty() {
            bail!("HMM START state (0) must be non-emitting (emissions[0] must be empty)");
        }
        if !self.emissions[finish].is_empty() {
            bail!(
                "HMM FINISH state ({finish}) must be non-emitting (emissions[last] must be empty)"
            );
        }

        // Check transitions in range
        for (i, outs) in self.transitions.iter().enumerate() {
            for &(next, lp) in outs {
                if next >= s {
                    bail!("transition out of range: {i} -> {next}, but S={s}");
                }
                if !lp.is_finite() {
                    bail!("non-finite transition log_prob at {i} -> {next}: {lp}");
                }
            }
        }

        // Infer M and check emitting row lengths
        let m = self
            .num_monomers()
            .ok_or_else(|| anyhow!("HMM has no emitting states (all emissions rows are empty)"))?;

        for (i, row) in self.emissions.iter().enumerate() {
            if row.is_empty() {
                continue; // non-emitting ok
            }
            if row.len() != m {
                bail!(
                    "emission row length mismatch at state {i}: expected M={m}, got {}",
                    row.len()
                );
            }
            for (j, &lp) in row.iter().enumerate() {
                if !lp.is_finite() {
                    bail!("non-finite emission log_prob at state {i}, monomer {j}: {lp}");
                }
            }
        }

        Ok(())
    }
}
