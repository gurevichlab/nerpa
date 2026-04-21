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

#[derive(Debug, Clone, Display, Hash, Default, PartialEq, Eq, Deserialize, Serialize)]
#[display("{bgc_id}:{variant_idx}")]
pub struct BGC_Variant_ID {
    pub bgc_id: BGC_ID,
    pub variant_idx: usize,
}


// HMM data structure (names preserved)
#[derive(Debug, Clone, Deserialize)]
pub struct HMM {
    pub bgc_variant_id: BGC_Variant_ID,

    // transitions[i]: vector of (next_state, log_prob)
    #[serde(deserialize_with = "crate::data_types::json_helpers::de_transitions_null_lp_as_neg_inf")]
    pub transitions: Vec<Vec<(StateIdx, LogProb)>>,
    // emissions[i][j]: log probability of emitting monomer j from state i
    #[serde(deserialize_with = "crate::data_types::json_helpers::de_vec_vec_logprob_null_as_neg_inf")]
    pub emissions: Vec<Vec<LogProb>>,
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
