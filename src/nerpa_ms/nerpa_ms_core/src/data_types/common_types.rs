use serde::{Deserialize, Serialize};

// rBAN index of a monomer within the monomer graph
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize)]
pub struct MonomerIdx(pub u32);

// Encoding of a monomer as an integer code --
// same as used in the HMM emissions and the DAG vertex labels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Deserialize, Serialize)]
pub struct MonomerCode(pub u32);

impl MonomerCode {
	pub fn as_usize(&self) -> usize {
		self.0 as usize
	}
}

pub type LogProb = f64;
