use serde::Deserialize;

// rBAN index of a monomer within the monomer graph
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct MonomerIdx(pub u32);

// Encoding of a monomer as an integer code --
// same as used in the HMM emissions and the DAG vertex labels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Deserialize)]
pub struct MonomerCode(pub u32);

pub type LogProb = f64;
