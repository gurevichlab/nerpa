use crate::data_types::common_types::{LogProb, MonomerIdx};
use crate::data_types::hmm::StateIdx;
use serde::{Deserialize, Deserializer};

impl<'de> Deserialize<'de> for MonomerIdx {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        #[derive(Deserialize)]
        #[serde(untagged)]
        enum Repr {
            Num(u32),
            Str(String),
        }

        match Repr::deserialize(deserializer)? {
            Repr::Num(n) => Ok(MonomerIdx(n)),
            Repr::Str(s) => {
                let n = s.parse::<u32>().map_err(serde::de::Error::custom)?;
                Ok(MonomerIdx(n))
            }
        }
    }
}

// transitions JSON: Vec<Vec<[next_state, (number|null)]>>
// null encodes -inf *log-prob*
pub fn de_transitions_null_lp_as_neg_inf<'de, D>(
    deserializer: D,
) -> Result<Vec<Vec<(StateIdx, LogProb)>>, D::Error>
where
    D: Deserializer<'de>,
{
    // Parse the second tuple element as Option<f64> so `null` is accepted.
    let raw: Vec<Vec<(StateIdx, Option<f64>)>> = Vec::deserialize(deserializer)?;

    Ok(raw
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|(next, lp)| (next, lp.unwrap_or(f64::NEG_INFINITY)))
                .collect()
        })
        .collect())
}

// As JSON doesn't support -inf, we allow `null` to represent -inf in the emissions matrix.
pub fn de_vec_vec_logprob_null_as_neg_inf<'de, D>(
    deserializer: D,
) -> Result<Vec<Vec<LogProb>>, D::Error>
where
    D: Deserializer<'de>,
{
    // First parse as Option so `null` is accepted.
    let raw: Vec<Vec<Option<f64>>> = Vec::deserialize(deserializer)?;

    // Then map `null` -> -inf.
    let mapped = raw
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|x| x.unwrap_or(f64::NEG_INFINITY))
                .collect::<Vec<f64>>()
        })
        .collect();

    Ok(mapped)
}

use crate::data_types::parsed_rban_record::{AtomId, MonomerEdge, MonomerEdgeInfo, MonomerInfo};

impl<'de> Deserialize<'de> for AtomId {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        #[derive(Deserialize)]
        #[serde(untagged)]
        enum Repr {
            Num(u32),
            Str(String),
        }

        match Repr::deserialize(deserializer)? {
            Repr::Num(n) => Ok(AtomId(n)),
            Repr::Str(s) => {
                let n = s.parse::<u32>().map_err(serde::de::Error::custom)?;
                Ok(AtomId(n))
            }
        }
    }
}

use crate::data_types::parsed_rban_record::{BondType, MonomerEdgeInfoSingle};

// In the JSON, `all_edges` entries are tuples like:
//   [ { "4": 6, "5": 7 }, 1, "AMINO" ]
// (not objects with named fields). Accept both representations.
impl<'de> Deserialize<'de> for MonomerEdgeInfoSingle {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        #[serde(untagged)]
        enum ArityRepr {
            Num(f64),
            Str(String),
        }

        #[derive(Deserialize)]
        #[serde(untagged)]
        enum Repr {
            Obj {
                monomer_to_atom: HashMap<MonomerIdx, AtomId>,
                arity: ArityRepr,
                bond_type: BondType,
            },
            Tup(HashMap<MonomerIdx, AtomId>, ArityRepr, BondType),
        }

        match Repr::deserialize(deserializer)? {
            Repr::Obj {
                monomer_to_atom,
                arity,
                bond_type,
            } => {
                let arity = match arity {
                    ArityRepr::Num(n) => n.to_string(),
                    ArityRepr::Str(s) => s,
                };
                Ok(MonomerEdgeInfoSingle {
                    monomer_to_atom,
                    arity,
                    bond_type,
                })
            }
            Repr::Tup(monomer_to_atom, arity, bond_type) => {
                let arity = match arity {
                    ArityRepr::Num(n) => n.to_string(),
                    ArityRepr::Str(s) => s,
                };
                Ok(MonomerEdgeInfoSingle {
                    monomer_to_atom,
                    arity,
                    bond_type,
                })
            }
        }
    }
}

use std::collections::HashMap;
use std::hash::Hash;

// JSON shape:
//   "atomic_bonds": [
//     [ [0, 1], { ... } ],
//     [ [1, 2], { ... } ]
//   ]
pub fn de_vec_pairs_to_hashmap<'de, D, K, V>(deserializer: D) -> Result<HashMap<K, V>, D::Error>
where
    D: Deserializer<'de>,
    K: Deserialize<'de> + Eq + Hash,
    V: Deserialize<'de>,
{
    let pairs: Vec<(K, V)> = Vec::deserialize(deserializer)?;
    Ok(pairs.into_iter().collect())
}

use serde::ser::SerializeSeq;
use serde::{Serialize, Serializer};

pub fn ser_hashmap_as_vec_pairs<S, K, V>(
    map: &HashMap<K, V>,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: Serializer,
    K: Serialize + Eq + Hash,
    V: Serialize,
{
    let mut seq = serializer.serialize_seq(Some(map.len()))?;
    for (k, v) in map {
        // serializes each entry as a 2-tuple: [key, value]
        seq.serialize_element(&(k, v))?;
    }
    seq.end()
}

// q: a function that takes a json value which is either a string or a number, and returns a string
pub fn de_str_or_num_to_str<'de, D>(deserializer: D) -> Result<String, D::Error>
where
    D: Deserializer<'de>,
{
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum Repr {
        Str(String),
        Num(f64),
    }

    match Repr::deserialize(deserializer)? {
        Repr::Str(s) => Ok(s),
        Repr::Num(n) => Ok(n.to_string()),
    }
}
