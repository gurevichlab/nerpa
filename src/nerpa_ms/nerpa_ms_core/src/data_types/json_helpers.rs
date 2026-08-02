use crate::algo::gen_new_variants::Altered_rBAN_Record;
use crate::data_types::common_types::{LogOdds, MonomerIdx};
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
) -> Result<Vec<Vec<(StateIdx, LogOdds)>>, D::Error>
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

// serialize transitions: Vec<Vec<(next_state, logprob)>>
// -inf *log-prob* is encoded as `null`
pub fn ser_transitions_neg_inf_lp_as_null<S>(
    transitions: &Vec<Vec<(StateIdx, LogOdds)>>,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    let mapped: Vec<Vec<(StateIdx, Option<f64>)>> = transitions
        .iter()
        .map(|row| {
            row.iter()
                .map(|&(next, lp)| {
                    let lp_opt = if lp.is_infinite() && lp.is_sign_negative() {
                        None
                    } else {
                        Some(lp)
                    };
                    (next, lp_opt)
                })
                .collect()
        })
        .collect();

    mapped.serialize(serializer)
}

// As JSON doesn't support -inf, we allow `null` to represent -inf in the emissions matrix.
pub fn de_vec_vec_logprob_null_as_neg_inf<'de, D>(
    deserializer: D,
) -> Result<Vec<Vec<LogOdds>>, D::Error>
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

// As JSON doesn't support -inf, we serialize -inf as `null` in the emissions matrix.
pub fn ser_vec_vec_logprob_neg_inf_as_null<S>(
    emissions: &Vec<Vec<LogOdds>>,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    let mapped: Vec<Vec<Option<f64>>> = emissions
        .iter()
        .map(|row| {
            row.iter()
                .map(|&lp| {
                    if lp.is_infinite() && lp.is_sign_negative() {
                        None
                    } else {
                        Some(lp)
                    }
                })
                .collect()
        })
        .collect();

    mapped.serialize(serializer)
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

use serde::ser::SerializeMap;

pub fn serialize_new_variants_as_id_map<S>(
    variants: &Vec<Altered_rBAN_Record>,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let mut map = serializer.serialize_map(Some(variants.len()))?;
    for v in variants {
        map.serialize_entry(&v.id(), v)?;
    }
    map.end()
}
