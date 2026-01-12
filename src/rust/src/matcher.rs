use crate::data_types::*;
use crate::viterbi_algorithm::{get_hmm_score, get_hmm_score_with_path};
use std::collections::HashMap;
use ordered_float::OrderedFloat;
use std::cmp::Reverse;

pub fn linearization_score(hmm: &HMM, linearization: &NRP_Linearization) -> LogProb {
    get_hmm_score(hmm, &linearization.mon_codes)
}

pub fn get_best_linearization<'a>(
    hmm: &HMM,
    nrp_linearizations: &'a [NRP_Linearization],
) -> Option<(LogProb, &'a NRP_Linearization)> {
    // Determine the best non-iterative linearization
    nrp_linearizations
        .iter()
        .map(|lin| (linearization_score(hmm, lin), lin))
        .max_by_key(|(score, _lin)| *score)
}

pub fn get_best_split_linearizations<'a>(
    hmm: &HMM,
    groups_split: &'a [Vec<NRP_Linearization>],
) -> (LogProb, Vec<&'a NRP_Linearization>) {
    // Each group is a vector of linearizations
    // For each group, determine the best linearization. Then combine them
    let best_linearizations_with_scores: Vec<(LogProb, &NRP_Linearization)> =
        groups_split.iter()
            .map(|group| get_best_linearization(hmm, group)
                .expect("A linearization should be always present"))
            .collect();

    let combined_score: LogProb = best_linearizations_with_scores.iter()
        .map(|(score, _)| score)
        .sum();

    let best_linearizations: Vec<&NRP_Linearization> =
        best_linearizations_with_scores.iter()
            .map(|(_score, lin)| *lin)
            .collect();

    (combined_score, best_linearizations)
}

pub fn get_best_linearizations_for_nrp<'a>(
    hmm: &HMM,
    nrp_linearizations: &'a NRP_Linearizations_Info,
) -> (LogProb, Vec<&'a NRP_Linearization>) {
    // return score and best linearizations (one for non-iterative or many for iterative)

    // Determine the best non-iterative linearization
    let (best_non_iterative_score, best_non_iterative_linearization) =
        get_best_linearization(hmm, &nrp_linearizations.non_iterative)
            .expect("A non-iterative linearization should be always present");

    // Determine the best iterative linearizations
    let (best_iterative_score, best_iterative_linearizations) =
        nrp_linearizations.iterative.iter()
            .map(|groups_split| get_best_split_linearizations(hmm, groups_split))
            .max_by_key(|(score, _lin)| *score)
            .unwrap_or((NEG_INFINITY, Vec::new()));

    if best_non_iterative_score > best_iterative_score {
        (best_non_iterative_score, vec![best_non_iterative_linearization])
    } else {
        (best_iterative_score, best_iterative_linearizations)
    }
    
}


pub fn build_match<'a>(hmm: &'a HMM,
                       nrp_linearizations_info: &'a NRP_Linearizations_Info) -> MatchInfoLight <'a>{
    let (raw_score, linearizations) = get_best_linearizations_for_nrp(hmm, nrp_linearizations_info);
    let score = raw_score - nrp_linearizations_info.score_vs_avg_bgc;

    MatchInfoLight {
        raw_score,
        score,
        bgc_variant_id: &hmm.bgc_variant_id,
        nrp_id: &nrp_linearizations_info.nrp_id,
        linearizations,
    }
}

pub fn get_best_matches_for_bgc_variant<'a>(
    hmm: &'a HMM,
    nrps_info: &'a [NRP_Linearizations_Info],
    max_num_matches_per_bgc: usize,
) -> Vec<MatchInfoLight<'a>> {

    let mut matches: Vec<MatchInfoLight> =
        nrps_info.iter()
        .map(|nrp_linearizations_info| build_match(hmm, nrp_linearizations_info))
        .collect();

    let max_num_matches =
        if max_num_matches_per_bgc > 0 { max_num_matches_per_bgc }
        else { matches.len() };

    matches.sort_by_key(|m| m.score);
    matches.truncate(max_num_matches);

    matches
}

pub fn get_best_matches_for_nrp<'a>(
    nrp_info: &'a NRP_Linearizations_Info,
    hmms: &'a [HMM],
    max_num_matches_per_nrp: usize,
) -> Vec<MatchInfoLight<'a>> {

    let mut matches: Vec<MatchInfoLight> =
        hmms.iter()
            .map(|hmm| build_match(hmm, nrp_info))
            .collect();

    let max_num_matches =
        if max_num_matches_per_nrp > 0 { max_num_matches_per_nrp }
        else { matches.len() };

    matches.sort_by_key(|m| m.score);
    matches.truncate(max_num_matches);

    matches
}

pub fn filter_and_sort_matches<'a, 'b>(
    matches: &'b [MatchInfoLight<'a>],
    config: &MatchingConfig,
) -> Vec<&'b MatchInfoLight <'a>> {
    // Step 0: keep only the best match for each (BGC, NRP) pair
    let mut unique_matches_map: HashMap<(&BGC_ID, &NRP_ID), &MatchInfoLight> = HashMap::new();

    for m in matches.iter() {
        let bgc_id = &m.bgc_variant_id.bgc_id;
        let nrp_id  = &m.nrp_id;

        unique_matches_map.entry((bgc_id, nrp_id))
            .and_modify(|e| { if m.score > e.score { *e = &m; } })
            .or_insert(m);
    }

    let mut unique_matches: Vec<&MatchInfoLight> = unique_matches_map.values().cloned().collect();

    // Step 1: Sort all matches by descending score
    unique_matches.sort_by_key(|m| Reverse(m.score));

    // Step 2: Keep only top matches for each BGC and NRP
    let mut bgc_counts: HashMap<&BGC_ID, usize> = HashMap::new();
    let mut nrp_counts: HashMap<&NRP_ID, usize> = HashMap::new();
    let mut best_matches = Vec::new();

    for m in unique_matches.iter() {
        let bgc_id = &m.bgc_variant_id.bgc_id;
        let nrp_id  = &m.nrp_id;

        let bgc_counts_entry = bgc_counts.entry(bgc_id).or_insert(0);
        let nrp_counts_entry = nrp_counts.entry(nrp_id).or_insert(0);

        let match_rank_bgc = *bgc_counts_entry;
        let match_rank_nrp = *nrp_counts_entry;

        if match_rank_bgc < config.min_num_matches_per_bgc ||
            match_rank_nrp < config.min_num_matches_per_nrp ||
            (match_rank_bgc < config.max_num_matches_per_bgc && match_rank_nrp < config.max_num_matches_per_nrp)
        {
            best_matches.push(m.clone());
            *bgc_counts_entry += 1;
            *nrp_counts_entry += 1;
        }
    }

    // Step 3: Return only the top max_num_matches overall
    best_matches.truncate(config.max_num_matches);

    best_matches
}

pub fn get_full_match_info<'a>(match_light: &MatchInfoLight<'a>, hmm: &'a HMM) -> MatchInfo<'a> {
    let opt_paths: Vec<Vec<StateIdx>> =
        match_light.linearizations.iter()
        .map(|linearization| get_hmm_score_with_path(hmm, &linearization.mon_codes).1)
        .collect();

    MatchInfo {
        raw_score: match_light.raw_score,
        score: match_light.score,
        bgc_variant_id: match_light.bgc_variant_id,
        nrp_id: match_light.nrp_id,
        linearizations: match_light.linearizations.clone(),
        optimal_paths: opt_paths,
    }
}

pub fn get_matches<'a>(
    hmms: &'a [HMM],
    nrps_info: &'a [NRP_Linearizations_Info],
    config: &MatchingConfig,
    _num_threads: usize, // ignored as requested
) -> Vec<MatchInfo<'a>> {
    let mut matches_light = Vec::new();

    // Process each BGC variant
    for hmm in hmms.iter() {
        matches_light.extend(
            get_best_matches_for_bgc_variant(hmm, nrps_info,
                                             config.max_num_matches_per_bgc)
        );
    }

    // Process each NRP
    for nrp_info in nrps_info.iter() {
        matches_light.extend(
            get_best_matches_for_nrp(nrp_info, hmms,
                                     config.max_num_matches_per_nrp)
        );
    }

    // Filter and sort the matches
    let matches_light_filtered = filter_and_sort_matches(&matches_light, config);

    // Reconstruct the full match info
    let hmms_by_id: HashMap<&BGC_Variant_ID, &HMM> =
        hmms.iter()
            .map(|hmm| (&hmm.bgc_variant_id, hmm))
            .collect();

    matches_light_filtered.iter()
        .map(|match_light| get_full_match_info(match_light, hmms_by_id[&match_light.bgc_variant_id]))
        .collect()
}
