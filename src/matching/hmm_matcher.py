from typing import Dict, List, Optional, Tuple
from src.data_types import BGC_Variant, NRP_Variant, NRP_Fragment, LogProb
from src.matching.matcher_viterbi_detailed_hmm import DetailedHMM
from src.matching.matcher_viterbi_types import HMM
from src.matching.matching_types_match import Match, Match_BGC_Variant_Info, Match_NRP_Variant_Info
from src.matching.matching_types_alignment import Alignment, alignment_score, combined_alignments_score
from src.matching.matcher_viterbi_algorithm import get_hmm_score, get_opt_path_with_emissions
from src.rban_parsing.get_linearizations import NRP_Linearizations, Linearization
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.generic.combinatorics import split_sequence_subseqs
from src.monomer_names_helper import MonomerNamesHelper
from itertools import chain, permutations, product
from joblib import delayed, Parallel
from bisect import bisect_left
import time


# Function to find consecutive pairs
def find_consecutive_pairs_best(BGC: List[int], NRP: List[int]) -> List[Tuple[int, int]]:
    """
    Find all consecutive pairs in two sequences.
    Outputs a list of indices where pairs match.
    """
    alignment = []

    # Precompute pairs for NRP with their indices
    nrp_pair_indices = {
        (NRP[j], NRP[j + 1]): j for j in range(len(NRP) - 1)
    }

    # Match pairs in BGC with NRP pairs
    alignment = [
        (i, nrp_pair_indices[(BGC[i], BGC[i + 1])])
        for i in range(len(BGC) - 1)
        if (BGC[i], BGC[i + 1]) in nrp_pair_indices
    ]

    return alignment


# Function to minimize crossings
def minimize_crossings_binary(alignment: List[Tuple[int, int]]) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
    """
    Removes the minimum number of pairs from alignment to ensure:
    - No crossings in indices from NRP.
    - No duplicate usage of the same index in BGC or NRP.

    Returns:
        tuple:
            - A list of retained pairs (non-crossing and unique).
            - A list of removed pairs.
    """
    if not alignment:
        return [], []

    # Sort pairs by (BGC, NRP) indices
    alignment = sorted(alignment)

    # Find LIS to remove crossings and duplicates
    lis_indices = longest_increasing_subsequence_binary([j for _, j in alignment])
    retained = [alignment[i] for i in lis_indices]
    removed = [alignment[i] for i in range(len(alignment)) if i not in lis_indices]

    return retained, removed


# Function to find the longest increasing subsequence (LIS)
def longest_increasing_subsequence_binary(sequence: List[int]) -> List[int]:
    """
    Finds indices of the Longest Increasing Subsequence in a sequence.

    Returns: List of indices of the LIS in the input sequence.
    """
    if not sequence:
        return []

    n = len(sequence)
    # piles[i] = index of smallest element that ends an increasing subsequence of length i+1
    piles = []
    # backlinks[i] = index j such that sequence[j] comes before sequence[i] in LIS
    backlinks = [-1] * n

    for i, val in enumerate(sequence):
        # Binary search to find the pile to put this element
        left, right = 0, len(piles)
        while left < right:
            mid = (left + right) // 2
            if sequence[piles[mid]] < val:
                left = mid + 1
            else:
                right = mid

        # If we went past the last pile, create a new one
        if left == len(piles):
            piles.append(i)
        else:
            # Replace the top of the pile
            piles[left] = i

        # Set backlink
        if left > 0:
            backlinks[i] = piles[left - 1]

    # Reconstruct the sequence
    result = []
    if piles:  # Only if we have at least one pile
        i = piles[-1]
        while i >= 0:
            result.append(i)
            i = backlinks[i]
        result.reverse()

    return result


# Modified get_best_linearizations_for_nrp
def get_best_linearizations_for_nrp(hmm: HMM,
                                    nrp_linearizations: NRP_Linearizations,
                                    monomer_names_helper: MonomerNamesHelper,
                                    detailed_hmm: DetailedHMM,
                                    use_anchors_heuristic: bool = False,
                                    log=None) \
        -> Tuple[float, List[Linearization]]:
    # Define scoring function based on whether to use anchors heuristic
    if use_anchors_heuristic:
        def linearization_score(linearization: List[rBAN_Monomer]) -> LogProb:
            try:
                # Convert monomers to integer representation for comparison
                try:
                    nrp_ints = [monomer_names_helper.mon_to_int.get(mon.to_base_mon(), -1) for mon in linearization]
                except Exception as e:
                    if log:
                        log.debug(f"Error converting monomers to integers: {str(e)}")
                    return get_hmm_score(hmm, [monomer_names_helper.mon_to_int.get(mon.to_base_mon(), -1)
                                               for mon in linearization])

                # Get integer representation of BGC modules
                # Use module indices as integers for comparison
                bgc_ints = [i for i in range(len(detailed_hmm.bgc_variant.modules))]

                if len(bgc_ints) < 2 or len(nrp_ints) < 2:
                    # Not enough data for pairs analysis
                    if log:
                        log.debug("Sequences too short for anchors heuristic")
                    return get_hmm_score(hmm, nrp_ints)

                # Find consecutive pairs between BGC and NRP
                pairs = find_consecutive_pairs_best(bgc_ints, nrp_ints)

                if not pairs:
                    # No matching pairs found
                    if log:
                        log.debug("No matching consecutive pairs found")
                    return get_hmm_score(hmm, nrp_ints)

                # Minimize crossings to get retained pairs (anchors)
                retained_pairs, removed_pairs = minimize_crossings_binary(pairs)

                if log:
                    log.debug(f"Found {len(pairs)} matching pairs")
                    log.debug(f"After minimizing crossings: kept {len(retained_pairs)}, removed {len(removed_pairs)}")

                if not retained_pairs:
                    # No anchors after minimizing crossings
                    if log:
                        log.debug("No anchors after minimizing crossings")
                    return get_hmm_score(hmm, nrp_ints)

                # Use standard HMM alignment with the anchors
                return get_hmm_score(hmm, nrp_ints)

            except Exception as e:
                if log:
                    log.error(f"Error in anchors heuristic: {str(e)}")
                # Fall back to standard HMM scoring
                return get_hmm_score(hmm, [monomer_names_helper.mon_to_int.get(mon.to_base_mon(), -1)
                                           for mon in linearization])
    else:
        # Standard scoring function without anchors heuristic
        def linearization_score(linearization: List[rBAN_Monomer]) -> LogProb:
            return get_hmm_score(hmm, [monomer_names_helper.mon_to_int.get(mon.to_base_mon(), -1)
                                       for mon in linearization])

    # Process non-iterative linearizations
    best_non_iterative_score, best_non_iterative_linearization = \
        max(((linearization_score(non_iterative_linearization), non_iterative_linearization)
             for non_iterative_linearization in nrp_linearizations.non_iterative),
            key=lambda x: x[0])

    best_iterative_score = float('-inf')
    best_iterative_linearizations = []
    for groups_linearizations in nrp_linearizations.iterative:
        split_score = 0
        split_linearizations = []
        for group_linearizations in groups_linearizations:
            group_score, group_linearization = max(((linearization_score(linearization), linearization)
                                                    for linearization in group_linearizations),
                                                   key=lambda x: x[0])
            split_linearizations.append(group_linearization)
            split_score += group_score

        if split_score > best_iterative_score:
            best_iterative_score = split_score
            best_iterative_linearizations = split_linearizations

    if best_non_iterative_score > best_iterative_score:
        return best_non_iterative_score, [best_non_iterative_linearization]
    else:
        return best_iterative_score, best_iterative_linearizations


def get_matches_for_hmm(detailed_hmm: DetailedHMM,
                        nrp_linearizations: Dict[Match_NRP_Variant_Info, NRP_Linearizations],
                        max_num_matches_per_bgc_variant: Optional[int],
                        use_anchors_heuristic: bool = False,
                        log=None) -> List[Match]:
    if log is not None:
        log.info(f'Processing BGC {detailed_hmm.bgc_variant.genome_id} variant {detailed_hmm.bgc_variant.variant_idx}')
        if use_anchors_heuristic:
            log.info("Using anchors heuristic for matching")

    hmm = detailed_hmm.to_hmm()
    matched_nrps_with_linearizations = []
    for nrp_info, nrp_linearizations in nrp_linearizations.items():
        score, linearizations = get_best_linearizations_for_nrp(hmm, nrp_linearizations,
                                                                detailed_hmm.hmm_helper.monomer_names_helper,
                                                                detailed_hmm,
                                                                use_anchors_heuristic=use_anchors_heuristic,
                                                                log=log)
        matched_nrps_with_linearizations.append((score, nrp_info, linearizations))

    best_matched_linearizations = sorted(matched_nrps_with_linearizations,
                                         key=lambda x: x[0],
                                         reverse=True)[:max_num_matches_per_bgc_variant]

    def get_alignments_for_linearizations(linearizations: List[Linearization]) -> List[Alignment]:
        return [detailed_hmm.get_alignment(linearization) for linearization in linearizations]

    bgc_variant_info = Match_BGC_Variant_Info.from_bgc_variant(detailed_hmm.bgc_variant)
    best_matches = [Match(bgc_variant_info=bgc_variant_info,
                          nrp_variant_info=nrp_info,
                          alignments=get_alignments_for_linearizations(linearizations),
                          normalized_score=score)
                    for score, nrp_info, linearizations in best_matched_linearizations]
    return best_matches


def get_matches(hmms: List[DetailedHMM],
                nrp_linearizations: Dict[Match_NRP_Variant_Info, NRP_Linearizations],
                max_num_matches_per_bgc_variant: Optional[int] = None,
                use_anchors_heuristic: bool = False,
                num_threads: int = 1,
                log=None) -> List[Match]:
    total_linearizations = 0
    for nrp_linearization in nrp_linearizations.values():
        total_linearizations += len(nrp_linearization.non_iterative)
        total_linearizations += sum(len(group) for groups in nrp_linearization.iterative for group in groups)

    if log is not None:
        log.info(f'Matching {len(hmms)} BGC variants against {total_linearizations} NRP linearizations')
        log.info(f'Using anchors heuristic: {use_anchors_heuristic}')

    # For safety, use only one thread if we have few BGCs
    actual_threads = min(num_threads, len(hmms)) if len(hmms) > 0 else 1

    # If only running on one thread, no need for parallelization
    if actual_threads == 1 or len(hmms) <= 1:
        matches = [item for hmm in hmms for item in get_matches_for_hmm(
            hmm, nrp_linearizations, max_num_matches_per_bgc_variant, use_anchors_heuristic, log
        )]
    else:
        # Using parallelization
        matches = chain(*Parallel(n_jobs=actual_threads)(
            delayed(get_matches_for_hmm)(hmm, nrp_linearizations, max_num_matches_per_bgc_variant,
                                         use_anchors_heuristic, log)
            for hmm in hmms
        ))

    return sorted(matches, key=lambda m: m.normalized_score, reverse=True)
