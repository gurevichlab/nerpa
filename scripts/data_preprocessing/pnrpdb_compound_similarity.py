from __future__ import annotations
import csv
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import (
    Callable,
    Dict,
    Hashable,
    List,
    NamedTuple,
    Tuple,
    Optional,
)

import yaml
import copy

from src.monomer_names_helper import Chirality
from src.rban_parsing.nrp_variant_types import NRP_Variant
from src.rban_parsing.rban_monomer import (
    rBAN_Monomer,
    rBAN_idx,
)
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.rban_parsing.retrieve_nrp_variants import (
    rban_records_to_nrp_variants,
)
from itertools import combinations
from enum import Enum
import argparse

def node_key_ignore_chr(mon: rBAN_Monomer) -> tuple:
    return (mon.residue, mon.methylated, mon.is_pks_hybrid)

class GraphType(Enum):
    RBAN = "rban"
    NERPA = "nerpa"

    def strictness_level(self) -> int:
        match self:
            case GraphType.NERPA:
                return 0
            case GraphType.RBAN:
                return 1

class ChiralityHandling(Enum):
    STRICT = "strict"

    # if both compounds have identified chiralities, STRICT is applied;
    # if at least one has only unknown chiralities, IGNORE is applied
    # motivation: not to discard matches when one compound was in isomeric SMILES format
    # and the other in non-isomeric SMILES format
    WEAK = "weak"  

    IGNORE = "ignore"

    def strictness_level(self) -> int:
        match self:
            case ChiralityHandling.IGNORE:
                return 0
            case ChiralityHandling.WEAK:
                return 1
            case ChiralityHandling.STRICT:
                return 2





class ComparisonMode(NamedTuple):
    graph_type: GraphType
    max_distance: int
    chirality_handling: ChiralityHandling

    @classmethod
    def list_all_modes(
            cls,
            graph_type: Optional[GraphType]=None
    ) -> List[ComparisonMode]:
        max_distances = (0, 1)
        graph_types = (
            [graph_type]
            if graph_type is not None
            else list(GraphType)
        )
        return [
            cls(graph_type=gt, max_distance=md, chirality_handling=ch)
            for gt in graph_types
            for md in max_distances
            for ch in ChiralityHandling
        ]

    def __str__(self) -> str:
        return f"gt-{self.graph_type.value}__d<={self.max_distance}__chr-{self.chirality_handling.value}"

    @classmethod
    def from_str(cls, s: str) -> ComparisonMode:
        parts = s.split("__")
        if len(parts) != 3:
            raise ValueError(f"Invalid ComparisonMode string: {s}")
        graph_type_str, max_distance_str, chirality_handling_str = parts
        graph_type = GraphType(graph_type_str.split("-")[1])
        max_distance = int(max_distance_str.split("<=")[1])
        chirality_handling = ChiralityHandling(chirality_handling_str.split("-")[1])
        return cls(
            graph_type=graph_type,
            max_distance=max_distance,
            chirality_handling=chirality_handling
        )

    def _strictness_tuple(self) -> Tuple[int, int, int]:
        return (
            self.graph_type.strictness_level(),
            self.max_distance,
            self.chirality_handling.strictness_level()
        )

    def is_stricter_than(self, other: ComparisonMode) -> bool:
        self_strictness_levels = self._strictness_tuple()
        other_strictness_levels = other._strictness_tuple()

        all_geq = all(
            self_level >= other_level
            for self_level, other_level in zip(self_strictness_levels,
                                               other_strictness_levels)
        )

        any_gr = any(
            self_level > other_level
            for self_level, other_level in zip(self_strictness_levels,
                                               other_strictness_levels)
        )
        return all_geq and any_gr


def all_unknown_chr(nrp: NRP_Variant) -> bool:
    return all(
        mon.chirality == Chirality.UNKNOWN
        for fragment in nrp.fragments
        for mon in fragment.monomers
    )

def tentative_subs_dist1(
        source: NRP_Variant,
        target: NRP_Variant,
        node_key: Optional[Callable[[rBAN_Monomer], Hashable]] = None,
) -> List[Tuple[rBAN_idx, rBAN_Monomer]]:
    """
    Returns the list of substitutions to be applied to source
    to make its monomer composition match that of target
    """
    if node_key is None:
        node_key = lambda mon: mon.to_base_mon()

    source_mons = [
        mon
        for fragment in source.fragments
        for mon in fragment.monomers
    ]
    target_mons = [
        mon
        for fragment in target.fragments
        for mon in fragment.monomers
    ]
    
    source_mon_keys_counter = Counter(
        node_key(mon)
        for mon in source_mons
    )
    target_mon_keys_counter = Counter(
        node_key(mon)
        for mon in target_mons
    )

    needs = target_mon_keys_counter - source_mon_keys_counter
    excess = source_mon_keys_counter - target_mon_keys_counter
    if sum(needs.values()) != 1 or sum(excess.values()) != 1:
        return []  # not one substitution away

    needed_key = next(iter(needs))
    excess_key = next(iter(excess))

    target_mon = next(
        mon
        for mon in target_mons
        if node_key(mon) == needed_key
    )

    source_indices = [
        mon.rban_idx
        for mon in source_mons
        if node_key(mon) == excess_key
    ]

    return [(idx, target_mon) for idx in source_indices]
    

def make_sub(
        nrp_variant: NRP_Variant,
        idx: rBAN_idx,
        new_mon: rBAN_Monomer
) -> NRP_Variant:
    """
    Returns a new NRP_Variant with the monomer at idx replaced by new_mon
    """
    new_variant = copy.deepcopy(nrp_variant)
    for fragment in new_variant.fragments:
        for i, mon in enumerate(fragment.monomers):
            if mon.rban_idx == idx:
                # don't change the rban_idx of the new monomer
                fragment.monomers[i] = new_mon._replace(rban_idx=idx)
                return new_variant

    raise ValueError(f"Monomer with rban_idx {idx} not found in NRP_Variant")
             

def one_substitution_away(
        source: NRP_Variant,
        target: NRP_Variant,
        node_key: Optional[Callable[[rBAN_Monomer], Hashable]] = None,
) -> bool:
    tentative_subs = tentative_subs_dist1(source, target, node_key=node_key)
    for idx, new_mon in tentative_subs:
        new_variant = make_sub(source, idx, new_mon)
        if new_variant.is_isomorphic_to(target, node_key=node_key):
            return True

    return False


def compare_nrp_variants(
        nrp1: NRP_Variant,
        nrp2: NRP_Variant
) -> Dict[ComparisonMode, bool]:
    # 1. max_distance = 0: check isomorphism
    isomorphic = nrp1.is_isomorphic_to(nrp2)
    isomorphic_ignore_chr = nrp1.is_isomorphic_to(
        nrp2,
        node_key=node_key_ignore_chr
    )
    isomorphic_weak_chr = (
        isomorphic_ignore_chr
        if (all_unknown_chr(nrp1) or all_unknown_chr(nrp2))
        else isomorphic
    )

    # 2. max_distance = 1: check if they are one substitution away
    one_sub_away = one_substitution_away(nrp1, nrp2)
    at_most_one_sub_away = one_sub_away or isomorphic

    one_sub_away_ignore_chr = one_substitution_away(
        nrp1, nrp2,
        node_key=node_key_ignore_chr
    )
    at_most_one_sub_away_ignore_chr = one_sub_away_ignore_chr or isomorphic_ignore_chr

    one_sub_away_weak_chr = (
        one_sub_away_ignore_chr
        if (all_unknown_chr(nrp1) or all_unknown_chr(nrp2))
        else one_sub_away
    )
    at_most_one_sub_away_weak_chr = one_sub_away_weak_chr or isomorphic_weak_chr

    return {
        ComparisonMode(GraphType.NERPA, 0, ChiralityHandling.STRICT): isomorphic,
        ComparisonMode(GraphType.NERPA, 0, ChiralityHandling.IGNORE): isomorphic_ignore_chr,
        ComparisonMode(GraphType.NERPA, 0, ChiralityHandling.WEAK): isomorphic_weak_chr,
        ComparisonMode(GraphType.NERPA, 1, ChiralityHandling.STRICT): at_most_one_sub_away,
        ComparisonMode(GraphType.NERPA, 1, ChiralityHandling.IGNORE): at_most_one_sub_away_ignore_chr,
        ComparisonMode(GraphType.NERPA, 1, ChiralityHandling.WEAK): at_most_one_sub_away_weak_chr,
    }


def parse_args(nerpa_dir: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="Compute compound similarity for pnrpdb2")
    default_preprocessed_nrps = (
        nerpa_dir
        / "data"
        / "input"
        / "preprocessed"
        / "pnrpdb2_parsed_rban_records.yaml"
    )
    parser.add_argument(
        "--preprocessed-nrps",
        dest="preprocessed_nrps",
        type=Path,
        default=default_preprocessed_nrps,
        help=f"Path to preprocessed rBAN records (default: {default_preprocessed_nrps})"
    )

    default_output_path = (
        nerpa_dir
        / "data"
        / "for_training_and_testing"
        / "pnrpdb2_compound_similarity.tsv"
    )
    parser.add_argument(
        "-o", "--out",
        dest="out",
        type=Path,
        default=default_output_path,
        help=f"Output path (default: {default_output_path})"
    )

    return parser.parse_args()


def sanity_check(
        nrp1: NRP_Variant,
        nrp2: NRP_Variant,
        sim_info: Dict[ComparisonMode, bool]
) -> None:
    """
    Check that if a stricter mode is True, then all less strict modes are also True
    """
    nrp1_id = nrp1.nrp_variant_id.nrp_id
    nrp2_id = nrp2.nrp_variant_id.nrp_id
    for mode1, mode2 in combinations(sim_info.keys(), 2):
        if mode1.is_stricter_than(mode2):
            if sim_info[mode1] and not sim_info[mode2]:
                raise ValueError(
                    f"Sanity check failed for {nrp1_id} vs {nrp2_id}: {mode1} is stricter than {mode2}, "
                    f"but {mode1} is True and {mode2} is False\n"
                    f"{nrp1_id}:\n"
                    f"{nrp1.to_dict()}\n\n"
                    f"{nrp2_id}:\n"
                    f"{nrp2.to_dict()}\n\n"
                )

def main():
    nerpa_dir = Path(__file__).resolve().parent.parent
    args = parse_args(nerpa_dir)

    print(f'Loading rBAN records variants from {args.preprocessed_nrps}...')
    rban_records = [
         Parsed_rBAN_Record.from_dict(record)
         for record in yaml.safe_load(args.preprocessed_nrps.open())
    ]
    nrp_variants = rban_records_to_nrp_variants(rban_records)
    nrp_variants.sort(key=lambda nrp: nrp.nrp_variant_id.nrp_id)

    print(f'Comparing {len(nrp_variants)} compounds...')

    cmp_modes = ComparisonMode.list_all_modes(GraphType.NERPA)

    fieldnames = (
        ["nrp1_id", "nrp2_id",]
        +
        [str(mode) for mode in cmp_modes]
    )

    rows = []
    for (nrp1, nrp2) in combinations(nrp_variants, 2):
        sim_info = compare_nrp_variants(nrp1, nrp2)

        # check that if a stricter mode is True, then all less strict modes are also True
        sanity_check(nrp1, nrp2, sim_info)

        if any(sim_info.values()):  # record only if there is at least one similarity
            rows.append({
                "nrp1_id": nrp1.nrp_variant_id.nrp_id,
                "nrp2_id": nrp2.nrp_variant_id.nrp_id,
                **{
                    str(mode): sim_info[mode]
                    for mode in cmp_modes
                }
            })

    with args.out.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with args.out.with_suffix(".tsv.readme").open("w") as f:
        f.write(
            "This file contains pairwise similarity information for "
            f"compounds in PNRPDB2 (compounds taken from {args.preprocessed_nrps}).\n"

            "For memory saving purposes the output is NOT symmetric "
            "(i.e. only pairs (nrp1_id, nrp2_id) with nrp1_id < nrp2_id are included) "
            "and NOT reflective (i.e. nrp1 is always different from nrp2).\n"

            "Comparison modes are encoded as follows:\n"
            "- gt-{graph_type}__d<={max_distance}__chr-{chirality_handling}\n"
            "where:\n"
            "-- graph_type: 'rban' or 'nerpa' --- "
            "which molecule representation was used for the comparison: "
            "rBAN monomer graph or Nerpa NRP_Variant\n"

            "-- max_distance: 0 or 1 --- "
            "maximum number of monomer substitutions allowed for two compounds to be considered similar\n"
            "-- chirality_handling: 'strict', 'ignore', or 'weak' --- "
            "how chiralities of monomers are treated in the comparison:\n"
            "  --- 'strict': chiralities must match\n"
            "  --- 'ignore': chiralities are ignored\n"
            "  --- 'weak': if both compounds have identified chiralities, "
            "'strict' is applied; if at least one has only unknown chiralities, "
            "'ignore' is applied.\n "
            "The motivation for 'weak' is to not discard matches when "
            "one compound was in isomeric SMILES format and the other "
            " in non-isomeric SMILES format.\n"
        )
    
    print(f"Wrote comparison results to {args.out}")


if __name__ == "__main__":
    main()
