from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List

from src.antismash_parsing.antismash_parser_types import GeneId
from src.antismash_parsing.bgc_variant_types import A_Domain_ID
from src.monomer_names_helper import NorineMonomerName
from src.testing.simplified_alignment import (
    SimplifiedAlignment,
    simplified_alignment_from_match,
)
from src.matching.match_type import Match


@dataclass
class SimplifiedMatch:
    bgc_id: str
    nrp_id: str
    alignment: SimplifiedAlignment


def nerpa2_match_to_simplified_match(match: Match,
                                     index_offset: int = 0) -> SimplifiedMatch:
    """
    Convert a Nerpa2 Match object to a SimplifiedMatch object.

    Args:
        match: The Nerpa2 Match object to convert.

    Returns:
        A SimplifiedMatch object representing the Nerpa2 match.
    """
    bgc_id = match.genome_id
    if bgc_id.endswith('.gbk'):
        bgc_id = bgc_id[:-4]
    if bgc_id == 'converted_antiSMASH_v5_outputs':
        bgc_id = match.nrp_variant_id.nrp_id.split('.')[0]

    return SimplifiedMatch(
        bgc_id=bgc_id,
        nrp_id=match.nrp_variant_id.nrp_id,
        alignment=simplified_alignment_from_match(match, index_offset=index_offset),
    )
