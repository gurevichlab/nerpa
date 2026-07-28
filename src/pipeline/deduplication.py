from typing import List, Dict
from src.antismash_parsing.bgc_variant_types import (
    BGC_Variant_ID,
    BGC_Variant,
)
from src.rban_parsing.nrp_variant_types import (
    NRP_ID,
    NRP_Variant_ID,
    NRP_Variant,
    is_mibig_nrp,
    is_norine_nrp,
)
from collections import defaultdict


def reassign_representative_nrps(
        nrp_id_to_repr_id: Dict[NRP_ID, NRP_ID]
) -> Dict[NRP_ID, NRP_ID]:
    """
    Try to make MIBiG or NORINE variants the representatives of their isomorphic groups
    """
    repr_id_to_members: Dict[NRP_ID, List[NRP_ID]] = defaultdict(list)
    for nrp_id, repr_id in nrp_id_to_repr_id.items():
        repr_id_to_members[repr_id].append(nrp_id)

    reassigned: Dict[NRP_ID, NRP_ID] = {}

    for members in repr_id_to_members.values():
        members = sorted(members)
        mibig_members = [nrp_id for nrp_id in members if is_mibig_nrp(nrp_id)]
        norine_members = [nrp_id for nrp_id in members if is_norine_nrp(nrp_id)]

        if mibig_members:
            new_repr_id = mibig_members[0]
        elif norine_members:
            new_repr_id = norine_members[0]
        else:
            new_repr_id = members[0]

        for m in members:
            reassigned[m] = new_repr_id
        reassigned[new_repr_id] = new_repr_id

    return reassigned


def cluster_isomorphic_nrp_variants(nrp_variants: List[NRP_Variant]) -> Dict[NRP_Variant_ID, NRP_Variant_ID]:
    variant_id_to_variant = {
        variant.nrp_variant_id.nrp_id: variant
        for variant in nrp_variants
    }
    variant_id_to_repr_id: Dict[NRP_ID, NRP_ID] = {}

    for nrp_variant in nrp_variants:
        nrp_id = nrp_variant.nrp_variant_id.nrp_id
        repr_id = next((repr_id
                        for repr_id in variant_id_to_repr_id.values()
                        if nrp_variant.is_isomorphic_to(variant_id_to_variant[repr_id])),
                       None)
        variant_id_to_repr_id[nrp_id] = (
            nrp_id
            if repr_id is None
            else repr_id
        )

    # Reassign representatives: if possible, make MIBiG or NORINE variants the representatives of their isomorphic groups
    variant_id_to_repr_id = reassign_representative_nrps(variant_id_to_repr_id)

    # I still have NRP_Variant_ID and NRP_ID even though they are the same, so I convert back to NRP_Variant_ID for backwards compatibility
    return {
        NRP_Variant_ID(nrp_id): NRP_Variant_ID(repr_id)
        for nrp_id, repr_id in variant_id_to_repr_id.items()
    }


def cluster_isomorphic_bgc_variants(bgc_variants: List[BGC_Variant]) -> Dict[BGC_Variant_ID, BGC_Variant_ID]:
    variant_id_to_variant = {variant.bgc_variant_id: variant
                             for variant in bgc_variants}
    variant_id_to_repr_id: Dict[BGC_Variant_ID, BGC_Variant_ID] = {}

    for bgc_variant in bgc_variants:
        repr_id = next((repr_id
                        for repr_id in variant_id_to_repr_id.values()
                        if bgc_variant == variant_id_to_variant[repr_id]),
                       None)
        variant_id_to_repr_id[bgc_variant.bgc_variant_id] = bgc_variant.bgc_variant_id \
            if repr_id is None else repr_id

    return variant_id_to_repr_id
