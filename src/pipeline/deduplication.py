from typing import List, Dict
from src.data_types import (
    BGC_Variant_ID,
    BGC_Variant,
    NRP_Variant,
    NRP_Variant_ID,
)


def cluster_isomorphic_nrp_variants(nrp_variants: List[NRP_Variant]) -> Dict[NRP_Variant_ID, NRP_Variant_ID]:
    variant_id_to_variant = {variant.nrp_variant_id: variant
                             for variant in nrp_variants}
    variant_id_to_repr_id: Dict[NRP_Variant_ID, NRP_Variant_ID] = {}

    for nrp_variant in nrp_variants:
        repr_id = next((repr_id
                        for repr_id in variant_id_to_repr_id.values()
                        if nrp_variant == variant_id_to_variant[repr_id]),
                       None)
        variant_id_to_repr_id[nrp_variant.nrp_variant_id] = nrp_variant.nrp_variant_id \
            if repr_id is None else repr_id

    return variant_id_to_repr_id


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
