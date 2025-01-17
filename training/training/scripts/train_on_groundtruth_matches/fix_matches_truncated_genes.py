from typing import List, Dict, Tuple
from src.matching.matching_types_match import Match
from src.data_types import BGC_Variant, BGC_Module
from fix_match import is_subsequence, filter_unique, bgc_variant_match_compatible


def bgc_module_loose_key(module: BGC_Module) -> Tuple[str, int]:
    return module.aa34_code, module.a_domain_idx
def find_relaxed_compatible_bgc_variants(
        match: Match, bgc_variants: List[BGC_Variant]
) -> List[BGC_Variant]:
    """
    Finds BGC variants compatible with the match, using aa34_code for alignment
    and ignoring gene names. Ensures all alignments fit in order.
    """
    compatible_variants = []
    for bgc_variant in bgc_variants:
        bgc_codes = [bgc_module_loose_key(module)
                     for module in bgc_variant.modules]

        # Check that all alignments in the match fit into the BGC variant
        all_alignments_fit = all(
            is_subsequence(
                list(filter_unique((step.bgc_module.aa34_code, step.bgc_module.a_domain_idx)
                                    for step in alignment if step.bgc_module)),
                bgc_codes
            )
            for alignment in match.alignments
        )

        if all_alignments_fit:
            compatible_variants.append(bgc_variant)

    return compatible_variants


def update_gene_names_using_aa34(
        match: Match, bgc_variant: BGC_Variant
) -> Match:
    """
    Updates gene names in the match's alignments using aa34_code to identify corresponding modules.
    Raises an exception if any aa34_code in the match is not found in the BGC variant.
    """
    # Map aa34_code to gene_id from the BGC variant
    code_to_gene_map = {
        bgc_module_loose_key(module): module.gene_id for module in bgc_variant.modules
    }

    # Update gene names in the match alignments
    for alignment in match.alignments:
        for step in alignment:
            if step.bgc_module:
                module_code = (step.bgc_module.aa34_code, step.bgc_module.a_domain_idx)
                if module_code not in code_to_gene_map:
                    raise ValueError(
                        f"module_code '{module_code}' in match is not present in the BGC variant."
                    )
                step.bgc_module = step.bgc_module._replace(
                    gene_id=code_to_gene_map[module_code]
                )
    return match


def resolve_ambiguity(variants: List[BGC_Variant], match: Match) -> BGC_Variant:
    """
    Resolves ambiguity when multiple BGC variants are compatible.
    """
    print(f"Ambiguity detected for match {match.nrp_variant_info.nrp_id}:")
    for variant in variants:
        print(f" - Compatible BGC Variant ID: {variant.variant_idx}")
    return variants[0]  # Select the first variant as a fallback


def fix_approved_matches(
        approved_matches: List[Match], bgc_variants_dict: Dict[str, List[BGC_Variant]]
) -> List[Match]:
    """
    Fixes gene names in approved matches using compatible BGC variants.
    """
    fixed_matches = []

    for match in approved_matches:
        nrp_id = match.nrp_variant_info.nrp_id
        bgc_variants = bgc_variants_dict.get(nrp_id, [])

        # Find compatible variants using relaxed compatibility (all alignments must fit)
        if any(bgc_variant_match_compatible(bgc_variant, match)
               for bgc_variant in bgc_variants):
            fixed_matches.append(match)
            continue

        loosely_compatible_variants = find_relaxed_compatible_bgc_variants(match, bgc_variants)

        if not loosely_compatible_variants:
            print(f"No loosely compatible variant found for match: {nrp_id}")
            raise

        # Handle ambiguities
        if len(loosely_compatible_variants) > 1:
            chosen_variant = resolve_ambiguity(loosely_compatible_variants, match)
        else:
            chosen_variant = loosely_compatible_variants[0]

        # Update gene names using the chosen BGC variant
        fixed_match = update_gene_names_using_aa34(match, chosen_variant)
        fixed_matches.append(fixed_match)

    return fixed_matches
