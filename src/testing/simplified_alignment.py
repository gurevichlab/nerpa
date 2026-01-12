from __future__ import annotations

import yaml
from itertools import chain, product, permutations
from typing import Optional, Tuple, List, NamedTuple, Iterable

from more_itertools import peekable

from src.antismash_parsing.bgc_variant_types import A_Domain_ID, BGC_Variant
from src.generic.combinatorics import split_into_iterations, is_subsequence, subseq_occurences, filter_unique
from src.matching.alignment_step_type import AlignmentStepLight
from src.matching.alignment_type import Alignment, AlignmentLight
from src.hmm.detailed_hmm import DetailedHMM
#from src.matching.detailed_hmm import DetailedHMM
from src.matching.match_type import Match
from src.monomer_names_helper import NorineMonomerName
from src.rban_parsing.get_linearizations import generate_linearizations_uniq
from src.rban_parsing.nrp_variant_types import NRP_Variant
from src.rban_parsing.rban_monomer import rBAN_Monomer, rban_name_no_unk_idx


# SimplifiedAlignmentStep = Tuple[Optional[A_Domain_ID], Optional[NorineMonomerName]]
class SimplifiedAlignmentStep(NamedTuple):
    a_domain_id: Optional[A_Domain_ID]
    rban_name: Optional[NorineMonomerName]

SimplifiedAlignment = List[SimplifiedAlignmentStep]


# ---- PyYAML helper to force flow style on selected sequences ----
class FlowSeq(list):
    """Marker for sequences that should be dumped in flow style (inline)."""
    pass


def _repr_flowseq(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)


yaml.add_representer(FlowSeq, _repr_flowseq)


def _wrap_alignment(al: SimplifiedAlignment):
    wrapped = []
    for ad_id, monomer in al:
        first = FlowSeq(list(ad_id)) if ad_id is not None else None
        wrapped.append(FlowSeq([first, monomer]))
    return wrapped


def simplified_alignment_to_str(al: SimplifiedAlignment) -> str:
    return yaml.dump(_wrap_alignment(al), default_flow_style=False, sort_keys=False)


def simplify_alignment(al: Alignment, index_offset: int = 0) -> SimplifiedAlignment:
    simplified_alignment = []
    for step in al:
        a_domain_id = (
            A_Domain_ID(step.bgc_module.gene_id,
                        step.bgc_module.a_domain_idx - index_offset)
            if step.bgc_module is not None
            else None
        )

        rban_monomer_name = (
            step.nrp_monomer.rban_name
            if step.nrp_monomer is not None
            else None
        )
        if a_domain_id is not None or rban_monomer_name is not None:
            simplified_alignment.append(SimplifiedAlignmentStep(a_domain_id, rban_monomer_name))

    return simplified_alignment


def simplified_alignment_from_match(match: Match, index_offset: int = 0) -> SimplifiedAlignment:
    return list(chain.from_iterable(
        simplify_alignment(alignment, index_offset=index_offset)
        for alignment in match.alignments)
    )


def simplified_alignment_from_str(s: str) -> SimplifiedAlignment:
    match = Match.from_str(s)
    return simplified_alignment_from_match(match)


def match_steps_coincide(step1: SimplifiedAlignmentStep,
                         step2: SimplifiedAlignmentStep,
                         ignore_chiralities: bool = True) -> bool:
    a_domain_id_1, rban_name_1 = step1
    a_domain_id_2, rban_name_2 = step2
    assert all([a_domain_id_1 is not None,
                a_domain_id_2 is not None,
                rban_name_1 is not None,
                rban_name_2 is not None]), \
        "Cannot compare steps with missing A-domain ID or rBAN name"

    if ignore_chiralities:
        if rban_name_1.startswith('D-'):
            rban_name_1 = rban_name_1[2:]
        if rban_name_2.startswith('D-'):
            rban_name_2 = rban_name_2[2:]

    rban_names_coinside = (rban_name_1 == rban_name_2) \
                          or (rban_name_1[0] == 'X' and rban_name_2[0] == 'X')
    return a_domain_id_1 == a_domain_id_2 and rban_names_coinside


# def fst_mismatched_step(al1: SimplifiedAlignment,
#                         al2: SimplifiedAlignment) -> Optional[
#     Tuple[Optional[SimplifiedAlignmentStep], Optional[SimplifiedAlignmentStep]]]:
#     only_matches1 = [(a_domain_id, rban_name)
#                      for (a_domain_id, rban_name) in al1
#                      if a_domain_id is not None and rban_name is not None]
#     only_matches2 = [(a_domain_id, rban_name)
#                      for (a_domain_id, rban_name) in al2
#                      if a_domain_id is not None and rban_name is not None]
#
#     fst_discrepancy = next(((step1, step2)
#                             for step1, step2 in zip(only_matches1, only_matches2)
#                             if not match_steps_coincide(step1, step2)),
#                            None)
#     if fst_discrepancy is not None:
#         return fst_discrepancy
#
#     if len(only_matches1) > len(only_matches2):
#         return (only_matches1[len(only_matches2)], None)
#     if len(only_matches2) > len(only_matches1):
#         return (None, only_matches2[len(only_matches1)])
#
#     return None

def fst_mismatched_step(al1: SimplifiedAlignment,
                        al2: SimplifiedAlignment) -> Optional[
    Tuple[Optional[SimplifiedAlignmentStep], Optional[SimplifiedAlignmentStep]]]:

    al1_matched_steps = [step for step in al1
                         if step.a_domain_id is not None and step.rban_name is not None]
    al2_matched_steps = [step for step in al2
                         if step.a_domain_id is not None and step.rban_name is not None]

    fst_discrepancy = next(((step1, step2)
                            for step1, step2 in zip(al1_matched_steps, al2_matched_steps)
                            if not match_steps_coincide(step1, step2)),
                           None)
    if fst_discrepancy is not None:
        return fst_discrepancy

    if len(al1_matched_steps) > len(al2_matched_steps):
        return (al1_matched_steps[len(al2_matched_steps)], None)
    if len(al2_matched_steps) > len(al1_matched_steps):
        return (None, al2_matched_steps[len(al1_matched_steps)])

    return None

def get_mismatched_steps(
        true_alignment: SimplifiedAlignment,
        other_alignment: SimplifiedAlignment
) -> List[Tuple[Optional[SimplifiedAlignmentStep], Optional[SimplifiedAlignmentStep]]]:
    true_dict = {
        a_domain_id: monomer_name
        for a_domain_id, monomer_name in true_alignment
        if a_domain_id is not None and monomer_name is not None
    }
    other_dict = {
        a_domain_id: monomer_name
        for a_domain_id, monomer_name in other_alignment
        if a_domain_id is not None and monomer_name is not None
    }

    all_a_domain_ids = set(true_dict.keys()) | set(other_dict.keys())
    discrepancies = []

    for a_domain_id in all_a_domain_ids:
        true_monomer = true_dict.get(a_domain_id)
        other_monomer = other_dict.get(a_domain_id)

        if true_monomer is not None and other_monomer is not None:
            # Both have mappings; check if they coincide
            monomers_coincide = (true_monomer == other_monomer) or (true_monomer[0] == 'X' and other_monomer[0] == 'X')
            if not monomers_coincide:
                discrepancies.append(((a_domain_id, true_monomer), (a_domain_id, other_monomer)))
        elif true_monomer is not None:
            # Only in true
            discrepancies.append(((a_domain_id, true_monomer), None))
        else:
            # Only in other
            discrepancies.append((None, (a_domain_id, other_monomer)))

    return discrepancies


def check_simplified_alignments_equal(al1: SimplifiedAlignment,
                                      al2: SimplifiedAlignment) -> bool:
    return fst_mismatched_step(al1, al2) is None

def split_into_bgc_iterations(
        alignment: SimplifiedAlignment
) -> List[SimplifiedAlignment]:
    '''
    TODO: currently all nrp skips are tossed
    because I don't know how to handle them to be consistent with NRP fragmentations.
    !!!!! That's bad -- the function removes information from the alignment !!!!!
    '''
    # print(f'Reconstructing alignment:\n{simplified_alignment_to_str(alignment)}')
    gene_ids = set(step.a_domain_id.gene_id
                   for step in alignment
                   if step.a_domain_id is not None)
    if len(gene_ids) == 1:
        key_fun = lambda step: step.a_domain_id.a_domain_idx
    else:
        key_fun = lambda step: step.a_domain_id.gene_id

    alignment_wo_none = [step for step in alignment
                         if step.a_domain_id is not None]

    bgc_iterations_wo_none = list(split_into_iterations(alignment_wo_none, key=key_fun))
    # print(f'Split into {len(bgc_iterations_wo_none)} BGC iterations (without NRP skips).')
    # for bgc_iteration in bgc_iterations_wo_none:
    #     print(simplified_alignment_to_str(bgc_iteration))
    #     print('---')

    bgc_iterations = []
    alignment_iter = peekable(alignment)

    # remove nrp skips at the very beginning, so that each bgc iteration starts with a_domain_id not None
    bgc_iteration = []
    while alignment_iter.peek().a_domain_id is None:
        step = next(alignment_iter)
        bgc_iteration.append(step)

    # reinsert nrp skips into bgc iterations
    for bgc_iteration_wo_none in bgc_iterations_wo_none:
        for step in bgc_iteration_wo_none:
            assert step == alignment_iter.peek(), \
                ("Alignment steps mismatch during BGC iteration splitting."
                 f" Expected: {step}, got: {alignment_iter.peek()}")
            bgc_iteration.append(next(alignment_iter))
            while True:
                step = alignment_iter.peek(None)
                if step is None or step.a_domain_id is not None:
                    break
                step = next(alignment_iter)
                bgc_iteration.append(step)

        bgc_iterations.append(bgc_iteration)
        bgc_iteration = []

    assert (alignment_iter.peek(None) is None), "Not all alignment steps were processed during BGC iteration splitting."

    # remove skips in between iterations to avoid ambiguity (I don't know which iteration they belong to)
    # by design, they are assigned to the end of the previous iteration
    for i in range(len(bgc_iterations) - 1):  # all but last
        while bgc_iterations[i][-1].a_domain_id is None:
            bgc_iterations[i].pop()

    return bgc_iterations


class StepWithIndexes(NamedTuple):
    module_idx: Optional[int]
    monomer_idx: Optional[int]

def make_alignment_with_indexes(
        simplified_alignment: SimplifiedAlignment,
        module_idxs_in_bgc: dict[A_Domain_ID, int],
        matched_mon_idxs: List[int]
) -> List[StepWithIndexes]:
    alignment_with_indexes = []
    mon_idx_iter = iter(matched_mon_idxs)
    for a_domain_id, rban_name in simplified_alignment:
        mon_idx = next(mon_idx_iter) if rban_name is not None else None
        module_idx = module_idxs_in_bgc[a_domain_id] if a_domain_id is not None else None
        alignment_with_indexes.append(StepWithIndexes(module_idx, mon_idx))

    return alignment_with_indexes


def simplified_alignment_to_light_alignment_one_iter(
        simplified_alignment: SimplifiedAlignment,
        bgc_variant: BGC_Variant,
        rban_monomers: List[rBAN_Monomer],
) -> Iterable[AlignmentLight]:
    assert len(split_into_bgc_iterations(simplified_alignment)) == 1, \
        "The provided simplified alignment contains multiple BGC iterations."

    bgc_aids = [A_Domain_ID.from_module(module) for module in bgc_variant.modules]
    nrp_monomer_names = [rban_name_no_unk_idx(monomer.rban_name)
                         for monomer in rban_monomers]

    alignment_aids = [a_domain_id
                      for a_domain_id, _ in simplified_alignment
                      if a_domain_id is not None]
    alignment_monomer_names = [rban_name_no_unk_idx(rban_name)
                               for _, rban_name in simplified_alignment
                               if rban_name is not None]

    if not is_subsequence(list(filter_unique(alignment_aids)),
                          bgc_aids):
        return
    if not is_subsequence(list(alignment_monomer_names),
                          nrp_monomer_names):
        return

    module_idxs_in_bgc = {a_id: i for i, a_id in enumerate(bgc_aids)}

    for matched_mon_idxs in subseq_occurences(alignment_monomer_names,
                                              nrp_monomer_names):
        alignment_with_indexes = make_alignment_with_indexes(simplified_alignment,
                                                             module_idxs_in_bgc,
                                                             matched_mon_idxs)

        light_alignment = []
        next_module_idx = 0
        next_monomer_idx = 0
        for step_mod_idx, step_mon_idx in alignment_with_indexes + [StepWithIndexes(len(bgc_aids), len(rban_monomers))]:
            # artificially adding end step to flush remaining modules/monomers

            # prepending modules missing in alignment
            if step_mod_idx is not None:
                light_alignment.extend(
                    AlignmentStepLight(a_domain_id=bgc_aids[module_idx],
                                       nrp_monomer=None)
                    for module_idx in range(next_module_idx, step_mod_idx)
                )
                next_module_idx = step_mod_idx + 1  # note: the step_mod_idx itself is not added yet

            # prepending monomers missing in alignment
            if step_mon_idx is not None:
                light_alignment.extend(
                    AlignmentStepLight(a_domain_id=None,
                                       nrp_monomer=rban_monomers[mon_idx])
                    for mon_idx in range(next_monomer_idx, step_mon_idx)
                )
                next_monomer_idx = step_mon_idx + 1  # note: the step_mon_idx itself is not added yet

            # adding alignment step
            if step_mod_idx != len(bgc_aids):  # skip the artificial end step
                light_alignment.append(
                    AlignmentStepLight(a_domain_id=(bgc_aids[step_mod_idx]
                                                    if step_mod_idx is not None else None),
                                       nrp_monomer=(rban_monomers[step_mon_idx]
                                                    if step_mon_idx is not None else None))
                )

        yield light_alignment


def _simplified_alignment_to_light_alignments(
        simplified_alignment: SimplifiedAlignment,
        hmm: DetailedHMM,
        nrp_variant: NRP_Variant,
        #ambiguity_handling: Literal['any', 'max', 'crash'] = 'crash',
) -> Iterable[Tuple[AlignmentLight, ...]]:

    if nrp_variant.nrp_variant_id.nrp_id == 'BGC0000985.4':
        pass
    alignments = split_into_bgc_iterations(simplified_alignment)
    nrp_linearizations_iter = generate_linearizations_uniq(nrp_variant,
                                                           num_fragment_groups=len(alignments))

    for linearizations in nrp_linearizations_iter:
        # match each alignment to corresponding linearization
        for perm_linearizations in permutations(linearizations):
            alignments_per_linearizations: Tuple[Iterable[AlignmentLight], ...] = tuple(
                simplified_alignment_to_light_alignment_one_iter(
                    alignment,
                    hmm.bgc_variant,
                    linearization
                )
                for alignment, linearization in zip(alignments, perm_linearizations)
            )

            yield from product(*alignments_per_linearizations)


def simplified_alignment_to_light_alignments(
        simplified_alignment: SimplifiedAlignment,
        hmm: DetailedHMM,
        nrp_variant: NRP_Variant,
        #ambiguity_handling: Literal['any', 'max', 'crash'] = 'crash',
) -> Iterable[Tuple[AlignmentLight, ...]]:
    def al_step_key(step: AlignmentStepLight) -> tuple:
        return (
            step.a_domain_id,
            step.nrp_monomer.to_base_mon() if step.nrp_monomer is not None else None
        )

    def alignments_key(alignments: Tuple[AlignmentLight, ...]) -> tuple:
        return tuple(
            sorted(
                tuple(al_step_key(step) for step in alignment)
                for alignment in alignments
            )
        )

    return filter_unique(
        _simplified_alignment_to_light_alignments(
            simplified_alignment,
            hmm,
            nrp_variant
        ),
        key=alignments_key
    )




