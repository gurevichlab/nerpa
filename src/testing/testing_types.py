from __future__ import annotations

from dataclasses import dataclass, asdict

import yaml
from enum import Enum
from itertools import chain
from pathlib import Path
from typing import Optional, Tuple, List

from src.antismash_parsing.antismash_parser_types import BGC_ID
from src.data_types import A_Domain_ID
from src.matching.alignment_type import Alignment
from src.matching.match_type import Match
from src.monomer_names_helper import NorineMonomerName
from src.write_results import write_yaml

SimplifiedAlignment = List[Tuple[Optional[A_Domain_ID], Optional[NorineMonomerName]]]

# ---- PyYAML helper to force flow style on selected sequences ----
class FlowSeq(list):
    """Marker for sequences that should be dumped in flow style (inline)."""
    pass

def _repr_flowseq(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)

yaml.add_representer(FlowSeq, _repr_flowseq)


def simplify_alignment(al: Alignment) -> SimplifiedAlignment:
    simplified_alignment = []
    for step in al:
        a_domain_id = (
            A_Domain_ID(step.bgc_module.gene_id, 
                        step.bgc_module.a_domain_idx)
            if step.bgc_module is not None
            else None
        )
            
        rban_monomer_name = (
            step.nrp_monomer.rban_name
            if step.nrp_monomer is not None
            else None
        )
        simplified_alignment.append((a_domain_id, rban_monomer_name))
        
    return simplified_alignment


def simplified_alignment_from_match(match: Match) -> SimplifiedAlignment:
    return list(chain.from_iterable(
        simplify_alignment(alignment)
        for alignment in match.alignments)
    )


def simplified_alignment_from_str(s: str) -> SimplifiedAlignment:
    match = Match.from_str(s)
    return simplified_alignment_from_match(match)


def check_simplified_alignments_equal(al1: SimplifiedAlignment,
                                      al2: SimplifiedAlignment) -> bool:
    only_matches1 = [(a_domain_id, rban_name) 
                     for (a_domain_id, rban_name) in al1 
                     if a_domain_id is not None and rban_name is not None]
    only_matches2 = [(a_domain_id, rban_name)
                     for (a_domain_id, rban_name) in al2 
                     if a_domain_id is not None and rban_name is not None]
    
    if len(only_matches1) != len(only_matches2):
        return False
    
    for (a_domain_id_1, rban_name_1), (a_domain_id_2, rban_name_2) in zip(al1, al2):
        rban_names_coinside = (rban_name_1 == rban_name_2) \
                or (rban_name_1[0] == 'X' and rban_name_2[0] == 'X')
        if a_domain_id_1 != a_domain_id_2 or not rban_names_coinside:
            return False
        
    return True


class TestResult(Enum):
    OK = 'OK'
    ACCEPTABLE_ALTERNATIVE_ALIGNMENT = 'ACCEPTABLE_ALTERNATIVE_ALIGNMENT'
    WRONG_ALIGNMENT = 'WRONG_ALIGNMENT'


@dataclass
class TestMatch:
    bgc_id: str
    nrp_id: str
    true_alignment: SimplifiedAlignment
    acceptable_alternative_alignments: List[SimplifiedAlignment]
    
    def test(self, match: Match) -> TestResult:
        match_bgc_id = match.bgc_variant_id.bgc_id.genome_id
        if match_bgc_id.endswith('.gbk'):
            match_bgc_id = match_bgc_id[:-4]
        if (match_bgc_id != self.bgc_id
                or match.nrp_variant_id.nrp_id != self.nrp_id):
            raise ValueError('BGC ID or NRP ID do not match the test case')
        
        simplified_alignment = simplified_alignment_from_match(match)
        if check_simplified_alignments_equal(simplified_alignment, self.true_alignment):
            return TestResult.OK
        
        for alt_al in self.acceptable_alternative_alignments:
            if check_simplified_alignments_equal(simplified_alignment, alt_al):
                return TestResult.ACCEPTABLE_ALTERNATIVE_ALIGNMENT
            
        return TestResult.WRONG_ALIGNMENT

    @classmethod
    def from_match(cls, match: Match) -> TestMatch:
        match_nrp_id = match.nrp_variant_id.nrp_id
        match_bgc_id = match.bgc_variant_id.bgc_id.genome_id
        if match_bgc_id.endswith('.gbk'):
            match_bgc_id = match_bgc_id[:-4]
        if match_bgc_id == 'converted_antiSMASH_v5_outputs':
            match_bgc_id = match_nrp_id.split('.')[0]

        return cls(bgc_id=match_bgc_id,
                   nrp_id=match.nrp_variant_id.nrp_id,
                   true_alignment=simplified_alignment_from_match(match),
                   acceptable_alternative_alignments=[])

    def _wrap_alignment(self, al: SimplifiedAlignment):
        wrapped = []
        for ad_id, monomer in al:
            first = FlowSeq(list(ad_id)) if ad_id is not None else None
            wrapped.append(FlowSeq([first, monomer]))
        return wrapped

    def to_yaml_obj(self):
        """Return a dict representation with flow sequences, ready for YAML dump."""
        data = asdict(self)
        data["true_alignment"] = self._wrap_alignment(self.true_alignment)
        data["acceptable_alternative_alignments"] = [
            self._wrap_alignment(al) for al in self.acceptable_alternative_alignments
        ]
        return data

    def to_yaml_str(self) -> str:
        """Dump a single TestMatch as a YAML string."""
        return yaml.dump(self.to_yaml_obj(), default_flow_style=False, sort_keys=False)

    # ---- Dump a list of TestMatch instances ----
    @classmethod
    def dump_list_to_str(cls, matches: List[TestMatch]) -> str:
        return yaml.dump(
            [m.to_yaml_obj() for m in matches],
            default_flow_style=False,
            sort_keys=False
        )

    @classmethod
    def from_yaml_dict(cls, data: dict) -> TestMatch:
        """Reconstruct a TestMatch instance from a YAML-decoded dictionary."""

        def unwrap_alignment(al):
            out = []
            for item in al:
                ad_id, monomer = item
                ad_id = tuple(ad_id) if ad_id is not None else None
                out.append((ad_id, monomer))
            return out

        true_alignment = unwrap_alignment(data["true_alignment"])
        acceptable_alternative_alignments = [
            unwrap_alignment(al) for al in data["acceptable_alternative_alignments"]
        ]

        return cls(
            bgc_id=data["bgc_id"],
            nrp_id=data["nrp_id"],
            true_alignment=true_alignment,
            acceptable_alternative_alignments=acceptable_alternative_alignments,
        )


def load_matches_from_txt(matches_txt: Path) -> List[Match]:
    matches_strs = matches_txt.read_text().split('\n\n')
    matches_strs = [match_str for match_str in matches_strs
                    if match_str.strip()]
    return [Match.from_str(matches_str)
            for matches_str in matches_strs]


if __name__ == '__main__':
    approved_matches_txt = Path('/home/ilianolhin/git/nerpa2/matches_inspection/approved_matches.txt')
    approved_matches = load_matches_from_txt(approved_matches_txt)
    tests = sorted([TestMatch.from_match(match)
                    for match in approved_matches],
                   key=lambda test: (test.bgc_id, test.nrp_id))
    test_ids = [(test.bgc_id, test.nrp_id) for test in tests]
    assert len(test_ids) == len(set(test_ids)), "Duplicate test cases found!"
    approved_matches_yaml = Path('/home/ilianolhin/git/nerpa2/data/for_training_and_testing/approved_matches.yaml')
    approved_matches_yaml.write_text(TestMatch.dump_list_to_str(tests))


