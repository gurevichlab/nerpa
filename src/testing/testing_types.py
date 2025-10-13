from __future__ import annotations

from dataclasses import dataclass, asdict

import yaml
from enum import Enum
from itertools import chain
from pathlib import Path
from typing import Optional, Tuple, List

from src.antismash_parsing.bgc_variant_types import A_Domain_ID
from src.matching.alignment_type import Alignment
from src.matching.match_type import Match
from src.monomer_names_helper import NorineMonomerName
from src.testing.simplified_alignment import (
    SimplifiedAlignment,
    simplified_alignment_from_match,
    check_simplified_alignments_equal,
    _wrap_alignment,
)
class TestResult(Enum):
    CORRECT = 'CORRECT'
    ACCEPTABLE_ALTERNATIVE = 'ACCEPTABLE_ALTERNATIVE'
    WRONG = 'WRONG'


@dataclass
class TestMatch:
    bgc_id: str
    nrp_id: str
    true_alignment: SimplifiedAlignment
    acceptable_alternative_alignments: List[SimplifiedAlignment]
    
    def test(self, match: Match) -> TestResult:
        match_bgc_id = match.genome_id
        if match_bgc_id.endswith('.gbk'):
            match_bgc_id = match_bgc_id[:-4]
        if (match_bgc_id != self.bgc_id
                or match.nrp_variant_id.nrp_id != self.nrp_id):
            raise ValueError('BGC ID or NRP ID do not match the test case')
        
        simplified_alignment = simplified_alignment_from_match(match)
        if check_simplified_alignments_equal(simplified_alignment, self.true_alignment):
            return TestResult.CORRECT
        
        for alt_al in self.acceptable_alternative_alignments:
            if check_simplified_alignments_equal(simplified_alignment, alt_al):
                return TestResult.ACCEPTABLE_ALTERNATIVE
            
        return TestResult.WRONG

    @classmethod
    def from_match(cls, match: Match) -> TestMatch:
        match_nrp_id = match.nrp_variant_id.nrp_id
        match_bgc_id = match.genome_id
        if match_bgc_id.endswith('.gbk'):
            match_bgc_id = match_bgc_id[:-4]
        if match_bgc_id == 'converted_antiSMASH_v5_outputs':
            match_bgc_id = match_nrp_id.split('.')[0]

        if min(
            step.bgc_module.a_domain_idx
            for alignment in match.alignments
            for step in alignment
            if step.bgc_module is not None
        ) != 0:
            raise ValueError('The match alignment does not start with A-domain index 0; '
                             'cannot create a TestMatch from it.')

        return cls(bgc_id=match_bgc_id,
                   nrp_id=match.nrp_variant_id.nrp_id,
                   true_alignment=simplified_alignment_from_match(match),
                   acceptable_alternative_alignments=[])


    def to_yaml_obj(self):
        """Return a dict representation with flow sequences, ready for YAML dump."""
        data = asdict(self)
        data["true_alignment"] = _wrap_alignment(self.true_alignment)
        data["acceptable_alternative_alignments"] = [
            _wrap_alignment(al) for al in self.acceptable_alternative_alignments
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