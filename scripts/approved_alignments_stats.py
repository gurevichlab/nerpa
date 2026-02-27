from pathlib import Path
from typing import List, Set

import yaml

from src.antismash_parsing.bgc_variant_types import A_Domain_ID
from src.testing.simplified_alignment import SimplifiedAlignmentStep
from src.testing.testing_types import TestMatch

def main():
    nerpa_dir = Path(__file__).resolve().parent.parent

    approved_matches_yaml = nerpa_dir / 'data/for_training_and_testing/approved_matches.yaml'
    with approved_matches_yaml.open() as f:
        approved_matches = [
            TestMatch.from_yaml_dict(match_dict)
            for match_dict in yaml.safe_load(f)
        ]

    def get_match_steps(match: TestMatch) -> List[SimplifiedAlignmentStep]:
        return [
            step
            for step in match.true_alignment
            if step.a_domain_id is not None and step.rban_name is not None
        ]

    def get_a_domain_ids(match: TestMatch) -> Set[A_Domain_ID]:
        return {
            step.a_domain_id
            for step in match.true_alignment
            if step.a_domain_id is not None
        }

    total_match_steps = sum(len(get_match_steps(match))
                            for match in approved_matches)
    total_a_domain_ids = sum(len(get_a_domain_ids(match))
                            for match in approved_matches)
    print(f"Total approved alignments: {len(approved_matches)}")
    print(f"Total match steps: {total_match_steps}")
    print(f"Total unique A-domain IDs: {total_a_domain_ids}")

    # BGC0000395.0 -- 1 skip in the beginning
    # BGC0000459.1 -- ins at start, 2 skips at end
    # 985.1
    # 1035.0
    # 1214.4
    # 1230.3
    # 1462.0
    # 1953.1 !
    # 2170.0
if __name__ == '__main__':
    main()