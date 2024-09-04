from typing import (
    Dict,
    Iterable,
    List,
    NamedTuple,
    Set,
)
from src.monomer_names_helper import antiSMASH_MonomerName, MonomerResidue
from src.generic.string import hamming_distance
from training.specificity_prediction.scripts.training_table_generator.load_data import (
    AA_code,
    AA10_code,
    AA34_code,
    TrainingTableInputData
)
from dataclasses import dataclass


# note: the scores below are different only for "validated" substrates
class AA_scores(NamedTuple):
    score1: float  # keep true answers in the codes dictionary, i.e., both aa10 and aa34 scores are always 1.0
    score2: float  # exclude one copy of true answers but keep the duplicates (if any, i.e., num_sources > 1)
    score3: float  # exclude duplicates of the true answers but keep the equivalent codes from other entries (if any)
    score4: float  # exclude the equivalent codes, i.e., compute scores based on the closest analogue


class SVM_scores(NamedTuple):
    single_amino: float
    small_cluster: float
    large_cluster: float
    physicochemical_class: float


@dataclass
class NerpaTableEntry:
    substrate: MonomerResidue
    validated: bool
    num_sources: int
    aa10_scores: AA_scores
    aa34_scores: AA_scores
    svm_scores: SVM_scores

    def to_str(self, sep: str) -> str:
        return sep.join(map(str, [self.substrate, int(self.validated),
                                  *self.aa10_scores, *self.aa34_scores, *self.svm_scores,
                                  self.num_sources]))

    @classmethod
    def header(cls) -> List[str]:
        return ["substrate", "is_correct",
                "aa10_score1", "aa10_score2", "aa10_score3", "aa10_score4",
                "aa34_score1", "aa34_score2", "aa34_score3", "aa34_score4",
                "svm_single_amino_score", "svm_small_cluster_score",
                "svm_large_cluster_score", "svm_class_score",
                "num_sources"]


def get_aa_score(code: AA_code, all_codes: Iterable[AA_code]) -> float:
    def similarity_score(code1: AA_code, code2: AA_code) -> float:
        assert len(code1) == len(code2), "Codes should have the same length"
        return 1 - hamming_distance(code1, code2) / len(code1)
    return max(similarity_score(code, other_code) for other_code in all_codes)


def get_aa_scores(code: AA_code,
                  codes_counter: Dict[AA_code, int],
                  num_sources: int) -> AA_scores:
    known_codes = set(codes_counter.keys())
    known_codes_versions = [
        known_codes,
        known_codes - {code} if codes_counter[code] <= 1 else known_codes,
        known_codes - {code} if codes_counter[code] <= num_sources else known_codes,
        known_codes - {code}
    ]
    return AA_scores(*(get_aa_score(code, known_codes_version)
                       for known_codes_version in known_codes_versions))


def get_max_svm_scores(svm_scores: Dict[antiSMASH_MonomerName, SVM_scores],
                       core_svm_short_names: List[antiSMASH_MonomerName]) -> SVM_scores:
    missing_svm_name = next((svm_short_name for svm_short_name in core_svm_short_names
                             if svm_short_name not in svm_scores), None)
    assert missing_svm_name is None, f"{missing_svm_name} not in 'svm_scores'! Should be there by design"

    if not any(core_svm_name in svm_scores
               for core_svm_name in core_svm_short_names):
        return SVM_scores(-1.0, -1.0, -1.0, -1.0)  # default scores if no match  TODO: put -1 in config
    else:
        return SVM_scores(*(max(svm_scores[svm_short_name][scores_idx]
                                for svm_short_name in core_svm_short_names)
                            for scores_idx in range(4)))


def create_nerpa_table_entry(aa10_code: AA10_code,
                             aa34_code: AA34_code,
                             num_domains: int,
                             nerpa_monomer: MonomerResidue,
                             validated: bool,
                             svm_scores: Dict[antiSMASH_MonomerName, SVM_scores],
                             svm_short_names: List[antiSMASH_MonomerName],
                             nerpa_monomer_aa10_codes: Dict[AA10_code, int],
                             nerpa_monomer_aa34_codes: Dict[AA34_code, int]) -> NerpaTableEntry:
    aa10_scores = get_aa_scores(aa10_code, nerpa_monomer_aa10_codes, num_domains)
    aa34_scores = get_aa_scores(aa34_code, nerpa_monomer_aa34_codes, num_domains)

    if validated:
        assert aa10_scores.score1 == 1.0, "score1 for validated substrates should be 1.0 by design (aa10)"
        assert aa34_scores.score1 == 1.0, "score1 for validated substrates should be 1.0 by design (aa34)"
    else:
        aa10_scores = AA_scores(*((aa10_scores.score1,)*4))  # TODO: won't that happen automatically?
        aa34_scores = AA_scores(*((aa34_scores.score1,)*4))

    return NerpaTableEntry(nerpa_monomer,
                           validated,
                           num_domains,
                           aa10_scores,
                           aa34_scores,
                           get_max_svm_scores(svm_scores, svm_short_names))


def generate_nerpa_table_entries(input_data: TrainingTableInputData) -> List[NerpaTableEntry]:
    # Initialize the list of NerpaTableEntry objects
    nerpa_table_entries = []

    # Iterate through each row and generate NerpaTableEntry objects (exactly len(nerpa_supported_monomers) + 1 per row)
    for _, row in input_data.extended_signatures_df.iterrows():
        validated_nerpa_monomers: Set[MonomerResidue] = {input_data.short_name_to_nerpa_monomer[substrate]
                                                         for substrate in eval(row['substrates'])}

        # Collect SVM scores (ignoring the first four columns)
        svm_monomer_to_its_scores: Dict[antiSMASH_MonomerName, SVM_scores] = {col: eval(row[col])
                                                                              for col in input_data.extended_signatures_df.columns[4:]}

        missing_nerpa_monomer = next((nerpa_monomer
                                      for nerpa_monomer in input_data.nerpa_supported_monomers
                                      if nerpa_monomer not in input_data.nerpa_monomer_to_SVM_short_names),
                                     None)
        # assert missing_nerpa_monomer is None, f"{missing_nerpa_monomer} not in 'nerpa_monomer_to_SVM_short_names'! Should be there by design"
        # it's okay, nerpa_monomer_to_SVM_short_names is defaultdict(list) and will return an empty list for missing keys

        nerpa_table_entries.extend(
            create_nerpa_table_entry(row['aa10'], row['aa34'], row['num_domains'], nerpa_monomer,
                                     validated=nerpa_monomer in validated_nerpa_monomers,
                                     svm_scores=svm_monomer_to_its_scores,
                                     svm_short_names=input_data.nerpa_monomer_to_SVM_short_names[nerpa_monomer],
                                     nerpa_monomer_aa10_codes=input_data.substrate_to_aa10_codes[nerpa_monomer],
                                     nerpa_monomer_aa34_codes=input_data.substrate_to_aa34_codes[nerpa_monomer])
            for nerpa_monomer in input_data.nerpa_supported_monomers + [input_data.nerpa_unknown_monomer]
        )

    return nerpa_table_entries
