import argparse
import yaml
import pandas as pd
from pathlib import Path
from typing import (
    Dict,
    Iterable,
    List,
    Literal,
    NamedTuple,
    Optional,
    Tuple,
    Set,
    Union
)
from collections import Counter, defaultdict, OrderedDict
from src.monomer_names_helper import antiSMASH_MonomerName, MonomerResidue
from src.generic.string import hamming_distance
from dataclasses import dataclass

# Get the absolute path of the directory where the script is located
script_dir = Path(__file__).resolve().parent
SEPARATOR = '\t'
nerpa_table_column_headers = ["substrate", "is_correct",
                              "aa10_score1", "aa10_score2", "aa10_score3", "aa10_score4",
                              "aa34_score1", "aa34_score2", "aa34_score3", "aa34_score4",
                              "svm_single_amino_score", "svm_small_cluster_score",
                              "svm_large_cluster_score", "svm_class_score",
                              "num_sources"]  # see the fields in NerpaTableEntry


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Process monomer signature data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Input group
    input_group = parser.add_argument_group("input", "Input options")
    input_group.add_argument("--num-monomers", type=int, default=None,
                             help="Number of most frequent monomers to consider. By default all monomers are included.")
    input_group.add_argument("--unknown-monomer-name", type=str, default="Unknown",
                             help="Special name designating all rare (unsupported) monomers.")
    input_group.add_argument("--monomer-frequencies", type=str, default="data/core_frequency.tsv",
                             help="Path to a TSV file with all core monomers and their frequencies sorted from the most common to the rarest.")
    input_group.add_argument("--extended-signatures", type=str, default="data/extended_signatures.tsv",
                             help="Path to the extended signatures file.")
    input_group.add_argument("--monomers-table", type=str, default="data/monomers_unique.tsv",
                             help="Path to the monomer names conversion table.")

    # Output group
    output_group = parser.add_argument_group("output", "Output options")
    output_group.add_argument("--monomer-signatures-dict", type=str, default="data/monomer_signatures.yaml",
                              help="Path to a YAML file to store the dictionary with all aa10 and aa34 signatures per every supported core monomer.")
    output_group.add_argument("--nerpa-training-table", type=str, default="data/nerpa_scoring_table.tsv",
                              help="Path to a TSV file to store the main specificity training table.")

    return parser.parse_args()


def get_nerpa_supported_monomers(monomer_frequencies_tsv: Path,
                                 num_monomers: Optional[int]) -> List[MonomerResidue]:
    df = pd.read_csv(monomer_frequencies_tsv, sep=SEPARATOR)
    selected_monomers = df['core'].head(num_monomers) if num_monomers is not None else df['core']
    return list(selected_monomers)


def get_short_name_to_nerpa_monomer_mapping(monomers_table_tsv: Path,
                                            nerpa_supported_monomers: List[MonomerResidue],
                                            nerpa_unknown_monomer: MonomerResidue) -> Dict[antiSMASH_MonomerName, MonomerResidue]:
    df = pd.read_csv(monomers_table_tsv, sep=SEPARATOR)

    short_name_to_nerpa_monomer = defaultdict(lambda: nerpa_unknown_monomer)
    for _, row in df.iterrows():
        if row['core'] in nerpa_supported_monomers:
            short_name_to_nerpa_monomer[row['as_short']] = row['core']

    return short_name_to_nerpa_monomer


def get_nerpa_monomer_to_SVM_short_names_mapping(as_short_names_svm: List[antiSMASH_MonomerName],
                                                 short_name_to_nerpa_monomer: Dict[antiSMASH_MonomerName, MonomerResidue],
                                                 nerpa_unknown_monomer: MonomerResidue) -> Dict[MonomerResidue, List[antiSMASH_MonomerName]]:
    nerpa_monomer_to_SVM_short_names: Dict[MonomerResidue, List[antiSMASH_MonomerName]] = defaultdict(list)

    # Populate the dictionary
    for as_short_name in as_short_names_svm:
        # Map the SVM short name to the corresponding nerpa monomer
        nerpa_monomer = short_name_to_nerpa_monomer.get(as_short_name, nerpa_unknown_monomer)
        if nerpa_monomer != nerpa_unknown_monomer:  # !! see below
            nerpa_monomer_to_SVM_short_names[nerpa_monomer].append(as_short_name)

    # !! Think about me !!
    # !! Some SVM-supported monomers are rare (e.g. 'Pgl' is present only in 3 A-domains in the signatures file).
    # !! As the result, they could be not in the Nerpa core monomers and be converted to "Unknown" (special_monomer_name)
    # !! Without the "if" check above, this special monomer will map to only these few (e.g., 1) SVM rare monomers.
    # !! This could bias the downstream results -- we will treat this as if SVM predict some scores for all "Unknown"s,
    # !! not just a few supported by SVM and ignored by Nerpa. Let's be safe and just ignore SVM predictions for them.

    return nerpa_monomer_to_SVM_short_names


AA10_code = str
AA34_code = str
AA_code = Union[AA10_code, AA34_code]
SubstrateCodesDict = Dict[MonomerResidue, Dict[AA_code, int]]


def read_substrate_to_aa_codes(extended_signatures_tsv: Path,
                               short_name_to_nerpa_monomer: Dict[antiSMASH_MonomerName, MonomerResidue]) \
        -> Tuple[SubstrateCodesDict, SubstrateCodesDict]:
    df = pd.read_csv(extended_signatures_tsv, sep=SEPARATOR)

    # Initialize the defaultdicts for aa10 and aa34 codes
    substrate_to_aa10_codes = defaultdict(Counter)
    substrate_to_aa34_codes = defaultdict(Counter)

    # Iterate over each row in the DataFrame
    for _, row in df.iterrows():
        substrates: List[antiSMASH_MonomerName] = eval(row['substrates'])  # Evaluates the list from string to Python list

        nerpa_monomers = [short_name_to_nerpa_monomer[substrate] for substrate in substrates]

        aa10_code = row['aa10']
        aa34_code = row['aa34']

        num_domains = row['num_domains']

        # Update the dictionaries for each nerpa_monomer
        for nerpa_monomer in nerpa_monomers:
            substrate_to_aa10_codes[nerpa_monomer][aa10_code] += num_domains
            substrate_to_aa34_codes[nerpa_monomer][aa34_code] += num_domains

    return substrate_to_aa10_codes, substrate_to_aa34_codes


def get_sorted_nerpa_monomer_names(nerpa_supported_monomers: List[MonomerResidue],
                                   nerpa_unknown_monomer: MonomerResidue) -> List[MonomerResidue]:
    return nerpa_supported_monomers + [nerpa_unknown_monomer]


def write_monomer_signatures(monomer_signatures_yaml: Path,
                             substrate_to_aa10_codes: SubstrateCodesDict,
                             substrate_to_aa34_codes: SubstrateCodesDict,
                             nerpa_supported_monomers: List[MonomerResidue],
                             nerpa_unknown_monomer: MonomerResidue):
    # Reformatting the two internal dicts into one codes_dict with lexicographical order of monomers
    codes_dict: Dict[MonomerResidue, Tuple[List[AA10_code], List[AA34_code]]] = OrderedDict()

    for substrate in nerpa_supported_monomers + [nerpa_unknown_monomer]:
        aa10_list = list(substrate_to_aa10_codes[substrate].keys())
        aa34_list = list(substrate_to_aa34_codes[substrate].keys())
        codes_dict[substrate] = (aa10_list, aa34_list)

    # Dumping codes_dict to YAML file
    with open(monomer_signatures_yaml, 'w') as yaml_file:
        yaml.dump(codes_dict, yaml_file, default_flow_style=False)

    print(f"Saved monomer aa10 and aa34 codes in {monomer_signatures_yaml}")


def ensure_output_directories_exist(output_paths: List[Path]):
    for path in output_paths:
        # Convert the file path to a Path object
        parent_dir = Path(path).parent

        # Create the directory (and any missing parents) if it doesn't exist
        if not parent_dir.exists():
            parent_dir.mkdir(parents=True, exist_ok=True)
            print(f"Created directory: {parent_dir}")
        # else:
        #     print(f"Directory already exists: {parent_dir}")


# note: the scores below are different only for "validated" substrates
# score1 -- keep true answers in the codes dictionary, i.e., both aa10 and aa34 scores are always 1.0
# score2 -- exclude one copy of true answers but keep the duplicates (if any, i.e., num_sources > 1)
# score3 -- exclude duplicates of the true answers but keep the equivalent codes from other entries (if any)
# score4 -- exclude the equivalent codes, i.e., compute scores based on the closest analogue
class AA_scores(NamedTuple):
    score1: float
    score2: float
    score3: float
    score4: float


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
    aa10_scores: AA_scores = AA_scores(0.0, 0.0, 0.0, 0.0)
    aa34_scores: AA_scores = AA_scores(0.0, 0.0, 0.0, 0.0)
    svm_scores: SVM_scores = SVM_scores(-1.0, -1.0, -1.0, -1.0)

    def __str__(self):
        return SEPARATOR.join(map(str, [self.substrate, int(self.validated),
                                        *self.aa10_scores, *self.aa34_scores, *self.svm_scores,
                                        self.num_sources]))


def get_aa_score(code: AA_code, all_codes: Iterable[AA_code]) -> float:
    def similarity_score(code: AA_code, other_code: AA_code) -> float:
        assert len(code) == len(other_code), "Codes should have the same length"
        return 1 - hamming_distance(code, other_code) / len(code)
    return max(similarity_score(code, other_code) for other_code in all_codes)


def get_aa_scores(code: AA_code,
                  codes_counter: Dict[AA_code, int],
                  num_sources: int) -> AA_scores:
    known_codes = set(codes_counter.keys())
    score1 = get_aa_score(code, known_codes)
    score2 = get_aa_score(code, known_codes if codes_counter[code] <= 1 else known_codes - {code})
    score3 = get_aa_score(code, known_codes if codes_counter[code] <= num_sources else known_codes - {code})
    score4 = get_aa_score(code, known_codes - {code})
    return AA_scores(score1, score2, score3, score4)


def get_svm_max_scores(svm_scores: Dict[antiSMASH_MonomerName, SVM_scores],
                       core_svm_short_names: List[antiSMASH_MonomerName]) -> SVM_scores:
    missing_svm_name = next((svm_short_name for svm_short_name in core_svm_short_names
                             if svm_short_name not in svm_scores), None)
    assert missing_svm_name is None, f"{missing_svm_name} not in 'svm_scores'! Should be there by design"

    if not any(core_svm_name in svm_scores
               for core_svm_name in core_svm_short_names):
        return SVM_scores(-1.0, -1.0, -1.0, -1.0)  # default scores if no match
    else:
        return SVM_scores(*(max(svm_scores[svm_short_name][scores_idx]
                                for svm_short_name in core_svm_short_names)
                            for scores_idx in range(4)))


def generate_nerpa_table_entries(extended_signatures_file: Path,
                                 nerpa_supported_monomers: List[MonomerResidue],
                                 nerpa_unknown_monomer: MonomerResidue,
                                 short_name_to_nerpa_monomer: Dict[antiSMASH_MonomerName, MonomerResidue],
                                 nerpa_monomer_to_SVM_short_names: Dict[MonomerResidue, List[antiSMASH_MonomerName]],
                                 substrate_to_aa10_codes: SubstrateCodesDict,
                                 substrate_to_aa34_codes: SubstrateCodesDict) -> List[NerpaTableEntry]:


    def _create_nerpa_table_entry(row,
                                  nerpa_monomer: MonomerResidue,
                                  validated: bool,
                                  svm_scores: Dict[antiSMASH_MonomerName, SVM_scores]) -> NerpaTableEntry:
        entry = NerpaTableEntry(nerpa_monomer, validated, row['num_domains'])

        # Compute and set aa10 and aa34 scores
        aa10_scores = get_aa_scores(row['aa10'], substrate_to_aa10_codes[nerpa_monomer])
        aa34_score1 = get_aa_scores(row['aa34'], substrate_to_aa34_codes[nerpa_monomer])

        if validated:
            assert aa10_score1 == 1.0, "score1 for validated substrates should be 1.0 by design (aa10)"
            assert aa34_score1 == 1.0, "score1 for validated substrates should be 1.0 by design (aa34)"
            entry.aa10_scores.score1 = aa10_score1
            entry.aa34_scores.score1 = aa34_score1
            _compute_and_set_scores(row['aa10'], substrate_to_aa10_codes[nerpa_monomer], entry, 'aa10')
            _compute_and_set_scores(row['aa34'], substrate_to_aa34_codes[nerpa_monomer], entry, 'aa34')
        else:
            entry.aa10_scores =
            _set_all_scores(entry, aa10_score1, 'aa10')
            _set_all_scores(entry, aa34_score1, 'aa34')

        # Set SVM scores
        assert nerpa_monomer in nerpa_monomer_to_SVM_short_names, f"{nerpa_monomer} is not in nerpa_monomer_to_SVM_short_names!"
        svm_short_names = nerpa_monomer_to_SVM_short_names[nerpa_monomer]
        max_scores = _get_max_scores(svm_scores, svm_short_names)
        (entry.svm_single_amino_score, entry.svm_small_cluster_score,
         entry.svm_large_cluster_score, entry.svm_class_score) = max_scores

        return entry


    # Load the extended signatures
    df = pd.read_csv(extended_signatures_file, sep=SEPARATOR)

    # Initialize the list of NerpaTableEntry objects
    nerpa_table_entries = []

    # Iterate through each row and generate NerpaTableEntry objects (exactly len(sorted_nerpa_monomer_names) per row)
    for _, row in df.iterrows():
        validated_nerpa_monomers: Set[MonomerResidue] = {short_name_to_nerpa_monomer[substrate]
                                                         for substrate in eval(row['substrates'])}

        # Collect SVM scores (ignoring the first four columns)
        svm_monomer_to_its_scores: Dict[antiSMASH_MonomerName, SVM_scores] = {col: eval(row[col])
                                                                              for col in df.columns[4:]}

        nerpa_table_entries.extend(
            _create_nerpa_table_entry(row, nerpa_monomer,
                                      validated=nerpa_monomer in validated_nerpa_monomers,
                                      svm_scores=svm_monomer_to_its_scores)
            for nerpa_monomer in nerpa_supported_monomers + [nerpa_unknown_monomer]  # ? why not just nerpa_supported_monomers?
        )

    return nerpa_table_entries


def write_nerpa_table(nerpa_training_table_path, nerpa_table_entries):
    with open(nerpa_training_table_path, 'w') as f:
        f.write(SEPARATOR.join(nerpa_table_column_headers) + "\n")
        f.write("\n".join(map(str, nerpa_table_entries)))
    print(f"Saved the Nerpa training table in {nerpa_training_table_path}")


def main():
    # Parse arguments
    args = parse_args()

    # reading
    nerpa_supported_monomers = get_nerpa_supported_monomers(args.monomer_frequencies, args.num_monomers)
    short_name_to_nerpa_monomer = get_short_name_to_nerpa_monomer_mapping(args.monomers_table, nerpa_supported_monomers,
                                                                          args.special_monomer_name)
    nerpa_monomer_to_SVM_short_names = get_nerpa_monomer_to_SVM_short_names_mapping(args.extended_signatures,
                                                                                    short_name_to_nerpa_monomer,
                                                                                    nerpa_supported_monomers,
                                                                                    args.special_monomer_name)

    # reading and processing
    sorted_nerpa_monomer_names = get_sorted_nerpa_monomer_names(nerpa_supported_monomers, args.special_monomer_name)
    substrate_to_aa10_codes, substrate_to_aa34_codes = read_substrate_to_aa_codes(args.extended_signatures,
                                                                                  short_name_to_nerpa_monomer)
    nerpa_table_entries = generate_nerpa_table_entries(args.extended_signatures, sorted_nerpa_monomer_names,
                                                       short_name_to_nerpa_monomer, nerpa_monomer_to_SVM_short_names,
                                                       substrate_to_aa10_codes, substrate_to_aa34_codes)

    # writing
    ensure_output_directories_exist([args.monomer_signatures_dict, args.nerpa_training_table])
    write_monomer_signatures(args.monomer_signatures_dict, substrate_to_aa10_codes, substrate_to_aa34_codes,
                             sorted_nerpa_monomer_names)
    write_nerpa_table(args.nerpa_training_table, nerpa_table_entries)


if __name__ == "__main__":
    main()
