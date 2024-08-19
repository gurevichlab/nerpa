import argparse
import yaml
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple
from collections import Counter, defaultdict, OrderedDict


# Get the absolute path of the directory where the script is located
script_dir = Path(__file__).resolve().parent
SEPARATOR = '\t'
nerpa_table_column_headers = ["substrate", "is_correct",
                              "aa10_score1", "aa10_score2", "aa10_score3", "aa10_score4",
                              "aa34_score1", "aa34_score2", "aa34_score3", "aa34_score4",
                              "svm_single_amino_score", "svm_small_cluster_score",
                              "svm_large_cluster_score", "svm_class_score",
                              "num_sources"]  # see the fields in NerpaTableEntry


def parse_args():
    parser = argparse.ArgumentParser(
        description="Process monomer signature data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Input group
    input_group = parser.add_argument_group("input", "Input options")
    input_group.add_argument("--num-monomers", type=int, default=-1,
                             help="Number of most frequent monomers to consider, '-1' for all.")
    input_group.add_argument("--special-monomer-name", type=str, default="Unknown",
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


def get_nerpa_supported_monomers(monomer_frequencies_path, num_monomers):
    # Read the TSV file into a DataFrame
    df = pd.read_csv(monomer_frequencies_path, sep=SEPARATOR)

    # Select the first 'num_monomers' values from the first column
    if num_monomers == -1:
        # If num_monomers is -1, select all values
        selected_monomers = df['core']
    else:
        selected_monomers = df['core'].head(num_monomers)

    # Convert the selected values to a set
    return set(selected_monomers)


def get_short_name_to_nerpa_monomer_mapping(monomers_table_path, nerpa_supported_monomers, special_monomer_name):
    # Read the TSV file into a DataFrame
    df = pd.read_csv(monomers_table_path, sep=SEPARATOR)

    # Create the dictionary to store the mappings
    short_name_to_nerpa_monomer = {}

    # Iterate over each row in the DataFrame
    for _, row in df.iterrows():
        short_name = row['as_short']
        core_monomer = row['core']

        # Map to the core monomer if it's in nerpa_supported_monomers, otherwise use special_monomer_name
        if core_monomer in nerpa_supported_monomers:
            short_name_to_nerpa_monomer[short_name] = core_monomer
        else:
            short_name_to_nerpa_monomer[short_name] = special_monomer_name

    return short_name_to_nerpa_monomer


def get_nerpa_monomer_to_SVM_short_names_mapping(extended_signatures_path, short_name_to_nerpa_monomer,
                                                 nerpa_supported_monomers, special_monomer_name):
    # Read the TSV file into a DataFrame
    df = pd.read_csv(extended_signatures_path, sep=SEPARATOR)

    # Extract the relevant column names (ignoring the first 4 columns)
    svm_short_names = df.columns[4:]

    # Initialize the dictionary with an empty list for each nerpa_supported_monomer
    nerpa_monomer_to_SVM_short_names = {monomer: [] for monomer in nerpa_supported_monomers}
    nerpa_monomer_to_SVM_short_names[special_monomer_name] = []

    # Populate the dictionary
    for short_name in svm_short_names:
        # Map the SVM short name to the corresponding nerpa monomer
        nerpa_monomer = short_name_to_nerpa_monomer.get(short_name, special_monomer_name)
        nerpa_monomer_to_SVM_short_names[nerpa_monomer].append(short_name)

    # !! Think about me !!
    # !! Some SVM-supported monomers are rare (e.g. 'Pgl' is present only in 3 A-domains in the signatures file).
    # !! As the result, they could be not in the Nerpa core monomers and be converted to "Unknown" (special_monomer_name)
    # !! Without the command below, this special monomer will map to only these few (e.g., 1) SVM rare monomers.
    # !! This could bias the downstream results -- we will treat this as if SVM predict some scores for all "Unknown"s,
    # !! not just a few supported by SVM and ignored by Nerpa. Let's be safe and just ignore SVM predictions for them.
    nerpa_monomer_to_SVM_short_names[special_monomer_name] = []

    return nerpa_monomer_to_SVM_short_names


def read_substrate_to_aa_codes(extended_signatures_path, short_name_to_nerpa_monomer):
    # Read the TSV file into a DataFrame
    df = pd.read_csv(extended_signatures_path, sep=SEPARATOR)

    # Initialize the defaultdicts for aa10 and aa34 codes
    substrate_to_aa10_codes = defaultdict(Counter)
    substrate_to_aa34_codes = defaultdict(Counter)

    # Iterate over each row in the DataFrame
    for _, row in df.iterrows():
        # Convert substrates column to a list of nerpa_monomer names
        substrates = eval(row['substrates'])  # Evaluates the list from string to Python list

        # Map substrates to nerpa_monomer names
        nerpa_monomers = [short_name_to_nerpa_monomer.get(substrate, substrate) for substrate in substrates]

        # Get the aa10 and aa34 codes
        aa10_code = row['aa10']
        aa34_code = row['aa34']

        # Get the num_domains count
        num_domains = row['num_domains']

        # Update the dictionaries for each nerpa_monomer
        for nerpa_monomer in nerpa_monomers:
            substrate_to_aa10_codes[nerpa_monomer][aa10_code] += num_domains
            substrate_to_aa34_codes[nerpa_monomer][aa34_code] += num_domains

    return substrate_to_aa10_codes, substrate_to_aa34_codes


def get_sorted_nerpa_monomer_names(nerpa_supported_monomers, special_monomer_name):
    def _custom_sort(substrate):
        return substrate == special_monomer_name, substrate  # to place `special_monomer_name` the last in the list

    return sorted(nerpa_supported_monomers.union([special_monomer_name]), key=_custom_sort)


def write_monomer_signatures(monomer_signatures_dict_path, substrate_to_aa10_codes, substrate_to_aa34_codes,
                             sorted_nerpa_monomer_names):
    # Reformatting the two internal dicts into one codes_dict with lexicographical order of monomers
    codes_dict: Dict[str, Tuple[List[str], List[str]]] = OrderedDict()

    for substrate in sorted_nerpa_monomer_names:
        aa10_list = list(substrate_to_aa10_codes[substrate].keys())
        aa34_list = list(substrate_to_aa34_codes[substrate].keys())
        codes_dict[substrate] = (aa10_list, aa34_list)

    # Dumping codes_dict to YAML file
    with open(monomer_signatures_dict_path, 'w') as yaml_file:
        yaml.dump(codes_dict, yaml_file, default_flow_style=False)

    print(f"Saved monomer aa10 and aa34 codes in {monomer_signatures_dict_path}")


def ensure_output_directories_exist(output_paths):
    for path in output_paths:
        # Convert the file path to a Path object
        parent_dir = Path(path).parent

        # Create the directory (and any missing parents) if it doesn't exist
        if not parent_dir.exists():
            parent_dir.mkdir(parents=True, exist_ok=True)
            print(f"Created directory: {parent_dir}")
        # else:
        #     print(f"Directory already exists: {parent_dir}")


class NerpaTableEntry:
    def __init__(self, substrate, validated, num_sources):
        self.substrate = substrate
        self.validated = validated
        # note: the scores below are different only for "validated" substrates
        # score1 -- keep true answers in the codes dictionary, i.e., both aa10 and aa34 scores are always 1.0
        # score2 -- exclude one copy of true answers but keep the duplicates (if any, i.e., num_sources > 1)
        # score3 -- exclude duplicates of the true answers but keep the equivalent codes from other entries (if any)
        # score4 -- exclude the equivalent codes, i.e., compute scores based on the closest analogue
        self.aa10_score1 = 0.0
        self.aa10_score2 = 0.0
        self.aa10_score3 = 0.0
        self.aa10_score4 = 0.0
        self.aa34_score1 = 0.0
        self.aa34_score2 = 0.0
        self.aa34_score3 = 0.0
        self.aa34_score4 = 0.0

        # for substrates not supported in SVM
        self.svm_single_amino_score = -1.0
        self.svm_small_cluster_score = -1.0
        self.svm_large_cluster_score = -1.0
        self.svm_class_score = -1.0
        self.num_sources = num_sources

    def __str__(self):
        return SEPARATOR.join(map(str, [self.substrate, int(self.validated),
                                        self.aa10_score1, self.aa10_score2, self.aa10_score3, self.aa10_score4,
                                        self.aa34_score1, self.aa34_score2, self.aa34_score3, self.aa34_score4,
                                        self.svm_single_amino_score, self.svm_small_cluster_score,
                                        self.svm_large_cluster_score, self.svm_class_score,
                                        self.num_sources]))


def generate_nerpa_table_entries(extended_signatures_file, sorted_nerpa_monomer_names, short_name_to_nerpa_monomer,
                                 nerpa_monomer_to_SVM_short_names, substrate_to_aa10_codes, substrate_to_aa34_codes):

    def _get_max_scores(svm_scores, svm_short_names):
        max_scores = [-1.0, -1.0, -1.0, -1.0]  # default scores if no match
        for svm_short_name in svm_short_names:
            assert svm_short_name in svm_scores, f"{svm_short_name} not in 'svm_scores'! Should be there by design"
            max_scores = [
                max(current_max, new_score)
                for current_max, new_score in zip(max_scores, svm_scores[svm_short_name])
            ]
        return max_scores

    def _create_nerpa_table_entry(row, nerpa_monomer, validated, svm_scores):
        entry = NerpaTableEntry(nerpa_monomer, validated, row['num_domains'])

        # Compute and set aa10 and aa34 scores
        aa10_score1 = _get_aa_score(row['aa10'], substrate_to_aa10_codes[nerpa_monomer])
        aa34_score1 = _get_aa_score(row['aa34'], substrate_to_aa34_codes[nerpa_monomer])

        if validated:
            assert aa10_score1 == 1.0, "score1 for validated substrates should be 1.0 by design (aa34)"
            assert aa34_score1 == 1.0, "score1 for validated substrates should be 1.0 by design (aa34)"
            entry.aa10_score1 = aa10_score1
            entry.aa34_score1 = aa34_score1
            _compute_and_set_scores(row['aa10'], substrate_to_aa10_codes[nerpa_monomer], entry, 'aa10')
            _compute_and_set_scores(row['aa34'], substrate_to_aa34_codes[nerpa_monomer], entry, 'aa34')
        else:
            _set_all_scores(entry, aa10_score1, 'aa10')
            _set_all_scores(entry, aa34_score1, 'aa34')

        # Set SVM scores
        assert nerpa_monomer in nerpa_monomer_to_SVM_short_names, f"{nerpa_monomer} is not in nerpa_monomer_to_SVM_short_names!"
        svm_short_names = nerpa_monomer_to_SVM_short_names[nerpa_monomer]
        max_scores = _get_max_scores(svm_scores, svm_short_names)
        (entry.svm_single_amino_score, entry.svm_small_cluster_score,
         entry.svm_large_cluster_score, entry.svm_class_score) = max_scores

        return entry

    def _get_aa_score(code, all_codes):
        max_matched_letters = 0
        for c in all_codes:
            max_matched_letters = max(max_matched_letters, sum(code[i] == c[i] for i in range(len(c))))
        return float(max_matched_letters) / len(code)

    def _compute_and_set_scores(code, all_codes, entry, aa_type):
        codes_copy = all_codes.copy()
        if code in codes_copy:
            codes_copy[code] -= 1  # removing the true answer
            if codes_copy[code] <= 0:
                del codes_copy[code]
        setattr(entry, f'{aa_type}_score2', _get_aa_score(code, codes_copy))

        if code in codes_copy:
            codes_copy[code] -= (entry.num_sources - 1)  # removing duplicates (from the same signature entry)
            if codes_copy[code] <= 0:
                del codes_copy[code]
        setattr(entry, f'{aa_type}_score3', _get_aa_score(code, codes_copy))

        if code in codes_copy:
            del codes_copy[code]  # removing all equivalent codes
        setattr(entry, f'{aa_type}_score4', _get_aa_score(code, codes_copy))

    def _set_all_scores(entry, score_value, aa_type):
        for idx in range(1, 5):
            setattr(entry, f'{aa_type}_score{idx}', score_value)

    # Load the extended signatures
    df = pd.read_csv(extended_signatures_file, sep=SEPARATOR)

    # Initialize the list of NerpaTableEntry objects
    nerpa_table_entries = []

    # Iterate through each row and generate NerpaTableEntry objects (exactly len(sorted_nerpa_monomer_names) per row)
    for _, row in df.iterrows():
        validated_nerpa_monomers = {short_name_to_nerpa_monomer[substrate] for substrate in eval(row['substrates'])}

        # Collect SVM scores (ignoring the first four columns)
        svm_monomer_to_its_scores = {col: eval(row[col]) for col in df.columns[4:]}

        for nerpa_monomer in sorted_nerpa_monomer_names:
            nerpa_table_entries.append(
                _create_nerpa_table_entry(row, nerpa_monomer,
                                          validated=nerpa_monomer in validated_nerpa_monomers,
                                          svm_scores=svm_monomer_to_its_scores)
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
