#!/usr/bin/env python3
"""
Generator of the table for the new Nerpa scoring
"""
from copy import copy
from pathlib import Path
import nrpys
import antismash.modules.nrps_pks.nrpys as antismash_nrpys
from antismash.modules.nrps_pks.nrpys import PredictorSVMResult
import os
import pandas as pd
from typing import Dict, List, Tuple
from collections import Counter, defaultdict, OrderedDict
import yaml

# TODO: make it options
DEPENDENCIES_ROOT = "/Users/agu22/biolab/nerpa/NEW_SCORING/useful_files/"
SIGNATURES_FPATH = "./signatures.tsv"
# SIGNATURES_FPATH = "./aux_files/signatures_tmp.tsv"
NERPA_MONOMERS_FPATH = "./aux_files/monomersLogP.tsv"
OUTPUT_DIR = "nerpa_tables"
OUTPUT_FNAME = "nerpa_scoring_table.tsv"
OUTPUT_SEPARATOR = "\t"

# TODO: make it less ugly
NERPA_SUBSTRATES = pd.read_csv(NERPA_MONOMERS_FPATH, sep='\t')["NameID"].tolist()
NERPA_SUBSTRATES[-1] = "unknown"
IGNORED_MODS = ["me", "oh", "branched", "dh", "cl"]
parsed_fully = 0
parsed_with_mods = 0
not_parsed = 0

# TODO
AS_TO_NERPA = {
    # "easy" cases
    "Cya": "cysa",
    "Glyca": "glycolic-acid",
    "Aol": "alaninol",
    "bAla": "b-ala",
    "bLys": "b-lys",
    "D-Lya": "lyserg",
    "Pgl": "phg",
    "Piz": "piperazic",
    "3clLeu": "tcl",
    "Valol": "vol",
    "pPro": "4ppro",
    "2,3-dohBza": "dhb",
    # nontrivial ones: CHECKME once again!
    "Kiv": "2-oxo-isovaleric-acid",
    "pyrAla": "33p-l-ala",
    "R-ohTyr": "bht",
    "ohTyr": "bht",
    "dhAbu": "dht",
    "3-ohVal": "hyv",
    "me-clHic": "indole-3-carboxylic",
    "D-pheLac": "phe-ac",
    "n-epox-oxoDec": "aeo",
    "2-oh-4-mePen": "HICA",

    # CHECKME: nontrivial adjustment (but we ignore modifications anyway)
    "ohOrn": "orn",
    "aThr": "thr",
    "aIle": "ile"

    '''
    TODO: not matched so far 
    ???": "HICA --> 2-oh-4-mePen (according to the Google Doc this matched to rBAN's "4Me-Hva")
    ???": "LDAP --> ignore! Based on his "old" code (DAQDLAVVNK -- used in old Nerpa), it is just Dpr but we have such substrate already 
    ???": "aeo  --> n-epox-oxoDec
    '''
}

SVM_SUBSTRATES_ORIGINAL = ["Arg", "Asp", "Glu", "Asn", "Lys", "Gln", "Orn", "ohOrn", "Aad",
                           "Ala", "Gly", "Val", "Leu", "Ile", "Abu", "Iva", "Ser", "Thr", "Hpg", "dHpg", "Cys", "Pro", "Pip",
                           "Phe", "Tyr", "2,3-dohBza", "Pgl", "R-ohTyr"]

SVM_SUBSTRATES = [AS_TO_NERPA.get(substrate, substrate).lower() for substrate in SVM_SUBSTRATES_ORIGINAL]


class SignatureEntry:
    def __init__(self, raw_entry_line):
        assert len(raw_entry_line.split()) == 5, "Error: corrupted line in the signatures file"
        self.aa10, self.aa34, raw_substrates, _, sources = raw_entry_line.split()
        self.substrates = set(map(antismash_substrate_to_nerpa_substrate, raw_substrates.split('|')))
        self.num_sources = sources.count('|') + 1

    def __str__(self):
        return " ".join([self.aa10, self.aa34, ",".join(self.substrates), str(self.num_sources)])


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
        return OUTPUT_SEPARATOR.join(map(str, [self.substrate, int(self.validated),
                                               self.aa10_score1, self.aa10_score2, self.aa10_score3, self.aa10_score4,
                                               self.aa34_score1, self.aa34_score2, self.aa34_score3, self.aa34_score4,
                                               self.svm_single_amino_score, self.svm_small_cluster_score,
                                               self.svm_large_cluster_score, self.svm_class_score,
                                               self.num_sources]))


def antismash_substrate_to_nerpa_substrate(substrate):
    lc_name = AS_TO_NERPA.get(substrate, substrate.lower())
    if lc_name in NERPA_SUBSTRATES:
        global parsed_fully
        parsed_fully += 1
        return lc_name

    # ignoring modifications and clarifications
    '''
        examples: 
        omeTyr --> tyr
        4R-ohPro --> pro
        2S-Hiv --> hiv
    '''
    unmodified_substrate = substrate.split('-')[-1] if IGNORED_MODS else substrate
    for mod in IGNORED_MODS:
        unmodified_substrate = unmodified_substrate.split(mod)[-1]
    if unmodified_substrate != substrate:
        lc_name = AS_TO_NERPA.get(unmodified_substrate, unmodified_substrate).lower()
        if lc_name in NERPA_SUBSTRATES:
            global parsed_with_mods
            parsed_with_mods += 1
            return lc_name

    print(f"Failed to parse '{substrate}', tried '{lc_name}'")
    global not_parsed
    not_parsed += 1
    return NERPA_SUBSTRATES[-1]  # special case: the last entry is "unknown" or "none" or something like this


def run_nrpys(signature_entries):
    config = nrpys.Config()
    config.model_dir = os.path.join(DEPENDENCIES_ROOT, "models")
    config.stachelhaus_signatures = os.path.join(DEPENDENCIES_ROOT, "signatures.tsv")
    config.skip_v1 = True
    config.skip_v3 = True

    names: list[str] = [f"{idx}" for idx in range(len(signature_entries))]
    signatures: list[str] = [signature_entry.aa34 for signature_entry in signature_entries]

    return [antismash_nrpys.PredictorSVMResult.from_nrpys(nrpys_res)
            for nrpys_res in nrpys.run(config, names, signatures)]


def get_nerpa_table_entry(signature_entry: SignatureEntry, svm_result: PredictorSVMResult,
                          substrate: str, aa10_codes: Counter, aa34_codes: Counter):
    def __to_nerpa_substrates(svm_substrates):
        return [antismash_substrate_to_nerpa_substrate(s.short) for s in svm_substrates]

    def __get_aa_score(code, all_codes):
        if len(all_codes):
            assert len(code) == len(next(iter(all_codes)))
        max_matched_letters = 0
        for c in all_codes:
            max_matched_letters = max(max_matched_letters, sum(code[i] == c[i] for i in range(len(c))))
        return float(max_matched_letters) / len(code)

    nerpa_table_entry = NerpaTableEntry(substrate, substrate in signature_entry.substrates, signature_entry.num_sources)
    nerpa_table_entry.aa10_score1 = __get_aa_score(signature_entry.aa10, aa10_codes)
    nerpa_table_entry.aa34_score1 = __get_aa_score(signature_entry.aa34, aa34_codes)
    if nerpa_table_entry.validated:
        code = signature_entry.aa10
        all_codes = copy(aa10_codes)
        assert code in all_codes
        all_codes[code] -= 1  # removing the true answer
        if all_codes[code] <= 0:
            del all_codes[code]
        nerpa_table_entry.aa10_score2 = __get_aa_score(code, all_codes)
        all_codes[code] -= (nerpa_table_entry.num_sources - 1)  # removing duplicates (from the same signature entry)
        if all_codes[code] <= 0:
            del all_codes[code]
        nerpa_table_entry.aa10_score3 = __get_aa_score(code, all_codes)
        del all_codes[code]  # removing all equivalent codes
        nerpa_table_entry.aa10_score4 = __get_aa_score(code, all_codes)

        # the same for aa34  FIXME: not so ugly please
        code = signature_entry.aa34
        all_codes = copy(aa34_codes)
        assert code in all_codes
        all_codes[code] -= 1  # removing the true answer
        if all_codes[code] <= 0:
            del all_codes[code]
        nerpa_table_entry.aa34_score2 = __get_aa_score(code, all_codes)
        all_codes[code] -= (nerpa_table_entry.num_sources - 1)  # removing duplicates (from the same signature entry)
        if all_codes[code] <= 0:
            del all_codes[code]
        nerpa_table_entry.aa34_score3 = __get_aa_score(code, all_codes)
        del all_codes[code]  # removing all equivalent codes
        nerpa_table_entry.aa34_score4 = __get_aa_score(code, all_codes)
    else:
        nerpa_table_entry.aa10_score2, nerpa_table_entry.aa10_score3, nerpa_table_entry.aa10_score4 = (
                [nerpa_table_entry.aa10_score1] * 3)
        nerpa_table_entry.aa34_score2, nerpa_table_entry.aa34_score3, nerpa_table_entry.aa34_score4 = (
                [nerpa_table_entry.aa34_score1] * 3)


    if substrate in SVM_SUBSTRATES:
        nerpa_table_entry.svm_single_amino_score = svm_result.single_amino.score \
            if substrate in __to_nerpa_substrates(svm_result.single_amino.substrates) else 0.0
        nerpa_table_entry.svm_small_cluster_score = svm_result.small_cluster.score \
            if substrate in __to_nerpa_substrates(svm_result.small_cluster.substrates) else 0.0
        nerpa_table_entry.svm_large_cluster_score = svm_result.large_cluster.score \
            if substrate in __to_nerpa_substrates(svm_result.large_cluster.substrates) else 0.0
        nerpa_table_entry.svm_class_score = svm_result.physicochemical_class.score \
            if substrate in __to_nerpa_substrates(svm_result.physicochemical_class.substrates) else 0.0

    return nerpa_table_entry


def main():
    signature_entries = list()
    substrate_to_aa10_codes = defaultdict(Counter)
    substrate_to_aa34_codes = defaultdict(Counter)
    nerpa_table_entries = list()

    # for debug
    MAX_LINES_TO_PROCESS = -1  # -1 to process everything

    with open(SIGNATURES_FPATH) as f:
        for idx, raw_entry in enumerate(f):
            if MAX_LINES_TO_PROCESS > 0 and idx > MAX_LINES_TO_PROCESS:
                break
            signature_entry = SignatureEntry(raw_entry)
            # CHECKME: we might want to skip entries with not parsed substrates
            # right now, let's treat them as just another type of substrate ("unknown")
            for substrate in signature_entry.substrates:
                substrate_to_aa10_codes[substrate][signature_entry.aa10] += signature_entry.num_sources
                substrate_to_aa34_codes[substrate][signature_entry.aa34] += signature_entry.num_sources
            signature_entries.append(signature_entry)

    global parsed_fully
    global parsed_with_mods
    global not_parsed
    print(f"Stats (fully, with mods, not parsed): {parsed_fully}, {parsed_with_mods}, {not_parsed}")

    # Optimized code for reformatting into codes_dict with lexicographical order
    codes_dict: Dict[str, Tuple[List[str], List[str]]] = OrderedDict()

    def custom_sort(substrate):
        return (substrate == 'unknown', substrate)

    for substrate in sorted(substrate_to_aa10_codes.keys(), key=custom_sort):
        aa10_list = list(substrate_to_aa10_codes[substrate].keys())
        aa34_list = list(substrate_to_aa34_codes[substrate].keys())
        codes_dict[substrate] = (aa10_list, aa34_list)

    # Dumping codes_dict to YAML file
    yaml_fpath = os.path.join(OUTPUT_DIR, 'residue_signatures.yaml')
    with open(yaml_fpath, 'w') as yaml_file:
        yaml.dump(codes_dict, yaml_file, default_flow_style=False)

    print(f"Saved codes in {yaml_fpath}")
    exit(0)

    # filling in the Nerpa table
    svm_results = run_nrpys(signature_entries)

    nerpa_table_column_headers = ["substrate", "is_correct",
                                  "aa10_score1", "aa10_score2", "aa10_score3", "aa10_score4",
                                  "aa34_score1", "aa34_score2", "aa34_score3", "aa34_score4",
                                  "svm_single_amino_score", "svm_small_cluster_score",
                                  "svm_large_cluster_score", "svm_class_score",
                                  "num_sources"]  # see the fields in NerpaTableEntry

    for signature_entry, svm_result in zip(signature_entries, svm_results):
        for substrate in NERPA_SUBSTRATES:
            nerpa_table_entries.append(get_nerpa_table_entry(signature_entry, svm_result, substrate,
                                                             substrate_to_aa10_codes[substrate],
                                                             substrate_to_aa34_codes[substrate]))

    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    with open(output_fname := os.path.join(OUTPUT_DIR, OUTPUT_FNAME), 'w') as f:
        f.write(OUTPUT_SEPARATOR.join(nerpa_table_column_headers) + "\n")
        f.write("\n".join(map(str, nerpa_table_entries)))
    print(f"Done! See results in {output_fname}")


if __name__ == "__main__":
    main()