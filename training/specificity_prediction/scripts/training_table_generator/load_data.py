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
from itertools import islice


def get_nerpa_supported_monomers(as_residue_frequencies_tsv: Path,
                                 norine_residue_frequncies_tsv: Path,
                                 num_monomers: Optional[int],
                                 SEPARATOR: str) -> List[MonomerResidue]:
    as_table = pd.read_csv(as_residue_frequencies_tsv, sep=SEPARATOR)
    norine_table = pd.read_csv(norine_residue_frequncies_tsv, sep=SEPARATOR)
    selected_monomers = list(islice((res for res in as_table['core'] if res in norine_table['core'].values),
                                    num_monomers))
    return sorted(selected_monomers)


def get_short_name_to_nerpa_monomer_mapping(monomers_table_tsv: Path,
                                            nerpa_supported_monomers: List[MonomerResidue],
                                            nerpa_unknown_monomer: MonomerResidue,
                                            SEPARATOR: str) -> Dict[antiSMASH_MonomerName, MonomerResidue]:
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
        nerpa_monomer = short_name_to_nerpa_monomer[as_short_name]
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
                               short_name_to_nerpa_monomer: Dict[antiSMASH_MonomerName, MonomerResidue],
                               SEPARATOR: str) \
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


@dataclass
class TrainingTableInputData:
    nerpa_supported_monomers: List[MonomerResidue]
    nerpa_unknown_monomer: MonomerResidue
    short_name_to_nerpa_monomer: Dict[antiSMASH_MonomerName, MonomerResidue]
    nerpa_monomer_to_SVM_short_names: Dict[MonomerResidue, List[antiSMASH_MonomerName]]
    substrate_to_aa10_codes: SubstrateCodesDict
    substrate_to_aa34_codes: SubstrateCodesDict
    extended_signatures_df: pd.DataFrame


def load_data(as_monomer_frequencies_tsv: Path,
              norine_monomer_frequencies_tsv: Path,
              monomers_table_tsv: Path,
              extended_signatures_tsv: Path,
              num_monomers: Optional[int],
              nerpa_unknown_monomer: MonomerResidue,
              SEPARATOR: str) -> TrainingTableInputData:
    # Get the list of supported monomers
    nerpa_supported_monomers = get_nerpa_supported_monomers(as_monomer_frequencies_tsv,
                                                            norine_monomer_frequencies_tsv,
                                                            num_monomers, SEPARATOR)

    # Get the mapping from short names to Nerpa monomers
    short_name_to_nerpa_monomer = get_short_name_to_nerpa_monomer_mapping(monomers_table_tsv,
                                                                          nerpa_supported_monomers,
                                                                          nerpa_unknown_monomer,
                                                                          SEPARATOR)

    # Get the mapping from Nerpa monomers to SVM short names
    extended_signatures_df = pd.read_csv(extended_signatures_tsv, sep=SEPARATOR)
    as_short_names_svm = list(extended_signatures_df.columns[4:])
    nerpa_monomer_to_SVM_short_names = get_nerpa_monomer_to_SVM_short_names_mapping(as_short_names_svm,
                                                                                    short_name_to_nerpa_monomer,
                                                                                    nerpa_unknown_monomer)

    # Get the dictionaries for aa10 and aa34 codes
    substrate_to_aa10_codes, substrate_to_aa34_codes = read_substrate_to_aa_codes(extended_signatures_tsv,
                                                                                  short_name_to_nerpa_monomer,
                                                                                  SEPARATOR)

    return TrainingTableInputData(nerpa_supported_monomers,
                                  nerpa_unknown_monomer,
                                  short_name_to_nerpa_monomer,
                                  nerpa_monomer_to_SVM_short_names,
                                  substrate_to_aa10_codes,
                                  substrate_to_aa34_codes,
                                  pd.read_csv(extended_signatures_tsv, sep=SEPARATOR))
