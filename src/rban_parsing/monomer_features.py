from typing import List, Tuple
from enum import Enum, auto
from src.generic.graphs import BackboneSequence
from src.rban_parsing.rban_parser import MonomerInfo
from src.monomer_names_helper import UNKNOWN_RESIDUE, MonomerNamesHelper
from src.rban_parsing.monomer_features_types import MonomerFeature, MonomerFeatures


def get_monomer_features(monomer_idx: int,
                         monomer_info: MonomerInfo,
                         backbone_sequence: BackboneSequence,
                         names_helper: MonomerNamesHelper) -> MonomerFeatures:
    if backbone_sequence.is_cyclic:
        return ()

    monomer_pos = backbone_sequence.node_idxs.index(monomer_idx)
    pairs = [
        (MonomerFeature.START_OF_FRAGMENT, monomer_pos == 0),
        (MonomerFeature.END_OF_FRAGMENT, monomer_pos == len(backbone_sequence) - 1),
        (MonomerFeature.UNKNOWN_RESIDUE, names_helper.parsed_name(monomer_info.name, 'norine').residue == UNKNOWN_RESIDUE),
        (MonomerFeature.PKS_HYBRID, monomer_info.is_hybrid),
        (MonomerFeature.PKS, monomer_info.name in names_helper.pks_names)
    ]
    return tuple(feature for feature, is_present in pairs if is_present)