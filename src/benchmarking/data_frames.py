from __future__ import annotations
import polars as pl
from pathlib import Path
from typing import List

class PNRPDB_Compound_Similarity(pl.DataFrame):
    FST_COMPOUND_ID = 'nrp1_id'
    SND_COMPOUND_ID = 'nrp2_id'

    NERPA_ISO_WEAK_CHR = 'gt-nerpa__d<=0__chr-weak'
    NERPA_NO_MORE_ONE_SUB_WEAK_CHR = 'gt-nerpa__d<=1__chr-weak'

    NERPA_ISO = 'gt-nerpa__d<=0__chr-strict'
    NERPA_NO_MORE_ONE_SUB = 'gt-nerpa__d<=1__chr-strict'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def from_csv(cls, filepath: Path) -> PNRPDB_Compound_Similarity:
        """Load MIBiG BGCs info from a CSV file."""
        df = pl.read_csv(filepath, separator='\t')
        return df

NRP_CMP_METHOD = str

NRP_CMP_METHODS: List[NRP_CMP_METHOD] = [
    PNRPDB_Compound_Similarity.NERPA_ISO,
    PNRPDB_Compound_Similarity.NERPA_NO_MORE_ONE_SUB,
    PNRPDB_Compound_Similarity.NERPA_ISO_WEAK_CHR,
    PNRPDB_Compound_Similarity.NERPA_NO_MORE_ONE_SUB_WEAK_CHR,
]
    

class MIBiG_BGCs_Info(pl.DataFrame):
    """DataFrame for MIBiG BGC information."""
    BGC_ID = 'bgc_id'
    NRP_IDS = 'nrp_ids'
    BIGSCAPE_FAMILIES = 'bigscape_families'
    NUM_A_DOMAINS = 'num_a_domains'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def from_csv(cls, filepath: Path) -> MIBiG_BGCs_Info:
        """Load MIBiG BGCs info from a CSV file."""
        df = pl.read_csv(filepath, separator='\t')
        return df


class PNRPDB_Info(pl.DataFrame):
    """DataFrame for PNRPDB information."""
    COMPOUND_ID = 'compound_id'
    REPR_ID = 'nrp_variant_iso_class_representative'
    NUM_NRPS_MONOMERS = 'num_nrps_monomers'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def from_csv(cls, filepath: Path) -> PNRPDB_Info:
        """Load PNRPDB info from a CSV file."""
        df = pl.read_csv(filepath, separator='\t')
        return df

