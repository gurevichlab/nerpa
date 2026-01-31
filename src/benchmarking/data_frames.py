from __future__ import annotations
import polars as pl
from pathlib import Path

class PNRPDB_Compound_Similarity(pl.DataFrame):
    FST_COMPOUND_ID = 'nrp1_id'
    SND_COMPOUND_ID = 'nrp2_id'
    # rBAN_ISO_ALLOW_UNK_CHR = 'rban_isomorphic_unknown_chr_equal_known_cmp'
    NERPA_EQUAL_ALLOW_UNK_CHR = 'nerpa_equal_unknown_chr_equal_known_cmp'
    # rBAN_ONE_SUB_ALLOW_UNK_CHR = 'rban_one_sub_away_unknown_chr_equal_known_cmp'
    NERPA_ONE_SUB_ALLOW_UNK_CHR = 'nerpa_one_sub_away_unknown_chr_equal_known_cmp'
    NERPA_NO_MORE_ONE_SUB_ALLOW_UNK_CHR = 'nerpa_no_more_one_sub_away_unknown_chr_equal_known_cmp'
    NERPA_EQUAL = 'nerpa_equal_nerpa_mon_cmp'
    NERPA_ONE_SUB = 'nerpa_one_sub_away_nerpa_mon_cmp'
    NERPA_NO_MORE_ONE_SUB = 'nerpa_no_more_one_sub_away_nerpa_mon_cmp'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def from_csv(cls, filepath: Path) -> PNRPDB_Compound_Similarity:
        """Load MIBiG BGCs info from a CSV file."""
        df = pl.read_csv(filepath, separator='\t')
        return df

class MIBiG_BGCs_Info(pl.DataFrame):
    """DataFrame for MIBiG BGC information."""
    BGC_ID = 'bgc_id'
    NRP_ID = 'compound_id'
    IN_APPROVED_MATCHES = 'in_approved_matches'
    BIGSCAPE_FAMILIES = 'bigscape_families'
    NUM_A_DOMAINS = 'num_a_domains'
    _DEPRECATED = 'deprecated in mibig4'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def from_csv(cls, filepath: Path) -> MIBiG_BGCs_Info:
        """Load MIBiG BGCs info from a CSV file."""
        df = pl.read_csv(filepath, separator='\t')
        df = df.filter(~pl.col(cls._DEPRECATED))
        return df


class PNRPDB_Info(pl.DataFrame):
    """DataFrame for PNRPDB information."""
    COMPOUND_ID = 'compound_id'
    ISO_CLASS_ID = 'nrp_variant_iso_class_representative'
    NUM_NRPS_MONOMERS = 'num_nrps_monomers'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def from_csv(cls, filepath: Path) -> PNRPDB_Info:
        """Load PNRPDB info from a CSV file."""
        df = pl.read_csv(filepath, separator='\t')
        return df

