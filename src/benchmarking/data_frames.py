from __future__ import annotations
import polars as pl
from pathlib import Path
from collections import defaultdict
from typing import Dict

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

    @classmethod
    def load_from_nerpa_root(cls, nerpa_root: Path) -> MIBiG_BGCs_Info:
        """Load MIBiG BGCs info from the nerpa repo's data folder."""
        assert (nerpa_root / "nerpa.py").exists(), f"Invalid nerpa_root: {nerpa_root}"
        return cls.from_csv(
            nerpa_root
            / 'data'
            / 'for_training_and_testing'
            / 'mibig_bgcs_info.tsv'
        )


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

    @classmethod
    def load_from_nerpa_root(cls, nerpa_root: Path) -> PNRPDB_Info:
        """Load PNRPDB info from the nerpa repo's data folder."""
        assert (nerpa_root / "nerpa.py").exists(), f"Invalid nerpa_root: {nerpa_root}"
        return cls.from_csv(
            nerpa_root
            / 'data'
            / 'for_training_and_testing'
            / 'pnrpdb2_info.tsv'
        )

    @classmethod
    def is_mibig_nrp(cls) -> pl.Expr:
        return (
            pl.col(PNRPDB_Info.COMPOUND_ID)
            .str
            .starts_with('BGC')
        )

def get_good_mibig_bgs_to_nrps(
        mibig_bgcs_info: MIBiG_BGCs_Info,
        pnrpdb_info: PNRPDB_Info,
        min_num_a_domains: int = 3,
        min_num_nrps_monomers: int = 3,
) -> Dict[str, list[str]]:
    '''
    Get a mapping from MIBiG BGC IDs to their NRPs.
    Only include BGCs that have at least min_num_a_domains A domains and NRPs that have at least min_num_nrps_monomers monomers
    '''
    good_nrps = set(
        pnrpdb_info
        .filter(
            PNRPDB_Info.is_mibig_nrp() &
            (pl.col(PNRPDB_Info.NUM_NRPS_MONOMERS) >= min_num_nrps_monomers)
        )
        [PNRPDB_Info.COMPOUND_ID]
    )

    bgc_to_nrps = defaultdict(list)
    for row in mibig_bgcs_info.iter_rows(named=True):
        bgc_id = row[MIBiG_BGCs_Info.BGC_ID]
        nrp_id = row[MIBiG_BGCs_Info.NRP_ID]
        num_a_domains = row[MIBiG_BGCs_Info.NUM_A_DOMAINS]

        if num_a_domains >= min_num_a_domains and nrp_id in good_nrps:
            bgc_to_nrps[bgc_id].append(nrp_id)

    return bgc_to_nrps
    
