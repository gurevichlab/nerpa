from __future__ import annotations
import polars as pl
from pathlib import Path


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

