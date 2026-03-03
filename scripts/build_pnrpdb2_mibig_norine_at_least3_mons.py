from pathlib import Path
import polars as pl

from src.benchmarking.data_frames import PNRPDB_Info

nerpa_dir = Path(__file__).parent.parent

# PNRPDB_Info is inherited from polars.Dataframe
pnrpdb_info = PNRPDB_Info.from_csv(
    nerpa_dir / 'data/for_training_and_testing/pnrpdb2_info.tsv'
)

def is_mibig_norine(nrp_id: str) -> bool:
    return nrp_id.startswith('NOR') or nrp_id.startswith('BGC')

# q: filter to only MIBiG/Norine NRPs with at least 3 monomers (use 'num_nrps_monomers' column)
filtered = pnrpdb_info.filter(
    (pl.col(PNRPDB_Info.NUM_NRPS_MONOMERS) >= 3) &
    (pl.col(PNRPDB_Info.COMPOUND_ID).apply(is_mibig_norine))
)[PNRPDB_Info.COMPOUND_ID].to_list()

pnrpdb = pl.read_csv(nerpa_dir / 'data/input/pnrpdb2_expanded.tsv', separator='\t')

filtered_pnrpdb = pnrpdb.filter(pl.col('ID').is_in(filtered))

filtered_pnrpdb.write_csv(nerpa_dir / 'data/input/pnrpdb2_mibig_norine_at_least3_mons.tsv', separator='\t')


