import subprocess
import os
from pathlib import Path
from src.benchmarking.data_frames import (
    MIBiG_BGCs_Info,
    PNRPDB_Info,
    get_good_mibig_bgs_to_nrps,
)

os.environ["RUST_BACKTRACE"] = "1"
os.environ["RUSTFLAGS"] = "-Awarnings"
nerpa_root = (
    Path(__file__).resolve()
    .parent.parent.parent.parent.parent
)

mibig_bgcs_info = MIBiG_BGCs_Info.load_from_nerpa_root(nerpa_root)
pnrpdb_info = PNRPDB_Info.load_from_nerpa_root(nerpa_root)
