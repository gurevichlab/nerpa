from collections import defaultdict
from typing import (
    Dict,
    List,
    Tuple,
    Union
)
from src.data_types import SMILES, rBAN_Residue_Name
from command_line_args_helper import CommandLineArgs
from src.nerpa_pipeline import nerpa_utils
from src.nerpa_pipeline.logger import NerpaLogger
from config import rBAN_Config


import csv
import json
import os
from pathlib import Path
from zipfile import ZipFile
from itertools import chain
from functools import reduce
from dataclasses import dataclass

MonomerData = dict  # I.O.: I have no idea what is actually stored there
from src.nerpa_pipeline.rban_parser import (
    MonomerIdx,
    Raw_rBAN_Record,
    trim_unknown_monomers
)

@dataclass
class antiSMASH_Helper:
    antismash_exec: Path

    def run_antismash(self,
                      seqs: Path,  # fasta or gbk file with sequences
                      threads: int,
                      log: NerpaLogger,
                      output_dir: Path) -> Union[None, Path]:
        if self.antismash_exec is None:
            log.error("Can't find antismash 5.x executable. Please make sure that you have antismash 5.x installed "
                      "in your system or provide path to antismash source directory via --antismash-path option.")
            return None

        command = [self.antismash_exec,
                   '--genefinding-tool', 'prodigal',
                   '--output-dir', output_dir,
                   '--minimal', '--skip-zip-file', '--enable-nrps-pks',
                   '--cpus', str(threads), seqs]
        nerpa_utils.sys_call(command, log, cwd=output_dir)
        return output_dir
