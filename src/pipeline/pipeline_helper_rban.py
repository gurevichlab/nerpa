from typing import Dict, List
from src.pipeline.command_line_args_helper import CommandLineArgs
from src.pipeline.logger import NerpaLogger
from src.config import Config, rBAN_Processing_Config
from src.data_types import Chirality, NRP_Variant

from src.rban_parsing.retrieve_nrp_variants import retrieve_nrp_variants
from src.rban_parsing.rban_parser import (
    Parsed_rBAN_Record,
    parsed_chiralities,
    Raw_rBAN_Record
)
from src.monomer_names_helper import MonomerNamesHelper
from src.rban_parsing.handle_monomers import get_monomers_chirality

from src.pipeline.rban_helper import rBAN_Helper
from pathlib import Path
import csv
import json
import yaml
from dataclasses import dataclass
from itertools import chain


@dataclass
class PipelineHelper_rBAN:
    config: Config
    args: CommandLineArgs
    log: NerpaLogger
    monomer_names_helper: MonomerNamesHelper = None

    def __post_init__(self):
        self.set_rban_helper()

    def preprocessed_nrp_variants(self) -> bool:
        return self.args.structures is not None

    def load_nrp_variants(self) -> List[NRP_Variant]:
        self.log.info('Loading preprocessed NRP variants')
        nrp_variants = []
        for file_with_nrp_variants in filter(lambda f: f.suffix in ('.yml', '.yaml'),
                                             self.args.structures.iterdir()):
            nrp_variants.extend(NRP_Variant.from_yaml_dict(yaml_record)
                                for yaml_record in yaml.safe_load(file_with_nrp_variants.read_text()))
        return nrp_variants

    def set_rban_helper(self):
        custom_monomers = list(chain(*(json.load(path_to_monomers.open('r'))
                                       for path_to_monomers in (self.config.rban_config.nerpa_monomers,
                                                                self.args.rban_monomers)
                                       if path_to_monomers is not None)))
        self.rban_helper = rBAN_Helper(self.config.rban_config,
                                       self.config.output_config.rban_output_config,
                                       custom_monomers)

    def get_input_for_rban(self) -> List[dict]:
        def default_id(i: int) -> str:
            return f'compound_{i:06d}'

        if self.args.smiles_tsv:
            reader = csv.DictReader(self.args.smiles_tsv.open('r'),
                                    delimiter=self.args.sep, quoting=csv.QUOTE_NONE)
            return [{'id': row[self.args.col_id] if self.args.col_id else default_id(i),
                     'smiles': row[self.args.col_smiles]}
                    for i, row in enumerate(reader)]
        elif self.args.smiles:
            return [{'id': default_id(i), 'smiles': smiles}
                    for i, smiles in enumerate(self.args.smiles)]

    def get_raw_rban_results(self) -> List[Raw_rBAN_Record]:
        if self.args.rban_output:
            rban_records = json.load(self.args.rban_output.open('r'))
            return rban_records if isinstance(rban_records, list) else [rban_records]

        smiles_with_ids = self.get_input_for_rban()

        self.log.info('\n======= Structures preprocessing with rBAN')
        return self.rban_helper.run_rban(smiles_with_ids,
                                         self.log, report_not_processed=True)

    def get_chiralities_per_record(self, rban_records: List[Raw_rBAN_Record]) -> Dict[int, Dict[int, Chirality]]:
        chiralities_dict = {}
        for i, rban_record in enumerate(rban_records):
            try:
                chiralities_dict[i] = parsed_chiralities(get_monomers_chirality(rban_record))
            except:
                self.log.warning(
                    f'Structure "{rban_record["id"]}": unexpected error while determining stereoisomeric configuration '
                    f'of monomers. Stereoisomeric configuration will be omitted.')
        return chiralities_dict

    def get_rban_results(self) -> List[Parsed_rBAN_Record]:
        rban_records = self.get_raw_rban_results()
        hybrid_monomers_per_record = self.rban_helper.get_hybrid_monomers_per_record(rban_records, self.log)
        chiralities_per_record = self.get_chiralities_per_record(rban_records)

        self.log.info("\n======= Done with Structures preprocessing with rBAN")
        return [Parsed_rBAN_Record(rban_record,
                                   hybrid_monomers_per_record[i],
                                   chiralities_per_record[i])
                for i, rban_record in enumerate(rban_records)]

    def get_nrp_variants(self,
                         parsed_rban_records: List[Parsed_rBAN_Record]) -> List[NRP_Variant]:
        self.log.info('\n======= Processing rBAN output')
        self.log.info(f'results will be in {self.config.output_config.main_out_dir / Path("NRP_Variants")}')
        nrp_variants = retrieve_nrp_variants(parsed_rban_records,
                                             self.monomer_names_helper,
                                             self.config.rban_processing_config,
                                             self.log)
        self.log.info('\n======= Done with Processing rBAN output')
        return nrp_variants
