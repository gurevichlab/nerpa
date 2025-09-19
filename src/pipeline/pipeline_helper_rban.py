from collections import defaultdict
from typing import Dict, List, Tuple
from src.pipeline.command_line_args_helper import CommandLineArgs
from src.pipeline.logging.logger import NerpaLogger
from src.config import Config
from src.data_types import Chirality, NRP_Variant

from src.rban_parsing.retrieve_nrp_variants import retrieve_nrp_variants
from src.rban_parsing.rban_parser import (
    Parsed_rBAN_Record,
    parsed_chiralities,
    Raw_rBAN_Record, NRP_metadata
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
        return self.args.nrp_variants is not None

    def load_nrp_variants(self) -> List[NRP_Variant]:
        self.log.info('Loading preprocessed NRP variants')

        if self.args.nrp_variants.is_dir():
            files_with_nrp_variants = filter(lambda f: f.suffix in ('.yml', '.yaml'),
                                             self.args.nrp_variants.iterdir())
        else:
            files_with_nrp_variants = [self.args.nrp_variants]

        return [NRP_Variant.from_yaml_dict(yaml_record)
                for file_path in files_with_nrp_variants
                for yaml_record in yaml.safe_load(file_path.read_text())]

    def set_rban_helper(self):
        custom_monomers = list(chain(*(json.load(path_to_monomers.open('r'))
                                       for path_to_monomers in (self.config.rban_config.nerpa_monomers,
                                                                self.args.rban_monomers)
                                       if path_to_monomers is not None)))
        self.rban_helper = rBAN_Helper(self.config.rban_config,
                                       self.config.output_config.rban_output_config,
                                       custom_monomers)

    def get_input_for_rban_with_metadata(self) -> Tuple[List[dict], Dict[str, NRP_metadata]]:
        def default_id(i: int) -> str:
            return f'compound_{i:06d}'

        rows = []
        metadata = {}
        if self.args.smiles_tsv:
            reader = csv.DictReader(self.args.smiles_tsv.open('r'),
                                    delimiter=self.args.sep, quoting=csv.QUOTE_NONE)
            for i, row in enumerate(reader):
                compound_id = row[self.args.col_id] if self.args.col_id else default_id(i)
                metadata[compound_id] = NRP_metadata.from_dict(row)
                rows.append({'id': compound_id,
                             'smiles': row[self.args.col_smiles]})
        elif self.args.smiles:
            for i, smiles in enumerate(self.args.smiles):
                compound_id = default_id(i)
                metadata[compound_id] = NRP_metadata.from_dict({'smiles': smiles})
                rows.append({'id': compound_id,
                             'smiles': smiles})

        return rows, metadata

    def get_raw_rban_results(self) -> Tuple[List[Raw_rBAN_Record], Dict[str, NRP_metadata]]:
        if self.args.rban_output:
            metadata = defaultdict(lambda : NRP_metadata(None, None, None, None, None))
            rban_records = json.load(self.args.rban_output.open('r'))
            return (rban_records, metadata) if isinstance(rban_records, list) else ([rban_records], metadata)

        smiles_with_ids, metadata = self.get_input_for_rban_with_metadata()

        self.log.info('\n======= Structures preprocessing with rBAN')
        rban_output =  self.rban_helper.run_rban(smiles_with_ids,
                                                 self.log,
                                                 report_not_processed=True)

        return rban_output, metadata


    def get_chiralities_per_record(self, rban_records: List[Raw_rBAN_Record]) -> Dict[int, Dict[int, Chirality]]:
        chiralities_dict = {}
        for i, rban_record in enumerate(rban_records):
            try:
                chiralities_dict[i] = parsed_chiralities(get_monomers_chirality(rban_record))
            except:
                self.log.warning(
                    f'Structure "{rban_record["id"]}": unexpected error while determining stereoisomeric configuration '
                    f'of monomers. Stereoisomeric configuration will be omitted.')
                chiralities_dict[i] = defaultdict(lambda: Chirality.UNKNOWN)
        return chiralities_dict

    def get_rban_results(self) -> List[Parsed_rBAN_Record]:
        rban_records, metadata = self.get_raw_rban_results()
        hybrid_monomers_per_record = self.rban_helper.get_hybrid_monomers_per_record(rban_records, self.log)
        chiralities_per_record = self.get_chiralities_per_record(rban_records)

        self.log.info("\n======= Done with Structures preprocessing with rBAN")
        return [Parsed_rBAN_Record(rban_record,
                                   hybrid_monomers_per_record[i],
                                   chiralities_per_record[i],
                                   metadata)
                for i, rban_record in enumerate(rban_records)]

    def get_nrp_variants(self,
                         parsed_rban_records: List[Parsed_rBAN_Record]) -> List[NRP_Variant]:
        self.log.info('\n======= Processing rBAN output')
        self.log.info(f'results will be in {self.config.output_config.nrp_variants}')
        nrp_variants = retrieve_nrp_variants(parsed_rban_records,
                                             self.monomer_names_helper,
                                             self.config.rban_processing_config,
                                             self.log)
        self.log.info('\n======= Done with Processing rBAN output')
        return nrp_variants
