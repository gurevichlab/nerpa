from typing import List
from command_line_args_helper import CommandLineArgs
from src.nerpa_pipeline.logger import NerpaLogger
from config import Config
from src.data_types import NRP_Variant

from src.nerpa_pipeline import handle_rban
from src.nerpa_pipeline.rban_parser import (
    Parsed_rBAN_Record,
    parsed_chiralities,
)
from src.nerpa_pipeline.rban_names_helper import rBAN_Names_Helper
from src.nerpa_pipeline.handle_monomers import get_monomers_chirality

from rban_helper import rBAN_Helper
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
    rban_helper: rBAN_Helper = None

    def __post_init__(self):
        self.set_rban_helper()

    def preprocessed_nrp_variants(self) -> bool:
        return self.args.structures is not None

    def load_nrp_variants(self) -> List[NRP_Variant]:
        self.log.info('Loading preprocessed NRP variants')
        nrp_variants = []
        for file_with_nrp_variants in filter(lambda f: f.suffix in ('.yml', '.yaml'),
                                             Path(self.args.structures).iterdir()):
            nrp_variants.extend(NRP_Variant.from_yaml_dict(yaml_record)
                                for yaml_record in yaml.safe_load(file_with_nrp_variants.read_text()))
        return nrp_variants

    def set_rban_helper(self):
        custom_monomers = list(chain(*(json.load(path_to_monomers.open('r'))
                                       for path_to_monomers in (self.config.paths.nerpa_monomers,
                                                                self.args.rban_monomers)
                                       if path_to_monomers is not None)))
        self.rban_helper = rBAN_Helper(self.config.rban_config, custom_monomers)

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

    def get_rban_results(self) -> List[Parsed_rBAN_Record]:
        smiles_with_ids = self.get_input_for_rban()
        self.log.info('\n======= Structures preprocessing with rBAN')
        rban_records = self.rban_helper.run_rban(smiles_with_ids,
                                                 self.log, report_not_processed=True)
        hybrid_monomers_dict = self.rban_helper.get_hybrid_monomers(rban_records, self.log)
        chiralities_dict = {}
        for i, rban_record in enumerate(rban_records):
            try:
                chiralities_dict[i] = parsed_chiralities(get_monomers_chirality(rban_record))
            except:
                self.log.warning(
                    f'Structure "{rban_record["id"]}": unexpected error while determining stereoisomeric configuration '
                    f'of monomers. Stereoisomeric configuration will be omitted.')

        self.log.info("\n======= Done with Structures preprocessing with rBAN")
        return [Parsed_rBAN_Record(rban_record,
                                   hybrid_monomers_dict[i],
                                   chiralities_dict[i])
                for i, rban_record in enumerate(rban_records)]

    def get_nrp_variants(self, parsed_rban_records: List[Parsed_rBAN_Record]) -> List[NRP_Variant]:
        self.log.info('\n======= Processing rBAN output')
        self.log.info(f'results will be in {self.config.paths.main_out_dir / Path("NRP_Variants")}')
        recognized_monomers = [row['Code']
                               for row in csv.DictReader(self.config.paths.nerpa_monomers_info.open('r'),
                                                         delimiter='\t')]
        rban_names_helper = rBAN_Names_Helper(self.config.paths.nerpa_monomers_info)
        nrp_variants = handle_rban.retrieve_nrp_variants(parsed_rban_records,
                                                         recognized_monomers, rban_names_helper,
                                                         self.log)
        self.log.info('\n======= Done with Processing rBAN output')
        return nrp_variants
