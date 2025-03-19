from collections import defaultdict
from typing import (
    Dict,
    List,
    Optional,
    Tuple,
    Union
)
from src.monomer_names_helper import NorineMonomerName
from src.pipeline import nerpa_utils
from src.pipeline.logger import NerpaLogger
from src.config import rBAN_Config, rBAN_Output_Config

import json
from pathlib import Path
from zipfile import ZipFile
from dataclasses import dataclass

MonomerData = dict  # I.O.: I have no idea what is actually stored there
from src.rban_parsing.rban_parser import (
    MonomerIdx,
    Raw_rBAN_Record,
    get_hybrid_monomers_smiles
)

@dataclass
class rBAN_Helper:
    config: rBAN_Config
    output_cfg: rBAN_Output_Config
    monomers_db: List[MonomerData]

    def __init__(self,
                 config: rBAN_Config, output_cfg: rBAN_Output_Config,
                 custom_monomers: Union[List[MonomerData], None] = None):
        self.config = config
        self.output_cfg = output_cfg
        with ZipFile(self.config.rban_jar) as zf:
            self.monomers_db = json.loads(zf.read('molecules/monomer/nrproMonomers.json'))  # default db

        if custom_monomers is not None:
            start_id = 1 + max(m['id'] for m in self.monomers_db)
            self.monomers_db.extend([{**monomer_data, 'id': i}
                                     for i, monomer_data in enumerate(custom_monomers, start=start_id)])

    def run_rban(self,
                 smiles_with_ids: List[dict],  # [{'id': id, 'smiles': smiles}]
                 log: NerpaLogger,
                 input_file: Optional[Path] = None,
                 output_file_name: Optional[str] = None,
                 report_not_processed: bool = False) -> List[Raw_rBAN_Record]:
        if input_file is None:
            input_file = self.output_cfg.default_input_file
        if output_file_name is None:
            output_file_name = self.output_cfg.default_output_file_name

        json.dump(smiles_with_ids, input_file.open('w'))
        json.dump(self.monomers_db, self.output_cfg.default_monomers_file.open('w'))
        command = ['java', '-jar', str(self.config.rban_jar),
                   '-inputFile', str(input_file),
                   '-outputFolder', str(self.output_cfg.rban_output_dir) + '/',
                   '-outputFileName', str(output_file_name),
                   '-monomersDB', str(self.output_cfg.default_monomers_file)]  # TODO: refactor

        nerpa_utils.sys_call(command, log)

        rban_records = json.load((self.output_cfg.rban_output_dir / output_file_name).open('r'))

        if report_not_processed:
            ids_in = set(record['id'] for record in smiles_with_ids)
            ids_out = set(record['id'] for record in rban_records)
            for idx in ids_in - ids_out:
                log.warning(f'No rBAN output for structure "{idx}"')

        return rban_records

    def get_hybrid_monomers_per_record(self,
                                       rban_records: List[Raw_rBAN_Record],
                                       log: NerpaLogger) -> Dict[int, Dict[MonomerIdx, NorineMonomerName]]:
        # check all unrecognized monomers for PK involvement
        def joined_id(struct_id: int, monomer_id: int) -> str:
            return f'{struct_id}_{monomer_id}'

        def parse_joined_id(joined_id: str) -> Tuple[int, int]:
            return tuple(map(int, joined_id.rsplit('_', 1)))

        trimmed_monomers_per_struct = {}
        for struct_idx, rban_record in enumerate(rban_records):
            try:
                trimmed_monomers_per_struct[struct_idx] = get_hybrid_monomers_smiles(rban_record)
            except:
                log.warning(
                    f'\nStructure "{struct_idx}": unexpected error while resolving NRP-PK hybrid monomer '
                    f'candidate. Skipping.')

        trimmed_monomers = [{'id': joined_id(struct_idx, monomer_idx), 'smiles': trimmed_smiles}
                            for struct_idx in trimmed_monomers_per_struct
                            for monomer_idx, trimmed_smiles in trimmed_monomers_per_struct[struct_idx]]

        if not trimmed_monomers:
            return defaultdict(dict)

        # run rBAN on the trimmed monomers
        log.info('\n=== Resolving NRP-PK hybrid monomers candidates')
        self.run_rban(trimmed_monomers,
                      input_file=self.output_cfg.putative_hybrids_input_file,
                      output_file_name=self.output_cfg.putative_hybrids_output_file_name,
                      log=log)

        # Add newly recognized monomers to the dictionary
        hybrid_monomers_dict = defaultdict(dict)
        new_monomers_processed = json.load((self.output_cfg.rban_output_dir /
                                            Path(self.output_cfg.putative_hybrids_output_file_name)).open('r'))
        for rban_record in new_monomers_processed:
            struct_id, monomer_id = parse_joined_id(rban_record['id'])
            aa_smi = rban_record["isomericSmiles"]
            aa_code = rban_record["monomericGraph"]["monomericGraph"]['monomers'][0]['monomer']['monomer'][
                'monomer']
            if not aa_code.startswith('X'):  # rBAN managed to recognize the trimmed monomer
                hybrid_monomers_dict[struct_id][monomer_id] = NorineMonomerName(aa_code)

        return hybrid_monomers_dict
