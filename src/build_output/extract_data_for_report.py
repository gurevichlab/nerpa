from typing import NamedTuple
import graphviz
import json
from pathlib import Path
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.build_output.draw_graph import (
    GraphDrawingHelper,
    MoleculeDrawingHelper
)
from src.monomer_names_helper import MonomerNamesHelper


class GeneratedNRPs_DataForReport:
    raw_data_filtered: dict[str, dict]
    generated_nrps_graph_data: dict[str, graphviz.Digraph]
    generated_nrps_molecule_data: dict[str, dict]

    def __init__(
            self,
            generated_nrps_path: Path,
            monomer_names_helper: MonomerNamesHelper,
            ids_to_keep: set[str] | None = None,
    ):
        with generated_nrps_path.open(encoding="utf-8") as input_file:
            output_items = json.load(input_file)

        # Preserve the same structure as output_items,
        # but filter the variants based on ids_to_keep
        self.raw_data_filtered = {}
        for item_id, item in output_items.items():
            filtered_nrps_for_match = {
                variant_id: variant_data
                for variant_id, variant_data in item["new_variants"].items()
                if ids_to_keep is None or variant_id in ids_to_keep
            }
            new_item = item.copy()
            new_item["new_variants"] = filtered_nrps_for_match
            if new_item["new_variants"]:
                self.raw_data_filtered[item_id] = new_item
                
                
        new_records: list[Parsed_rBAN_Record] = [
            Parsed_rBAN_Record.from_dict(variant["new_record"])
            for item in self.raw_data_filtered.values()
            for variant in item["new_variants"].values()
        ]

        self.generated_nrps_graph_data: dict[str, graphviz.DiGraph] = {
            record.compound_id:
            GraphDrawingHelper(
                record,
                monomer_names_helper=monomer_names_helper
            ).render()
            for record in new_records
        }

        self.generated_nrps_molecule_data: dict[str, dict] = {
            record.compound_id:
            MoleculeDrawingHelper(
                record,
                monomer_names_helper=monomer_names_helper
            ).get_drawing_data()
            for record in new_records
        }
            

            
