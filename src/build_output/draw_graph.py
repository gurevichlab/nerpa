def draw_rban_record(rban_record: Parsed_rBAN_Record, out_file: Path):
    draw_graph(MonomerGraph.from_rban_record(rban_record), out_dir / Path(f'{rban_record.compound_id}.png'))