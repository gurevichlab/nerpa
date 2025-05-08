from pathlib import Path
import argparse
from argparse import Namespace as CommandLineArgs
from src.config import Config
from io import StringIO
import csv

def add_genomic_arguments(parser: argparse.ArgumentParser):
    genomic_group = parser.add_argument_group('Genomic input', 'Genomes of NRP-producing organisms or BGC predictions')
    genomic_group.add_argument("--antismash", "-a", dest="antismash", action='append', type=Path,
                               metavar='DIR',
                               help="single antiSMASH output directory or directory with many antiSMASH outputs inside")
    genomic_group.add_argument("--antismash-paths-file", dest="antismash_outpaths_file", metavar='FILE',
                               help="file with list of paths to antiSMASH output directories", type=Path)
    genomic_group.add_argument("--antismash-job-ids", dest="antismash_job_ids", nargs='*', metavar='STR',
                               help="space-separated job IDs to download from the antiSMASH webserver", type=str)
    genomic_group.add_argument("--genome", dest="seqs", action='append', metavar='FILE',
                               help="genome sequence in the GenBank/EMBL/FASTA format", type=Path)
    genomic_group.add_argument("--paras-results",
                               help="path to a directory with PARAS results", type=Path)


def add_struct_arguments(parser: argparse.ArgumentParser):
    struct_group = parser.add_argument_group('Chemical input', 'Structures of NRP molecules')
    struct_input_group = struct_group.add_mutually_exclusive_group()
    struct_input_group.add_argument("--smiles", dest="smiles", nargs='*', metavar='STR',
                                    help="space-separated structures in the SMILES format", type=str)
    struct_input_group.add_argument("--smiles-tsv", dest="smiles_tsv", metavar='FILE',
                                    help="multi-column file with structures in the SMILES format and metadata", type=Path)
    struct_group.add_argument("--col-smiles", dest="col_smiles", metavar='STR',
                              help="column name in smiles-tsv for structures in the SMILES format [default: '%(default)s']",
                              type=str, default='SMILES')
    struct_group.add_argument("--col-id", dest="col_id", metavar='STR',
                              help="column name in smiles-tsv for structure identifier [if not provided, row index will be used]",
                              type=str)
    struct_group.add_argument("--sep", dest="sep", metavar='CHAR',
                              help="column separator in smiles-tsv [default: '\\t']", type=str, default='\t')
    struct_input_group.add_argument("--rban-json", dest="rban_output", metavar='FILE',
                                    help="rBAN-preprocessed NRP structures in the JSON file", type=Path)

def add_advanced_arguments(parser: argparse.ArgumentParser):
    advanced_input_group = parser.add_argument_group('Advanced input',
                                                     'Preprocessed data '
                                                     'in custom Nerpa-compliant formats')
    advanced_input_group.add_argument("--bgc-variants", dest="bgc_variants", metavar='DIR',
                                      help="directory with predicted BGC variants (yaml files)", type=Path)
    advanced_input_group.add_argument("--nrp-variants", dest="nrp_variants", metavar='DIR',
                                      help="directory with predicted NRP variants (yaml files)", type=Path)
    advanced_input_group.add_argument("--configs-dir", help="custom directory with Nerpa configs",
                                      metavar='DIR', action="store", type=Path)
    advanced_input_group.add_argument('--rban-monomers-db', dest='rban_monomers', type=Path, default=None,
                                      metavar='FILE', help='file with custom monomers in rBAN compatible format')

def add_debug(parser: argparse.ArgumentParser):
    debug_input_group = parser.add_argument_group('Debug options',
                                                     'Tweak Nerpa behavior for debugging purposes ')
    debug_input_group.add_argument('--disable-calibration', action='store_true',
                                   help='specificity predictions will be used as is, without calibration')
    debug_input_group.add_argument('--disable-dictionary-lookup', action='store_true',
                                   help='do not use dictionary of known A domain specificities')
    debug_input_group.add_argument('--draw-hmms', action='store_true',
                                   help='draw HMMs with optimal paths for all matches')

def add_pipeline_arguments(parser: argparse.ArgumentParser, default_cfg: Config):
    configs_group = parser.add_argument_group('Nerpa pipeline',
                                              'Nerpa running configuration')

    configs_group.add_argument("--output-dir", "-o", type=Path, metavar='DIR',
                               help="output directory "
                                    "[default: {CWD}/"
                                    f"{default_cfg.output_config.main_out_dir.parent.name}/" "{CURRENT_TIME}]")
    configs_group.add_argument("--force-output-dir", dest="output_dir_reuse", action="store_true",
                               help="do not crash if the output directory already exists and rewrite its content")

    configs_group.add_argument("--threads", "-t", default=1, type=int, metavar='INT',
                               help="number of threads for running Nerpa [default: %(default)s]", action="store")

    configs_group.add_argument("--process-hybrids", dest="process_hybrids", action="store_true", default=False,
                               help="process NRP-PK hybrid monomers (requires the use of rBAN)")

    configs_group.add_argument("--antismash-installation-dir", dest="antismash_path", type=Path,
                               default=None, metavar='DIR',
                               help="path to the antiSMASH installation directory, i.e., the one containing the "
                                    "'run_antismash.py' script")

    configs_group.add_argument("--max-num-matches-per-bgc", default=None, type=int, metavar='INT',
                               help="maximum number of matches to report per BGC; set 0 for unlimited "
                                    f"[default: {default_cfg.matching_config.max_num_matches_per_bgc}]",
                               action="store")
    configs_group.add_argument("--max-num-matches-per-nrp", default=None, type=int, metavar='INT',
                               help="maximum number of matches to report per NRP; set 0 for unlimited "
                                    f"[default: {default_cfg.matching_config.max_num_matches_per_nrp}]",
                               action="store")
    configs_group.add_argument("--max-num-matches", default=None, type=int, metavar='INT',
                               help="maximum number of matches to report in total; set 0 for unlimited "
                                    f"[default: {default_cfg.matching_config.max_num_matches}]",
                               action="store")
    #configs_group.add_argument("--heuristic-discard", default=False,
    #                           help="immediately discard bad matches based on heuristics", action="store_true")
    configs_group.add_argument("--skip-molecule-drawing",
                               action="store_true", default=False,
                               help="do not draw NRP molecules and monomer graphs "
                                    "(faster and saves space but they will be missing in the HTML report)")

    configs_group.add_argument("--fast-matching",
                               action="store_true", default=False,
                               help="use C++ executable to perform matching (requires compilation)")

    # configs_group.add_argument("--only-preprocessing", action="store_true", default=False,
    #                            help="only generate NRP and BGC variants, do not perform matching (useful for debugging)")
    configs_group.add_argument("--debug", action="store_true", default=False,
                               help="run in the debug mode and keep all intermediate files")


def build_cmdline_args_parser(default_cfg: Config) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    add_genomic_arguments(parser)
    add_struct_arguments(parser)
    add_advanced_arguments(parser)
    add_pipeline_arguments(parser, default_cfg)
    add_debug(parser)
    return parser


def get_command_line_args(default_cfg: Config) -> CommandLineArgs:
    parser = build_cmdline_args_parser(default_cfg)
    parsed_args = parser.parse_args()
    try:
        validate_arguments(parsed_args)
    except ValidationError as e:
        help_message = StringIO()
        parser.print_help(help_message)
        error_message = str(e) if str(e) else 'Options validation failed!'
        raise ValidationError(error_message + '\n' + help_message.getvalue())
    return parsed_args


def check_tsv_ids_duplicates(reader, col_id):
    from collections import defaultdict
    counts = defaultdict(int)
    for row in reader:
        counts[row[col_id]] += 1
    duplicates = [(k, v) for k,v in counts.items() if v > 1]
    return duplicates


class ValidationError(Exception):
    pass


def validate(expr, msg=''):
    if not expr:
        raise ValidationError(msg)


def validate_arguments(args):  # TODO: I think it all could be done with built-in argparse features
    if not any([args.bgc_variants,
                args.antismash,
                args.antismash_outpaths_file,
                args.antismash_job_ids,
                args.seqs]):
        raise ValidationError(f'at least one genome/BGC input is required')
    if args.bgc_variants and (args.antismash or args.antismash_outpaths_file or args.antismash_job_ids or args.seqs):
        # TODO: what's wrong with having both?
        raise ValidationError(f'argument --bgc-variants is not compatible with other genome/BGC input options')
    if not any([args.nrp_variants,
                args.smiles,
                args.smiles_tsv,
                args.rban_output]):
        raise ValidationError(f'at least one NRP structure input is required')
    if args.nrp_variants and (args.smiles or args.smiles_tsv or args.rban_output):
        # TODO: what's wrong with having both?
        raise ValidationError('argument --nrp-variants is not compatible with other NRP input options')
    if args.smiles_tsv:
        try:
            with open(args.smiles_tsv, newline='') as f_in:
                reader = csv.DictReader(f_in, delimiter=args.sep, quoting=csv.QUOTE_NONE)
                validate(args.col_smiles in reader.fieldnames,
                         f'Column "{args.col_smiles}" was specified but does not exist in {args.smiles_tsv}.')
                if args.col_id:
                    validate(args.col_id in reader.fieldnames,
                             f'Column "{args.col_id}" was specified but does not exist in {args.smiles_tsv}.')
                    duplicates = check_tsv_ids_duplicates(reader, args.col_id)
                    validate(len(duplicates) == 0, f'Duplicate IDs are found: {duplicates}')
        except FileNotFoundError:
            raise ValidationError(f'No such file: "{args.smiles_tsv}".')
        except csv.Error as e:
            raise ValidationError(f'Cannot parse "{args.smiles_tsv}": {e}.')
        except Exception as e:
            raise ValidationError(f'Invalid input file "{args.smiles_tsv}": {e}.')


