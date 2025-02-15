from pathlib import Path
import argparse
from argparse import Namespace as CommandLineArgs
from src.config import Config
from io import StringIO
import csv

def add_genomic_arguments(parser: argparse.ArgumentParser):
    genomic_group = parser.add_argument_group('Genomic input', 'Genomes of NRP-producing organisms (i.e. BGC predictions)')
    genomic_group.add_argument("--antismash_outpaths_file", dest="antismash_outpaths_file",
                               help="file with list of paths to antiSMASH output directories", type=Path)
    genomic_group.add_argument("--antismash", "-a", dest="antismash", action='append', type=Path,
                               help="single antiSMASH output directory or directory with many antiSMASH outputs")
    genomic_group.add_argument("--antismash-job-ids", dest="antismash_job_ids", nargs='*',
                               help="job IDs for antiSMASH results to download", type=str)
    genomic_group.add_argument("--sequences", dest="seqs",
                               help="GenBank/EMBL/FASTA file containing DNA sequences", type=Path)


def add_struct_arguments(parser: argparse.ArgumentParser):
    struct_group = parser.add_argument_group('Chemical input', 'Structures of NRP molecules')
    struct_input_group = struct_group.add_mutually_exclusive_group()
    struct_input_group.add_argument("--rban-json", dest="rban_output",
                                    help="json file with rBAN-preprocessed NRP structures", type=Path)
    struct_input_group.add_argument("--smiles", dest="smiles", nargs='*',
                                    help="string (or several strings) with structures in the SMILES format", type=str)
    struct_input_group.add_argument("--smiles-tsv", dest="smiles_tsv",
                                    help="multi-column file containing structures in the SMILES format", type=Path)
    struct_group.add_argument("--col-smiles", dest="col_smiles",
                              help="column name in smiles-tsv for structures in the SMILES format [default: 'SMILES']",
                              type=str, default='SMILES')
    struct_group.add_argument("--col-id", dest="col_id",
                              help="column name in smiles-tsv for structure identifier (if not provided, row index will be used)",
                              type=str)
    struct_group.add_argument("--sep", dest="sep",
                              help="column separator in smiles-tsv", type=str, default='\t')


def add_advanced_arguments(parser: argparse.ArgumentParser):
    advanced_input_group = parser.add_argument_group('Advanced input',
                                                     'Preprocessed BGC predictions and NRP structures '
                                                     'in custom Nerpa-compliant formats')
    advanced_input_group.add_argument("--bgc-variants", "-b", dest="bgc_variants",
                                      help="Folder with predicted BGC variants (yaml files)", type=Path)
    advanced_input_group.add_argument("--nrp-variants", "-s", dest="nrp_variants",
                                      help="Folder with predicted NRP variants (yaml files)", type=Path)
    advanced_input_group.add_argument("--configs_dir", help="custom directory with adjusted Nerpa configs", action="store",
                                      type=Path)
    advanced_input_group.add_argument("--force-existing-outdir", dest="output_dir_reuse", action="store_true",
                                      default=False,
                                      help="don't crash if the output dir already exists")
    #configs_group.add_argument("--only-preprocessing", action="store_true", default=False,
    #                           help="only generate NRP and BGC variants, do not perform matching (useful for debugging)")
    advanced_input_group.add_argument("--debug", action="store_true", default=False,
                                      help="run in the debug mode (keep intermediate files)")


def add_config_arguments(parser: argparse.ArgumentParser, default_cfg: Config):
    configs_group = parser.add_argument_group('Nerpa config',
                                              'Nerpa running configuration')
    configs_group.add_argument('--rban-monomers-db', dest='rban_monomers', type=Path, default=None,
                        help='file with custom monomers in rBAN compatible format')
    configs_group.add_argument("--process-hybrids", dest="process_hybrids", action="store_true", default=False,
                        help="process NRP-PK hybrid monomers (requires use of rBAN)")
    configs_group.add_argument("--threads", default=1, type=int,
                               help="number of threads for running Nerpa", action="store")
    configs_group.add_argument("--max-num-matches-per-bgc", default=None, type=int,
                               help="maximum number of matches to report per BGC. "
                                    f"By default {default_cfg.matching_config.max_num_matches_per_bgc}",
                               action="store")
    configs_group.add_argument("--max-num-matches-per-nrp", default=None, type=int,
                               help="maximum number of matches to report per NRP. "
                                    f"By default {default_cfg.matching_config.max_num_matches_per_nrp}",
                               action="store")
    configs_group.add_argument("--max-num-matches", default=None, type=int,
                               help="maximum number of matches to report in total. "
                                    f"By default {default_cfg.matching_config.max_num_matches}",
                               action="store")
    #configs_group.add_argument("--heuristic-discard", default=False,
    #                           help="immediately discard bad matches based on heuristics", action="store_true")
    configs_group.add_argument("--dont-draw-molecules", action="store_true", default=False,
                               help="do not draw NRP molecule and monomer graps (they will not appear in the report)")


def add_external_tools_args(parser: argparse.ArgumentParser):
    parser.add_argument('--antismash-path', dest='antismash_path', type=Path, default=None,
                        help='path to antismash source directory')


def build_cmdline_args_parser(default_cfg: Config) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    add_genomic_arguments(parser)
    add_struct_arguments(parser)
    add_external_tools_args(parser)
    add_advanced_arguments(parser)
    add_config_arguments(parser, default_cfg)

    default_out_dir = default_cfg.output_config.main_out_dir.parent
    # q: default_out_dir has format '.../results_%Y-%m-%d_%H-%M-%S', but we want '.../results'
    parser.add_argument("--output_dir", "-o",
                        help="output directory. "
                             f"If not provided, results will be saved in {default_out_dir}/" "{CURRENT_TIME}",
                        type=Path)
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
                args.seqs]):
        raise ValidationError(f'one of the arguments --predictions --antismash/-a --antismash_output_list '
                              f'--sequences is required')
    if args.bgc_variants and (args.antismash or args.antismash_outpaths_file or args.seqs):
        # TODO: what's wrong with having both?
        raise ValidationError(f'argument --predictions: not allowed with argument --antismash/-a '
                              f'or --antismash_output_list or --sequences')
    if not any([args.nrp_variants,
                args.smiles,
                args.smiles_tsv,
                args.rban_output]):
        raise ValidationError(f'one of the arguments --rban-json --smiles-tsv --smiles --structures/-s'
                              f'is required')
    if args.nrp_variants and (args.smiles or args.smiles_tsv or args.rban_output):
        # TODO: what's wrong with having both?
        raise ValidationError('argument --structures/-s: not allowed with argument --rban-json or --smiles '
                              'or --smiles-tsv')
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


