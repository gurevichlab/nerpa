#!/usr/bin/env python3

import sys
import os
import argparse
import csv
import site
import shutil

import nerpa_init
nerpa_init.init()

site.addsitedir(nerpa_init.python_modules_dir)

import nerpa_utils
import handle_TE
import handle_rban
import logger


def parse_args(log):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    genomic_group = parser.add_argument_group('Genomic input', 'description')
    genomic_group.add_argument("--predictions", "-p", nargs=1, dest="predictions",
                                     help="path to file with paths to prediction files", type=str)
    genomic_group.add_argument("--antismash_output_list", dest="antismash_out",
                               help="path to file with list of paths to antiSMASH output folders", type=str)
    genomic_group.add_argument("--antismash", "-a", dest="antismash", help="path to antiSMASH output or to directory with many antiSMASH outputs", type=str)

    struct_group = parser.add_argument_group('Chemical input', 'description')
    struct_input_group = struct_group.add_mutually_exclusive_group(required=True)
    struct_input_group.add_argument("--smiles", dest="smiles", nargs='*',
                        help="string (or several strings) with structures in the SMILES format", type=str)
    struct_input_group.add_argument("--smiles-tsv", dest="smiles_tsv",
                        help="multi-column file containing structures in the SMILES format", type=str)
    struct_input_group.add_argument("--graphs", dest="graphs", help="", type=str)
    struct_group.add_argument("--col_smiles", dest="col_smiles",
                        help="name of column with structures in the SMILES format [default: 'SMILES']", type=str, default='SMILES')
    struct_group.add_argument("--col_id", dest="col_id",
                        help="name of col with ids (if not provided, id=[number of row])", type=str)
    struct_group.add_argument("--sep", dest="sep",
                        help="separator in smiles-tsv", type=str, default='\t')

    # parser.add_argument("--insertion", help="insertion score [default=-2.8]", default=-2.8, action="store")
    # parser.add_argument("--deletion", help="deletion score [default=-5]", default=-5, action="store")
    parser.add_argument("--configs_dir", help="custom directory with configs", action="store", type=str)
    parser.add_argument("--threads", default=1, type=int, help="number of threads for running Nerpa", action="store")
    parser.add_argument("--output_dir", "-o", help="output dir [default: nerpa_results/results_<datetime>]",
                        type=str, default=None)

    parsed_args = parser.parse_args()

    validate_arguments(parsed_args, parser, log)
    return parsed_args


def check_tsv_ids_duplicates(reader, col_id):
    ids = [row[col_id] for row in reader]
    return len(ids) == len(set(ids))


class ValidationError(Exception):
    pass


def validate(expr, msg=''):
    if not expr:
        raise ValidationError(msg)


def validate_arguments(args, parser, log):
    try:
        if not (args.predictions or args.antismash or args.antismash_out):
            raise ValidationError(f'one of the arguments --predictions --antismash/-a --antismash_output_list is required')
        if args.predictions and (args.antismash or args.antismash_out):
            raise ValidationError(f'argument --predictions: not allowed with argument --antismash/-a or --antismash_output_list')

        if args.smiles_tsv:
            try:
                with open(args.smiles_tsv, newline='') as f_in:
                    reader = csv.DictReader(f_in, delimiter=args.sep, quoting=csv.QUOTE_NONE)
                    validate(args.col_smiles in reader.fieldnames,
                        f'Column "{args.col_smiles}" was specified but does not exist in {args.smiles_tsv}.')
                    if args.col_id:
                        validate(args.col_id in reader.fieldnames,
                            f'Column "{args.col_id}" was specified but does not exist in {args.smiles_tsv}.')
                        validate(check_tsv_ids_duplicates(reader, args.col_id),
                            f'Duplicate IDs are found.')
            except FileNotFoundError:
                raise ValidationError(f'No such file: "{args.smiles_tsv}".')
            except csv.Error as e:
                raise ValidationError(f'Cannot parse "{args.smiles_tsv}": {e}.')
            except Exception as e:
                raise ValidationError(f'Invalid input file "{args.smiles_tsv}": {e}.')

    except ValidationError as e:
        parser.print_help()
        error_msg = f'{e}\n' if str(e) else 'Options validation failed!'
        log.error(error_msg, to_stderr=True)


def gen_graphs_from_smiles_tsv(args, main_out_dir,
                               path_to_monomers_tsv, path_to_graphs, path_to_rban_jar, path_to_monomers,
                               log):
    path_to_rban_input = os.path.join(main_out_dir, 'rban.input.json')
    if args.smiles_tsv:
        handle_rban.generate_rban_input_from_smiles_tsv(
            args.smiles_tsv, path_to_rban_input, sep=args.sep, id_col_name=args.col_id, smi_col_name=args.col_smiles)
    else:
        handle_rban.generate_rban_input_from_smiles_string(args.smiles, path_to_rban_input)

    path_to_rban_output = os.path.join(main_out_dir, 'rban.output.json')
    log.info('\n======= Structures preprocessing with rBAN')
    command = ['java', '-jar', path_to_rban_jar,
               '-monomersDB', path_to_monomers,
               '-inputFile', path_to_rban_input,
               '-outputFolder', main_out_dir + '/',  # rBAN specifics
               '-outputFileName', os.path.basename(path_to_rban_output)]
    nerpa_utils.sys_call(command, log)

    handle_rban.generate_graphs_from_rban_output(path_to_rban_output, path_to_monomers_tsv, path_to_graphs,
                                                 main_out_dir, path_to_rban_jar, log)


def copy_prediction_list(args, main_out_dir):
    new_prediction_path = os.path.join(main_out_dir, "prediction.info")
    with open(new_prediction_path, 'w') as f:
        with open(args.predictions[0]) as fr:
            for line in fr:
                line_parts = line.split()
                file = line_parts[0]
                if (file[0] != '/'):
                    file = os.path.join(os.path.dirname(os.path.abspath(args.predictions[0])), file)
                f.write(file + "\n")
    return new_prediction_path


def get_antismash_list(args):
    aouts = []
    if args.antismash_out is not None:
        with open(args.antismash_out) as f:
            for dir in f:
                aouts.append(dir)

    is_one = False
    if args.antismash is not None:
        for filename in os.listdir(args.antismash):
            if filename.endswith(".json"):
                is_one = True
            if filename == "nrpspks_predictions_txt":
                is_one = True
        if is_one:
            aouts.append(args.antismash)
        else:
            for filename in os.listdir(args.antismash):
                aouts.append(os.path.join(args.antismash, filename))

    return aouts


def run(args, log):
    output_dir = nerpa_utils.set_up_output_dir(output_dirpath=args.output_dir)
    log.set_up_file_handler(output_dir)
    log.start()

    if args.predictions is not None:
        path_predictions = copy_prediction_list(args, output_dir)
    else:
        path_predictions = handle_TE.create_predictions_by_antiSMASHout(get_antismash_list(args), output_dir)

    input_configs_dir = args.configs_dir if args.configs_dir else nerpa_init.configs_dir
    current_configs_dir = os.path.join(output_dir, "configs")
    # remember shutil.copytree caveat (compared to dir_util.copy_tree):
    # directory metadata will be copied that may cause potential problems
    # if src dir is too old and there is a cluster cronjob which
    # automatically remove old files from the temporary workspace
    shutil.copytree(input_configs_dir, current_configs_dir, copy_function=shutil.copy)

    path_to_graphs = os.path.join(output_dir, 'path_to_graphs')
    local_monomers_cfg = os.path.join(current_configs_dir, "monomers.tsv")
    if args.graphs:
        shutil.copyfile(args.graphs, path_to_graphs)
    else:
        path_to_rban = os.path.join(nerpa_init.external_tools_dir, 'rBAN', 'rBAN-1.0.jar')
        path_to_monomers = os.path.join(nerpa_init.external_tools_dir, 'rBAN', 'nrproMonomers_nerpa.json')
        gen_graphs_from_smiles_tsv(args, output_dir,
                                   local_monomers_cfg, path_to_graphs, path_to_rban, path_to_monomers,
                                   log)

    details_mol_dir = os.path.join(output_dir, 'details_mols')
    if not os.path.exists(details_mol_dir):
        os.makedirs(details_mol_dir)

    command = [os.path.join(nerpa_init.bin_dir, "NRPsMatcher"),
               path_predictions, path_to_graphs, '--configs_dir', current_configs_dir]
    log.info("\n======= Nerpa matching")
    nerpa_utils.sys_call(command, log, cwd=output_dir)
    log.info("RESULTS:")
    log.info("Main report is saved to " + os.path.join(output_dir, 'report.csv'), indent=1)
    log.info("Detailed reports are saved to " + output_dir, indent=1)
    log.finish()


if __name__ == "__main__":
    log = logger.NerpaLogger()
    try:
        args = parse_args(log)
        run(args, log)
    except Exception:
        _, exc_value, _ = sys.exc_info()
        log.exception(exc_value)
    finally:
        # TODO: clean up: remove all intermediate files
        pass
