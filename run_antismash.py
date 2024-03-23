import os
from src.nerpa_pipeline.NRPSPredictor_utils.json_handler import get_main_json_fpath
from src.nerpa_pipeline.NRPSPredictor_utils.main import main as convert_antiSMASH_v5

def get_antismash_v3_compatible_input_paths(listing_fpath, list_of_paths, output_dir, log):
    '''
    Parses all antiSMASH-related options,
    detects all relevant output dirs (either with .json [aS v.5] or with ./txt/ & ./nrpspks_predictions_txt [aS v.3],
    converts aS v.5 to aS v.3-compliants if needed,
    returns list of paths to each v.3-compliant directory
    :param args:
    :return:
    '''

    def _get_input_antiSMASH_paths(lookup_paths):
        def _is_antiSMASHv3_path(path):
            if os.path.isdir(path) and \
                    os.path.isdir(os.path.join(path, 'txt')) and \
                    os.path.isdir(os.path.join(path, 'nrpspks_predictions_txt')):
                return True
            return False

        def _is_antiSMASHv5_path(path):
            if os.path.isfile(path) and path.endswith('.json'):
                return True
            if os.path.isdir(path) and get_main_json_fpath(dirpath=path) is not None:
                return True
            return False

        antiSMASHv3_paths = []
        antiSMASHv5_paths = []
        for entry in lookup_paths:
            if _is_antiSMASHv3_path(entry):
                antiSMASHv3_paths.append(entry)
            elif _is_antiSMASHv5_path(entry):
                antiSMASHv5_paths.append(entry)
            elif os.path.isdir(entry):
                # excluding dirs in runtime in os.walk: https://stackoverflow.com/questions/19859840/excluding-directories-in-os-walk
                for root, dirs, files in os.walk(entry, topdown=True):
                    # always ignore files since a single json in a dir should be caught one step before when path was the dir
                    # (see _is_antiSMASHv5_path() ) while multiple jsons in a dir probably means a false positive
                    dirs_to_keep_walking = []
                    for dir in dirs:
                        full_dir_path = os.path.join(root, dir)
                        if _is_antiSMASHv3_path(full_dir_path):
                            antiSMASHv3_paths.append(full_dir_path)
                        elif _is_antiSMASHv5_path(full_dir_path):
                            antiSMASHv5_paths.append(full_dir_path)
                        else:
                            dirs_to_keep_walking.append(dir)
                    dirs[:] = dirs_to_keep_walking
        return antiSMASHv3_paths, antiSMASHv5_paths

    lookup_locations = []
    if listing_fpath is not None:
        with open(listing_fpath) as f:
            for path in f:
                lookup_locations.append(path.strip())

    if list_of_paths:
        lookup_locations += list_of_paths

    antiSMASHv3_paths, antiSMASHv5_paths = _get_input_antiSMASH_paths(lookup_locations)
    log.info("\n=== Genome predictions found: %d antiSMASH v3 inputs; %d antiSMASH v5 inputs" %
             (len(antiSMASHv3_paths), len(antiSMASHv5_paths)))
    if antiSMASHv5_paths:
        log.info("\n======= Preprocessing antiSMASH v5 inputs")
        converted_antiSMASH_v5_outputs_dir = os.path.join(output_dir, "converted_antiSMASH_v5_outputs")
        log.info(f'results will be in {converted_antiSMASH_v5_outputs_dir}', indent=1)
        converted_antiSMASH_v5_paths = convert_antiSMASH_v5(antiSMASHv5_paths +
                                                            ['-o', converted_antiSMASH_v5_outputs_dir, '-m', 'hybrid', "-n", "v3" ])
        antiSMASHv3_paths += converted_antiSMASH_v5_paths
        log.info("\n======= Done with Preprocessing antiSMASH v5 inputs")

    return antiSMASHv3_paths
