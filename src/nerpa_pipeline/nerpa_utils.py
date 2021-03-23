import os
import datetime
import shlex
import subprocess

import nerpa_config
import logger


def set_up_output_dir(output_dirpath):
    make_latest_symlink = False

    if not output_dirpath:  # 'output dir was not specified with -o option'
        output_dirpath = os.path.join(os.path.abspath(
            nerpa_config.default_results_root_dirname),
            nerpa_config.default_results_dirname_prefix +
            datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S'))
        make_latest_symlink = True

        # in case of starting two jobs in the same second
        if os.path.isdir(output_dirpath):
            i = 2
            base_dirpath = output_dirpath
            while os.path.isdir(output_dirpath):
                output_dirpath = str(base_dirpath) + '__' + str(i)
                i += 1

    if not os.path.isdir(output_dirpath):
        os.makedirs(output_dirpath)

    # 'latest' symlink
    if make_latest_symlink:
        prev_dirpath = os.getcwd()
        os.chdir(nerpa_config.default_results_root_dirname)

        latest_symlink = 'latest'
        if os.path.islink(latest_symlink):
            os.remove(latest_symlink)
        os.symlink(os.path.basename(output_dirpath), latest_symlink)

        os.chdir(prev_dirpath)

    return os.path.abspath(output_dirpath)


def sys_call(cmd, indent='  ', cwd=None, verbose=True):
    def _process_readline(line):
        return str(line, "utf-8").rstrip()

    if isinstance(cmd, list):
        cmd_list = cmd
    else:
        cmd_list = shlex.split(cmd)

    if verbose:
        logger.info("\n== Running: %s\n" % (' '.join(cmd_list)))
    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=cwd)

    # output = ""
    while not proc.poll():
        line = _process_readline(proc.stdout.readline())
        if line:
            logger.info(indent + line)
            # if log:
            #     log.info(line)
            # else:
            #     output += line + "\n"
        if proc.returncode is not None:
            break

    for line in proc.stdout.readlines():
        line = _process_readline(line)
        if line:
            logger.info(indent + line)
            # if log:
            #     log.info(line)
            # else:
            #     output += line + "\n"

    if proc.returncode:
        logger.error("system call for: \"%s\" finished abnormally, OS return value: %d" % (cmd, proc.returncode))

    if verbose:
        logger.info("\n== Done\n")
    # return output
