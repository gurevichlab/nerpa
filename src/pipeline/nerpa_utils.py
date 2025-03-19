import os
import re
import shutil
import datetime
import shlex
import subprocess
from src.config import OutputConfig
from src.pipeline.logger import NerpaLogger
from pathlib import Path


def set_up_output_dir(output_cfg: OutputConfig,
                      crash_if_exists: bool,
                      log: NerpaLogger,
                      make_sym_link_to_latest: bool = True):
    if output_cfg.main_out_dir.exists():
        if crash_if_exists:
            log.error(f"output directory ({output_cfg.main_out_dir}) already exists! "
                      f"Rerun with --force-output-dir if you still want to use it as the output dir "
                      f"OR specify another directory. Exiting now..", to_stderr=True)
            exit(1)
        if output_cfg.main_out_dir.is_dir():  # TODO: check whether we want to completely remove it always
            shutil.rmtree(output_cfg.main_out_dir)
        else:
            output_cfg.main_out_dir.unlink()

    output_cfg.main_out_dir.mkdir(parents=True)

    # 'latest' symlink
    if make_sym_link_to_latest:
        if output_cfg.symlink_to_latest.is_symlink():
            output_cfg.symlink_to_latest.unlink(missing_ok=True)
        output_cfg.symlink_to_latest.symlink_to(output_cfg.main_out_dir, target_is_directory=True)

    return output_cfg.main_out_dir


def sys_call(cmd, log, indent='  ', cwd=None, verbose=True):
    def _process_readline(line):
        return str(line, "utf-8").rstrip()

    if isinstance(cmd, list):
        cmd_list = cmd
    else:
        cmd_list = shlex.split(cmd)

    if verbose:
        log.info("\n== Running: %s\n" % (' '.join(cmd_list)))
    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=cwd)

    # output = ""
    while not proc.poll():
        line = _process_readline(proc.stdout.readline())
        if line:
            log.info(indent + line)
            # if log:
            #     log.info(line)
            # else:
            #     output += line + "\n"
        if proc.returncode is not None:
            break

    for line in proc.stdout.readlines():
        line = _process_readline(line)
        if line:
            log.info(indent + line)
            # if log:
            #     log.info(line)
            # else:
            #     output += line + "\n"

    if proc.returncode:
        log.error("system call for: \"%s\" finished abnormally, OS return value: %d" % (cmd, proc.returncode))

    if verbose:
        log.info("\n== Done\n")
    # return output


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)

    cur_path = os.path.dirname(os.path.abspath(__file__))
    if is_exe(os.path.join(cur_path, fname)):
        return os.path.join(cur_path, fname)

    if fpath:
        if is_exe(program):
            return program
    elif "PATH" in os.environ:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def get_path_to_program(program, dirpath=None, min_version=None):
    """
    returns the path to an executable or None if it can't be found
    """
    def is_exe(fpath):
        if os.path.isfile(fpath) and os.access(fpath, os.X_OK):
            if not min_version or check_version(fpath, min_version):
                return True

    def check_version(fpath, min_version):
        p = subprocess.Popen([fpath, '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = p.communicate()

        version_pattern = re.compile(r'(?P<major_version>\d+)\.(?P<minor_version>\d+)')
        searchstring = stdout.decode('utf8').strip()

        # ad hoc workaround to AS 5.2.0 printing FutureWarning to stdout
        searchstring = searchstring.split('\n')[-1]

        v = version_pattern.search(searchstring)
        if not v.group('major_version') or not v.group('minor_version'):
            return False
        version, minor_version = map(int, min_version.split('.'))
        if int(v.group('major_version')) == version and int(v.group('minor_version')) >= minor_version:
            return True

    if dirpath:
        exe_file = os.path.join(dirpath, program)
        if is_exe(exe_file):
            return exe_file
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None
