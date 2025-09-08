import sys
import logging
import datetime
import platform
import os
from dataclasses import dataclass
from logging import Logger
from pathlib import Path
from typing import Literal

UNRECOVERABLE_ERROR_MSG = ("An unrecoverable error has occurred, and Nerpa cannot continue.\n"
                           "In case you have troubles running our tool, "
                           "you can post an issue on https://github.com/gurevichlab/nerpa/issues "
                           "or write to alexey.gurevich@helmholtz-hips.de")

@dataclass
class LoggingConfig:
    log_file: Path
    warnings_file: Path
    stdout_log_level: Literal['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    nerpa_version: str

    def __init__(self,
                 cfg: dict,
                 nerpa_dir: Path):
        self.log_file = nerpa_dir / cfg['log_file']
        self.warnings_file = nerpa_dir / cfg['warnings_file']

        stdout_log_level = cfg.get('stdout_log_level', 'INFO')
        if stdout_log_level not in ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'):
            raise ValueError(f'Invalid log level: {stdout_log_level}')
        self.stdout_log_level = stdout_log_level
        self.nerpa_version = (nerpa_dir / 'VERSION.txt').read_text().strip()

class PreliminaryLogger(Logger):
    """
    A simple logger used before the main NerpaLogger is set up.
    """
    def __init__(self):
        super().__init__(name='preliminary_nerpa_logger', level=logging.DEBUG)
        console_handler = logging.StreamHandler(sys.stderr)
        console_handler.setLevel(logging.ERROR)
        self.addHandler(console_handler)

    def exception(self, msg='', *args, **kwargs):
        if msg:
            msg += '\n\n'
        msg += UNRECOVERABLE_ERROR_MSG
        super().error(msg, *args, **kwargs)

        for handler in list(self.handlers):  # iterate over a copy
            self.removeHandler(handler)
            try:
                handler.close()
            except Exception:
                pass


class NerpaLogger(Logger):
    cfg: LoggingConfig
    _start_time: datetime.datetime
    _num_warnings: int
    _num_errors: int
    _indent = '    '


    def __init__(self, cfg: LoggingConfig):
        super().__init__(name='nerpa', level=logging.DEBUG)
        self.cfg = cfg
        self._num_warnings = 0
        self._num_errors = 0

        # Set up console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(cfg.stdout_log_level)
        self.addHandler(console_handler)

        # Set up file handlers:
        # verbose log (DEBUG+) and warnings log (WARNING+)
        self.cfg.log_file.parent.mkdir(parents=True, exist_ok=True)
        self.cfg.warnings_file.parent.mkdir(parents=True, exist_ok=True)

        verbose_handler = logging.FileHandler(self.cfg.log_file, mode='w', encoding='utf-8')
        verbose_handler.setLevel(logging.DEBUG)
        verbose_fmt = logging.Formatter(
            fmt='%(asctime)s [%(levelname)s] %(filename)s:%(lineno)d (%(funcName)s) - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        verbose_handler.setFormatter(verbose_fmt)
        self.addHandler(verbose_handler)

        warn_handler = logging.FileHandler(self.cfg.warnings_file, mode='w', encoding='utf-8')
        warn_handler.setLevel(logging.WARNING)
        warn_fmt = logging.Formatter(
            fmt='%(asctime)s [%(levelname)s] %(filename)s:%(lineno)d - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        warn_handler.setFormatter(warn_fmt)
        self.addHandler(warn_handler)

    def info(self, msg: str, *args, **kwargs):
        indent = kwargs.pop('indent', 0)
        if indent > 0:
            msg = self._indent * indent + msg.replace('\n', '\n' + self._indent * indent)
        super().info(msg, *args, **kwargs)

    def warning(self, *args, **kwargs):
        super().warning(*args, **kwargs)
        self._num_warnings += 1

    def error(self, *args, **kwargs):
        super().error(*args, **kwargs)
        self._num_errors += 1

    def log_command_line_args(self):
        args = sys.argv[:]
        for i, arg in enumerate(args):
            if ' ' in arg or '\t' in arg:
                args[i] = "'" + arg + "'"
        self.info('')
        self.info('Started with command: ' + ' '.join(args))

    def log_system_info(self):
        self.info('')
        self.info("System information:")
        self.info("OS: " + platform.platform(), indent=1)
        self.info("Python version: " + '.'.join(map(str, sys.version_info[:2])), indent=1)
        try:
            import multiprocessing
            self.info("CPUs number: " + str(multiprocessing.cpu_count()), indent=1)
        except ImportError:
            self.info("Problem occurred when getting CPUs number information", indent=1)

    def unrecoverable_error(self, msg='', *args, **kwargs):
        if msg:
            msg += '\n\n'
        msg += UNRECOVERABLE_ERROR_MSG
        self.error(msg, *args, **kwargs)

    def log_tool_version(self):
        self.info('\nNerpa version: ' + self.cfg.nerpa_version)

    def start(self):
        self.log_command_line_args()
        self.log_tool_version()
        self.log_system_info()

        self._start_time = datetime.datetime.now()
        self.info('Started: ' + self._start_time.strftime("%Y-%m-%d %H:%M:%S"))
        self.info(f'Logging to {self.cfg.log_file} (verbose) and {self.cfg.warnings_file} (warnings)')

    def finish(self):
        self.info(f'Verbose log is saved to {self.cfg.log_file}')
        self.info(f'Warnings log is saved to {self.cfg.warnings_file}')

        finish_time = datetime.datetime.now()
        self.info('Finished: ' + finish_time.strftime("%Y-%m-%d %H:%M:%S"))
        self.info(f'Elapsed time: {finish_time - self._start_time}')
        if self._num_warnings:
            self.info(f'WARNINGs: {self._num_warnings}')
        if self._num_errors:
            self.info(f'ERRORs: {self._num_errors}')

        self.info('\nThank you for using Nerpa!')

        for handler in list(self.handlers):  # iterate over a copy
            self.removeHandler(handler)
            try:
                handler.close()
            except Exception:
                pass
