import logging
import sys
from pathlib import Path


def configure_logging(log_file: Path | None = None) -> logging.Logger:


    logger = logging.getLogger('train_nerpa')
    logger.setLevel(logging.DEBUG)  # overall logger level: DEBUG

    # clear any existing handlers to avoid duplicate messages
    if logger.handlers:
        logger.handlers.clear()

    formatter = logging.Formatter(
        "%(levelname)s\n%(filename)s:%(lineno)d\n%(funcName)s\n\n%(message)s\n"
        "----------------------------------------------------------\n"
    )

    # stdout handler set to INFO
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    # file handler set to DEBUG
    # ensure output dir exists
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        if log_file.exists():
            log_file.unlink()

        file_handler = logging.FileHandler(str(log_file))
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        logger.addHandler(file_handler)

    return logger
