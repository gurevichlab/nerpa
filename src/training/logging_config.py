import logging
import sys
from pathlib import Path


def configure_logging(output_dir: Path) -> logging.Logger:
    # ensure output dir exists
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger('train_nerpa')
    logger.setLevel(logging.DEBUG)  # overall logger level: DEBUG

    # clear any existing handlers to avoid duplicate messages
    if logger.handlers:
        logger.handlers.clear()

    formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s:%(lineno)d: %(message)s')

    # stdout handler set to INFO
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)

    # file handler set to DEBUG
    file_handler = logging.FileHandler(str(output_dir / 'train_nerpa.log'))
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    return logger
