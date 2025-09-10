from src.pipeline.logging.logger import NerpaLogger, PreliminaryLogger
from src.pipeline.logging.dummy_logger import DummyLogger
from typing import Union

AnyLogger = Union[NerpaLogger, PreliminaryLogger, DummyLogger]