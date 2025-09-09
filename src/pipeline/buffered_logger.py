# buffered_logger.py  (Python 3.11+)
from __future__ import annotations

import inspect
import logging
import time
from typing import Any, Dict, Iterable, NamedTuple


class LogCall(NamedTuple):
    """One recorded log call."""
    timestamp: float
    method: str                 # 'debug' | 'info' | 'warning' | 'error' | 'critical' | 'exception' | 'log'
    level: int | None           # only set for method == 'log'
    message: str
    kwargs: Dict[str, Any]
    caller_filename: str
    caller_lineno: int
    caller_func: str


_LEVEL_BY_NAME: dict[str, int] = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warning": logging.WARNING,
    "error": logging.ERROR,
    "critical": logging.CRITICAL,
    "exception": logging.ERROR,  # .exception logs at ERROR + exc_info=True
}


class _CallsiteRewriteFilter(logging.Filter):
    """
    If a LogRecord has caller_* fields, rewrite standard fields so a normal
    formatter using %(filename)s, %(lineno)d, %(funcName)s shows the *original* call site.
    """
    def filter(self, record: logging.LogRecord) -> bool:
        d = record.__dict__
        fn = d.get("caller_filename")
        ln = d.get("caller_lineno")
        fu = d.get("caller_func")
        if fn is not None:
            record.pathname = fn
            record.filename = fn.rsplit("/", 1)[-1].rsplit("\\", 1)[-1]
        if ln is not None:
            record.lineno = ln
        if fu is not None:
            record.funcName = fu
        return True


class TemporaryCallsiteRewrite:
    """Context manager to install/remove _CallsiteRewriteFilter during replay."""
    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self._filter = _CallsiteRewriteFilter()

    def __enter__(self):
        for h in self.logger.handlers:
            h.addFilter(self._filter)
        return self

    def __exit__(self, exc_type, exc, tb):
        for h in self.logger.handlers:
            h.removeFilter(self._filter)
        return False


class _NullContext:
    def __enter__(self): return self
    def __exit__(self, *exc): return False


class BufferedLogger:
    """A logger that records log calls (no I/O) for later replay."""
    level: int
    calls: list[LogCall]

    def __init__(self, level: int = logging.INFO):
        self.level = level
        self.calls = []

    def setLevel(self, level: int) -> None:
        self.level = level

    def getEffectiveLevel(self) -> int:
        return self.level

    def isEnabledFor(self, level: int) -> bool:
        return level >= self.level

    def log(self, level: int, msg: str, **kwargs: Any) -> None:
        """Generic log with explicit level."""
        if self.isEnabledFor(level):
            self._record("log", level, msg, kwargs, skip=2)

    def __getattr__(self, name: str):
        """Provide .debug/.info/.warning/.error/.critical/.exception dynamically."""
        if name not in _LEVEL_BY_NAME:
            raise AttributeError(name)
        lvl = _LEVEL_BY_NAME[name]

        def wrapper(msg: str, **kwargs: Any):
            if name == "exception":
                kwargs.setdefault("exc_info", True)
            if self.isEnabledFor(lvl):
                self._record(name, None, msg, kwargs, skip=2)

        return wrapper


    def _record(self, method: str, level: int | None, msg: str,
                kwargs: Dict[str, Any], *, skip: int) -> None:
        frame = inspect.currentframe()
        for _ in range(skip):
            if frame is not None:
                frame = frame.f_back

        filename = frame.f_code.co_filename if frame else "<unknown>"
        lineno = frame.f_lineno if frame else 0
        funcname = frame.f_code.co_name if frame else "<unknown>"

        self.calls.append(
            LogCall(
                timestamp=time.time(),
                method=method,
                level=level,
                message=msg,
                kwargs=dict(kwargs) if kwargs else {},
                caller_filename=filename,
                caller_lineno=lineno,
                caller_func=funcname,
            )
        )

    @classmethod
    def replay(
        cls,
        real_logger: logging.Logger,
        buffers: Iterable[BufferedLogger],
        *,
        sort_by_time: bool = True,
        rewrite_callsite: bool = True,
    ) -> None:
        """Flush all recorded calls to `real_logger`."""
        calls = [c for b in buffers for c in b.calls]
        if sort_by_time:
            calls.sort(key=lambda c: c.timestamp)

        ctx = TemporaryCallsiteRewrite(real_logger) if rewrite_callsite else _NullContext()
        with ctx:
            for c in calls:
                kwargs = dict(c.kwargs) if c.kwargs else {}
                extra = kwargs.pop("extra", {}) or {}
                # Provide callsite for the filter to rewrite std fields
                extra.setdefault("caller_filename", c.caller_filename)
                extra.setdefault("caller_lineno", c.caller_lineno)
                extra.setdefault("caller_func", c.caller_func)
                kwargs["extra"] = extra

                if c.method == "log":
                    level = c.level if c.level is not None else logging.INFO
                    real_logger.log(level, c.message, **kwargs)
                else:
                    getattr(real_logger, c.method)(c.message, **kwargs)


# --- Example (remove or keep as a quick test) ---
if __name__ == "__main__":
    def worker_task(i: int) -> tuple[int, BufferedLogger]:
        blog = BufferedLogger()
        blog.info(f"Processing {i}")
        if i == 42:
            blog.warning(f"Found the answer: {i}")
        blog.log(logging.DEBUG, f"Debug detail for {i}")
        return i * 2, blog

    results_and_logs = [worker_task(x) for x in (1, 42, 3)]
    results = [r for r, _ in results_and_logs]
    buffers = [b for _, b in results_and_logs]

    logger = logging.getLogger("main")
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(filename)s:%(lineno)d in %(funcName)s - %(message)s"
    ))
    logger.handlers[:] = [handler]

    BufferedLogger.replay(logger, buffers, sort_by_time=True, rewrite_callsite=True)
    print(f"Results: {results}")
