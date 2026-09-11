"""File logging for individual workflow steps.

Example::

    with Logger("results/progress.log") as log:
        log.info("Starting analysis")
        log.step("Sequence search", "Search completed")

Each instance owns its handlers and leaves the root logger unchanged.

For a workflow spanning modules, call ``RunLogger.set_logger(log)`` at startup and
``RunLogger.set_logger(None)`` before closing it. ``RunLogger.message`` and ``RunLogger.stage`` use
that registered logger only in its owning process. Without a registered logger,
decorated functions execute normally without writing stage messages. Worker
errors must be forwarded to the parent process, which logs the traceback once.
"""

import logging
import os
from functools import wraps
from os import PathLike
from pathlib import Path


class Logger:
    """Write timestamped messages to a custom file, optionally also to stderr.

    Parent directories are created automatically. Files use UTF-8 and are
    appended to unless ``mode="w"`` is explicitly selected. ``level`` accepts
    standard logging levels, such as ``"INFO"`` or ``logging.DEBUG``.
    Use a context manager or call :meth:`close` to release the file handle.
    Instances are intended for use within one process.
    """

    def __init__(
        self,
        path: str | PathLike[str],
        *,
        level: int | str = "INFO",
        console: bool = False,
        mode: str = "a",
        fmt: str = "%(asctime)s [%(levelname)s] %(message)s",
    ) -> None:
        if mode not in {"a", "w"}:
            raise ValueError("mode must be 'a' (append) or 'w' (overwrite)")
        self.path = Path(path).expanduser()
        self._logger = logging.Logger(__name__, level=level)
        self._logger.propagate = False
        self._closed = False
        formatter = logging.Formatter(fmt, datefmt="%Y-%m-%d %H:%M:%S")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        handler = logging.FileHandler(self.path, mode=mode, encoding="utf-8")
        handler.setFormatter(formatter)
        self._logger.addHandler(handler)
        if console:
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(formatter)
            self._logger.addHandler(console_handler)

    def log(
        self,
        text: str,
        *args: object,
        level: int | str = "INFO",
        step: str | None = None,
        exc_info: bool = False,
    ) -> None:
        """Log text with optional ``%`` formatting, step name, and traceback.

        Set ``exc_info=True`` inside an exception handler to include the
        current exception. Messages below the configured level are omitted.
        """
        if self._closed:
            raise ValueError("Cannot write to a closed logger")
        if isinstance(level, str):
            levels = logging.getLevelNamesMapping()
            if level not in levels:
                raise ValueError(f"Unknown log level: {level!r}")
            level = levels[level]
        # Format the text before adding the step so '%' in a step is literal.
        if step is not None:
            text = f"[{step}] {text % args if args else text}"
            args = ()
        self._logger.log(level, text, *args, exc_info=exc_info)

    def step(self, name: str, text: str, *args: object, level: int | str = "INFO") -> None:
        """Record a message associated with a named workflow step."""
        self.log(text, *args, level=level, step=name)

    def debug(self, text: str, *args: object) -> None:
        self.log(text, *args, level=logging.DEBUG)

    def info(self, text: str, *args: object) -> None:
        self.log(text, *args, level=logging.INFO)

    def warning(self, text: str, *args: object) -> None:
        self.log(text, *args, level=logging.WARNING)

    def error(self, text: str, *args: object) -> None:
        self.log(text, *args, level=logging.ERROR)

    def exception(self, text: str, *args: object) -> None:
        """Log an error with its traceback from inside an exception handler."""
        self.log(text, *args, level=logging.ERROR, exc_info=True)

    def close(self) -> None:
        """Flush and close this instance's handlers; safe to call repeatedly."""
        for handler in self._logger.handlers[:]:
            self._logger.removeHandler(handler)
            handler.close()
        self._closed = True

    def __enter__(self) -> "Logger":
        if self._closed:
            raise ValueError("Cannot reuse a closed logger")
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()


class RunLogger:
    """Share one run's Logger across modules in its owning process.

    Register at startup and detach before closing the Logger. Class methods
    allow decorators to be declared before a run starts; the active logger is
    looked up when the decorated function executes. Worker errors must still
    be forwarded to the parent for logging.
    """

    _logger: Logger | None = None
    _pid: int | None = None
    _failed_stage: str | None = None

    @classmethod
    def set_logger(cls, logger: Logger | None) -> None:
        """Register a logger, or detach with None; reset the failed stage."""
        cls._logger = logger
        cls._pid = os.getpid() if logger is not None else None
        cls._failed_stage = None

    @classmethod
    def _current(cls) -> Logger | None:
        return cls._logger if cls._pid == os.getpid() else None

    @classmethod
    def failed_stage(cls, default: str) -> str:
        """Return the innermost failed stage in this process, or a fallback."""
        return (cls._failed_stage or default) if cls._current() is not None else default

    @classmethod
    def stage(cls, name: str, require_result: bool = False):
        """Record stage boundaries without swallowing exceptions.

        With require_result=True, a None return indicates an early stop.
        The main error handler logs the traceback once.
        """
        def decorate(function):
            @wraps(function)
            def run(*args, **kwargs):
                logger = cls._current()
                if logger is not None:
                    logger.step(name, "Started")
                try:
                    result = function(*args, **kwargs)
                except BaseException:
                    if logger is not None and cls._failed_stage is None:
                        cls._failed_stage = name
                    raise
                if logger is not None:
                    status = "Stopped before completion" if require_result and result is None else "Completed"
                    logger.step(name, status)
                return result
            return run
        return decorate

    @classmethod
    def message(cls, message: str, level: int | str = "INFO") -> None:
        """Write a milestone without Rich markup in the owning process."""
        logger = cls._current()
        if logger is not None:
            from rich.text import Text
            logger.log(Text.from_markup(str(message)).plain.strip(), level=level)
