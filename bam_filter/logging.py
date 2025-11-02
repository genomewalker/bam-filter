"""Lightweight logging facade with level-aware routing to the nogil logger.

The module keeps Python and C loggers aligned with the same verbosity value,
adds per-tag defaults so the CLI stays concise by default, and still exposes
small helper functions (`info`, `stats`, `graph`, `log`, `verbose`, …) used
throughout the project.
"""

from __future__ import annotations

import sys
from contextlib import contextmanager
from enum import IntEnum
from time import perf_counter
from typing import Any, Dict, Optional, Union

try:
    from bam_filter import _c_logging as _c_logging_ext  # type: ignore

    _HAS_CLOG = True
except Exception:
    _c_logging_ext = None  # type: ignore
    _HAS_CLOG = False

LogLevelInput = Union[int, str, "LogLevel"]


class LogLevel(IntEnum):
    QUIET = -1
    SUMMARY = 0
    INFO = 1
    DEBUG = 2
    TRACE = 3


_INT_TO_LEVEL = {level.value: level for level in LogLevel}
_NAME_TO_LEVEL: Dict[str, LogLevel] = {
    "QUIET": LogLevel.QUIET,
    "SILENT": LogLevel.QUIET,
    "SUMMARY": LogLevel.SUMMARY,
    "INFO": LogLevel.INFO,
    "INFORMATION": LogLevel.INFO,
    "DEBUG": LogLevel.DEBUG,
    "TRACE": LogLevel.TRACE,
    "VERBOSE": LogLevel.TRACE,
}

# Default level for frequently used tags (case-insensitive).
_TAG_LEVEL_DEFAULTS: Dict[str, LogLevel] = {
    "": LogLevel.SUMMARY,
    "CLI": LogLevel.SUMMARY,
    "REASSIGN": LogLevel.SUMMARY,
    "BUILD-TAXONOMY": LogLevel.SUMMARY,
    "PROCESS": LogLevel.INFO,
    "GRAPH": LogLevel.INFO,
    "BAM-WRITER": LogLevel.SUMMARY,
    "NETWORK-AWARE ANALYSIS": LogLevel.INFO,
    "BATCH": LogLevel.DEBUG,
    "BROKEN-STICK": LogLevel.DEBUG,
    "IGRAPH OPS": LogLevel.DEBUG,
    "LEIDEN": LogLevel.DEBUG,
    "COMMUNITY": LogLevel.TRACE,
    "ADAPTIVE": LogLevel.TRACE,
    "PHASE 1": LogLevel.DEBUG,
    "PHASE 2": LogLevel.DEBUG,
    "PHASE 3": LogLevel.DEBUG,
    "PHASE 4": LogLevel.DEBUG,
    "PHASE 5": LogLevel.DEBUG,
    "PHASE 6": LogLevel.DEBUG,
    "EMERGENCY REGULARIZATION": LogLevel.INFO,
    "SQUAREM": LogLevel.INFO,
}


def _coerce_level(level: LogLevelInput) -> LogLevel:
    if isinstance(level, LogLevel):
        return level
    if isinstance(level, str):
        key = level.strip().replace("-", "_").replace(" ", "_").upper()
        if key in _NAME_TO_LEVEL:
            return _NAME_TO_LEVEL[key]
        raise ValueError(f"Unknown log level name: {level!r}")
    try:
        value = int(level)
    except Exception as exc:  # pragma: no cover - defensive path
        raise TypeError(f"Unsupported level type: {level!r}") from exc
    if value in _INT_TO_LEVEL:
        return _INT_TO_LEVEL[value]
    if value < LogLevel.QUIET:
        return LogLevel.QUIET
    if value > LogLevel.TRACE:
        return LogLevel.TRACE
    return LogLevel(value)


def _normalize_tag(tag: Optional[str]) -> str:
    if not tag:
        return ""
    return str(tag).strip()


def _normalize_tag_key(tag: Optional[str]) -> str:
    return _normalize_tag(tag).upper()


class _LoggerState:
    __slots__ = ("level", "tag_levels")

    def __init__(self) -> None:
        self.level: LogLevel = LogLevel.SUMMARY
        self.tag_levels: Dict[str, LogLevel] = dict(_TAG_LEVEL_DEFAULTS)
        if _HAS_CLOG:
            try:
                _c_logging_ext.set_verbosity(int(self.level))
            except Exception:
                pass


_state = _LoggerState()


def _format(fmt: str, *args: Any) -> str:
    return fmt % args if args else fmt


def _resolve_level(tag: Optional[str], level: Optional[LogLevelInput]) -> LogLevel:
    if level is not None:
        return _coerce_level(level)
    key = _normalize_tag_key(tag)
    if key in _state.tag_levels:
        return _state.tag_levels[key]
    if key:
        return LogLevel.INFO
    return LogLevel.SUMMARY


def _should_emit(level: LogLevel) -> bool:
    try:
        if _HAS_CLOG:
            return bool(_c_logging_ext.should_log(int(level)))
    except Exception:
        pass
    return level <= _state.level


def _format_line(tag: str, message: str, level: LogLevel) -> str:
    text = message.rstrip("\n")
    show_tag = bool(tag) and _state.level >= LogLevel.INFO
    if not show_tag:
        return f"{text}\n"
    return f"{tag:<12} | {text}\n"


def _emit_to_stderr(tag: str, message: str, level: LogLevel) -> None:
    sys.stderr.write(_format_line(tag, message, level))


def _emit(tag: Optional[str], message: str, level: LogLevel) -> None:
    output_tag = _normalize_tag(tag)
    if _HAS_CLOG:
        try:
            _c_logging_ext.nogil_log(output_tag, message.rstrip("\n"), int(level))
            return
        except Exception:
            pass
    _emit_to_stderr(output_tag, message, level)


def set_tag_level(tag: str, level: LogLevelInput) -> None:
    """Override the default level for a specific tag (case-insensitive)."""
    _state.tag_levels[_normalize_tag_key(tag)] = _coerce_level(level)


def reset_tag_level(tag: str) -> None:
    """Remove a tag override so the default tier applies again."""
    _state.tag_levels.pop(_normalize_tag_key(tag), None)


def set_level(level: LogLevelInput) -> LogLevel:
    """Set the global logging level."""
    coerced = _coerce_level(level)
    _state.level = coerced
    if _HAS_CLOG:
        try:
            _c_logging_ext.set_verbosity(int(coerced))
        except Exception:
            pass
    return coerced


def get_level() -> LogLevel:
    if _HAS_CLOG:
        try:
            return _coerce_level(_c_logging_ext.get_verbosity())
        except Exception:
            pass
    return _state.level


def set_verbosity(level: LogLevelInput) -> None:
    set_level(level)


def get_verbosity() -> int:
    return int(get_level())


def log(tag: str, fmt: str, *args: Any, level: Optional[LogLevelInput] = None) -> None:
    resolved_level = _resolve_level(tag, level)
    if not _should_emit(resolved_level):
        return
    _emit(tag, _format(fmt, *args), resolved_level)


def info(fmt: str, *args: Any) -> None:
    log("INFO", fmt, *args, level=LogLevel.INFO)


def debug(tag: str, fmt: str, *args: Any) -> None:
    log(tag, fmt, *args, level=LogLevel.DEBUG)


def trace(tag: str, fmt: str, *args: Any) -> None:
    log(tag, fmt, *args, level=LogLevel.TRACE)


def summary(fmt: str, *args: Any) -> None:
    log("", fmt, *args, level=LogLevel.SUMMARY)


def stats(fmt: str, *args: Any) -> None:
    log("", fmt, *args, level=LogLevel.SUMMARY)


def graph(fmt: str, *args: Any) -> None:
    log("GRAPH", fmt, *args, level=LogLevel.INFO)


def warn(fmt: str, *args: Any) -> None:
    log("WARN", fmt, *args, level=LogLevel.SUMMARY)


def error(fmt: str, *args: Any) -> None:
    log("ERROR", fmt, *args, level=LogLevel.QUIET)


def verbose(level: int, tag: str, fmt: str, *args: Any) -> None:
    desired = _coerce_level(level)
    if not _should_emit(desired):
        return
    resolved = max(desired, _resolve_level(tag, None))
    _emit(tag, _format(fmt, *args), resolved)


def should_log(level: int) -> bool:
    return _should_emit(_coerce_level(level))


@contextmanager
def time_block(tag: str, label: str, level: int = 0):
    start = perf_counter()
    try:
        yield
    finally:
        duration = perf_counter() - start
        verbose(level, tag, "%s completed in %.3fs", label, duration)


def start_timer() -> float:
    return perf_counter()


def log_duration(tag: str, label: str, start: float, level: int = 0) -> float:
    duration = perf_counter() - start
    verbose(level, tag, "%s completed in %.3fs", label, duration)
    return duration


__all__ = [
    "LogLevel",
    "debug",
    "error",
    "get_level",
    "get_verbosity",
    "graph",
    "info",
    "log",
    "log_duration",
    "reset_tag_level",
    "set_level",
    "set_tag_level",
    "set_verbosity",
    "should_log",
    "start_timer",
    "stats",
    "summary",
    "time_block",
    "trace",
    "verbose",
    "warn",
]
