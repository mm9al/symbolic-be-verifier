from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass
class ProfileCounters:
    simplify_sec: float = 0.0
    combine_sec: float = 0.0


_current: Optional[ProfileCounters] = None


def current() -> Optional[ProfileCounters]:
    return _current


def set_current(counters: ProfileCounters) -> None:
    global _current
    _current = counters


def clear_current() -> None:
    global _current
    _current = None


def add_simplify(seconds: float) -> None:
    if _current is not None:
        _current.simplify_sec += seconds


def add_combine(seconds: float) -> None:
    if _current is not None:
        _current.combine_sec += seconds
