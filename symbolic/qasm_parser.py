from __future__ import annotations

from pathlib import Path
import re
from typing import Iterable, List, Optional

import sympy as sp

from .branch_state import Gate
from .scalar import parse_scalar


_GATE_RE = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)\s*(?:\(([^)]*)\))?\s+(.+);$")
_QUBIT_RE = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)\[(\d+)\]$")
_QREG2_RE = re.compile(r"^qreg\s+([A-Za-z_][A-Za-z0-9_]*)\[(\d+)\];$")
_QREG3_RE = re.compile(r"^qubit\[(\d+)\]\s+([A-Za-z_][A-Za-z0-9_]*);$")


class QASMParseError(ValueError):
    pass


def parse_qasm_file(path: str | Path) -> List[Gate]:
    return parse_qasm_text(Path(path).read_text())


def parse_qasm_text(text: str) -> List[Gate]:
    gates: List[Gate] = []
    known_register: Optional[str] = None

    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = _strip_comment(raw_line).strip()
        if not line:
            continue
        if _is_ignored_statement(line):
            continue

        qreg_name = _parse_register_declaration(line)
        if qreg_name is not None:
            known_register = qreg_name
            continue

        match = _GATE_RE.match(line)
        if match is None:
            raise QASMParseError(f"Unsupported QASM statement on line {line_number}: {raw_line.strip()}")

        name, parameter_text, operands_text = match.groups()
        qubits = _parse_operands(operands_text, known_register, line_number)
        parameter = _parse_parameter(parameter_text) if parameter_text is not None else None
        gates.append(Gate(name=name.lower(), qubits=tuple(qubits), parameter=parameter, raw=line))

    return gates


def _strip_comment(line: str) -> str:
    return line.split("//", 1)[0]


def _is_ignored_statement(line: str) -> bool:
    return (
        line.startswith("OPENQASM")
        or line.startswith("include")
        or line.startswith("opaque ")
        or line.startswith("creg ")
        or line.startswith("barrier ")
    )


def _parse_register_declaration(line: str) -> Optional[str]:
    match2 = _QREG2_RE.match(line)
    if match2:
        return match2.group(1)

    match3 = _QREG3_RE.match(line)
    if match3:
        return match3.group(2)

    return None


def _parse_operands(text: str, known_register: Optional[str], line_number: int) -> Iterable[int]:
    operands = []
    for raw_operand in text.split(","):
        operand = raw_operand.strip()
        match = _QUBIT_RE.match(operand)
        if match is None:
            raise QASMParseError(f"Expected qubit operand on line {line_number}: {operand!r}")
        register, index_text = match.groups()
        if known_register is not None and register != known_register:
            raise QASMParseError(
                f"Only one quantum register is supported; found {register!r}, expected {known_register!r}"
            )
        operands.append(int(index_text))
    return operands


def _parse_parameter(text: str) -> sp.Expr:
    return parse_scalar(text)
