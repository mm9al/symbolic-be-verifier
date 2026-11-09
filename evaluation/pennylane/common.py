"""Portable Ising inputs and artifact helpers; no PennyLane/verifier imports."""
from __future__ import annotations

import hashlib
import json
import re
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SIZES = "2,3,4,5,6,7,8,9,10,16,17,32,33,64,65,128,256"
BENCHMARKS = ROOT / "benchmarks" / "pennylane"
RESULTS = Path(__file__).resolve().parent / "results"


def parse_sizes(text: str) -> list[int]:
    sizes = set()
    for piece in text.split(","):
        if ".." in piece:
            start, end = map(int, piece.split(".."))
            if end < start:
                raise ValueError("Size ranges must be increasing")
            sizes.update(range(start, end + 1))
        else:
            sizes.add(int(piece))
    if not sizes or min(sizes) < 2:
        raise ValueError("Periodic Ising sizes must be at least 2")
    return sorted(sizes)


def ising_terms(n: int) -> list[dict]:
    if n < 2:
        raise ValueError("Periodic Ising sizes must be at least 2")
    terms = [dict(coefficient="1", paulis={str(i): "Z", str((i + 1) % n): "Z"})
             for i in range(n)]
    terms += [dict(coefficient="1", paulis={str(i): "X"}) for i in range(n)]
    return terms


def target_expression(n: int) -> str:
    return " + ".join(" ".join(f"{op}_{index}" for index, op in term["paulis"].items())
                      for term in ising_terms(n))


def gate(name: str, wires, angle: str | None = None) -> str:
    parameter = "" if angle is None else f"({angle})"
    return f"{name}{parameter} " + ", ".join(f"q[{w}]" for w in wires) + ";"


def inverse(lines: list[str]) -> list[str]:
    result = []
    for line in reversed(lines):
        match = re.fullmatch(r"(rx|ry|rz)\(([^)]+)\)( .+)", line)
        if match:
            name, angle, operands = match.groups()
            opposite = angle[1:] if angle.startswith("-") else "-" + angle
            result.append(f"{name}({opposite}){operands}")
        else:
            name, operands = line.split(" ", 1)
            result.append({"s": "sdg", "sdg": "s", "t": "tdg", "tdg": "t"}.get(name, name)
                          + " " + operands)
    return result


def signatures(lines: list[str]) -> list[tuple]:
    """Compare exact decimal values, allowing harmless spellings such as -0."""
    result = []
    for line in lines:
        match = re.fullmatch(r"(\w+)(?:\(([^)]+)\))? (.+);", line)
        if match is None:
            raise ValueError(f"Invalid gate text: {line}")
        name, angle, operands = match.groups()
        result.append((name, None if angle is None else Fraction(angle), operands))
    return result


def qasm_body(path: Path) -> list[str]:
    return [line.split("//", 1)[0].strip() for line in path.read_text().splitlines()
            if line.split("//", 1)[0].strip()
            and not line.strip().startswith(("OPENQASM", "include", "qreg", "opaque"))]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
