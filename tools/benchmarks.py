from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from symbolic.expr import OpExpr

DEFAULT_MIN_SIZE = 2
DEFAULT_MAX_SIZE = 128
DEFAULT_BLOCK_ENCODING_SIZES = tuple(range(DEFAULT_MIN_SIZE, DEFAULT_MAX_SIZE))


@dataclass(frozen=True)
class PauliTerm:
    coefficient: Fraction
    paulis: tuple[tuple[int, str], ...]

    @staticmethod
    def from_mapping(coefficient: int | float | Fraction, paulis: dict[int, str]) -> "PauliTerm":
        return PauliTerm(_fraction(coefficient), tuple(sorted(paulis.items())))

    @property
    def locality(self) -> int:
        return len(self.paulis)


@dataclass(frozen=True)
class HamiltonianBenchmark:
    benchmark_id: str
    model: str
    n_system: int
    terms: tuple[PauliTerm, ...]
    graph: str = ""

    @property
    def alpha(self) -> Fraction:
        return sum((abs(term.coefficient) for term in self.terms), Fraction(0, 1))

    @property
    def selector_ancillas(self) -> int:
        return math.ceil(math.log2(len(self.terms)))

    @property
    def n_ancilla(self) -> int:
        return self.selector_ancillas + 1

    @property
    def max_locality(self) -> int:
        return max((term.locality for term in self.terms), default=0)

    @property
    def ancillas(self) -> tuple[int, ...]:
        return tuple(range(self.n_ancilla))

    @property
    def systems(self) -> tuple[int, ...]:
        offset = self.n_ancilla
        return tuple(range(offset, offset + self.n_system))

    def expected_operator(self) -> OpExpr:
        if self.alpha == 0:
            raise ValueError("Hamiltonian normalization alpha must be nonzero")

        terms: dict[tuple[str, ...], sp.Expr] = {}
        for term in self.terms:
            key = ["I"] * self.n_system
            for index, op in term.paulis:
                key[index] = op
            pauli_string = tuple(key)
            terms[pauli_string] = terms.get(pauli_string, sp.Integer(0)) + _sympy_fraction(term.coefficient / self.alpha)
        return OpExpr(terms, num_qubits=self.n_system)

    def manifest_row(self, qasm_path: Path) -> dict[str, str]:
        return {
            "benchmark_id": self.benchmark_id,
            "rq": "RQ1",
            "stage": "block_encoding",
            "model": self.model,
            "graph": self.graph,
            "n_system": str(self.n_system),
            "n_terms": str(len(self.terms)),
            "alpha": _format_fraction(self.alpha),
            "selector_ancillas": str(self.selector_ancillas),
            "n_ancilla": str(self.n_ancilla),
            "max_locality": str(self.max_locality),
            "ancillas": " ".join(str(index) for index in self.ancillas),
            "systems": " ".join(str(index) for index in self.systems),
            "qasm_path": _format_path(qasm_path),
            "expected": str(self.expected_operator()),
        }


def make_hamiltonian(model: str, n: int, *, graph: str | None = None) -> HamiltonianBenchmark:
    model = model.lower()
    if model == "ising":
        return HamiltonianBenchmark(
            benchmark_id=f"ising_periodic_n{n}",
            model=model,
            graph="periodic",
            n_system=n,
            terms=tuple(ising_terms(n)),
        )
    if model == "maxcut":
        graph = graph or "cycle"
        if graph != "cycle":
            raise ValueError(f"Unsupported deterministic MaxCut graph family: {graph}")
        edges = cycle_edges(n)
        return HamiltonianBenchmark(
            benchmark_id=f"maxcut_cycle_n{n}",
            model=model,
            graph=graph,
            n_system=n,
            terms=tuple(maxcut_terms_from_edges(edges)),
        )
    if model == "heisenberg":
        return HamiltonianBenchmark(
            benchmark_id=f"heisenberg_n{n}",
            model=model,
            n_system=n,
            terms=tuple(heisenberg_terms(n)),
        )
    raise ValueError(f"Unsupported benchmark model: {model}")


def block_encoding_suite(sizes: Sequence[int] = DEFAULT_BLOCK_ENCODING_SIZES) -> list[HamiltonianBenchmark]:
    benchmarks = []
    for model in ("ising", "maxcut", "heisenberg"):
        for n in sizes:
            benchmarks.append(make_hamiltonian(model, n))
    return benchmarks


def normalize_terms(terms: Sequence[PauliTerm]) -> tuple[list[PauliTerm], Fraction]:
    alpha = sum((abs(term.coefficient) for term in terms), Fraction(0, 1))
    if alpha == 0:
        raise ValueError("Hamiltonian normalization alpha must be nonzero")
    return [
        PauliTerm(coefficient=term.coefficient / alpha, paulis=term.paulis)
        for term in terms
    ], alpha


def pauli_sum_to_expected_string(terms: Sequence[PauliTerm], *, n_system: int) -> str:
    normalized, _alpha = normalize_terms(terms)
    op_terms: dict[tuple[str, ...], sp.Expr] = {}
    for term in normalized:
        key = ["I"] * n_system
        for index, op in term.paulis:
            key[index] = op
        pauli_string = tuple(key)
        op_terms[pauli_string] = op_terms.get(pauli_string, sp.Integer(0)) + _sympy_fraction(term.coefficient)
    return str(OpExpr(op_terms, num_qubits=n_system))


def ising_terms(n: int, J: int | float | Fraction = 1, g: int | float | Fraction = 1) -> list[PauliTerm]:
    terms = []
    for i in range(n):
        terms.append(PauliTerm.from_mapping(J, {i: "Z", (i + 1) % n: "Z"}))
    for i in range(n):
        terms.append(PauliTerm.from_mapping(g, {i: "X"}))
    return terms


def cycle_edges(n: int) -> list[tuple[int, int]]:
    return [(i, (i + 1) % n) for i in range(n)]


def maxcut_terms_from_edges(edges: Iterable[tuple[int, int]]) -> list[PauliTerm]:
    terms = []
    for i, j in edges:
        terms.append(PauliTerm.from_mapping(Fraction(1, 2), {}))
        terms.append(PauliTerm.from_mapping(Fraction(-1, 2), {i: "Z", j: "Z"}))
    return terms


def heisenberg_terms(
    n: int,
    Jx: int | float | Fraction = 1,
    Jy: int | float | Fraction = 1,
    Jz: int | float | Fraction = 1,
) -> list[PauliTerm]:
    terms = []
    for i in range(n - 1):
        terms.append(PauliTerm.from_mapping(Jx, {i: "X", i + 1: "X"}))
        terms.append(PauliTerm.from_mapping(Jy, {i: "Y", i + 1: "Y"}))
        terms.append(PauliTerm.from_mapping(Jz, {i: "Z", i + 1: "Z"}))
    return terms


def write_block_encoding_qasm(benchmark: HamiltonianBenchmark, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = _block_encoding_qasm_lines(benchmark)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_manifest(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "benchmark_id",
        "rq",
        "stage",
        "model",
        "graph",
        "n_system",
        "n_terms",
        "alpha",
        "selector_ancillas",
        "n_ancilla",
        "max_locality",
        "ancillas",
        "systems",
        "qasm_path",
        "expected",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _block_encoding_qasm_lines(benchmark: HamiltonianBenchmark) -> list[str]:
    selector_width = benchmark.selector_ancillas
    work = selector_width
    system_offset = benchmark.n_ancilla
    num_qubits = benchmark.n_ancilla + benchmark.n_system

    prep = _uniform_first_l_state_prep(selector_width, len(benchmark.terms), work)
    lines = [
        "OPENQASM 2.0;",
        'include "qelib1.inc";',
        "",
        f"qreg q[{num_qubits}];",
        "",
        f"// RQ1 block-encoding benchmark: {benchmark.benchmark_id}",
        f"// q[0]..q[{selector_width - 1}] = selector ancillas",
        f"// q[{work}] = work ancilla",
        f"// q[{system_offset}]..q[{num_qubits - 1}] = system qubits",
        f"// top-left block = H / alpha, alpha = {_format_fraction(benchmark.alpha)}",
        "",
        "// PREP: uniform state over Hamiltonian terms",
    ]
    lines.extend(_format_gate(gate) for gate in prep)
    lines.extend(["", "// SELECT: branch-conditioned signed Pauli terms"])

    for index, term in enumerate(benchmark.terms):
        select = _select_term_gates(
            index=index,
            selector_width=selector_width,
            work=work,
            system_offset=system_offset,
            term=term,
        )
        if not select:
            lines.append(f"// term {index}: +I")
            continue
        lines.append(f"// term {index}: {_format_term_comment(term)}")
        lines.extend(_format_gate(gate) for gate in select)

    lines.extend(["", "// PREP dagger"])
    lines.extend(_format_gate(gate) for gate in _inverse_gates(prep))
    return lines


def _uniform_first_l_state_prep(selector_width: int, term_count: int, work: int) -> list[tuple[str, tuple[int, ...], float | None]]:
    if term_count < 1:
        raise ValueError("term_count must be positive")
    if term_count > 2**selector_width:
        raise ValueError("term_count does not fit in selector_width bits")

    gates: list[tuple[str, tuple[int, ...], float | None]] = []
    for depth in range(selector_width):
        for prefix_value in range(2**depth):
            prefix = _bits(prefix_value, depth)
            count = _prefix_count(prefix, selector_width, term_count)
            if count == 0:
                continue
            count0 = _prefix_count(prefix + (0,), selector_width, term_count)
            count1 = _prefix_count(prefix + (1,), selector_width, term_count)
            target = depth
            if count1 == 0:
                continue
            if count0 == 0:
                _append_controlled_x(gates, prefix, target, work)
                continue
            theta = 2.0 * math.asin(math.sqrt(count1 / count))
            _append_controlled_ry(gates, prefix, target, theta, work)
    return gates


def _select_term_gates(
    *,
    index: int,
    selector_width: int,
    work: int,
    system_offset: int,
    term: PauliTerm,
) -> list[tuple[str, tuple[int, ...], float | None]]:
    if term.coefficient > 0 and not term.paulis:
        return []

    gates: list[tuple[str, tuple[int, ...], float | None]] = []
    pattern = _bits(index, selector_width)
    _append_pattern_toggles(gates, pattern)
    _append_compute_flag(gates, selector_width, work)

    if term.coefficient < 0:
        gates.append(("z", (work,), None))

    for system_index, op in term.paulis:
        system = system_offset + system_index
        if op == "X":
            gates.append(("cx", (work, system), None))
        elif op == "Z":
            gates.append(("cz", (work, system), None))
        elif op == "Y":
            gates.append(("sdg", (work,), None))
            gates.append(("cx", (work, system), None))
            gates.append(("cz", (work, system), None))
        else:
            raise ValueError(f"Unsupported Pauli operator: {op}")

    _append_compute_flag(gates, selector_width, work)
    _append_pattern_toggles(gates, pattern)
    return gates


def _append_controlled_ry(
    gates: list[tuple[str, tuple[int, ...], float | None]],
    prefix: tuple[int, ...],
    target: int,
    theta: float,
    work: int,
) -> None:
    if not prefix:
        gates.append(("ry", (target,), theta))
        return

    _append_pattern_toggles(gates, prefix)
    if len(prefix) == 1:
        _append_cry(gates, 0, target, theta)
    else:
        controls = tuple(range(len(prefix)))
        gates.append(("mcx", controls + (work,), None))
        _append_cry(gates, work, target, theta)
        gates.append(("mcx", controls + (work,), None))
    _append_pattern_toggles(gates, prefix)


def _append_controlled_x(
    gates: list[tuple[str, tuple[int, ...], float | None]],
    prefix: tuple[int, ...],
    target: int,
    work: int,
) -> None:
    if not prefix:
        gates.append(("x", (target,), None))
        return

    _append_pattern_toggles(gates, prefix)
    if len(prefix) == 1:
        gates.append(("cx", (0, target), None))
    else:
        controls = tuple(range(len(prefix)))
        gates.append(("mcx", controls + (work,), None))
        gates.append(("cx", (work, target), None))
        gates.append(("mcx", controls + (work,), None))
    _append_pattern_toggles(gates, prefix)


def _append_cry(
    gates: list[tuple[str, tuple[int, ...], float | None]],
    control: int,
    target: int,
    theta: float,
) -> None:
    gates.append(("ry", (target,), theta / 2))
    gates.append(("cx", (control, target), None))
    gates.append(("ry", (target,), -theta / 2))
    gates.append(("cx", (control, target), None))


def _append_pattern_toggles(gates: list[tuple[str, tuple[int, ...], float | None]], pattern: tuple[int, ...]) -> None:
    for qubit, bit in enumerate(pattern):
        if bit == 0:
            gates.append(("x", (qubit,), None))


def _append_compute_flag(
    gates: list[tuple[str, tuple[int, ...], float | None]],
    selector_width: int,
    work: int,
) -> None:
    controls = tuple(range(selector_width))
    if selector_width == 1:
        gates.append(("cx", (0, work), None))
    else:
        gates.append(("mcx", controls + (work,), None))


def _inverse_gates(
    gates: Sequence[tuple[str, tuple[int, ...], float | None]]
) -> list[tuple[str, tuple[int, ...], float | None]]:
    inverse = []
    for name, qubits, parameter in reversed(gates):
        if name == "ry":
            inverse.append((name, qubits, -float(parameter)))
        else:
            inverse.append((name, qubits, parameter))
    return inverse


def _prefix_count(prefix: tuple[int, ...], width: int, term_count: int) -> int:
    return sum(1 for value in range(term_count) if _bits(value, width)[: len(prefix)] == prefix)


def _bits(value: int, width: int) -> tuple[int, ...]:
    return tuple((value >> shift) & 1 for shift in reversed(range(width)))


def _format_gate(gate: tuple[str, tuple[int, ...], float | None]) -> str:
    name, qubits, parameter = gate
    operands = ", ".join(f"q[{qubit}]" for qubit in qubits)
    if parameter is None:
        return f"{name} {operands};"
    return f"{name}({_format_angle(parameter)}) {operands};"


def _format_angle(value: float) -> str:
    if abs(value) < 1e-15:
        value = 0.0
    return f"{value:.17g}"


def _format_term_comment(term: PauliTerm) -> str:
    sign = "+" if term.coefficient > 0 else "-"
    if not term.paulis:
        return f"{sign}I"
    body = " ".join(f"{op}{index}" for index, op in term.paulis)
    return f"{sign}{body}"


def _fraction(value: int | float | Fraction) -> Fraction:
    if isinstance(value, Fraction):
        return value
    return Fraction(value)


def _sympy_fraction(value: Fraction) -> sp.Rational:
    return sp.Rational(value.numerator, value.denominator)


def _format_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _format_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def _parse_sizes(text: str) -> tuple[int, ...]:
    sizes: list[int] = []
    for raw_piece in text.split(","):
        piece = raw_piece.strip()
        if not piece:
            continue
        range_match = re.fullmatch(r"(\d+)\s*(?:\.\.|-)\s*(\d+)", piece)
        if range_match:
            start, stop = (int(value) for value in range_match.groups())
            if start > stop:
                raise ValueError(f"Invalid descending size range: {piece}")
            sizes.extend(range(start, stop + 1))
        else:
            sizes.append(int(piece))
    return tuple(dict.fromkeys(sizes))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate RQ1 block-encoding benchmark QASM files.")
    parser.add_argument(
        "--sizes",
        default=f"{DEFAULT_MIN_SIZE}..{DEFAULT_MAX_SIZE}",
        help='Comma-separated sizes or ranges, e.g. "2..128" or "2,3,5,8".',
    )
    parser.add_argument("--out-dir", type=Path, default=ROOT / "benchmarks" / "block_encoding")
    args = parser.parse_args(argv)

    out_dir = args.out_dir
    generated_dir = out_dir / "generated"
    manifest_path = out_dir / "manifest.csv"
    if generated_dir.exists():
        shutil.rmtree(generated_dir)

    rows = []
    for benchmark in block_encoding_suite(_parse_sizes(args.sizes)):
        qasm_path = generated_dir / benchmark.model / f"{benchmark.benchmark_id}.qasm"
        write_block_encoding_qasm(benchmark, qasm_path)
        rows.append(benchmark.manifest_row(qasm_path))
    write_manifest(rows, manifest_path)

    print(f"Wrote {len(rows)} block-encoding benchmarks to {_format_path(out_dir)}")
    print(f"Wrote manifest to {_format_path(manifest_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
