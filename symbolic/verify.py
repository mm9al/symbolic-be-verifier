from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable, List, Optional, Sequence, Tuple

import sympy as sp

from .branch_state import BranchState, Gate
from .expr import OpExpr, parse_operator_expression
from .qasm_parser import parse_qasm_file
from .scalar import scalar_simplify


PASS = "PASS"
PASS_UP_TO_GLOBAL_PHASE = "PASS_UP_TO_GLOBAL_PHASE"
FAIL = "FAIL"


@dataclass(frozen=True)
class TraceStep:
    index: int
    gate: Optional[Gate]
    state: BranchState


@dataclass(frozen=True)
class VerificationResult:
    final_state: BranchState
    trace: Sequence[TraceStep]
    expected: Optional[OpExpr] = None
    ancilla_qubits: Tuple[int, ...] = (0,)
    system_qubits: Tuple[int, ...] = (1,)

    @property
    def success(self) -> Optional[bool]:
        if self.expected is None:
            return None
        return self.status in {PASS, PASS_UP_TO_GLOBAL_PHASE}

    @property
    def status(self) -> Optional[str]:
        if self.expected is None:
            return None
        if self.final_state.top_left().equals(self.expected):
            return PASS
        is_proportional, phase = proportional_phase(self.final_state.top_left(), self.expected)
        if is_proportional:
            return PASS_UP_TO_GLOBAL_PHASE
        return FAIL

    @property
    def global_phase(self) -> Optional[sp.Expr]:
        if self.expected is None or self.final_state.top_left().equals(self.expected):
            return None
        is_proportional, phase = proportional_phase(self.final_state.top_left(), self.expected)
        if is_proportional:
            return phase
        return None


def run_circuit(
    gates: Iterable[Gate],
    *,
    ancilla: int = 0,
    ancillas: Optional[Sequence[int]] = None,
    system: int = 1,
    systems: Optional[Sequence[int]] = None,
    keep_trace: bool = False,
) -> VerificationResult:
    ancilla_qubits = _normalize_ancillas(ancilla=ancilla, ancillas=ancillas)
    system_qubits = _normalize_systems(system=system, systems=systems)
    state = BranchState.initial(num_system_qubits=len(system_qubits), num_ancillas=len(ancilla_qubits))
    trace: List[TraceStep] = [TraceStep(0, None, state)] if keep_trace else []

    for index, gate in enumerate(gates, start=1):
        state = state.apply(gate, ancillas=ancilla_qubits, systems=system_qubits)
        if keep_trace:
            trace.append(TraceStep(index, gate, state))

    return VerificationResult(final_state=state, trace=trace, ancilla_qubits=ancilla_qubits, system_qubits=system_qubits)


def verify_qasm_file(
    path: str | Path,
    *,
    expected: str | OpExpr | None = None,
    ancilla: int = 0,
    ancillas: Optional[Sequence[int]] = None,
    system: int = 1,
    systems: Optional[Sequence[int]] = None,
    keep_trace: bool = False,
) -> VerificationResult:
    gates = parse_qasm_file(path)
    ancilla_qubits = _normalize_ancillas(ancilla=ancilla, ancillas=ancillas)
    system_qubits = _normalize_systems(system=system, systems=systems)
    expected_expr = _normalize_expected(expected, num_system_qubits=len(system_qubits))
    result = run_circuit(gates, ancillas=ancilla_qubits, systems=system_qubits, keep_trace=keep_trace)
    return VerificationResult(
        final_state=result.final_state,
        trace=result.trace,
        expected=expected_expr,
        ancilla_qubits=ancilla_qubits,
        system_qubits=system_qubits,
    )


def format_result(result: VerificationResult, *, show_trace: bool = False) -> str:
    lines: List[str] = []

    lines.extend(_format_mapping(result))

    if show_trace:
        trace = result.trace or [TraceStep(0, None, result.final_state)]
        for step in trace:
            lines.extend(_format_trace_step(step))

    if result.final_state.num_ancillas == 1:
        lines.append(f"Final B0 = {result.final_state.b0}")
        lines.append(f"Final B1 = {result.final_state.b1}")
    else:
        key = (0,) * result.final_state.num_ancillas
        lines.append(f"Final B[{_format_key(key)}] = {result.final_state.top_left()}")

    if result.expected is not None:
        expected_label = "Expected B0" if result.final_state.num_ancillas == 1 else f"Expected B[{_format_key((0,) * result.final_state.num_ancillas)}]"
        lines.append(f"{expected_label} = {result.expected}")
        if result.status == PASS:
            lines.append(PASS)
        elif result.status == PASS_UP_TO_GLOBAL_PHASE:
            lines.append(PASS_UP_TO_GLOBAL_PHASE)
            lines.append(f"phase = {_format_scalar(result.global_phase)}")
        else:
            diff = result.final_state.top_left() - result.expected
            lines.append(FAIL)
            diff_label = "B0" if result.final_state.num_ancillas == 1 else f"B[{_format_key((0,) * result.final_state.num_ancillas)}]"
            lines.append(f"{diff_label} - expected = {diff}")

    return "\n".join(lines)


def _normalize_expected(expected: str | OpExpr | None, *, num_system_qubits: int) -> Optional[OpExpr]:
    if expected is None:
        return None
    if isinstance(expected, OpExpr):
        return expected.with_num_qubits(num_system_qubits)
    return parse_operator_expression(expected, num_qubits=num_system_qubits)


def _format_mapping(result: VerificationResult) -> List[str]:
    lines = ["Ancilla qubits:"]
    lines.extend(f"  q[{qasm_index}] -> ancilla[{index}]" for index, qasm_index in enumerate(result.ancilla_qubits))
    lines.append("")
    lines.append("System qubits:")
    lines.extend(f"  q[{qasm_index}] -> system[{index}]" for index, qasm_index in enumerate(result.system_qubits))
    lines.append("")
    return lines


def _format_trace_step(step: TraceStep) -> List[str]:
    label = "Initial" if step.gate is None else str(step.gate)
    lines = [f"{step.index} {label}"]
    for key, branch in sorted(step.state.branches.items()):
        lines.append(f"  B[{_format_key(key)}] = {branch}")
    if not step.state.branches:
        lines.append("  <all zero branches>")
    return lines


def _format_key(key: Tuple[int, ...]) -> str:
    return "".join(str(bit) for bit in key)


def proportional_phase(actual: OpExpr, expected: OpExpr) -> Tuple[bool, Optional[sp.Expr]]:
    actual, expected = _align_operator_qubits(actual, expected)
    actual_terms = _nonzero_terms(actual)
    expected_terms = _nonzero_terms(expected)

    if set(actual_terms) != set(expected_terms):
        return False, None
    if not expected_terms:
        return False, None

    p0 = next(iter(expected_terms))
    phase = scalar_simplify(actual_terms[p0] / expected_terms[p0])

    for pauli_string, expected_coeff in expected_terms.items():
        if scalar_simplify(actual_terms[pauli_string] - phase * expected_coeff) != 0:
            return False, None

    if scalar_simplify(phase * sp.conjugate(phase) - 1) != 0:
        return False, phase

    return True, phase


def _align_operator_qubits(actual: OpExpr, expected: OpExpr) -> Tuple[OpExpr, OpExpr]:
    num_qubits = max(actual.num_qubits, expected.num_qubits)
    return actual.with_num_qubits(num_qubits), expected.with_num_qubits(num_qubits)


def _nonzero_terms(expr: OpExpr) -> dict:
    return {pauli_string: coeff for pauli_string, coeff in expr.terms.items() if scalar_simplify(coeff) != 0}


def _format_scalar(value: Optional[sp.Expr]) -> str:
    if value is None:
        return "None"
    return sp.sstr(scalar_simplify(value)).replace("I", "i")


def _normalize_systems(*, system: int, systems: Optional[Sequence[int]]) -> Tuple[int, ...]:
    if systems is None:
        return (system,)
    normalized = tuple(systems)
    if not normalized:
        raise ValueError("At least one system qubit is required")
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"Duplicate system qubits are not allowed: {normalized}")
    return normalized


def _normalize_ancillas(*, ancilla: int, ancillas: Optional[Sequence[int]]) -> Tuple[int, ...]:
    if ancillas is None:
        return (ancilla,)
    normalized = tuple(ancillas)
    if not normalized:
        raise ValueError("At least one ancilla qubit is required")
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"Duplicate ancilla qubits are not allowed: {normalized}")
    return normalized


def _parse_qubit_arg(value: str | int) -> int:
    if isinstance(value, int):
        return value
    match = re.fullmatch(r"q\[(\d+)\]", value.strip())
    if match:
        return int(match.group(1))
    return int(value)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Symbolically verify a block encoding with one or more ancillas.")
    parser.add_argument("qasm_path", type=Path)
    parser.add_argument("--expected", help='Expected top-left block, for example "(X0 X1 + Z0 Z1)/2".')
    parser.add_argument("--ancilla", default="0", help="Single ancilla qubit index. Ignored when --ancillas is set.")
    parser.add_argument("--ancillas", nargs="+", help='Ancilla qubits in branch-key order, e.g. q[0] q[1].')
    parser.add_argument("--system", default="1", help="Single system qubit index. Ignored when --systems is set.")
    parser.add_argument("--systems", nargs="+", help='System qubits in operator-index order, e.g. q[2] q[3].')
    parser.add_argument("--trace", action="store_true", help="Print every intermediate B0/B1 update.")
    args = parser.parse_args(argv)

    ancillas = tuple(_parse_qubit_arg(value) for value in args.ancillas) if args.ancillas is not None else None
    systems = tuple(_parse_qubit_arg(value) for value in args.systems) if args.systems is not None else None

    result = verify_qasm_file(
        args.qasm_path,
        expected=args.expected,
        ancilla=_parse_qubit_arg(args.ancilla),
        ancillas=ancillas,
        system=_parse_qubit_arg(args.system),
        systems=systems,
        keep_trace=args.trace,
    )
    print(format_result(result, show_trace=args.trace))
    return 0 if result.success is not False else 1


if __name__ == "__main__":
    raise SystemExit(main())
