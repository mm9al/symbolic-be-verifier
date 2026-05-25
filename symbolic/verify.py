from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence

from .branch_state import BranchState, Gate
from .expr import OpExpr, parse_operator_expression
from .qasm_parser import parse_qasm_file


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

    @property
    def success(self) -> Optional[bool]:
        if self.expected is None:
            return None
        return self.final_state.b0.equals(self.expected)


def run_circuit(gates: Iterable[Gate], *, ancilla: int = 0, system: int = 1, keep_trace: bool = False) -> VerificationResult:
    state = BranchState.initial()
    trace: List[TraceStep] = [TraceStep(0, None, state)] if keep_trace else []

    for index, gate in enumerate(gates, start=1):
        state = state.apply(gate, ancilla=ancilla, system=system)
        if keep_trace:
            trace.append(TraceStep(index, gate, state))

    return VerificationResult(final_state=state, trace=trace)


def verify_qasm_file(
    path: str | Path,
    *,
    expected: str | OpExpr | None = None,
    ancilla: int = 0,
    system: int = 1,
    keep_trace: bool = False,
) -> VerificationResult:
    gates = parse_qasm_file(path)
    expected_expr = _normalize_expected(expected)
    result = run_circuit(gates, ancilla=ancilla, system=system, keep_trace=keep_trace)
    return VerificationResult(final_state=result.final_state, trace=result.trace, expected=expected_expr)


def format_result(result: VerificationResult, *, show_trace: bool = False) -> str:
    lines: List[str] = []

    if show_trace:
        trace = result.trace or [TraceStep(0, None, result.final_state)]
        for step in trace:
            label = "Initial" if step.gate is None else str(step.gate)
            b0 = str(step.state.b0)
            lines.append(f"{step.index:>3}  {label:<20}  B0 = {b0:<18}  B1 = {step.state.b1}")

    lines.append(f"Final B0 = {result.final_state.b0}")
    lines.append(f"Final B1 = {result.final_state.b1}")

    if result.expected is not None:
        lines.append(f"Expected B0 = {result.expected}")
        if result.success:
            lines.append("PASS")
        else:
            diff = result.final_state.b0 - result.expected
            lines.append("FAIL")
            lines.append(f"B0 - expected = {diff}")

    return "\n".join(lines)


def _normalize_expected(expected: str | OpExpr | None) -> Optional[OpExpr]:
    if expected is None:
        return None
    if isinstance(expected, OpExpr):
        return expected
    return parse_operator_expression(expected)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Symbolically verify a one-ancilla one-system block encoding.")
    parser.add_argument("qasm_path", type=Path)
    parser.add_argument("--expected", help='Expected top-left block, for example "(X + Z)/2".')
    parser.add_argument("--ancilla", type=int, default=0)
    parser.add_argument("--system", type=int, default=1)
    parser.add_argument("--trace", action="store_true", help="Print every intermediate B0/B1 update.")
    args = parser.parse_args(argv)

    result = verify_qasm_file(
        args.qasm_path,
        expected=args.expected,
        ancilla=args.ancilla,
        system=args.system,
        keep_trace=args.trace,
    )
    print(format_result(result, show_trace=args.trace))
    return 0 if result.success is not False else 1


if __name__ == "__main__":
    raise SystemExit(main())
