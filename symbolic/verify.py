from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import io
from pathlib import Path
import re
import time
from typing import Iterable, List, Optional, Sequence, Tuple

import sympy as sp

from .branch_state import BranchState, Gate
from .expr import OpExpr, parse_operator_expression
from . import profile
from .qasm_parser import parse_qasm_file
from .scalar import scalar_simplify
from .word_expr import GARBAGE_ATOMS, WordExpr, eliminate


PASS = "PASS"
PASS_UP_TO_SCALE = "PASS_UP_TO_SCALE"
FAIL = "FAIL"
FAIL_GARBAGE = "FAIL_GARBAGE"
DEFAULT_TOLERANCE = 1e-8


@dataclass(frozen=True)
class TraceStep:
    index: int
    gate: Optional[Gate]
    state: BranchState


@dataclass(frozen=True)
class GateProfileStep:
    gate_id: int
    gate_name: str
    gate: str
    num_nonzero_branches: int
    total_operator_terms: int
    max_terms_per_branch: int
    time_this_gate: float
    time_simplify: float
    time_combine: float


@dataclass(frozen=True)
class VerificationResult:
    final_state: BranchState
    trace: Sequence[TraceStep]
    gate_profiles: Sequence[GateProfileStep] = ()
    expected: Optional[OpExpr] = None
    ancilla_qubits: Tuple[int, ...] = (0,)
    system_qubits: Tuple[int, ...] = (1,)
    qsp_normalized: Optional[WordExpr] = None
    qsp_polynomial: Optional[sp.Expr] = None
    qsp_expected_polynomial: Optional[sp.Expr] = None
    qsp_actual: Optional[OpExpr] = None
    qsp_expected: Optional[OpExpr] = None
    qsp_status: Optional[str] = None
    qsp_polynomial_only: bool = False
    tolerance: float = DEFAULT_TOLERANCE

    @property
    def success(self) -> Optional[bool]:
        if self.qsp_status is not None:
            return self.qsp_status == PASS
        if self.expected is None:
            return None
        return self.status in {PASS, PASS_UP_TO_SCALE}

    @property
    def status(self) -> Optional[str]:
        if self.qsp_status is not None:
            return self.qsp_status
        if self.expected is None:
            return None
        top_left = self.final_state.top_left()
        if not isinstance(top_left, OpExpr):
            return FAIL
        if pauli_expr_close(top_left, self.expected, tol=self.tolerance):
            return PASS
        is_proportional, _scale = proportional_scale(top_left, self.expected, tol=self.tolerance)
        if is_proportional:
            return PASS_UP_TO_SCALE
        return FAIL

    @property
    def scale(self) -> Optional[sp.Expr]:
        top_left = self.final_state.top_left()
        if self.expected is None or not isinstance(top_left, OpExpr) or pauli_expr_close(top_left, self.expected, tol=self.tolerance):
            return None
        is_proportional, scale = proportional_scale(top_left, self.expected, tol=self.tolerance)
        if is_proportional:
            return scale
        return None


def run_circuit(
    gates: Iterable[Gate],
    *,
    ancilla: int = 0,
    ancillas: Optional[Sequence[int]] = None,
    system: int = 1,
    systems: Optional[Sequence[int]] = None,
    keep_trace: bool = False,
    profile_gates: bool = False,
    expression_kind: Optional[str] = None,
) -> VerificationResult:
    gates = list(gates)
    ancilla_qubits = _normalize_ancillas(ancilla=ancilla, ancillas=ancillas)
    system_qubits = _normalize_systems(system=system, systems=systems)
    if expression_kind is None:
        expression_kind = "word" if any(_is_qsp_word_gate(gate.name) for gate in gates) else "op"
    state = BranchState.initial(
        num_system_qubits=len(system_qubits),
        num_ancillas=len(ancilla_qubits),
        expression_kind=expression_kind,
    )
    trace: List[TraceStep] = [TraceStep(0, None, state)] if keep_trace else []
    gate_profiles: List[GateProfileStep] = []

    for index, gate in enumerate(gates, start=1):
        counters = profile.ProfileCounters()
        if profile_gates:
            profile.set_current(counters)
        started = time.perf_counter()
        try:
            state = state.apply(gate, ancillas=ancilla_qubits, systems=system_qubits)
        finally:
            if profile_gates:
                profile.clear_current()
        elapsed = time.perf_counter() - started
        if keep_trace:
            trace.append(TraceStep(index, gate, state))
        if profile_gates:
            num_branches, total_terms, max_terms = _branch_profile_stats(state)
            gate_profiles.append(
                GateProfileStep(
                    gate_id=index,
                    gate_name=gate.name,
                    gate=str(gate),
                    num_nonzero_branches=num_branches,
                    total_operator_terms=total_terms,
                    max_terms_per_branch=max_terms,
                    time_this_gate=elapsed,
                    time_simplify=counters.simplify_sec,
                    time_combine=counters.combine_sec,
                )
            )

    return VerificationResult(
        final_state=state,
        trace=trace,
        gate_profiles=gate_profiles,
        ancilla_qubits=ancilla_qubits,
        system_qubits=system_qubits,
    )


def verify_qasm_file(
    path: str | Path,
    *,
    expected: str | OpExpr | None = None,
    ancilla: int = 0,
    ancillas: Optional[Sequence[int]] = None,
    system: int = 1,
    systems: Optional[Sequence[int]] = None,
    keep_trace: bool = False,
    profile_gates: bool = False,
    base: str | OpExpr | None = None,
    expected_polynomial: str | sp.Expr | None = None,
    hermitian_base: bool = False,
    compare_polynomial_only: bool = False,
    tolerance: float = DEFAULT_TOLERANCE,
) -> VerificationResult:
    gates = parse_qasm_file(path)
    ancilla_qubits = _normalize_ancillas(ancilla=ancilla, ancillas=ancillas)
    system_qubits = _normalize_systems(system=system, systems=systems)
    expected_expr = _normalize_expected(expected, num_system_qubits=len(system_qubits))
    expression_kind = (
        "word"
        if expected_polynomial is not None or any(_is_qsp_word_gate(gate.name) for gate in gates)
        else "op"
    )
    result = run_circuit(
        gates,
        ancillas=ancilla_qubits,
        systems=system_qubits,
        keep_trace=keep_trace,
        profile_gates=profile_gates,
        expression_kind=expression_kind,
    )
    qsp_normalized = None
    qsp_polynomial = None
    qsp_expected_polynomial = None
    qsp_actual = None
    qsp_expected = None
    qsp_status = None

    if expression_kind == "word":
        top_left = result.final_state.top_left()
        if not isinstance(top_left, WordExpr):
            raise ValueError("QSP word mode requires word-mode branch values")
        qsp_normalized = normalize_qsp_expression(top_left, hermitian_base=hermitian_base)

    if expected_polynomial is not None:
        if base is None and not compare_polynomial_only:
            raise ValueError("--base is required when --expected-polynomial is set")
        if qsp_normalized is None:
            raise ValueError("--expected-polynomial requires QSP word-mode branch values")
        if not is_h_polynomial_word(qsp_normalized):
            qsp_status = FAIL_GARBAGE
        else:
            qsp_polynomial = convert_h_word_to_polynomial(qsp_normalized)
            qsp_expected_polynomial = parse_polynomial(expected_polynomial)
            if compare_polynomial_only:
                qsp_status = PASS if polynomial_close(qsp_polynomial, qsp_expected_polynomial, tol=tolerance) else FAIL
            else:
                base_expr = _normalize_base(base, num_system_qubits=len(system_qubits))
                qsp_actual = eval_polynomial_on_pauliop(qsp_polynomial, base_expr)
                qsp_expected = eval_polynomial_on_pauliop(qsp_expected_polynomial, base_expr)
                qsp_status = PASS if pauli_expr_close(qsp_actual, qsp_expected, tol=tolerance) else FAIL

    return VerificationResult(
        final_state=result.final_state,
        trace=result.trace,
        gate_profiles=result.gate_profiles,
        expected=expected_expr,
        ancilla_qubits=ancilla_qubits,
        system_qubits=system_qubits,
        qsp_normalized=qsp_normalized,
        qsp_polynomial=qsp_polynomial,
        qsp_expected_polynomial=qsp_expected_polynomial,
        qsp_actual=qsp_actual,
        qsp_expected=qsp_expected,
        qsp_status=qsp_status,
        qsp_polynomial_only=compare_polynomial_only,
        tolerance=tolerance,
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
        elif result.status == PASS_UP_TO_SCALE:
            lines.append(PASS_UP_TO_SCALE)
            lines.append(f"scale = {_format_scalar(result.scale)}")
        else:
            diff = result.final_state.top_left() - result.expected
            lines.append(FAIL)
            diff_label = "B0" if result.final_state.num_ancillas == 1 else f"B[{_format_key((0,) * result.final_state.num_ancillas)}]"
            lines.append(f"{diff_label} - expected = {diff}")

    if result.qsp_normalized is not None or result.qsp_status is not None:
        if result.qsp_normalized is not None:
            key = (0,) * result.final_state.num_ancillas
            lines.append(f"Normalized B[{_format_key(key)}] = {result.qsp_normalized}")
        if result.qsp_polynomial is not None:
            lines.append(f"Polynomial = {_format_polynomial(result.qsp_polynomial)}")
        if result.qsp_expected_polynomial is not None:
            lines.append(f"Expected polynomial = {_format_polynomial(result.qsp_expected_polynomial)}")
        if result.qsp_polynomial_only:
            lines.append("Comparison mode = polynomial-only")
        if result.qsp_actual is not None:
            lines.append(f"Actual evaluated = {result.qsp_actual}")
        if result.qsp_expected is not None:
            lines.append(f"Expected evaluated = {result.qsp_expected}")
        if result.qsp_status is not None:
            lines.append(result.qsp_status)

    return "\n".join(lines)


def format_status(result: VerificationResult) -> str:
    if result.qsp_status is not None:
        return result.qsp_status
    if result.status is not None:
        return result.status
    return "NO_EXPECTED"


def format_gate_profiles_csv(result: VerificationResult) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=gate_profile_fieldnames())
    writer.writeheader()
    for row in gate_profile_rows(result):
        writer.writerow(row)
    return output.getvalue()


def gate_profile_rows(result: VerificationResult) -> list[dict[str, str]]:
    return [
        {
            "gate_id": str(step.gate_id),
            "gate_name": step.gate_name,
            "gate": step.gate,
            "num_nonzero_branches": str(step.num_nonzero_branches),
            "total_operator_terms": str(step.total_operator_terms),
            "max_terms_per_branch": str(step.max_terms_per_branch),
            "time_this_gate": f"{step.time_this_gate:.9f}",
            "time_simplify": f"{step.time_simplify:.9f}",
            "time_combine": f"{step.time_combine:.9f}",
        }
        for step in result.gate_profiles
    ]


def gate_profile_fieldnames() -> list[str]:
    return [
        "gate_id",
        "gate_name",
        "gate",
        "num_nonzero_branches",
        "total_operator_terms",
        "max_terms_per_branch",
        "time_this_gate",
        "time_simplify",
        "time_combine",
    ]


def _branch_profile_stats(state: BranchState) -> tuple[int, int, int]:
    term_counts = [_branch_term_count(branch) for branch in state.branches.values()]
    if not term_counts:
        return 0, 0, 0
    return len(term_counts), sum(term_counts), max(term_counts)


def _branch_term_count(branch: OpExpr | WordExpr) -> int:
    return len(branch.terms)


def _normalize_expected(expected: str | OpExpr | None, *, num_system_qubits: int) -> Optional[OpExpr]:
    if expected is None:
        return None
    if isinstance(expected, OpExpr):
        return expected.with_num_qubits(num_system_qubits)
    return parse_operator_expression(expected, num_qubits=num_system_qubits)


def _normalize_base(base: str | OpExpr, *, num_system_qubits: int) -> OpExpr:
    if isinstance(base, OpExpr):
        return base.with_num_qubits(num_system_qubits)
    return parse_operator_expression(base, num_qubits=num_system_qubits)


def _is_qsp_word_gate(name: str) -> bool:
    lowered = name.lower()
    return lowered in {"uh", "uhdg"} or re.fullmatch(r"c+uh(?:dg)?", lowered) is not None


def scalar_close(a: sp.Expr, b: sp.Expr, *, tol: float = DEFAULT_TOLERANCE) -> bool:
    diff = scalar_simplify(sp.sympify(a) - sp.sympify(b))
    if diff == 0:
        return True

    try:
        return abs(complex(sp.N(diff, 50))) <= tol
    except Exception:
        return False


def pauli_expr_close(actual: OpExpr, expected: OpExpr, *, tol: float = DEFAULT_TOLERANCE) -> bool:
    actual, expected = _align_operator_qubits(actual, expected)
    diff = (actual - expected).simplify()
    return all(scalar_close(coeff, 0, tol=tol) for coeff in diff.terms.values())


def polynomial_close(actual: sp.Expr, expected: sp.Expr, *, tol: float = DEFAULT_TOLERANCE) -> bool:
    x = sp.Symbol("x")
    diff = sp.expand(parse_polynomial(actual) - parse_polynomial(expected))
    if diff == 0:
        return True
    poly = sp.Poly(diff, x)
    return all(scalar_close(coeff, 0, tol=tol) for _, coeff in poly.terms())


def normalize_qsp_expression(expr: WordExpr, *, hermitian_base: bool = False) -> WordExpr:
    normalized = eliminate(expr)
    if hermitian_base:
        normalized = normalized.replace_atom("Hd", "H").simplify()
    return normalized


def convert_h_word_to_polynomial(expr: WordExpr) -> sp.Expr:
    x = sp.Symbol("x")
    polynomial = sp.Integer(0)
    for word, coeff in expr.terms.items():
        if any(atom != "H" for atom in word):
            raise ValueError(f"Cannot convert non-H word to polynomial: {word!r}")
        polynomial += coeff * x ** len(word)
    return sp.expand(polynomial)


def is_h_polynomial_word(expr: WordExpr) -> bool:
    return not expr.has_any_atom(GARBAGE_ATOMS | {"Hd"})


def parse_polynomial(polynomial: str | sp.Expr) -> sp.Expr:
    x = sp.Symbol("x")
    if isinstance(polynomial, str):
        text = polynomial.replace("^", "**")
        parsed = sp.sympify(text, locals={"x": x, "pi": sp.pi, "I": sp.I, "i": sp.I})
    else:
        parsed = sp.sympify(polynomial)
    return sp.expand(parsed)


def eval_polynomial_on_pauliop(polynomial: str | sp.Expr, base: OpExpr) -> OpExpr:
    x = sp.Symbol("x")
    parsed = parse_polynomial(polynomial)
    poly = sp.Poly(parsed, x)
    powers = {0: OpExpr.identity(base.num_qubits)}
    result = OpExpr.zero(base.num_qubits)

    for (degree,), coeff in poly.terms():
        if degree not in powers:
            current = powers[max(powers)]
            for exponent in range(max(powers) + 1, degree + 1):
                current = current * base
                powers[exponent] = current
        result += powers[degree].scale(coeff)
    return result


def _format_polynomial(polynomial: sp.Expr) -> str:
    return sp.sstr(scalar_simplify(polynomial)).replace("**", "^")


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


def proportional_scale(actual: OpExpr, expected: OpExpr, *, tol: float = DEFAULT_TOLERANCE) -> Tuple[bool, Optional[sp.Expr]]:
    actual, expected = _align_operator_qubits(actual, expected)
    actual_terms = _nonzero_terms(actual)
    expected_terms = _nonzero_terms(expected)

    if not expected_terms:
        return False, None

    p0 = next(iter(expected_terms))
    scale = scalar_simplify(actual_terms.get(p0, sp.Integer(0)) / expected_terms[p0])

    for pauli_string in set(actual_terms) | set(expected_terms):
        actual_coeff = actual_terms.get(pauli_string, sp.Integer(0))
        expected_coeff = expected_terms.get(pauli_string, sp.Integer(0))
        if not scalar_close(actual_coeff, scale * expected_coeff, tol=tol):
            return False, None

    return True, scale


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
    parser.add_argument("--base", help='Base block for QSP/QSVT polynomial checks, for example "(X + Z)/2".')
    parser.add_argument("--expected-polynomial", help='Expected QSP polynomial in x, for example "4*x^3 - 3*x".')
    parser.add_argument("--hermitian-base", action="store_true", help="Rewrite Hd atoms to H before QSP polynomial evaluation.")
    parser.add_argument(
        "--compare-polynomial-only",
        action="store_true",
        help="Compare QSP polynomials directly instead of evaluating them on --base.",
    )
    parser.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE, help="Numeric tolerance for coefficient comparisons.")
    parser.add_argument("--ancilla", default="0", help="Single ancilla qubit index. Ignored when --ancillas is set.")
    parser.add_argument("--ancillas", nargs="+", help='Ancilla qubits in branch-key order, e.g. q[0] q[1].')
    parser.add_argument("--system", default="1", help="Single system qubit index. Ignored when --systems is set.")
    parser.add_argument("--systems", nargs="+", help='System qubits in operator-index order, e.g. q[2] q[3].')
    parser.add_argument("--trace", action="store_true", help="Print every intermediate B0/B1 update.")
    parser.add_argument("--result-only", action="store_true", help="Print only PASS/FAIL status.")
    parser.add_argument("--profile-gates", action="store_true", help="Collect per-gate branch/term/timing counters.")
    parser.add_argument("--profile-output", type=Path, help="Write --profile-gates CSV output to this path instead of stdout.")
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
        base=args.base,
        expected_polynomial=args.expected_polynomial,
        hermitian_base=args.hermitian_base,
        compare_polynomial_only=args.compare_polynomial_only,
        tolerance=args.tolerance,
        profile_gates=args.profile_gates,
    )
    if args.profile_gates:
        profile_csv = format_gate_profiles_csv(result)
        if args.profile_output is not None:
            args.profile_output.parent.mkdir(parents=True, exist_ok=True)
            args.profile_output.write_text(profile_csv, encoding="utf-8")
        else:
            print(profile_csv, end="")
    if args.result_only:
        print(format_status(result))
    else:
        print(format_result(result, show_trace=args.trace))
    return 0 if result.success is not False else 1


if __name__ == "__main__":
    raise SystemExit(main())
