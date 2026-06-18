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
from .word_expr import GARBAGE_ATOMS, WordExpr, eliminate


PASS = "PASS"
PASS_UP_TO_GLOBAL_PHASE = "PASS_UP_TO_GLOBAL_PHASE"
FAIL = "FAIL"
FAIL_GARBAGE = "FAIL_GARBAGE"
DEFAULT_TOLERANCE = 1e-8
POLYNOMIAL_COMPARE_PARTS = ("full", "real", "imag")


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
    qsp_normalized: Optional[WordExpr] = None
    qsp_polynomial: Optional[sp.Expr] = None
    qsp_expected_polynomial: Optional[sp.Expr] = None
    qsp_actual: Optional[OpExpr] = None
    qsp_expected: Optional[OpExpr] = None
    qsp_status: Optional[str] = None
    qsp_compare_part: str = "full"
    qsp_polynomial_only: bool = False
    tolerance: float = DEFAULT_TOLERANCE

    @property
    def success(self) -> Optional[bool]:
        if self.qsp_status is not None:
            return self.qsp_status == PASS
        if self.expected is None:
            return None
        return self.status in {PASS, PASS_UP_TO_GLOBAL_PHASE}

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
        is_proportional, phase = proportional_phase(top_left, self.expected)
        if is_proportional:
            return PASS_UP_TO_GLOBAL_PHASE
        return FAIL

    @property
    def global_phase(self) -> Optional[sp.Expr]:
        top_left = self.final_state.top_left()
        if self.expected is None or not isinstance(top_left, OpExpr) or pauli_expr_close(top_left, self.expected, tol=self.tolerance):
            return None
        is_proportional, phase = proportional_phase(top_left, self.expected)
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
    expression_kind: Optional[str] = None,
) -> VerificationResult:
    gates = list(gates)
    ancilla_qubits = _normalize_ancillas(ancilla=ancilla, ancillas=ancillas)
    system_qubits = _normalize_systems(system=system, systems=systems)
    if expression_kind is None:
        expression_kind = "word" if any(gate.name.lower() in {"uh", "uhdg"} for gate in gates) else "op"
    state = BranchState.initial(
        num_system_qubits=len(system_qubits),
        num_ancillas=len(ancilla_qubits),
        expression_kind=expression_kind,
    )
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
    base: str | OpExpr | None = None,
    expected_polynomial: str | sp.Expr | None = None,
    hermitian_base: bool = False,
    compare_polynomial_part: str = "full",
    compare_polynomial_only: bool = False,
    tolerance: float = DEFAULT_TOLERANCE,
) -> VerificationResult:
    if compare_polynomial_part not in POLYNOMIAL_COMPARE_PARTS:
        raise ValueError(f"compare_polynomial_part must be one of {POLYNOMIAL_COMPARE_PARTS}")

    gates = parse_qasm_file(path)
    ancilla_qubits = _normalize_ancillas(ancilla=ancilla, ancillas=ancillas)
    system_qubits = _normalize_systems(system=system, systems=systems)
    expected_expr = _normalize_expected(expected, num_system_qubits=len(system_qubits))
    expression_kind = "word" if expected_polynomial is not None or any(gate.name.lower() in {"uh", "uhdg"} for gate in gates) else "op"
    result = run_circuit(
        gates,
        ancillas=ancilla_qubits,
        systems=system_qubits,
        keep_trace=keep_trace,
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
        if qsp_normalized.has_any_atom(GARBAGE_ATOMS):
            qsp_status = FAIL_GARBAGE
        else:
            qsp_polynomial = convert_h_word_to_polynomial(qsp_normalized)
            qsp_expected_polynomial = parse_polynomial(expected_polynomial)
            actual_polynomial = project_polynomial_part(qsp_polynomial, compare_polynomial_part)
            if compare_polynomial_only:
                qsp_status = PASS if polynomial_close(actual_polynomial, qsp_expected_polynomial, tol=tolerance) else FAIL
            else:
                base_expr = _normalize_base(base, num_system_qubits=len(system_qubits))
                qsp_actual = eval_polynomial_on_pauliop(actual_polynomial, base_expr)
                qsp_expected = eval_polynomial_on_pauliop(qsp_expected_polynomial, base_expr)
                qsp_status = PASS if pauli_expr_close(qsp_actual, qsp_expected, tol=tolerance) else FAIL

    return VerificationResult(
        final_state=result.final_state,
        trace=result.trace,
        expected=expected_expr,
        ancilla_qubits=ancilla_qubits,
        system_qubits=system_qubits,
        qsp_normalized=qsp_normalized,
        qsp_polynomial=qsp_polynomial,
        qsp_expected_polynomial=qsp_expected_polynomial,
        qsp_actual=qsp_actual,
        qsp_expected=qsp_expected,
        qsp_status=qsp_status,
        qsp_compare_part=compare_polynomial_part,
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
        elif result.status == PASS_UP_TO_GLOBAL_PHASE:
            lines.append(PASS_UP_TO_GLOBAL_PHASE)
            lines.append(f"phase = {_format_scalar(result.global_phase)}")
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
        if result.qsp_compare_part != "full":
            lines.append(f"Compared actual polynomial part = {result.qsp_compare_part}")
        if result.qsp_polynomial_only:
            lines.append("Comparison mode = polynomial-only")
        if result.qsp_actual is not None:
            lines.append(f"Actual evaluated = {result.qsp_actual}")
        if result.qsp_expected is not None:
            lines.append(f"Expected evaluated = {result.qsp_expected}")
        if result.qsp_status is not None:
            lines.append(result.qsp_status)

    return "\n".join(lines)


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
    diff = scalar_simplify(sp.expand(parse_polynomial(actual) - parse_polynomial(expected)))
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
    return scalar_simplify(sp.expand(polynomial))


def parse_polynomial(polynomial: str | sp.Expr) -> sp.Expr:
    x = sp.Symbol("x")
    if isinstance(polynomial, str):
        text = polynomial.replace("^", "**")
        parsed = sp.sympify(text, locals={"x": x, "pi": sp.pi, "I": sp.I, "i": sp.I})
    else:
        parsed = sp.sympify(polynomial)
    return scalar_simplify(sp.expand(parsed))


def project_scalar_part(value: sp.Expr, part: str) -> sp.Expr:
    if part == "full":
        return scalar_simplify(value)
    if part == "real":
        return scalar_simplify(sp.re(value))
    if part == "imag":
        return scalar_simplify(sp.im(value))
    raise ValueError(f"part must be one of {POLYNOMIAL_COMPARE_PARTS}")


def project_polynomial_part(polynomial: str | sp.Expr, part: str) -> sp.Expr:
    parsed = parse_polynomial(polynomial)
    if part == "full":
        return parsed

    x = sp.Symbol("x")
    poly = sp.Poly(parsed, x)
    projected = sp.Integer(0)
    for (degree,), coeff in poly.terms():
        projected += project_scalar_part(coeff, part) * x ** degree
    return scalar_simplify(sp.expand(projected))


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
    parser.add_argument("--base", help='Base block for QSP/QSVT polynomial checks, for example "(X + Z)/2".')
    parser.add_argument("--expected-polynomial", help='Expected QSP polynomial in x, for example "4*x^3 - 3*x".')
    parser.add_argument("--hermitian-base", action="store_true", help="Rewrite Hd atoms to H before QSP polynomial evaluation.")
    parser.add_argument(
        "--compare-polynomial-part",
        choices=POLYNOMIAL_COMPARE_PARTS,
        default="full",
        help="Compare the full, real, or imaginary part of the actual polynomial coefficients against --expected-polynomial.",
    )
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
        compare_polynomial_part=args.compare_polynomial_part,
        compare_polynomial_only=args.compare_polynomial_only,
        tolerance=args.tolerance,
    )
    print(format_result(result, show_trace=args.trace))
    return 0 if result.success is not False else 1


if __name__ == "__main__":
    raise SystemExit(main())
