from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from functools import wraps
import io
import math
from pathlib import Path
import re
import time
from typing import Any, Iterable, List, Mapping, Optional, Sequence, Tuple

import sympy as sp

from .branch_state import BranchState, Gate, UnsupportedGateError
from .expr import OpExpr, parse_operator_expression
from .interval import acb, arb, ctx, format_ball, is_ball, is_exact_zero, to_acb, to_arb
from .interval import divide as scalar_divide
from .interval import multiply as scalar_multiply
from .interval import subtract as scalar_subtract
from . import profile
from .qasm_parser import parse_qasm_file
from .scalar import scalar_simplify
from .word_expr import GARBAGE_ATOMS, WordExpr, eliminate


PASS = "PASS"
PASS_EXACT = "PASS_EXACT"
PASS_UP_TO_SCALE = "PASS_UP_TO_SCALE"
FAIL = "FAIL"
FAIL_GARBAGE = "FAIL_GARBAGE"
DEFAULT_TOLERANCE = 1e-8
DEFAULT_BLOCK_ENCODING_CHECK_PRECISION = 60
DEFAULT_QSP_CALCULATION_PRECISION = 60
BALL_PRECISION_DIGITS_PER_SYSTEM_QUBIT = 2
DEFAULT_MAX_APPROXIMATION_GRID_POINTS = 10_000_000


@dataclass(frozen=True)
class BallPolynomial:
    """A polynomial whose coefficients are rigorous complex balls."""

    terms: Mapping[int, acb]

    @property
    def degree(self) -> int:
        return max(self.terms, default=0)


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
class ApproximationCheck:
    tau: float
    epsilon: float
    scale: sp.Expr
    beta: complex
    polynomial_degree: int
    approximation_degree: int
    angular_spacing: float
    num_grid_points: int
    max_grid_error: float
    worst_x: float
    epsilon_ball: arb
    beta_ball: acb
    max_grid_error_upper: arb

    @property
    def success(self) -> bool:
        return bool(self.max_grid_error_upper <= (self.epsilon_ball / 2).lower())

    @property
    def polynomial_lipschitz(self) -> float:
        return self.polynomial_derivative_bound

    @property
    def polynomial_derivative_bound(self) -> float:
        return float(self.approximation_degree)

    @property
    def target_lipschitz(self) -> float:
        return abs(complex(sp.N(self.scale, 50)) * self.tau)

    @property
    def spacing(self) -> float:
        return self.angular_spacing


@dataclass(frozen=True)
class BlockEncodingProportionalityCheck:
    epsilon: float
    alpha: float
    projection_length: float
    residual_norm: float
    threshold: float
    coefficient_norm: float
    target_norm: float
    alpha_ball: arb
    residual_norm_upper: arb
    threshold_ball: arb

    @property
    def success(self) -> bool:
        return bool(
            self.alpha_ball.lower() > 0
            and self.residual_norm_upper <= self.threshold_ball.lower()
        )


@dataclass(frozen=True)
class VerificationResult:
    final_state: BranchState
    trace: Sequence[TraceStep]
    gate_profiles: Sequence[GateProfileStep] = ()
    expected: Optional[OpExpr] = None
    ancilla_qubits: Tuple[int, ...] = (0,)
    system_qubits: Tuple[int, ...] = (1,)
    qsp_normalized: Optional[WordExpr] = None
    qsp_polynomial: Optional[sp.Expr | BallPolynomial] = None
    qsp_expected_polynomial: Optional[sp.Expr | BallPolynomial] = None
    qsp_actual: Optional[OpExpr] = None
    qsp_expected: Optional[OpExpr] = None
    qsp_approximation: Optional[ApproximationCheck] = None
    qsp_status: Optional[str] = None
    qsp_polynomial_only: bool = False
    block_encoding_epsilon: Optional[float] = None
    tolerance: float = DEFAULT_TOLERANCE

    @property
    def success(self) -> Optional[bool]:
        if self.qsp_status is not None:
            return self.qsp_status == PASS
        if self.expected is None:
            return None
        return self.status in {PASS, PASS_EXACT, PASS_UP_TO_SCALE}

    @property
    def status(self) -> Optional[str]:
        if self.qsp_status is not None:
            return self.qsp_status
        if self.expected is None:
            return None
        top_left = self.final_state.top_left()
        if not isinstance(top_left, OpExpr):
            return FAIL
        if self.block_encoding_epsilon is not None:
            check = block_encoding_proportionality_check(
                top_left,
                self.expected,
                epsilon=self.block_encoding_epsilon,
            )
            return PASS if check.success else FAIL
        if pauli_expr_close(top_left, self.expected, tol=self.tolerance):
            return PASS_EXACT
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

    @property
    def block_encoding_check(self) -> Optional[BlockEncodingProportionalityCheck]:
        top_left = self.final_state.top_left()
        if self.expected is None or self.block_encoding_epsilon is None or not isinstance(top_left, OpExpr):
            return None
        return block_encoding_proportionality_check(
            top_left,
            self.expected,
            epsilon=self.block_encoding_epsilon,
        )


def _with_ball_precision(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        requested = kwargs.get("qsp_calculation_precision", DEFAULT_QSP_CALCULATION_PRECISION)
        systems = kwargs.get("systems")
        num_system_qubits = len(systems) if systems is not None else 1
        precision = max(requested, BALL_PRECISION_DIGITS_PER_SYSTEM_QUBIT * num_system_qubits)
        kwargs["qsp_calculation_precision"] = precision
        with ctx.workdps(precision):
            return function(*args, **kwargs)

    return wrapped


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
    rewrite_after_each_step: bool = True,
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
            if (
                rewrite_after_each_step
                and expression_kind == "word"
                and _is_qsp_word_gate(gate.name)
            ):
                state = _rewrite_word_branches(state)
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


@_with_ball_precision
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
    rewrite_after_each_step: bool = True,
    qsp_calculation_precision: int = DEFAULT_QSP_CALCULATION_PRECISION,
    base: str | OpExpr | None = None,
    expected_polynomial: str | sp.Expr | None = None,
    extract_qsp_polynomial: bool = False,
    hermitian_base: bool = False,
    compare_polynomial_only: bool = False,
    target_exp_tau: float | None = None,
    target_exp_epsilon: float | None = None,
    target_exp_scale: str | sp.Expr = 1,
    max_approx_grid_points: int | None = DEFAULT_MAX_APPROXIMATION_GRID_POINTS,
    block_encoding_epsilon: float | None = None,
    tolerance: float = DEFAULT_TOLERANCE,
) -> VerificationResult:
    if qsp_calculation_precision < 1:
        raise ValueError("qsp_calculation_precision must be positive")
    if (target_exp_tau is None) != (target_exp_epsilon is None):
        raise ValueError("--hamsim-tau and --hamsim-epsilon must be set together")
    if target_exp_epsilon is not None and target_exp_epsilon <= 0:
        raise ValueError("--hamsim-epsilon must be positive")
    if block_encoding_epsilon is not None:
        if block_encoding_epsilon <= 0:
            raise ValueError("--block-encoding-epsilon must be positive")
        if expected is None:
            raise ValueError("--block-encoding-epsilon requires --expected")
    selected_modes = sum(
        (
            expected is not None,
            expected_polynomial is not None,
            target_exp_tau is not None,
        )
    )
    if selected_modes > 1:
        raise ValueError("Choose exactly one verification route: --expected, --expected-polynomial, or --hamsim-*")

    # Decimal QASM parameters enter the verifier as Arb intervals.  This keeps
    # their source-text meaning and rounding enclosure through every gate,
    # including non-QSP block-encoding circuits.
    gates = parse_qasm_file(path, numeric_precision=qsp_calculation_precision)
    ancilla_qubits = _normalize_ancillas(ancilla=ancilla, ancillas=ancillas)
    system_qubits = _normalize_systems(system=system, systems=systems)
    expected_expr = _normalize_expected(expected, num_system_qubits=len(system_qubits))
    expression_kind = (
        "word"
        if expected_polynomial is not None
        or extract_qsp_polynomial
        or target_exp_tau is not None
        or target_exp_epsilon is not None
        or any(_is_qsp_word_gate(gate.name) for gate in gates)
        else "op"
    )
    if expression_kind == "word":
        _validate_qsp_canonical_form(gates, ancillas=ancilla_qubits, systems=system_qubits)
    result = run_circuit(
        gates,
        ancillas=ancilla_qubits,
        systems=system_qubits,
        keep_trace=keep_trace,
        profile_gates=profile_gates,
        expression_kind=expression_kind,
        rewrite_after_each_step=rewrite_after_each_step,
    )
    qsp_normalized = None
    qsp_polynomial = None
    qsp_expected_polynomial = None
    qsp_actual = None
    qsp_expected = None
    qsp_approximation = None
    qsp_status = None

    if expression_kind == "word":
        top_left = result.final_state.top_left()
        if not isinstance(top_left, WordExpr):
            raise ValueError("QSP word mode requires word-mode branch values")
        qsp_normalized = normalize_qsp_expression(top_left, hermitian_base=hermitian_base)

    needs_qsp_polynomial = expected_polynomial is not None or target_exp_tau is not None or extract_qsp_polynomial
    if needs_qsp_polynomial:
        if expected_polynomial is not None and base is None and not compare_polynomial_only:
            raise ValueError("--base is required when --expected-polynomial is set")
        if qsp_normalized is None:
            raise ValueError("QSP polynomial checks require QSP word-mode branch values")
        if not is_h_polynomial_word(qsp_normalized):
            qsp_status = FAIL_GARBAGE
        else:
            qsp_polynomial = convert_h_word_to_polynomial(qsp_normalized)

    check_statuses: list[str] = []
    if expected_polynomial is not None and qsp_status != FAIL_GARBAGE:
        qsp_expected_polynomial = parse_polynomial(expected_polynomial)
        if compare_polynomial_only:
            check_statuses.append(
                PASS if polynomial_close(qsp_polynomial, qsp_expected_polynomial, tol=tolerance) else FAIL
            )
        else:
            base_expr = _normalize_base(base, num_system_qubits=len(system_qubits))
            qsp_actual = eval_polynomial_on_pauliop(qsp_polynomial, base_expr)
            qsp_expected = eval_polynomial_on_pauliop(qsp_expected_polynomial, base_expr)
            check_statuses.append(PASS if pauli_expr_close(qsp_actual, qsp_expected, tol=tolerance) else FAIL)

    if target_exp_tau is not None and qsp_status != FAIL_GARBAGE:
        qsp_approximation = verify_polynomial_approximates_exp(
            qsp_polynomial,
            tau=target_exp_tau,
            epsilon=target_exp_epsilon,
            precision=qsp_calculation_precision,
            max_grid_points=max_approx_grid_points,
        )
        check_statuses.append(PASS if qsp_approximation.success else FAIL)

    if qsp_status != FAIL_GARBAGE and check_statuses:
        qsp_status = PASS if all(status == PASS for status in check_statuses) else FAIL

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
        qsp_approximation=qsp_approximation,
        qsp_status=qsp_status,
        qsp_polynomial_only=compare_polynomial_only,
        block_encoding_epsilon=block_encoding_epsilon,
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
        if result.block_encoding_epsilon is None:
            expected_label = "Expected B0" if result.final_state.num_ancillas == 1 else f"Expected B[{_format_key((0,) * result.final_state.num_ancillas)}]"
        else:
            expected_label = "Target H"
        lines.append(f"{expected_label} = {result.expected}")
        if result.status == PASS_EXACT:
            lines.append(PASS_EXACT)
        elif result.status == PASS_UP_TO_SCALE:
            lines.append(PASS_UP_TO_SCALE)
            lines.append(f"scale = {_format_scalar(result.scale)}")
        elif result.status == PASS:
            check = result.block_encoding_check
            lines.append(PASS)
            if check is not None:
                lines.append(f"alpha = {check.alpha:.12g}")
                lines.append(f"Projection length = {check.projection_length:.12g}")
                lines.append(f"Coefficient residual upper bound = {check.residual_norm:.12g}")
                lines.append(f"Coefficient residual threshold = {check.threshold:.12g}")
                lines.append(
                    "Coefficient residual threshold ball = "
                    f"{format_ball(check.threshold_ball, 12)}"
                )
        else:
            diff = result.final_state.top_left() - result.expected
            lines.append(FAIL)
            check = result.block_encoding_check
            if check is not None:
                lines.append(f"alpha = {check.alpha:.12g}")
                lines.append(f"Projection length = {check.projection_length:.12g}")
                lines.append(f"Coefficient residual upper bound = {check.residual_norm:.12g}")
                lines.append(f"Coefficient residual threshold = {check.threshold:.12g}")
                lines.append(
                    "Coefficient residual threshold ball = "
                    f"{format_ball(check.threshold_ball, 12)}"
                )
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
        if result.qsp_approximation is not None:
            check = result.qsp_approximation
            if scalar_close(check.scale, 1):
                target_label = f"exp(-i*x*{check.tau:.12g})"
            else:
                target_label = f"{_format_scalar(check.scale)} * exp(-i*x*{check.tau:.12g})"
            lines.append(f"Approximation target = {target_label}")
            lines.append(f"Approximation beta = {_format_complex(check.beta)}")
            lines.append(f"Approximation beta ball = {format_ball(check.beta_ball, 12)}")
            lines.append(f"Approximation epsilon = {check.epsilon:.12g}")
            lines.append(f"Approximation polynomial degree = {check.polynomial_degree}")
            lines.append(f"Approximation Chebyshev degree D = {check.approximation_degree}")
            lines.append(f"Approximation grid angular spacing = {check.angular_spacing:.12g}")
            lines.append(f"Approximation grid roots = {check.num_grid_points}")
            lines.append(
                f"Approximation max grid error upper bound = {check.max_grid_error:.12g} "
                f"at x = {check.worst_x:.12g}"
            )
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


def _rewrite_word_branches(state: BranchState) -> BranchState:
    """Normalize every live QSP word branch after a signal-oracle step."""
    return BranchState(
        state.num_ancillas,
        state.num_system_qubits,
        {
            key: eliminate(branch) if isinstance(branch, WordExpr) else branch
            for key, branch in state.branches.items()
        },
        expression_kind=state.expression_kind,
    )


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


def _validate_qsp_canonical_form(
    gates: Sequence[Gate],
    *,
    ancillas: Sequence[int],
    systems: Sequence[int],
) -> None:
    """Reject QSP gates that do not preserve the zero/nonzero block abstraction."""

    ancilla_set = set(ancillas)
    system_set = set(systems)
    signal_layouts: list[tuple[tuple[int, ...], tuple[int, ...]]] = []

    for gate in gates:
        control_count = _qsp_signal_control_count(gate.name)
        if control_count is None:
            continue
        controls = gate.qubits[:control_count]
        operands = gate.qubits[control_count:]
        block_ancillas = tuple(qubit for qubit in operands if qubit in ancilla_set)
        signal_systems = tuple(qubit for qubit in operands if qubit in system_set)
        unknown = tuple(qubit for qubit in operands if qubit not in ancilla_set and qubit not in system_set)
        if unknown or not block_ancillas or not signal_systems:
            raise UnsupportedGateError(
                f"{gate.name} must contain block ancillas and system qubits from the declared mappings"
            )
        if any(control not in ancilla_set for control in controls):
            raise UnsupportedGateError(f"{gate.name} controls must be declared ancillas")
        if set(controls) & set(block_ancillas):
            raise UnsupportedGateError(f"{gate.name} selector controls cannot also be block ancillas")
        signal_layouts.append((block_ancillas, signal_systems))

    if not signal_layouts:
        return

    block_ancillas, signal_systems = signal_layouts[0]
    block_set = set(block_ancillas)
    for other_blocks, other_systems in signal_layouts[1:]:
        if set(other_blocks) != block_set or tuple(other_systems) != tuple(signal_systems):
            raise UnsupportedGateError(
                "QSP canonical form requires every UH/UHdg invocation to use the same "
                "block-ancilla and system registers"
            )

    # With one block ancilla, the non-principal subspace has only one branch,
    # so the ordinary branch semantics already represents it exactly.
    if len(block_set) == 1:
        return

    index = 0
    while index < len(gates):
        gate = gates[index]
        if _is_qsp_word_gate(gate.name) or not (set(gate.qubits) & block_set):
            index += 1
            continue

        prefix: list[int] = []
        cursor = index
        while cursor < len(gates):
            candidate = gates[cursor]
            if candidate.name.lower() != "x" or len(candidate.qubits) != 1 or candidate.qubits[0] not in block_set:
                break
            prefix.append(candidate.qubits[0])
            cursor += 1

        if len(prefix) != len(block_set) or set(prefix) != block_set:
            raise UnsupportedGateError(
                f"Gate {gate} distinguishes individual nonzero block-ancilla states; "
                "QSP input must use canonical all-zero projector controls"
            )
        if cursor >= len(gates) or not _is_canonical_block_projector_control(
            gates[cursor], block_set=block_set, ancilla_set=ancilla_set
        ):
            raise UnsupportedGateError(
                "Open-control X gates on block ancillas must be followed by one controlled-X "
                "whose controls include the complete block register"
            )

        controlled_x = gates[cursor]
        cursor += 1
        suffix = tuple(reversed(prefix))
        if _matches_block_x_suffix(gates, cursor, suffix):
            index = cursor + len(suffix)
            continue

        while cursor < len(gates) and not _same_gate(gates[cursor], controlled_x):
            if set(gates[cursor].qubits) & block_set:
                raise UnsupportedGateError(
                    "A canonical block-projector phase sandwich may not access a block ancilla "
                    "between its paired controlled-X gates"
                )
            cursor += 1
        if cursor >= len(gates):
            raise UnsupportedGateError("Unpaired controlled-X in block-projector phase sandwich")
        cursor += 1
        if not _matches_block_x_suffix(gates, cursor, suffix):
            raise UnsupportedGateError(
                "Block-projector phase sandwich must restore every open-control X in reverse order"
            )
        index = cursor + len(suffix)


def _qsp_signal_control_count(name: str) -> int | None:
    lowered = name.lower()
    if lowered in {"uh", "uhdg"}:
        return 0
    match = re.fullmatch(r"(c+)uh(?:dg)?", lowered)
    return len(match.group(1)) if match is not None else None


def _is_canonical_block_projector_control(
    gate: Gate,
    *,
    block_set: set[int],
    ancilla_set: set[int],
) -> bool:
    if gate.name.lower() not in {"cx", "mcx"} or len(gate.qubits) < 2:
        return False
    controls = gate.qubits[:-1]
    target = gate.qubits[-1]
    return (
        block_set.issubset(controls)
        and len(set(controls)) == len(controls)
        and all(control in ancilla_set for control in controls)
        and target in ancilla_set - block_set
        and target not in controls
    )


def _matches_block_x_suffix(gates: Sequence[Gate], start: int, suffix: tuple[int, ...]) -> bool:
    if start + len(suffix) > len(gates):
        return False
    return all(
        gates[start + offset].name.lower() == "x" and gates[start + offset].qubits == (qubit,)
        for offset, qubit in enumerate(suffix)
    )


def _same_gate(left: Gate, right: Gate) -> bool:
    return left.name.lower() == right.name.lower() and left.qubits == right.qubits and left.parameter == right.parameter


def scalar_close(a: sp.Expr, b: sp.Expr, *, tol: float = DEFAULT_TOLERANCE) -> bool:
    try:
        diff = scalar_simplify(scalar_subtract(a, b))
    except Exception:
        return False
    if is_exact_zero(diff):
        return True

    try:
        with ctx.workdps(DEFAULT_BLOCK_ENCODING_CHECK_PRECISION):
            return bool(to_acb(diff).abs_upper() <= to_arb(tol).lower())
    except Exception:
        return False


def pauli_expr_close(actual: OpExpr, expected: OpExpr, *, tol: float = DEFAULT_TOLERANCE) -> bool:
    actual, expected = _align_operator_qubits(actual, expected)
    diff = (actual - expected).simplify()
    return all(scalar_close(coeff, 0, tol=tol) for coeff in diff.terms.values())


def polynomial_close(
    actual: sp.Expr | BallPolynomial,
    expected: sp.Expr | BallPolynomial,
    *,
    tol: float = DEFAULT_TOLERANCE,
) -> bool:
    actual_terms = _polynomial_terms(actual)
    expected_terms = _polynomial_terms(expected)
    return all(
        scalar_close(actual_terms.get(degree, 0), expected_terms.get(degree, 0), tol=tol)
        for degree in set(actual_terms) | set(expected_terms)
    )


def verify_polynomial_approximates_exp(
    polynomial: str | sp.Expr | BallPolynomial,
    *,
    tau: float,
    epsilon: float,
    scale: str | sp.Expr = 1,
    precision: int = DEFAULT_QSP_CALCULATION_PRECISION,
    max_grid_points: int | None = DEFAULT_MAX_APPROXIMATION_GRID_POINTS,
) -> ApproximationCheck:
    if epsilon <= 0:
        raise ValueError("epsilon must be positive")
    if precision < 1:
        raise ValueError("precision must be positive")
    parsed = parse_polynomial(polynomial)
    scale_expr = parse_scalar_expression(scale)
    coeffs = _acb_polynomial_coefficients(parsed, precision=precision)
    polynomial_degree = len(coeffs) - 1
    approximation_degree = hamsim_chebyshev_approximation_degree(tau=tau, epsilon=epsilon)
    num_grid_points = hamsim_approximation_grid_size(
        tau=tau,
        epsilon=epsilon,
        polynomial_degree=polynomial_degree,
    )
    if max_grid_points is not None:
        if max_grid_points < 1:
            raise ValueError("max_grid_points must be at least 1, or None to disable the guard")
        if num_grid_points > max_grid_points:
            raise ValueError(
                "Approximation grid requires "
                f"{num_grid_points} Chebyshev roots by M = 4*max(D,d) "
                f"(D={approximation_degree}, degree={polynomial_degree}, tau={tau:.12g}, epsilon={epsilon:.12g}); "
                f"max_grid_points is {max_grid_points}"
            )
    with ctx.workdps(precision):
        pi = arb.pi()
        angular_spacing_ball = pi / num_grid_points
        tau_ball = to_arb(tau)
        epsilon_ball = to_arb(epsilon)
        scale_value = to_acb(scale_expr)
        grid_values: list[tuple[arb, acb, acb]] = []
        inner_product = acb(0)
        polynomial_norm_sq = arb(0)
        for index in range(1, num_grid_points + 1):
            point = ((2 * index - 1) * pi / (2 * num_grid_points)).cos()
            actual = _eval_acb_polynomial(coeffs, point)
            expected = scale_value * (-acb(0, 1) * point * tau_ball).exp()
            grid_values.append((point, actual, expected))
            inner_product += actual.conjugate() * expected
            polynomial_norm_sq += _abs_square_ball(actual)
        if polynomial_norm_sq.contains(0):
            beta_ball = acb(0)
        else:
            beta_ball = inner_product / polynomial_norm_sq

        max_error_upper = arb(0)
        worst_x_ball = arb(-1)
        for point, actual, expected in grid_values:
            error_upper = (beta_ball * actual - expected).abs_upper()
            if error_upper > max_error_upper:
                max_error_upper = error_upper
                worst_x_ball = point

        angular_spacing = float(angular_spacing_ball.mid())
        beta = complex(float(beta_ball.real.mid()), float(beta_ball.imag.mid()))
        max_error = _arb_upper_float(max_error_upper)
        worst_x = float(worst_x_ball.mid())

    return ApproximationCheck(
        tau=tau,
        epsilon=epsilon,
        scale=scale_expr,
        beta=beta,
        polynomial_degree=polynomial_degree,
        approximation_degree=approximation_degree,
        angular_spacing=angular_spacing,
        num_grid_points=num_grid_points,
        max_grid_error=max_error,
        worst_x=worst_x,
        epsilon_ball=epsilon_ball,
        beta_ball=beta_ball,
        max_grid_error_upper=max_error_upper,
    )


def hamsim_chebyshev_approximation_degree(*, tau: float, epsilon: float) -> int:
    if epsilon <= 0:
        raise ValueError("epsilon must be positive")
    with ctx.workdps(DEFAULT_QSP_CALCULATION_PRECISION):
        tau_ball = to_arb(tau)
        epsilon_ball = to_arb(epsilon)
        bound = arb(1).exp() * abs(tau_ball) / 2 + (arb(32) / epsilon_ball).log()
        degree = bound.upper().ceil().unique_fmpz()
        if degree is None:
            raise ValueError("Could not certify the Hamiltonian-simulation approximation degree")
        return max(0, int(degree))


def hamsim_approximation_grid_size(*, tau: float, epsilon: float, polynomial_degree: int) -> int:
    if polynomial_degree < 0:
        raise ValueError("polynomial_degree must be nonnegative")
    approximation_degree = hamsim_chebyshev_approximation_degree(tau=tau, epsilon=epsilon)
    return max(1, 4 * max(approximation_degree, polynomial_degree))


def normalize_qsp_expression(expr: WordExpr, *, hermitian_base: bool = False) -> WordExpr:
    normalized = eliminate(expr)
    if hermitian_base:
        normalized = normalized.replace_atom("Hd", "H").simplify()
    return normalized


def convert_h_word_to_polynomial(expr: WordExpr) -> sp.Expr | BallPolynomial:
    if any(is_ball(coeff) for coeff in expr.terms.values()):
        terms: dict[int, acb] = {}
        for word, coeff in expr.terms.items():
            if any(atom != "H" for atom in word):
                raise ValueError(f"Cannot convert non-H word to polynomial: {word!r}")
            degree = len(word)
            terms[degree] = terms.get(degree, acb(0)) + to_acb(coeff)
        return BallPolynomial({degree: coeff for degree, coeff in terms.items() if not coeff.is_zero()})

    x = sp.Symbol("x")
    polynomial = sp.Integer(0)
    for word, coeff in expr.terms.items():
        if any(atom != "H" for atom in word):
            raise ValueError(f"Cannot convert non-H word to polynomial: {word!r}")
        polynomial += coeff * x ** len(word)
    return sp.expand(polynomial)


def is_h_polynomial_word(expr: WordExpr) -> bool:
    return not expr.has_any_atom(GARBAGE_ATOMS | {"Hd"})


def parse_polynomial(polynomial: str | sp.Expr | BallPolynomial) -> sp.Expr | BallPolynomial:
    if isinstance(polynomial, BallPolynomial):
        return polynomial
    x = sp.Symbol("x")
    if isinstance(polynomial, str):
        text = polynomial.replace("^", "**")
        parsed = sp.sympify(
            text,
            locals={"x": x, "pi": sp.pi, "I": sp.I, "i": sp.I},
            rational=True,
        )
    else:
        parsed = sp.sympify(polynomial)
    return sp.expand(parsed)


def parse_scalar_expression(value: str | sp.Expr) -> sp.Expr:
    if isinstance(value, str):
        return sp.sympify(
            value.replace("^", "**"),
            locals={"pi": sp.pi, "I": sp.I, "i": sp.I},
            rational=True,
        )
    return sp.sympify(value)


def rescale_polynomial_for_target(
    polynomial: str | sp.Expr | BallPolynomial,
    scale: str | sp.Expr,
) -> sp.Expr | BallPolynomial:
    scale_expr = parse_scalar_expression(scale)
    scale_ball = to_acb(scale_expr)
    if scale_ball.abs_lower() <= 0:
        raise ValueError("target_exp_scale must be nonzero for rescaled approximation checking")
    parsed = parse_polynomial(polynomial)
    if isinstance(parsed, BallPolynomial):
        return BallPolynomial({degree: coeff / scale_ball for degree, coeff in parsed.terms.items()})
    return sp.expand(parsed / scale_expr)


def eval_polynomial_on_pauliop(polynomial: str | sp.Expr | BallPolynomial, base: OpExpr) -> OpExpr:
    parsed = parse_polynomial(polynomial)
    terms = _polynomial_terms(parsed)
    powers = {0: OpExpr.identity(base.num_qubits)}
    result = OpExpr.zero(base.num_qubits)

    for degree, coeff in sorted(terms.items(), reverse=True):
        if degree not in powers:
            current = powers[max(powers)]
            for exponent in range(max(powers) + 1, degree + 1):
                current = current * base
                powers[exponent] = current
        result += powers[degree].scale(coeff)
    return result


def _polynomial_terms(polynomial: str | sp.Expr | BallPolynomial) -> dict[int, Any]:
    parsed = parse_polynomial(polynomial)
    if isinstance(parsed, BallPolynomial):
        return dict(parsed.terms)
    x = sp.Symbol("x")
    poly = sp.Poly(parsed, x)
    if poly.is_zero:
        return {}
    return {degree: coeff for (degree,), coeff in poly.terms()}


def _acb_polynomial_coefficients(
    polynomial: str | sp.Expr | BallPolynomial,
    *,
    precision: int,
) -> list[acb]:
    with ctx.workdps(precision):
        terms = _polynomial_terms(polynomial)
        if not terms:
            return [acb(0)]
        coeffs = [acb(0) for _ in range(max(terms) + 1)]
        for degree, coeff in terms.items():
            coeffs[degree] = to_acb(coeff)
        return coeffs


def _eval_acb_polynomial(coeffs: Sequence[acb], point: arb) -> acb:
    value = acb(0)
    for coeff in reversed(coeffs):
        value = value * point + coeff
    return value


def _format_polynomial(polynomial: sp.Expr | BallPolynomial) -> str:
    if isinstance(polynomial, BallPolynomial):
        pieces = []
        for degree in sorted(polynomial.terms, reverse=True):
            coeff = format_ball(polynomial.terms[degree])
            monomial = "" if degree == 0 else ("*x" if degree == 1 else f"*x^{degree}")
            pieces.append(f"({coeff}){monomial}")
        return " + ".join(pieces) if pieces else "0"
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
    scale = scalar_simplify(
        scalar_divide(actual_terms.get(p0, sp.Integer(0)), expected_terms[p0])
    )

    for pauli_string in set(actual_terms) | set(expected_terms):
        actual_coeff = actual_terms.get(pauli_string, sp.Integer(0))
        expected_coeff = expected_terms.get(pauli_string, sp.Integer(0))
        if not scalar_close(actual_coeff, scalar_multiply(scale, expected_coeff), tol=tol):
            return False, None

    return True, scale


def block_encoding_proportionality_check(
    actual: OpExpr,
    target: OpExpr,
    *,
    epsilon: float,
) -> BlockEncodingProportionalityCheck:
    if epsilon <= 0:
        raise ValueError("block_encoding_epsilon must be positive")

    with ctx.workdps(DEFAULT_BLOCK_ENCODING_CHECK_PRECISION):
        return _block_encoding_proportionality_check_balls(
            actual,
            target,
            epsilon=epsilon,
        )


def _block_encoding_proportionality_check_balls(
    actual: OpExpr,
    target: OpExpr,
    *,
    epsilon: float,
) -> BlockEncodingProportionalityCheck:
    actual, target = _align_operator_qubits(actual, target)
    actual_terms = _nonzero_terms(actual)
    target_terms = _nonzero_terms(target)
    pauli_strings = set(actual_terms) | set(target_terms)
    threshold_ball = to_arb(epsilon) / arb(2**actual.num_qubits).sqrt()
    if not pauli_strings:
        return _make_block_encoding_check(
            epsilon=epsilon,
            alpha_ball=arb(0),
            projection_length_ball=arb(0),
            residual_norm_upper=arb.pos_inf(),
            threshold_ball=threshold_ball,
            coefficient_norm_ball=arb(0),
            target_norm_ball=arb(0),
        )

    actual_coeffs: dict[tuple[str, ...], acb] = {}
    target_coeffs: dict[tuple[str, ...], acb] = {}
    for pauli_string in pauli_strings:
        actual_coeffs[pauli_string] = to_acb(actual_terms.get(pauli_string, 0))
        target_coeffs[pauli_string] = to_acb(target_terms.get(pauli_string, 0))

    actual_norm_sq = sum((_abs_square_ball(coeff) for coeff in actual_coeffs.values()), arb(0))
    target_norm_sq = sum((_abs_square_ball(coeff) for coeff in target_coeffs.values()), arb(0))
    if actual_norm_sq.contains(0) or target_norm_sq.contains(0):
        return _make_block_encoding_check(
            epsilon=epsilon,
            alpha_ball=arb(0),
            projection_length_ball=arb(0),
            residual_norm_upper=arb.pos_inf(),
            threshold_ball=threshold_ball,
            coefficient_norm_ball=actual_norm_sq.sqrt(),
            target_norm_ball=target_norm_sq.sqrt(),
        )

    inner = sum(
        (actual_coeffs[key].conjugate() * target_coeffs[key] for key in pauli_strings),
        acb(0),
    )
    actual_norm = actual_norm_sq.sqrt()
    projection_length_ball = inner.real / actual_norm
    alpha_ball = projection_length_ball / actual_norm
    if alpha_ball.lower() <= 0:
        residual_norm_upper = arb.pos_inf()
    else:
        residual_norm_sq_upper = sum(
            (
                _square_upper(alpha_ball * actual_coeffs[key] - target_coeffs[key])
                for key in pauli_strings
            ),
            arb(0),
        )
        residual_norm_upper = residual_norm_sq_upper.sqrt().abs_upper()

    return _make_block_encoding_check(
        epsilon=epsilon,
        alpha_ball=alpha_ball,
        projection_length_ball=projection_length_ball,
        residual_norm_upper=residual_norm_upper,
        threshold_ball=threshold_ball,
        coefficient_norm_ball=actual_norm,
        target_norm_ball=target_norm_sq.sqrt(),
    )


def _abs_square_ball(value: acb) -> arb:
    return (value.real * value.real + value.imag * value.imag).nonnegative_part()


def _square_upper(value: acb) -> arb:
    upper = value.abs_upper()
    return upper * upper


def _make_block_encoding_check(
    *,
    epsilon: float,
    alpha_ball: arb,
    projection_length_ball: arb,
    residual_norm_upper: arb,
    threshold_ball: arb,
    coefficient_norm_ball: arb,
    target_norm_ball: arb,
) -> BlockEncodingProportionalityCheck:
    return BlockEncodingProportionalityCheck(
        epsilon=epsilon,
        alpha=float(alpha_ball.mid()),
        projection_length=float(projection_length_ball.mid()),
        residual_norm=_arb_upper_float(residual_norm_upper),
        threshold=float(threshold_ball.mid()),
        coefficient_norm=float(coefficient_norm_ball.mid()),
        target_norm=float(target_norm_ball.mid()),
        alpha_ball=alpha_ball,
        residual_norm_upper=residual_norm_upper,
        threshold_ball=threshold_ball,
    )


def _arb_upper_float(value: arb) -> float:
    """Convert an Arb upper bound to binary64 without rounding it downward."""

    if value.is_zero():
        return 0.0
    converted = float(value)
    if not math.isfinite(converted):
        return converted
    return math.nextafter(converted, math.inf)


def _align_operator_qubits(actual: OpExpr, expected: OpExpr) -> Tuple[OpExpr, OpExpr]:
    num_qubits = max(actual.num_qubits, expected.num_qubits)
    return actual.with_num_qubits(num_qubits), expected.with_num_qubits(num_qubits)


def _nonzero_terms(expr: OpExpr) -> dict:
    return {
        pauli_string: coeff
        for pauli_string, coeff in expr.terms.items()
        if not is_exact_zero(scalar_simplify(coeff))
    }


def _format_scalar(value: Optional[sp.Expr]) -> str:
    if value is None:
        return "None"
    if is_ball(value):
        return format_ball(value)
    return sp.sstr(scalar_simplify(value)).replace("I", "i")


def _format_complex(value: complex) -> str:
    expr = sp.N(value.real, 16) + sp.N(value.imag, 16) * sp.I
    return _format_scalar(expr)


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
    parser = argparse.ArgumentParser(
        description="Symbolically verify one route: block encoding, Hamiltonian simulation, or general QSP."
    )
    parser.add_argument("qasm_path", type=Path)
    parser.add_argument("--expected", help='Expected top-left block, for example "(X0 X1 + Z0 Z1)/2".')
    parser.add_argument("--base", help='Base block for QSP/QSVT polynomial checks, for example "(X + Z)/2".')
    parser.add_argument("--expected-polynomial", help='Expected QSP polynomial in x, for example "4*x^3 - 3*x".')
    parser.add_argument(
        "--target-exp-tau",
        type=float,
        dest="target_exp_tau",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--target-exp-epsilon",
        type=float,
        dest="target_exp_epsilon",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--target-exp-scale",
        default="1",
        dest="target_exp_scale",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--hamsim-tau",
        type=float,
        dest="target_exp_tau",
        default=argparse.SUPPRESS,
        help="Hamiltonian-simulation route: simulation time t.",
    )
    parser.add_argument(
        "--hamsim-epsilon",
        type=float,
        dest="target_exp_epsilon",
        default=argparse.SUPPRESS,
        help="Hamiltonian-simulation route: approximation tolerance.",
    )
    parser.add_argument(
        "--hamsim-scale",
        dest="target_exp_scale",
        default=argparse.SUPPRESS,
        help="Deprecated Hamiltonian-simulation route scale; the approximation check now infers beta from the grid.",
    )
    parser.add_argument(
        "--max-approx-grid-points",
        type=int,
        default=DEFAULT_MAX_APPROXIMATION_GRID_POINTS,
        help=(
            "Fail fast if the Hamiltonian-simulation approximation grid would exceed this many points. "
            "Use 0 to disable."
        ),
    )
    parser.add_argument(
        "--block-encoding-epsilon",
        type=float,
        help=(
            "Block-encoding route: accept when the final all-zero branch is proportional "
            "to the target Hamiltonian within this operator-norm tolerance."
        ),
    )
    parser.add_argument("--hermitian-base", action="store_true", help="Rewrite Hd atoms to H before QSP polynomial evaluation.")
    parser.add_argument(
        "--compare-polynomial-only",
        action="store_true",
        help="Compare QSP polynomials directly instead of evaluating them on --base.",
    )
    parser.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE, help="Numeric tolerance for coefficient comparisons.")
    parser.add_argument(
        "--qsp-calculation-precision",
        type=int,
        default=DEFAULT_QSP_CALCULATION_PRECISION,
        help=(
            "Minimum decimal digits used to parse numeric angles, propagate ball coefficients, "
            "and evaluate the approximation grid. Large block encodings automatically use at least "
            f"{BALL_PRECISION_DIGITS_PER_SYSTEM_QUBIT} digits per system qubit. "
            f"Default minimum: {DEFAULT_QSP_CALCULATION_PRECISION}."
        ),
    )
    parser.add_argument("--ancilla", default="0", help="Single ancilla qubit index. Ignored when --ancillas is set.")
    parser.add_argument("--ancillas", nargs="+", help='Ancilla qubits in branch-key order, e.g. q[0] q[1].')
    parser.add_argument("--system", default="1", help="Single system qubit index. Ignored when --systems is set.")
    parser.add_argument("--systems", nargs="+", help='System qubits in operator-index order, e.g. q[2] q[3].')
    parser.add_argument("--trace", action="store_true", help="Print every intermediate B0/B1 update.")
    parser.add_argument("--result-only", action="store_true", help="Print only PASS/FAIL status.")
    parser.add_argument("--profile-gates", action="store_true", help="Collect per-gate branch/term/timing counters.")
    parser.add_argument("--profile-output", type=Path, help="Write --profile-gates CSV output to this path instead of stdout.")
    parser.add_argument(
        "--defer-qsp-rewrite",
        action="store_true",
        help="Disable eager word rewriting after QSP signal steps (useful for inspecting unreduced traces).",
    )
    args = parser.parse_args(argv)

    ancillas = tuple(_parse_qubit_arg(value) for value in args.ancillas) if args.ancillas is not None else None
    systems = tuple(_parse_qubit_arg(value) for value in args.systems) if args.systems is not None else None
    max_approx_grid_points = None if args.max_approx_grid_points == 0 else args.max_approx_grid_points

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
        target_exp_tau=args.target_exp_tau,
        target_exp_epsilon=args.target_exp_epsilon,
        target_exp_scale=args.target_exp_scale,
        max_approx_grid_points=max_approx_grid_points,
        block_encoding_epsilon=args.block_encoding_epsilon,
        tolerance=args.tolerance,
        profile_gates=args.profile_gates,
        rewrite_after_each_step=not args.defer_qsp_rewrite,
        qsp_calculation_precision=args.qsp_calculation_precision,
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
