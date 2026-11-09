"""PennyLane-only sparse propagation with a certified truncation error budget.

For a branch column V = sum_a |a> B_a, discarding terms D changes V in
operator norm by at most sum_(a,P in D) |c_(a,P)|. Subsequent unitary gates
preserve this error, and the final ancilla projection cannot increase it.
The accumulated bound eta therefore certifies ||B - B_tilde|| <= eta for
the ORIGINAL decimal QASM, not just for the sparsified computation.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
import time

from evaluation.pennylane.local_summaries import SummaryGate, apply_summary, compile_summaries
from symbolic.branch_state import BranchState
from symbolic.expr import OpExpr
from symbolic.interval import arb, ctx, is_ball, to_acb
from symbolic.qasm_parser import parse_qasm_file
from symbolic.verify import VerificationResult


DEFAULT_THRESHOLD = "1e-24"
DEFAULT_BUDGET = "1e-14"


@dataclass
class ErrorLedger:
    threshold: str = DEFAULT_THRESHOLD
    budget: str = DEFAULT_BUDGET
    error: arb = field(default_factory=lambda: arb(0))
    dropped_terms: int = 0
    skipped_batches: int = 0
    recentered_coefficients: int = 0
    roundoff_error: arb = field(default_factory=lambda: arb(0))
    pruning_error: arb = field(default_factory=lambda: arb(0))
    skipped_recentering_batches: int = 0
    events: list = field(default_factory=list)

    def __post_init__(self):
        threshold, budget = Fraction(self.threshold), Fraction(self.budget)
        if threshold <= 0 or budget < 0:
            raise ValueError("Pruning threshold must be positive; budget must be nonnegative")
        # Lower bounds make the eligibility/budget comparisons conservative.
        self.threshold_ball = arb(str(threshold))
        self.budget_ball = arb(str(budget))

    def prune(self, state: BranchState, step: int) -> BranchState:
        if state.expression_kind != "op":
            raise ValueError("Bounded pruning requires concrete Pauli branches")
        if self.budget_ball.is_zero():
            return state
        discard = {}
        batch_error = arb(0)
        count = 0
        for key, operator in state.branches.items():
            words = []
            for word, coefficient in operator.terms.items():
                bound = to_acb(coefficient).abs_upper()
                if bound.is_finite() and bound <= self.threshold_ball.lower():
                    words.append(word)
                    batch_error = (batch_error+bound).upper()
            if words:
                discard[key] = set(words)
                count += len(words)
        if not count:
            return state
        new_error = (self.error+batch_error).upper()
        if not new_error <= self.budget_ball.lower():
            # Preserve every term if this batch would exhaust the budget.
            self.skipped_batches += 1
            return state
        branches = {
            key: (OpExpr({word: c for word, c in operator.terms.items() if word not in discard[key]},
                         num_qubits=state.num_system_qubits) if key in discard else operator)
            for key, operator in state.branches.items()
        }
        self.error = new_error
        self.pruning_error = (self.pruning_error+batch_error).upper()
        self.dropped_terms += count
        self.events.append(dict(step=step, kind="prune", dropped_terms=count, batch_error_upper=str(batch_error),
                                cumulative_error_upper=str(self.error)))
        return BranchState(state.num_ancillas, state.num_system_qubits, branches)

    def recenter(self, state: BranchState, step: int) -> BranchState:
        """Charge ball radii before replacing them with exact dyadic midpoints.

        This transports numerical uncertainty through unitary suffixes in
        operator norm instead of repeatedly amplifying uncorrelated boxes.
        It is never an unaccounted conversion of an interval to its midpoint.
        """
        if self.budget_ball.is_zero():
            return state
        batch_error, count = arb(0), 0
        replacements = {}
        for key, operator in state.branches.items():
            coefficients = {}
            for word, coefficient in operator.terms.items():
                if not is_ball(coefficient):
                    continue
                ball = to_acb(coefficient)
                if not ball.real.is_finite() or not ball.imag.is_finite():
                    raise ValueError("Nonfinite coefficient; increase precision")
                radius = (ball.real.rad()+ball.imag.rad()).upper()
                if not radius.is_finite():
                    raise ValueError("Nonfinite coefficient uncertainty; increase precision")
                if radius.is_zero():
                    continue
                coefficients[word] = ball.mid()
                batch_error = (batch_error+radius).upper()
                count += 1
            if coefficients:
                replacements[key] = coefficients
        if not count:
            return state
        new_error = (self.error+batch_error).upper()
        if not new_error <= self.budget_ball.lower():
            self.skipped_recentering_batches += 1
            return state
        branches = dict(state.branches)
        for key, coefficients in replacements.items():
            branches[key] = OpExpr({**state.branches[key].terms, **coefficients},
                                   num_qubits=state.num_system_qubits)
        self.error = new_error
        self.roundoff_error = (self.roundoff_error+batch_error).upper()
        self.recentered_coefficients += count
        self.events.append(dict(step=step, kind="recenter", recentered_coefficients=count,
            batch_error_upper=str(batch_error), cumulative_error_upper=str(self.error)))
        return BranchState(state.num_ancillas, state.num_system_qubits, branches)

    def report(self) -> dict:
        return dict(threshold=self.threshold, budget=self.budget,
            operator_error_upper=str(self.error), dropped_terms=self.dropped_terms,
            pruning_error_upper=str(self.pruning_error), roundoff_error_upper=str(self.roundoff_error),
            recentered_coefficients=self.recentered_coefficients,
            skipped_recentering_batches=self.skipped_recentering_batches,
            skipped_batches=self.skipped_batches, events=self.events,
            norm="Operator norm on the complete ancilla branch column, then projected block",
            rule="Charge sum abs_upper(c) for discarded terms and sum (Re radius + Im radius) for recentered terms, over ALL branches",
            propagation="Unitary suffixes preserve operator norm; final projection and restored global phase do not increase it",
            budget_policy="Keep the entire batch if its upper bound would exceed the budget",
            interval_roundoff="Recenter to exact dyadic midpoints only after charging the complete ball radii; preserve balls if the shared budget would be exceeded",
            original_qasm_verified=True)


def _check_unitary_steps(plan):
    """Do not apply a unitary-error proof to opaque or nonunitary operations."""
    allowed = {"h", "x", "y", "z", "s", "sdg", "t", "tdg", "cx", "cy", "cz", "rx", "ry", "rz"}
    for step in plan.steps:
        if isinstance(step, SummaryGate):
            continue  # compile_summaries requires an exact full-unitary certificate.
        if step.name not in allowed:
            raise ValueError(f"No unitary error-propagation proof for {step.name}")
        expected_arity = 2 if step.name in {"cx", "cy", "cz"} else 1
        if len(step.qubits) != expected_arity or len(set(step.qubits)) != expected_arity:
            raise ValueError("Unitary propagation requires valid arity and distinct gate operands")
        if (step.name in {"rx", "ry", "rz"}) != (step.parameter is not None):
            raise ValueError("Only RX/RY/RZ may have a rotation parameter, and require one")
        if step.parameter is not None:
            angle = to_acb(step.parameter)
            if not angle.imag.is_zero() or not angle.real.is_finite():
                raise ValueError("Unitary error propagation requires finite real rotation angles")


def verify_with_pruning(path: str | Path, *, ancillas, systems, precision: int,
                        threshold: str = DEFAULT_THRESHOLD, budget: str = DEFAULT_BUDGET):
    precision = max(precision, 2*len(systems))
    if len(set(ancillas) | set(systems)) != len(ancillas)+len(systems):
        raise ValueError("Unitary propagation requires disjoint registers without duplicate wires")
    with ctx.workdps(precision):
        gates = parse_qasm_file(path, numeric_precision=precision)
        plan = compile_summaries(gates, ancillas=ancillas, systems=systems)
        _check_unitary_steps(plan)
        ledger = ErrorLedger(threshold, budget)
        state = BranchState.initial(num_ancillas=len(ancillas), num_system_qubits=len(systems))
        peak_total = peak_per_branch = 1
        started = time.perf_counter()
        for index, step in enumerate(plan.steps, 1):
            if isinstance(step, SummaryGate):
                state = apply_summary(state, step, ancillas=ancillas, systems=systems)
            else:
                state = state.apply(step, ancillas=ancillas, systems=systems)
            sizes = [len(operator.terms) for operator in state.branches.values()]
            peak_total = max(peak_total, sum(sizes))
            peak_per_branch = max(peak_per_branch, max(sizes, default=0))
            state = ledger.prune(state, index)
            state = ledger.recenter(state, index)
        result = VerificationResult(final_state=state, trace=(), ancilla_qubits=tuple(ancillas),
                                    system_qubits=tuple(systems))
        summary_report = dict(plan.report, summary_execution_sec=time.perf_counter()-started,
                              final_pauli_terms=len(state.top_left().terms))
        pruning_report = dict(ledger.report(), peak_total_branch_terms=peak_total,
                              peak_terms_per_branch=peak_per_branch,
                              final_pauli_terms=len(state.top_left().terms))
    return result, summary_report, pruning_report


def residual_certificate(residual, alpha, dimension_factor, epsilon: str, error):
    """The fitted alpha is only a candidate; certify the original block."""
    threshold = arb(str(Fraction(epsilon)))/dimension_factor
    truncation = (alpha.abs_upper()*error).upper()
    total = (dimension_factor*residual+truncation).upper()
    return dict(passed=bool(alpha.lower() > 0 and total <= arb(str(Fraction(epsilon))).lower()),
                threshold=threshold, operator_error_upper=total,
                scaled_truncation_error_upper=truncation)


def classify_with_pruning(actual, structure_proof: dict, delta: str, error):
    # Import here to avoid changing the old helper or its exact-mode behavior.
    from evaluation.pennylane.run_block_encoding import classify_hermiticity

    if error.is_zero():
        return classify_hermiticity(actual, structure_proof, delta)
    delta_fraction = Fraction(delta)
    if delta_fraction < 0:
        raise ValueError("Hermiticity delta must be nonnegative")
    bounds = {"".join(word): to_acb(c).imag.abs_upper() for word, c in actual.terms.items()}
    worst = max(bounds, key=bounds.get) if bounds else ""
    maximum = (bounds.get(worst, arb(0))+error).upper()
    violations = [dict(pauli="".join(word), approximate_imaginary_coefficient=str(to_acb(c).imag))
                  for word, c in actual.terms.items() if to_acb(c).imag.abs_lower() > error]
    if violations:
        exact = dict(status="FAIL_NON_HERMITIAN", violations=violations,
            proof="Imaginary magnitude lower bound exceeds the complete truncation error bound")
    elif structure_proof["status"] == "PASS_STRUCTURE":
        exact = structure_proof
    else:
        exact = dict(status="UNPROVEN", reason="Pruned coefficients cannot prove exact zero imaginary parts")
    within = bool(delta_fraction > 0 and maximum < arb(str(delta_fraction)).lower())
    status = exact["status"] if exact["status"].startswith("PASS_") else (
        "PASS_WITHIN_DELTA" if within else "FAIL_DELTA" if delta_fraction > 0 else exact["status"])
    return dict(status=status, exact_status=exact["status"], exact=exact, structure=structure_proof,
        delta=delta, within_delta=within, truncation_error_upper=str(error),
        max_imaginary_coefficient_upper=str(maximum), max_imaginary_coefficient_pauli=worst,
        criterion="Exact structure OR max_P |Im(c_tilde_P)| + eta < delta; missing coefficients bounded by eta",
        approximate_acceptance=status == "PASS_WITHIN_DELTA")
