"""Opt-in exact summaries for native PennyLane Clifford+T subcircuits.

This module does not change the core parser, branch semantics, or rewrite rules.
Only complete, parameter-free templates are summarized. Decimal rotations are
retained, except rotations whose specified angle is exactly zero.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time
from typing import Sequence

import sympy as sp

from symbolic.branch_state import BranchState, Gate
from symbolic.expr import OpExpr
from symbolic.interval import ctx, is_exact_zero
from symbolic.qasm_parser import parse_qasm_file
from symbolic.verify import VerificationResult


Operation = tuple[str, tuple[int, ...]]


@dataclass(frozen=True)
class Template:
    name: str
    primitive: str
    operations: tuple[Operation, ...]


# The exact native decompositions in PennyLane 0.45.1. The role numbers are
# placeholders, not physical wire numbers. Both controls may be permuted.
PENNYLANE_CCX = (
    ("h", (2,)), ("cx", (1, 2)), ("tdg", (2,)), ("cx", (0, 2)),
    ("t", (2,)), ("cx", (1, 2)), ("tdg", (2,)), ("cx", (0, 2)),
    ("t", (2,)), ("t", (1,)), ("cx", (0, 1)), ("h", (2,)),
    ("t", (0,)), ("tdg", (1,)), ("cx", (0, 1)),
)
PENNYLANE_CCZ = PENNYLANE_CCX[1:] + (("h", (2,)),)
TEMPLATES = (
    Template("pennylane_ccx", "ccx", PENNYLANE_CCX),
    Template("pennylane_ccz", "ccz", PENNYLANE_CCZ),
)


@lru_cache(maxsize=None)
def certify_template(template: Template) -> dict:
    """Prove equality on every input column with exact algebraic amplitudes."""
    if template.primitive not in {"ccx", "ccz"}:
        raise ValueError("Unknown summary primitive")
    root_half = sp.sqrt(2) / 2
    phases = {"t": (1+sp.I)*root_half, "tdg": (1-sp.I)*root_half}
    for basis in range(8):
        initial = tuple((basis >> (2-i)) & 1 for i in range(3))
        amplitudes = {initial: sp.Integer(1)}
        for name, wires in template.operations:
            updated = {}

            def add(bits, amplitude):
                updated[bits] = sp.expand(updated.get(bits, 0)+amplitude)

            for bits, amplitude in amplitudes.items():
                if name == "cx":
                    control, target = wires
                    output = list(bits)
                    output[target] ^= bits[control]
                    add(tuple(output), amplitude)
                elif name == "h":
                    target = wires[0]
                    for output_bit in (0, 1):
                        output = list(bits)
                        output[target] = output_bit
                        sign = -1 if bits[target] and output_bit else 1
                        add(tuple(output), amplitude*root_half*sign)
                elif name in phases:
                    add(bits, amplitude*(phases[name] if bits[wires[0]] else 1))
                else:
                    raise ValueError(f"Unsupported exact template gate: {name}")
            amplitudes = {bits: amplitude for bits, amplitude in updated.items() if amplitude != 0}
        expected_bits = list(initial)
        expected_phase = sp.Integer(1)
        if initial[0] and initial[1]:
            if template.primitive == "ccx":
                expected_bits[2] ^= 1
            elif initial[2]:
                expected_phase = sp.Integer(-1)
        expected = {tuple(expected_bits): expected_phase}
        for bits in set(amplitudes) | set(expected):
            if sp.simplify(amplitudes.get(bits, 0)-expected.get(bits, 0)) != 0:
                raise ValueError(f"Exact equivalence proof failed: {template.name}, input={initial}")
    normalized = json.dumps(template.operations, separators=(",", ":"))
    return dict(status="PASS_EXACT", primitive=template.primitive, input_columns=8,
                absolute_phase_checked=True, coefficient_domain="Q(sqrt(2), i)",
                template_sha256=hashlib.sha256(normalized.encode()).hexdigest(),
                method="Exact computational-basis evaluation; no numerical tolerance")


@dataclass(frozen=True)
class SummaryGate:
    name: str
    qubits: tuple[int, int, int]
    template: str
    source_gate_start: int
    source_gate_stop: int


@dataclass(frozen=True)
class ExecutionPlan:
    steps: tuple[Gate | SummaryGate, ...]
    report: dict


def _match(window: Sequence[Gate], template: Template) -> tuple[int, int, int] | None:
    if len(window) != len(template.operations):
        return None
    mapping = {}
    for actual, (name, roles) in zip(window, template.operations):
        if actual.name != name or actual.parameter is not None or len(actual.qubits) != len(roles):
            return None
        for role, wire in zip(roles, actual.qubits):
            if role in mapping and mapping[role] != wire:
                return None
            mapping[role] = wire
    if set(mapping) != {0, 1, 2} or len(set(mapping.values())) != 3:
        return None
    return tuple(mapping[role] for role in range(3))


def compile_summaries(gates: Sequence[Gate], *, ancillas: Sequence[int],
                      systems: Sequence[int]) -> ExecutionPlan:
    """Match complete templates after removing only exact identity rotations."""
    started = time.perf_counter()
    ancilla_set, system_set = set(ancillas), set(systems)
    stream, removed = [], []
    for index, gate in enumerate(gates):
        if (gate.name in {"rx", "ry", "rz"} and len(gate.qubits) == 1
                and gate.qubits[0] in ancilla_set | system_set
                and gate.parameter is not None and is_exact_zero(gate.parameter)):
            removed.append(index+1)
        else:
            stream.append((index, gate))
    steps, occurrences, proofs = [], [], {}
    cursor = 0
    while cursor < len(stream):
        for template in TEMPLATES:
            width = len(template.operations)
            section = stream[cursor:cursor+width]
            wires = _match([gate for _, gate in section], template)
            if wires is None or not set(wires[:2]) <= ancilla_set or wires[2] not in ancilla_set | system_set:
                continue
            # Matching alone is never sufficient: prove the template before use.
            proofs[template.name] = certify_template(template)
            first, last = section[0][0]+1, section[-1][0]+1
            steps.append(SummaryGate(template.primitive, wires, template.name, first, last))
            occurrences.append(dict(template=template.name, primitive=template.primitive,
                                    qubits=list(wires), source_gate_start=first, source_gate_stop=last,
                                    parameter_free_gates=width))
            cursor += width
            break
        else:
            steps.append(stream[cursor][1])
            cursor += 1
    counts = Counter(item["primitive"] for item in occurrences)
    return ExecutionPlan(tuple(steps), dict(
        mode="exact_local_summaries", source_gate_count=len(gates), executed_steps=len(steps),
        summary_count=len(occurrences), summary_counts=dict(counts),
        summarized_source_gates=sum(item["parameter_free_gates"] for item in occurrences),
        removed_identity_gates=len(removed), removed_identity_source_gate_indices=removed,
        certificates=proofs, occurrences=occurrences,
        summary_compilation_sec=time.perf_counter()-started,
        decimal_rotations="Retained exactly as parsed; only exact-zero rotations removed",
        approximation_error="0", input_qasm_modified=False))


def apply_summary(state: BranchState, gate: SummaryGate, *, ancillas: Sequence[int],
                  systems: Sequence[int]) -> BranchState:
    """Apply the certified full unitary, including when a work ancilla is dirty."""
    if state.expression_kind != "op" or gate.name not in {"ccx", "ccz"}:
        raise ValueError("Local summaries require concrete Pauli execution")
    controls = tuple(ancillas.index(wire) for wire in gate.qubits[:2])
    target = gate.qubits[2]
    if target in ancillas:
        if gate.name == "ccx":
            return state.apply(Gate("mcx", gate.qubits), ancillas=ancillas, systems=systems)
        target_index = ancillas.index(target)
        branches = {key: -value if all(key[i] for i in controls) and key[target_index] else value
                    for key, value in state.branches.items()}
    elif target in systems:
        operator = OpExpr.pauli("X" if gate.name == "ccx" else "Z", index=systems.index(target),
                                num_qubits=state.num_system_qubits)
        branches = {key: operator*value if all(key[i] for i in controls) else value
                    for key, value in state.branches.items()}
    else:
        raise ValueError("Unknown summary target")
    return BranchState(state.num_ancillas, state.num_system_qubits, branches, expression_kind="op")


def execute_plan(plan: ExecutionPlan, *, ancillas: Sequence[int], systems: Sequence[int]) -> VerificationResult:
    state = BranchState.initial(num_system_qubits=len(systems), num_ancillas=len(ancillas))
    for step in plan.steps:
        if isinstance(step, SummaryGate):
            state = apply_summary(state, step, ancillas=ancillas, systems=systems)
        else:
            state = state.apply(step, ancillas=ancillas, systems=systems)
    return VerificationResult(final_state=state, trace=(), ancilla_qubits=tuple(ancillas),
                              system_qubits=tuple(systems))


def verify_with_summaries(path: str | Path, *, ancillas: Sequence[int], systems: Sequence[int],
                          precision: int) -> tuple[VerificationResult, dict]:
    precision = max(precision, 2*len(systems))
    with ctx.workdps(precision):
        gates = parse_qasm_file(path, numeric_precision=precision)
        plan = compile_summaries(gates, ancillas=ancillas, systems=systems)
        started = time.perf_counter()
        result = execute_plan(plan, ancillas=ancillas, systems=systems)
        report = dict(plan.report, summary_execution_sec=time.perf_counter()-started,
                      final_pauli_terms=len(result.final_state.top_left().terms))
    return result, report
