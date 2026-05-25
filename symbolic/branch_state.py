from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import sympy as sp

from .expr import OpExpr, identity, pauli, zero
from .scalar import cos_half, exp_minus_i_half, exp_plus_i_half, sin_half


class UnsupportedGateError(ValueError):
    pass


@dataclass(frozen=True)
class Gate:
    name: str
    qubits: Tuple[int, ...]
    parameter: Optional[sp.Expr] = None
    raw: str = ""

    def __str__(self) -> str:
        operands = ", ".join(f"q[{index}]" for index in self.qubits)
        if self.parameter is None:
            return f"{self.name} {operands}"
        return f"{self.name}({sp.sstr(self.parameter)}) {operands}"


@dataclass(frozen=True)
class BranchState:
    b0: OpExpr
    b1: OpExpr

    @staticmethod
    def initial(num_system_qubits: int = 1) -> "BranchState":
        return BranchState(identity(num_system_qubits), zero(num_system_qubits))

    def apply(
        self,
        gate: Gate,
        *,
        ancilla: int = 0,
        system: int = 1,
        systems: Optional[Sequence[int]] = None,
    ) -> "BranchState":
        name = gate.name.lower()
        if gate.parameter is not None and name not in {"rx", "ry", "rz"}:
            raise UnsupportedGateError(f"Parameterized gate is not supported yet: {gate}")
        if gate.parameter is None and name in {"rx", "ry", "rz"}:
            raise UnsupportedGateError(f"Rotation gate requires a parameter: {gate}")

        system_qubits = _normalize_systems(system=system, systems=systems)
        if len(gate.qubits) == 1:
            return self._apply_single_qubit(
                name,
                gate.qubits[0],
                parameter=gate.parameter,
                ancilla=ancilla,
                systems=system_qubits,
            )
        if len(gate.qubits) == 2:
            return self._apply_two_qubit(name, gate.qubits, ancilla=ancilla, systems=system_qubits)
        raise UnsupportedGateError(f"Unsupported gate arity: {gate}")

    def _apply_single_qubit(
        self,
        name: str,
        target: int,
        *,
        parameter: Optional[sp.Expr],
        ancilla: int,
        systems: Tuple[int, ...],
    ) -> "BranchState":
        if target == ancilla:
            if name == "x":
                return BranchState(self.b1, self.b0)
            if name == "z":
                return BranchState(self.b0, -self.b1)
            if name == "s":
                return BranchState(self.b0, self.b1.scale(sp.I))
            if name == "sdg":
                return BranchState(self.b0, self.b1.scale(-sp.I))
            if name == "h":
                inv_sqrt2 = sp.sqrt(2) / 2
                return BranchState((self.b0 + self.b1).scale(inv_sqrt2), (self.b0 - self.b1).scale(inv_sqrt2))
            if name == "rx":
                c = cos_half(parameter)
                minus_i_s = -sp.I * sin_half(parameter)
                return BranchState(
                    self.b0.scale(c) + self.b1.scale(minus_i_s),
                    self.b0.scale(minus_i_s) + self.b1.scale(c),
                )
            if name == "ry":
                c = cos_half(parameter)
                s = sin_half(parameter)
                return BranchState(
                    self.b0.scale(c) - self.b1.scale(s),
                    self.b0.scale(s) + self.b1.scale(c),
                )
            if name == "rz":
                return BranchState(
                    self.b0.scale(exp_minus_i_half(parameter)),
                    self.b1.scale(exp_plus_i_half(parameter)),
                )
            raise UnsupportedGateError(f"Unsupported ancilla gate: {name}")

        system_index = _system_index(target, systems)
        if system_index is not None:
            num_system_qubits = len(systems)
            if name == "x":
                gate_op = pauli("X", index=system_index, num_qubits=num_system_qubits)
            elif name == "z":
                gate_op = pauli("Z", index=system_index, num_qubits=num_system_qubits)
            elif name == "s":
                gate_op = _s_operator(system_index=system_index, num_system_qubits=num_system_qubits)
            elif name == "sdg":
                gate_op = _sdg_operator(system_index=system_index, num_system_qubits=num_system_qubits)
            elif name == "h":
                gate_op = (
                    pauli("X", index=system_index, num_qubits=num_system_qubits)
                    + pauli("Z", index=system_index, num_qubits=num_system_qubits)
                ).scale(sp.sqrt(2) / 2)
            elif name in {"rx", "ry", "rz"}:
                gate_op = _rotation_operator(name, parameter, system_index=system_index, num_system_qubits=num_system_qubits)
            else:
                raise UnsupportedGateError(f"Unsupported system gate: {name}")
            return BranchState(gate_op * self.b0, gate_op * self.b1)

        raise UnsupportedGateError(f"Gate target q[{target}] is neither ancilla q[{ancilla}] nor a system qubit {systems}")

    def _apply_two_qubit(self, name: str, qubits: Tuple[int, int], *, ancilla: int, systems: Tuple[int, ...]) -> "BranchState":
        left, right = qubits

        if name == "cx":
            system_index = _system_index(right, systems)
            if left == ancilla and system_index is not None:
                return BranchState(self.b0, pauli("X", index=system_index, num_qubits=len(systems)) * self.b1)
            raise UnsupportedGateError("Only cx with ancilla as control and a system qubit as target is supported")

        if name == "cz":
            other = right if left == ancilla else left if right == ancilla else None
            system_index = _system_index(other, systems) if other is not None else None
            if system_index is not None:
                return BranchState(self.b0, pauli("Z", index=system_index, num_qubits=len(systems)) * self.b1)
            raise UnsupportedGateError("Only cz between the ancilla and a system qubit is supported")

        raise UnsupportedGateError(f"Unsupported two-qubit gate: {name}")


def _normalize_systems(*, system: int, systems: Optional[Sequence[int]]) -> Tuple[int, ...]:
    if systems is None:
        return (system,)
    normalized = tuple(systems)
    if not normalized:
        raise ValueError("At least one system qubit is required")
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"Duplicate system qubits are not allowed: {normalized}")
    return normalized


def _system_index(qasm_qubit: Optional[int], systems: Tuple[int, ...]) -> Optional[int]:
    if qasm_qubit is None:
        return None
    try:
        return systems.index(qasm_qubit)
    except ValueError:
        return None


def _rotation_operator(
    name: str,
    theta: sp.Expr,
    *,
    system_index: int,
    num_system_qubits: int,
) -> OpExpr:
    c = cos_half(theta)
    minus_i_s = -sp.I * sin_half(theta)
    axis = {"rx": "X", "ry": "Y", "rz": "Z"}[name]
    identity_op = identity(num_system_qubits)
    axis_op = pauli(axis, index=system_index, num_qubits=num_system_qubits)
    return identity_op.scale(c) + axis_op.scale(minus_i_s)


def _s_operator(*, system_index: int, num_system_qubits: int) -> OpExpr:
    identity_op = identity(num_system_qubits)
    z_op = pauli("Z", index=system_index, num_qubits=num_system_qubits)
    return identity_op.scale((1 + sp.I) / 2) + z_op.scale((1 - sp.I) / 2)


def _sdg_operator(*, system_index: int, num_system_qubits: int) -> OpExpr:
    identity_op = identity(num_system_qubits)
    z_op = pauli("Z", index=system_index, num_qubits=num_system_qubits)
    return identity_op.scale((1 - sp.I) / 2) + z_op.scale((1 + sp.I) / 2)
