from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Mapping, Optional, Sequence, Tuple

import sympy as sp

from .expr import OpExpr, identity, pauli, zero
from .scalar import cos_half, exp_minus_i_half, exp_plus_i_half, sin_half


BranchKey = Tuple[int, ...]


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
    num_ancillas: int
    num_system_qubits: int
    branches: Mapping[BranchKey, OpExpr]

    def __post_init__(self) -> None:
        if self.num_ancillas < 1:
            raise ValueError("num_ancillas must be at least 1")
        if self.num_system_qubits < 1:
            raise ValueError("num_system_qubits must be at least 1")

        cleaned: Dict[BranchKey, OpExpr] = {}
        for key, branch in self.branches.items():
            normalized_key = tuple(key)
            if len(normalized_key) != self.num_ancillas:
                raise ValueError(
                    f"Branch key {normalized_key!r} has length {len(normalized_key)}, expected {self.num_ancillas}"
                )
            if any(bit not in (0, 1) for bit in normalized_key):
                raise ValueError(f"Branch key must be a bitstring: {normalized_key!r}")
            normalized_branch = branch.with_num_qubits(self.num_system_qubits)
            if not normalized_branch.is_zero():
                cleaned[normalized_key] = normalized_branch
        object.__setattr__(self, "branches", cleaned)

    @staticmethod
    def initial(num_system_qubits: int = 1, num_ancillas: int = 1) -> "BranchState":
        zero_key = (0,) * num_ancillas
        return BranchState(num_ancillas, num_system_qubits, {zero_key: identity(num_system_qubits)})

    @property
    def b0(self) -> OpExpr:
        return self.branch((0,))

    @property
    def b1(self) -> OpExpr:
        return self.branch((1,))

    def top_left(self) -> OpExpr:
        return self.branch((0,) * self.num_ancillas)

    def branch(self, key: BranchKey) -> OpExpr:
        if len(key) != self.num_ancillas:
            if self.num_ancillas == 1 and len(key) == 1:
                pass
            else:
                raise ValueError(f"Branch key {key!r} does not match {self.num_ancillas} ancillas")
        return self.branches.get(tuple(key), zero(self.num_system_qubits))

    def apply(
        self,
        gate: Gate,
        *,
        ancilla: int = 0,
        ancillas: Optional[Sequence[int]] = None,
        system: int = 1,
        systems: Optional[Sequence[int]] = None,
    ) -> "BranchState":
        name = gate.name.lower()
        if gate.parameter is not None and name not in {"rx", "ry", "rz"}:
            raise UnsupportedGateError(f"Parameterized gate is not supported yet: {gate}")
        if gate.parameter is None and name in {"rx", "ry", "rz"}:
            raise UnsupportedGateError(f"Rotation gate requires a parameter: {gate}")

        ancilla_qubits = _normalize_ancillas(ancilla=ancilla, ancillas=ancillas)
        system_qubits = _normalize_systems(system=system, systems=systems)
        if len(ancilla_qubits) != self.num_ancillas:
            raise ValueError(f"State has {self.num_ancillas} ancillas, got mapping {ancilla_qubits}")
        if len(system_qubits) != self.num_system_qubits:
            raise ValueError(f"State has {self.num_system_qubits} system qubits, got mapping {system_qubits}")

        if len(gate.qubits) == 1:
            return self._apply_single_qubit(
                name,
                gate.qubits[0],
                parameter=gate.parameter,
                ancillas=ancilla_qubits,
                systems=system_qubits,
            )
        if len(gate.qubits) == 2:
            return self._apply_two_qubit(name, gate.qubits, ancillas=ancilla_qubits, systems=system_qubits)
        raise UnsupportedGateError(f"Unsupported gate arity: {gate}")

    def _apply_single_qubit(
        self,
        name: str,
        target: int,
        *,
        parameter: Optional[sp.Expr],
        ancillas: Tuple[int, ...],
        systems: Tuple[int, ...],
    ) -> "BranchState":
        ancilla_index = _qubit_index(target, ancillas)
        if ancilla_index is not None:
            matrix = _ancilla_gate_matrix(name, parameter)
            return self._apply_ancilla_matrix(ancilla_index, matrix)

        system_index = _qubit_index(target, systems)
        if system_index is not None:
            gate_op = _system_gate_operator(
                name,
                parameter,
                system_index=system_index,
                num_system_qubits=self.num_system_qubits,
            )
            return self._apply_system_operator(gate_op)

        raise UnsupportedGateError(f"Gate target q[{target}] is neither an ancilla {ancillas} nor a system qubit {systems}")

    def _apply_two_qubit(self, name: str, qubits: Tuple[int, int], *, ancillas: Tuple[int, ...], systems: Tuple[int, ...]) -> "BranchState":
        left, right = qubits
        left_ancilla = _qubit_index(left, ancillas)
        right_ancilla = _qubit_index(right, ancillas)
        left_system = _qubit_index(left, systems)
        right_system = _qubit_index(right, systems)

        if name == "cx":
            if left_ancilla is not None and right_system is not None:
                gate_op = pauli("X", index=right_system, num_qubits=self.num_system_qubits)
                return self._apply_controlled_pauli_ancilla_to_system(left_ancilla, gate_op)
            if left_system is not None and right_ancilla is not None:
                return self._apply_cx_system_to_ancilla(left_system, right_ancilla)
            if left_ancilla is not None and right_ancilla is not None:
                return self._apply_cx_ancilla_to_ancilla(left_ancilla, right_ancilla)
            if left_system is not None and right_system is not None:
                raise UnsupportedGateError("system-system cx is not supported in v0.4")
            raise UnsupportedGateError("cx must act on known ancilla/system qubits")

        if name == "cz":
            if left_ancilla is not None and right_system is not None:
                gate_op = pauli("Z", index=right_system, num_qubits=self.num_system_qubits)
                return self._apply_controlled_pauli_ancilla_to_system(left_ancilla, gate_op)
            if right_ancilla is not None and left_system is not None:
                gate_op = pauli("Z", index=left_system, num_qubits=self.num_system_qubits)
                return self._apply_controlled_pauli_ancilla_to_system(right_ancilla, gate_op)
            if left_ancilla is not None and right_ancilla is not None:
                return self._apply_cz_ancilla_to_ancilla(left_ancilla, right_ancilla)
            if left_system is not None and right_system is not None:
                raise UnsupportedGateError("system-system cz is not supported in v0.4")
            raise UnsupportedGateError("cz must act on known ancilla/system qubits")

        raise UnsupportedGateError(f"Unsupported two-qubit gate: {name}")

    def _apply_ancilla_matrix(self, ancilla_index: int, matrix: Tuple[Tuple[sp.Expr, sp.Expr], Tuple[sp.Expr, sp.Expr]]) -> "BranchState":
        new: Dict[BranchKey, OpExpr] = {}
        processed = set()
        relevant_keys = set(self.branches)
        relevant_keys.update(_set_bit(key, ancilla_index, 1 - key[ancilla_index]) for key in self.branches)

        for key in sorted(relevant_keys):
            key0 = _set_bit(key, ancilla_index, 0)
            key1 = _set_bit(key, ancilla_index, 1)
            if key0 in processed:
                continue

            branch0 = self.branch(key0)
            branch1 = self.branch(key1)
            new[key0] = branch0.scale(matrix[0][0]) + branch1.scale(matrix[0][1])
            new[key1] = branch0.scale(matrix[1][0]) + branch1.scale(matrix[1][1])
            processed.add(key0)
            processed.add(key1)

        return self._replace(new)

    def _apply_system_operator(self, operator: OpExpr) -> "BranchState":
        return self._replace({key: operator * branch for key, branch in self.branches.items()})

    def _apply_controlled_pauli_ancilla_to_system(self, control_index: int, operator: OpExpr) -> "BranchState":
        new: Dict[BranchKey, OpExpr] = {}
        for key, branch in self.branches.items():
            new[key] = operator * branch if key[control_index] == 1 else branch
        return self._replace(new)

    def _apply_cx_system_to_ancilla(self, system_index: int, target_ancilla_index: int) -> "BranchState":
        z_op = pauli("Z", index=system_index, num_qubits=self.num_system_qubits)
        p0 = (identity(self.num_system_qubits) + z_op).scale(sp.Rational(1, 2))
        p1 = (identity(self.num_system_qubits) - z_op).scale(sp.Rational(1, 2))
        new: Dict[BranchKey, OpExpr] = {}

        for key, branch in self.branches.items():
            flipped = _flip_bit(key, target_ancilla_index)
            _add_branch(new, key, p0 * branch, self.num_system_qubits)
            _add_branch(new, flipped, p1 * branch, self.num_system_qubits)

        return self._replace(new)

    def _apply_cx_ancilla_to_ancilla(self, control_index: int, target_index: int) -> "BranchState":
        new: Dict[BranchKey, OpExpr] = {}
        for key, branch in self.branches.items():
            new_key = _flip_bit(key, target_index) if key[control_index] == 1 else key
            _add_branch(new, new_key, branch, self.num_system_qubits)
        return self._replace(new)

    def _apply_cz_ancilla_to_ancilla(self, left_index: int, right_index: int) -> "BranchState":
        new: Dict[BranchKey, OpExpr] = {}
        for key, branch in self.branches.items():
            new[key] = -branch if key[left_index] == 1 and key[right_index] == 1 else branch
        return self._replace(new)

    def _replace(self, branches: Mapping[BranchKey, OpExpr]) -> "BranchState":
        return BranchState(self.num_ancillas, self.num_system_qubits, branches)


def _ancilla_gate_matrix(name: str, theta: Optional[sp.Expr]) -> Tuple[Tuple[sp.Expr, sp.Expr], Tuple[sp.Expr, sp.Expr]]:
    if name == "x":
        return ((0, 1), (1, 0))
    if name == "z":
        return ((1, 0), (0, -1))
    if name == "s":
        return ((1, 0), (0, sp.I))
    if name == "sdg":
        return ((1, 0), (0, -sp.I))
    if name == "h":
        inv_sqrt2 = sp.sqrt(2) / 2
        return ((inv_sqrt2, inv_sqrt2), (inv_sqrt2, -inv_sqrt2))
    if name == "rx":
        c = cos_half(theta)
        minus_i_s = -sp.I * sin_half(theta)
        return ((c, minus_i_s), (minus_i_s, c))
    if name == "ry":
        c = cos_half(theta)
        s = sin_half(theta)
        return ((c, -s), (s, c))
    if name == "rz":
        return ((exp_minus_i_half(theta), 0), (0, exp_plus_i_half(theta)))
    raise UnsupportedGateError(f"Unsupported ancilla gate: {name}")


def _system_gate_operator(name: str, theta: Optional[sp.Expr], *, system_index: int, num_system_qubits: int) -> OpExpr:
    if name == "x":
        return pauli("X", index=system_index, num_qubits=num_system_qubits)
    if name == "z":
        return pauli("Z", index=system_index, num_qubits=num_system_qubits)
    if name == "s":
        return _s_operator(system_index=system_index, num_system_qubits=num_system_qubits)
    if name == "sdg":
        return _sdg_operator(system_index=system_index, num_system_qubits=num_system_qubits)
    if name == "h":
        return (
            pauli("X", index=system_index, num_qubits=num_system_qubits)
            + pauli("Z", index=system_index, num_qubits=num_system_qubits)
        ).scale(sp.sqrt(2) / 2)
    if name in {"rx", "ry", "rz"}:
        return _rotation_operator(name, theta, system_index=system_index, num_system_qubits=num_system_qubits)
    raise UnsupportedGateError(f"Unsupported system gate: {name}")


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


def _normalize_ancillas(*, ancilla: int, ancillas: Optional[Sequence[int]]) -> Tuple[int, ...]:
    if ancillas is None:
        return (ancilla,)
    normalized = tuple(ancillas)
    if not normalized:
        raise ValueError("At least one ancilla qubit is required")
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"Duplicate ancilla qubits are not allowed: {normalized}")
    return normalized


def _normalize_systems(*, system: int, systems: Optional[Sequence[int]]) -> Tuple[int, ...]:
    if systems is None:
        return (system,)
    normalized = tuple(systems)
    if not normalized:
        raise ValueError("At least one system qubit is required")
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"Duplicate system qubits are not allowed: {normalized}")
    return normalized


def _qubit_index(qasm_qubit: int, qubits: Tuple[int, ...]) -> Optional[int]:
    try:
        return qubits.index(qasm_qubit)
    except ValueError:
        return None


def _set_bit(key: BranchKey, index: int, value: int) -> BranchKey:
    bits = list(key)
    bits[index] = value
    return tuple(bits)


def _flip_bit(key: BranchKey, index: int) -> BranchKey:
    return _set_bit(key, index, 1 - key[index])


def _add_branch(branches: Dict[BranchKey, OpExpr], key: BranchKey, value: OpExpr, num_system_qubits: int) -> None:
    branches[key] = branches.get(key, zero(num_system_qubits)) + value
