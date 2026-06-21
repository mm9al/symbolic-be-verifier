from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Literal, Mapping, Optional, Sequence, Tuple, Union

import sympy as sp

from .expr import OpExpr, identity, pauli, zero
from .scalar import cos_half, exp_minus_i_half, exp_plus_i_half, sin_half
from .word_expr import WordExpr, atom as word_atom


BranchKey = Tuple[int, ...]
BranchValue = Union[OpExpr, WordExpr]
ExpressionKind = Literal["op", "word"]


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
    branches: Mapping[BranchKey, BranchValue]
    expression_kind: ExpressionKind = "op"

    def __post_init__(self) -> None:
        if self.num_ancillas < 1:
            raise ValueError("num_ancillas must be at least 1")
        if self.num_system_qubits < 1:
            raise ValueError("num_system_qubits must be at least 1")

        expression_kind = self.expression_kind
        for branch in self.branches.values():
            if isinstance(branch, WordExpr):
                expression_kind = "word"
                break

        cleaned: Dict[BranchKey, BranchValue] = {}
        for key, branch in self.branches.items():
            normalized_key = tuple(key)
            if len(normalized_key) != self.num_ancillas:
                raise ValueError(
                    f"Branch key {normalized_key!r} has length {len(normalized_key)}, expected {self.num_ancillas}"
                )
            if any(bit not in (0, 1) for bit in normalized_key):
                raise ValueError(f"Branch key must be a bitstring: {normalized_key!r}")
            normalized_branch = branch if isinstance(branch, WordExpr) else branch.with_num_qubits(self.num_system_qubits)
            if not normalized_branch.is_zero():
                cleaned[normalized_key] = normalized_branch
        object.__setattr__(self, "branches", cleaned)
        object.__setattr__(self, "expression_kind", expression_kind)

    @staticmethod
    def initial(
        num_system_qubits: int = 1,
        num_ancillas: int = 1,
        *,
        expression_kind: ExpressionKind = "op",
    ) -> "BranchState":
        zero_key = (0,) * num_ancillas
        identity_branch: BranchValue
        if expression_kind == "word":
            identity_branch = WordExpr.identity()
        else:
            identity_branch = identity(num_system_qubits)
        return BranchState(num_ancillas, num_system_qubits, {zero_key: identity_branch}, expression_kind=expression_kind)

    @property
    def b0(self) -> BranchValue:
        return self.branch((0,))

    @property
    def b1(self) -> BranchValue:
        return self.branch((1,))

    def top_left(self) -> BranchValue:
        return self.branch((0,) * self.num_ancillas)

    def branch(self, key: BranchKey) -> BranchValue:
        if len(key) != self.num_ancillas:
            if self.num_ancillas == 1 and len(key) == 1:
                pass
            else:
                raise ValueError(f"Branch key {key!r} does not match {self.num_ancillas} ancillas")
        return self.branches.get(tuple(key), self._zero())

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

        if name in {"uh", "uhdg"}:
            return self._apply_block_encoding_gate(name, gate.qubits, ancillas=ancilla_qubits, systems=system_qubits)

        if name == "mcx":
            return self._apply_mcx_gate(gate.qubits, ancillas=ancilla_qubits)

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

    def _apply_block_encoding_gate(
        self,
        name: str,
        qubits: Tuple[int, ...],
        *,
        ancillas: Tuple[int, ...],
        systems: Tuple[int, ...],
    ) -> "BranchState":
        block_indices = tuple(index for qubit in qubits if (index := _qubit_index(qubit, ancillas)) is not None)
        system_indices = tuple(index for qubit in qubits if (index := _qubit_index(qubit, systems)) is not None)
        unknown = tuple(qubit for qubit in qubits if _qubit_index(qubit, ancillas) is None and _qubit_index(qubit, systems) is None)

        if unknown:
            raise UnsupportedGateError(f"{name} operands must be known ancilla/system qubits: {qubits}")
        if not block_indices or not system_indices:
            raise UnsupportedGateError(f"{name} must include at least one block ancilla and one system qubit")
        if len(set(block_indices)) != len(block_indices):
            raise UnsupportedGateError(f"{name} has duplicate block ancilla operands")
        if len(block_indices) == 1:
            return self._apply_uh(block_indices[0], dagger=(name == "uhdg"))
        return self._apply_multi_uh(block_indices, dagger=(name == "uhdg"))

    def _apply_mcx_gate(self, qubits: Tuple[int, ...], *, ancillas: Tuple[int, ...]) -> "BranchState":
        if len(qubits) < 2:
            raise UnsupportedGateError("mcx requires at least one control and one target")
        control_indices = tuple(_qubit_index(qubit, ancillas) for qubit in qubits[:-1])
        target_index = _qubit_index(qubits[-1], ancillas)
        if target_index is None or any(index is None for index in control_indices):
            raise UnsupportedGateError("mcx currently supports ancilla controls and an ancilla target")
        controls = tuple(index for index in control_indices if index is not None)
        if target_index in controls:
            raise UnsupportedGateError("mcx target cannot also be a control")
        if len(set(controls)) != len(controls):
            raise UnsupportedGateError("mcx controls must be distinct")
        if len(controls) == 1:
            return self._apply_cx_ancilla_to_ancilla(controls[0], target_index)
        return self._apply_mcx_ancilla_to_ancilla(controls, target_index)

    def _apply_ancilla_matrix(self, ancilla_index: int, matrix: Tuple[Tuple[sp.Expr, sp.Expr], Tuple[sp.Expr, sp.Expr]]) -> "BranchState":
        new: Dict[BranchKey, BranchValue] = {}
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
        if self.expression_kind == "word":
            raise UnsupportedGateError("Explicit system gates are not supported in QSP word mode")
        return self._replace({key: operator * branch for key, branch in self.branches.items()})

    def _apply_controlled_pauli_ancilla_to_system(self, control_index: int, operator: OpExpr) -> "BranchState":
        if self.expression_kind == "word":
            raise UnsupportedGateError("Controlled system Pauli gates are not supported in QSP word mode")
        new: Dict[BranchKey, BranchValue] = {}
        for key, branch in self.branches.items():
            new[key] = operator * branch if key[control_index] == 1 else branch
        return self._replace(new)

    def _apply_cx_system_to_ancilla(self, system_index: int, target_ancilla_index: int) -> "BranchState":
        if self.expression_kind == "word":
            raise UnsupportedGateError("System-controlled cx is not supported in QSP word mode")
        z_op = pauli("Z", index=system_index, num_qubits=self.num_system_qubits)
        p0 = (identity(self.num_system_qubits) + z_op).scale(sp.Rational(1, 2))
        p1 = (identity(self.num_system_qubits) - z_op).scale(sp.Rational(1, 2))
        new: Dict[BranchKey, BranchValue] = {}

        for key, branch in self.branches.items():
            flipped = _flip_bit(key, target_ancilla_index)
            _add_branch(new, key, p0 * branch, self.num_system_qubits)
            _add_branch(new, flipped, p1 * branch, self.num_system_qubits)

        return self._replace(new)

    def _apply_cx_ancilla_to_ancilla(self, control_index: int, target_index: int) -> "BranchState":
        new: Dict[BranchKey, BranchValue] = {}
        for key, branch in self.branches.items():
            new_key = _flip_bit(key, target_index) if key[control_index] == 1 else key
            _add_branch(new, new_key, branch, self.num_system_qubits)
        return self._replace(new)

    def _apply_mcx_ancilla_to_ancilla(self, control_indices: Tuple[int, ...], target_index: int) -> "BranchState":
        new: Dict[BranchKey, BranchValue] = {}
        for key, branch in self.branches.items():
            new_key = _flip_bit(key, target_index) if all(key[index] == 1 for index in control_indices) else key
            _add_branch(new, new_key, branch, self.num_system_qubits)
        return self._replace(new)

    def _apply_cz_ancilla_to_ancilla(self, left_index: int, right_index: int) -> "BranchState":
        new: Dict[BranchKey, BranchValue] = {}
        for key, branch in self.branches.items():
            new[key] = -branch if key[left_index] == 1 and key[right_index] == 1 else branch
        return self._replace(new)

    def _apply_uh(self, ancilla_index: int, *, dagger: bool) -> "BranchState":
        if self.expression_kind != "word":
            raise UnsupportedGateError("UH/UHdg requires QSP word mode")

        if dagger:
            top_left, top_right, bottom_left, bottom_right = "Hd", "Ad", "Gd", "Cd"
        else:
            top_left, top_right, bottom_left, bottom_right = "H", "G", "A", "C"

        new: Dict[BranchKey, BranchValue] = {}
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
            new[key0] = word_atom(top_left) * branch0 + word_atom(top_right) * branch1
            new[key1] = word_atom(bottom_left) * branch0 + word_atom(bottom_right) * branch1
            processed.add(key0)
            processed.add(key1)

        return self._replace(new)

    def _apply_multi_uh(self, block_indices: Tuple[int, ...], *, dagger: bool) -> "BranchState":
        if self.expression_kind != "word":
            raise UnsupportedGateError("UH/UHdg requires QSP word mode")

        if dagger:
            top_left, top_right, bottom_left, bottom_right = "Hd", "Ad", "Gd", "Cd"
        else:
            top_left, top_right, bottom_left, bottom_right = "H", "G", "A", "C"

        new: Dict[BranchKey, BranchValue] = {}
        processed = set()
        group_keys = {_clear_bits(key, block_indices) for key in self.branches}

        for zero_key in sorted(group_keys):
            if zero_key in processed:
                continue

            complement_key = _canonical_complement_key(zero_key, block_indices)
            branch0 = self.branch(zero_key)
            branch1 = self._block_complement_branch(zero_key, block_indices)

            new[zero_key] = word_atom(top_left) * branch0 + word_atom(top_right) * branch1
            new[complement_key] = word_atom(bottom_left) * branch0 + word_atom(bottom_right) * branch1
            processed.add(zero_key)

        return self._replace(new)

    def _block_complement_branch(self, zero_key: BranchKey, block_indices: Tuple[int, ...]) -> BranchValue:
        total = self._zero()
        for block_bits in _nonzero_bit_patterns(len(block_indices)):
            key = list(zero_key)
            for block_index, bit in zip(block_indices, block_bits):
                key[block_index] = bit
            total = total + self.branch(tuple(key))
        return total

    def _replace(self, branches: Mapping[BranchKey, BranchValue]) -> "BranchState":
        return BranchState(self.num_ancillas, self.num_system_qubits, branches, expression_kind=self.expression_kind)

    def _zero(self) -> BranchValue:
        if self.expression_kind == "word":
            return WordExpr.zero()
        return zero(self.num_system_qubits)


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


def _clear_bits(key: BranchKey, indices: Tuple[int, ...]) -> BranchKey:
    bits = list(key)
    for index in indices:
        bits[index] = 0
    return tuple(bits)


def _canonical_complement_key(zero_key: BranchKey, block_indices: Tuple[int, ...]) -> BranchKey:
    bits = list(zero_key)
    bits[block_indices[0]] = 1
    return tuple(bits)


def _nonzero_bit_patterns(width: int) -> Tuple[Tuple[int, ...], ...]:
    return tuple(
        tuple((value >> shift) & 1 for shift in reversed(range(width)))
        for value in range(1, 2**width)
    )


def _add_branch(branches: Dict[BranchKey, BranchValue], key: BranchKey, value: BranchValue, num_system_qubits: int) -> None:
    fallback: BranchValue = WordExpr.zero() if isinstance(value, WordExpr) else zero(num_system_qubits)
    branches[key] = branches.get(key, fallback) + value
