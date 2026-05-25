from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import sympy as sp

from .expr import OpExpr, identity, pauli, zero


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
    def initial() -> "BranchState":
        return BranchState(identity(), zero())

    def apply(self, gate: Gate, *, ancilla: int = 0, system: int = 1) -> "BranchState":
        name = gate.name.lower()
        if gate.parameter is not None:
            raise UnsupportedGateError(f"Parameterized gate is not supported yet: {gate}")

        if len(gate.qubits) == 1:
            return self._apply_single_qubit(name, gate.qubits[0], ancilla=ancilla, system=system)
        if len(gate.qubits) == 2:
            return self._apply_two_qubit(name, gate.qubits, ancilla=ancilla, system=system)
        raise UnsupportedGateError(f"Unsupported gate arity: {gate}")

    def _apply_single_qubit(self, name: str, target: int, *, ancilla: int, system: int) -> "BranchState":
        if target == ancilla:
            if name == "x":
                return BranchState(self.b1, self.b0)
            if name == "z":
                return BranchState(self.b0, -self.b1)
            if name == "h":
                inv_sqrt2 = sp.sqrt(2) / 2
                return BranchState((self.b0 + self.b1).scale(inv_sqrt2), (self.b0 - self.b1).scale(inv_sqrt2))
            raise UnsupportedGateError(f"Unsupported ancilla gate: {name}")

        if target == system:
            if name == "x":
                gate_op = pauli("X")
            elif name == "z":
                gate_op = pauli("Z")
            elif name == "h":
                gate_op = (pauli("X") + pauli("Z")).scale(sp.sqrt(2) / 2)
            else:
                raise UnsupportedGateError(f"Unsupported system gate: {name}")
            return BranchState(gate_op * self.b0, gate_op * self.b1)

        raise UnsupportedGateError(f"Gate target q[{target}] is neither ancilla q[{ancilla}] nor system q[{system}]")

    def _apply_two_qubit(self, name: str, qubits: Tuple[int, int], *, ancilla: int, system: int) -> "BranchState":
        left, right = qubits

        if name == "cx":
            if left == ancilla and right == system:
                return BranchState(self.b0, pauli("X") * self.b1)
            raise UnsupportedGateError("Only cx with ancilla as control and system as target is supported")

        if name == "cz":
            if {left, right} == {ancilla, system}:
                return BranchState(self.b0, pauli("Z") * self.b1)
            raise UnsupportedGateError("Only cz between the ancilla and system qubit is supported")

        raise UnsupportedGateError(f"Unsupported two-qubit gate: {name}")
