from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from symbolic.qasm_parser import parse_qasm_text  # noqa: E402
from symbolic.verify import parse_polynomial  # noqa: E402


SOURCE_QASM = ROOT / "examples" / "qsp_hamsim_full_t05_eps01" / "qsp_hamsim_full_t05_eps01_deg3.qasm"
SOURCE_METADATA = ROOT / "examples" / "qsp_hamsim_full_t05_eps01" / "expected_polynomial.json"
DEFAULT_EXPANDED_QASM = ROOT / "examples" / "full_hamsim_embedded_lcu_x_minus_z_block_encoding_deg3.qasm"

ANCILLAS = (0, 1, 2, 3)
SYSTEM = 4
TOL = 1e-8

Matrix = tuple[tuple[complex, ...], ...]


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Dense-check full hamsim after replacing abstract UH with the concrete "
            "LCU block encoding from examples/lcu_x_minus_z.qasm."
        )
    )
    parser.add_argument(
        "--write-expanded-qasm",
        action="store_true",
        help=f"Write the expanded circuit to {DEFAULT_EXPANDED_QASM.relative_to(ROOT)} for inspection.",
    )
    args = parser.parse_args()

    metadata = json.loads(SOURCE_METADATA.read_text(encoding="utf-8"))
    expanded_qasm = expand_abstract_uh(SOURCE_QASM.read_text(encoding="utf-8"))
    if args.write_expanded_qasm:
        DEFAULT_EXPANDED_QASM.write_text(expanded_qasm, encoding="utf-8")

    actual = dense_top_left_block(expanded_qasm)
    base = matrix_scale(matrix_add(pauli_x_matrix(), matrix_scale(pauli_z_matrix(), -1)), 0.5)
    expected = polynomial_on_base(metadata["polynomial"], base)
    diff = max_abs_diff(actual, expected)

    print("full hamsim + embedded block encoding dense check")
    print(f"source: {SOURCE_QASM.relative_to(ROOT)}")
    print("embedded U_H: examples/lcu_x_minus_z.qasm, top-left block H = (X - Z)/2")
    print(f"expected polynomial: {metadata['polynomial']}")
    print(f"max |actual - expected|: {diff:.3e}")

    if diff > TOL:
        print_matrix("actual", actual)
        print_matrix("expected", expected)
        raise SystemExit(1)

    print("PASS")
    return 0


def expand_abstract_uh(qasm: str) -> str:
    lines: list[str] = []
    for line in qasm.splitlines():
        stripped = line.strip()
        if stripped.startswith(("opaque UH ", "opaque UHdg ", "opaque cUH ", "opaque cUHdg ")):
            continue
        if stripped == "UH q[3], q[4];":
            lines.extend(lcu_x_minus_z_uh_lines(control=None, dagger=False))
            continue
        if stripped == "UHdg q[3], q[4];":
            lines.extend(lcu_x_minus_z_uh_lines(control=None, dagger=True))
            continue
        if stripped == "cUH q[0], q[3], q[4];":
            lines.extend(lcu_x_minus_z_uh_lines(control="q[0]", dagger=False))
            continue
        if stripped == "cUHdg q[0], q[3], q[4];":
            lines.extend(lcu_x_minus_z_uh_lines(control="q[0]", dagger=True))
            continue
        lines.append(line)
    return "\n".join(lines) + "\n"


def lcu_x_minus_z_uh_lines(*, control: str | None, dagger: bool) -> list[str]:
    if control is None:
        forward = [
            "// expanded UH from examples/lcu_x_minus_z.qasm: top-left block H = (X - Z)/2",
            "// PREPARE",
            "x q[3];",
            "h q[3];",
            "// SELECT X on selection |0>",
            "x q[3];",
            "cx q[3], q[4];",
            "x q[3];",
            "// SELECT Z on selection |1>",
            "cz q[3], q[4];",
            "// PREPARE dagger",
            "h q[3];",
        ]
        daggered = [
            "// expanded UHdg from examples/lcu_x_minus_z.qasm",
            "// PREPARE dagger adjoint",
            "h q[3];",
            "// SELECT Z adjoint",
            "cz q[3], q[4];",
            "// SELECT X adjoint",
            "x q[3];",
            "cx q[3], q[4];",
            "x q[3];",
            "// PREPARE adjoint",
            "h q[3];",
            "x q[3];",
        ]
    else:
        forward = [
            "// expanded cUH from examples/lcu_x_minus_z.qasm: selector-controlled LCU block",
            "// controlled PREPARE",
            "cx q[0], q[3];",
            "ch q[0], q[3];",
            "// controlled SELECT X on selection |0>",
            "cx q[0], q[3];",
            "mcx q[0], q[3], q[4];",
            "cx q[0], q[3];",
            "// controlled SELECT Z on selection |1>",
            "h q[4];",
            "mcx q[0], q[3], q[4];",
            "h q[4];",
            "// controlled PREPARE dagger",
            "ch q[0], q[3];",
        ]
        daggered = [
            "// expanded cUHdg from examples/lcu_x_minus_z.qasm",
            "// controlled PREPARE dagger adjoint",
            "ch q[0], q[3];",
            "// controlled SELECT Z adjoint",
            "h q[4];",
            "mcx q[0], q[3], q[4];",
            "h q[4];",
            "// controlled SELECT X adjoint",
            "cx q[0], q[3];",
            "mcx q[0], q[3], q[4];",
            "cx q[0], q[3];",
            "// controlled PREPARE adjoint",
            "ch q[0], q[3];",
            "cx q[0], q[3];",
        ]
    return daggered if dagger else forward


def dense_top_left_block(qasm: str) -> Matrix:
    gates = parse_qasm_text(qasm)
    num_qubits = max(qubit for gate in gates for qubit in gate.qubits) + 1
    dim = 2**num_qubits
    actual: list[list[complex]] = [[0j for _ in range(2)] for _ in range(2)]

    for input_system_bit in (0, 1):
        state = [0j] * dim
        state[input_system_bit << SYSTEM] = 1
        for gate in gates:
            state = apply_gate(state, gate.name, gate.qubits, gate.parameter)
        for output_system_bit in (0, 1):
            actual[output_system_bit][input_system_bit] = state[output_system_bit << SYSTEM]

    return tuple(tuple(row) for row in actual)


def apply_gate(state: list[complex], name: str, qubits: tuple[int, ...], parameter: object) -> list[complex]:
    name = name.lower()
    if name == "x":
        return apply_controlled_x(state, (), qubits[0])
    if name == "cx":
        return apply_controlled_x(state, qubits[:-1], qubits[-1])
    if name == "mcx":
        return apply_controlled_x(state, qubits[:-1], qubits[-1])
    if name == "cz":
        return apply_controlled_z(state, qubits[:-1], qubits[-1])
    if name == "ch":
        return apply_controlled_single_qubit(state, qubits[:-1], qubits[-1], single_qubit_matrix("h", None))
    if name in {"h", "z", "s", "sdg", "rz"}:
        return apply_single_qubit(state, qubits[0], single_qubit_matrix(name, parameter))
    raise ValueError(f"Unsupported gate in dense check: {name}")


def single_qubit_matrix(name: str, parameter: object) -> Matrix:
    if name == "h":
        inv_sqrt2 = 1 / math.sqrt(2)
        return ((inv_sqrt2, inv_sqrt2), (inv_sqrt2, -inv_sqrt2))
    if name == "z":
        return ((1, 0), (0, -1))
    if name == "s":
        return ((1, 0), (0, 1j))
    if name == "sdg":
        return ((1, 0), (0, -1j))
    if name == "rz":
        angle = float(parameter)
        return (
            (complex(math.cos(angle / 2), -math.sin(angle / 2)), 0),
            (0, complex(math.cos(angle / 2), math.sin(angle / 2))),
        )
    raise ValueError(f"Unsupported single-qubit gate: {name}")


def apply_single_qubit(state: list[complex], target: int, matrix: Matrix) -> list[complex]:
    out = state[:]
    target_bit = 1 << target
    for zero in range(len(state)):
        if zero & target_bit:
            continue
        one = zero | target_bit
        value0 = state[zero]
        value1 = state[one]
        out[zero] = matrix[0][0] * value0 + matrix[0][1] * value1
        out[one] = matrix[1][0] * value0 + matrix[1][1] * value1
    return out


def apply_controlled_x(state: list[complex], controls: tuple[int, ...], target: int) -> list[complex]:
    out = state[:]
    target_bit = 1 << target
    control_mask = sum(1 << control for control in controls)
    for zero in range(len(state)):
        if zero & target_bit or (zero & control_mask) != control_mask:
            continue
        one = zero | target_bit
        out[zero] = state[one]
        out[one] = state[zero]
    return out


def apply_controlled_z(state: list[complex], controls: tuple[int, ...], target: int) -> list[complex]:
    out = state[:]
    target_bit = 1 << target
    control_mask = sum(1 << control for control in controls)
    for index, value in enumerate(state):
        if index & target_bit and (index & control_mask) == control_mask:
            out[index] = -value
    return out


def apply_controlled_single_qubit(
    state: list[complex], controls: tuple[int, ...], target: int, matrix: Matrix
) -> list[complex]:
    out = state[:]
    target_bit = 1 << target
    control_mask = sum(1 << control for control in controls)
    for zero in range(len(state)):
        if zero & target_bit or (zero & control_mask) != control_mask:
            continue
        one = zero | target_bit
        value0 = state[zero]
        value1 = state[one]
        out[zero] = matrix[0][0] * value0 + matrix[0][1] * value1
        out[one] = matrix[1][0] * value0 + matrix[1][1] * value1
    return out


def polynomial_on_base(polynomial: str, base: Matrix) -> Matrix:
    x = sp.Symbol("x")
    parsed = sp.Poly(parse_polynomial(polynomial), x)
    result = zero_matrix(2)
    powers = {0: identity_matrix(2)}
    for (degree,), coeff in parsed.terms():
        while degree not in powers:
            powers[len(powers)] = matrix_mul(powers[len(powers) - 1], base)
        result = matrix_add(result, matrix_scale(powers[degree], complex(sp.N(coeff, 30))))
    return result


def pauli_x_matrix() -> Matrix:
    return ((0, 1), (1, 0))


def pauli_z_matrix() -> Matrix:
    return ((1, 0), (0, -1))


def identity_matrix(size: int) -> Matrix:
    return tuple(tuple(1 if row == column else 0 for column in range(size)) for row in range(size))


def zero_matrix(size: int) -> Matrix:
    return tuple(tuple(0 for _ in range(size)) for _ in range(size))


def matrix_add(left: Matrix, right: Matrix) -> Matrix:
    return tuple(tuple(a + b for a, b in zip(left_row, right_row)) for left_row, right_row in zip(left, right))


def matrix_scale(matrix: Matrix, scalar: complex) -> Matrix:
    return tuple(tuple(scalar * value for value in row) for row in matrix)


def matrix_mul(left: Matrix, right: Matrix) -> Matrix:
    width = len(right[0])
    return tuple(
        tuple(sum(left[row][inner] * right[inner][column] for inner in range(len(right))) for column in range(width))
        for row in range(len(left))
    )


def max_abs_diff(left: Matrix, right: Matrix) -> float:
    return max(abs(a - b) for left_row, right_row in zip(left, right) for a, b in zip(left_row, right_row))


def print_matrix(label: str, matrix: Matrix) -> None:
    print(f"{label}:")
    for row in matrix:
        print("  " + "  ".join(f"{value.real:+.12g}{value.imag:+.12g}j" for value in row))


if __name__ == "__main__":
    raise SystemExit(main())
