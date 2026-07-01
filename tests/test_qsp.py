import json
import math
from pathlib import Path

import pytest
import sympy as sp

from symbolic.expr import parse_operator_expression
from symbolic.qasm_parser import parse_qasm_file
from symbolic.qsp import (
    QasmGate,
    adjoint_gate,
    build_raw_qsp_gate_list,
    daggerize_gate_list,
    full_hamsim_qasm_snippet,
    hamsim_exp_polynomial_expr,
    polynomial_expr,
    qasm_phase_data,
)
from symbolic.verify import PASS, eval_polynomial_on_pauliop, parse_polynomial, pauli_expr_close, polynomial_close, verify_qasm_file


ROOT = Path(__file__).parents[1]
QSP_EXAMPLE_DIR = ROOT / "examples" / "qsp_hamsim_t05_eps1e-4"
QSP_M2_EXAMPLE_DIR = ROOT / "examples" / "qsp_hamsim_t05_eps1e-4_m2"
QSP_FULL_EXAMPLE_DIR = ROOT / "examples" / "qsp_hamsim_full_t05_eps01"
QSP_FULL_XYZ_EXAMPLE_DIR = ROOT / "examples" / "qsp_hamsim_full_xyz_t05_eps01"

COS_POLYNOMIAL = "0.49999966355545339 - 0.062493938728279574*x^2 + 0.0012858918109143005*x^4"
SIN_POLYNOMIAL = "0.24999991579484243*x - 0.010415992523176125*x^3 + 0.00012885803586171963*x^5"
FULL_HAMSIM_DEG3_POLYNOMIAL = (
    "0.24991946353954456 - 0.12497982382931781*i*x "
    "- 0.030604023458682638*x^2 + 0.005127459989174488*i*x^3"
)
RAW_COS_COMPONENT_Y_RESPONSE = "exp(67374387816677681*i/10000000000000000)"


def test_phase_conversion_uses_qsvt_projector_convention():
    qsvt_phases, rz_angles = qasm_phase_data([0.0, 0.0, 0.0])

    assert _close_list(qsvt_phases, [math.pi / 4, math.pi / 2, math.pi / 4])
    assert _close_list(rz_angles, [-math.pi / 2, -math.pi, -math.pi / 2])


def test_polynomial_expr_formats_sparse_qsp_coefficients():
    assert polynomial_expr([0.5, 0.0, -0.125, 0.0, 0.01]) == "0.5 - 0.125*x^2 + 0.01*x^4"


def test_hamsim_exp_polynomial_expr_combines_cos_minus_i_sin_with_scale():
    assert hamsim_exp_polynomial_expr([1.0, 0.0, -0.25], [0.0, 0.5], scale=0.5) == "0.5 - 0.25*i*x - 0.125*x^2"


def test_qsp_dagger_gate_list_is_reverse_adjoint():
    gates = [
        QasmGate("rz", ("q[0]",), 1.25),
        QasmGate("UH", ("q[1]", "q[2]")),
        QasmGate("UHdg", ("q[1]", "q[2]")),
        QasmGate("sdg", ("q[0]",)),
        QasmGate("x", ("q[1]",)),
    ]

    dag = daggerize_gate_list(gates)

    assert dag == [adjoint_gate(gate) for gate in reversed(gates)]
    assert dag == [
        QasmGate("x", ("q[1]",)),
        QasmGate("s", ("q[0]",)),
        QasmGate("UH", ("q[1]", "q[2]")),
        QasmGate("UHdg", ("q[1]", "q[2]")),
        QasmGate("rz", ("q[0]",), -1.25),
    ]


def test_raw_qsp_dagger_gate_list_reverses_phase_and_signal_order():
    raw = build_raw_qsp_gate_list(
        [0.0, 0.0, 0.0, 0.0],
        [0.1, 0.2, 0.3, 0.4],
        [1.0, 2.0, 3.0, 4.0],
        phase_qubit="q[0]",
        block_ancillas=["q[1]"],
        system_qubits=["q[2]"],
        signal_gate="UH",
        signal_gate_dagger="UHdg",
    )

    raw_signals = [gate.name for gate in raw if gate.name in {"UH", "UHdg"}]
    dag_signals = [gate.name for gate in daggerize_gate_list(raw) if gate.name in {"UH", "UHdg"}]
    raw_angles = [gate.parameter for gate in raw if gate.name == "rz"]
    dag_angles = [gate.parameter for gate in daggerize_gate_list(raw) if gate.name == "rz"]

    assert raw_signals == ["UH", "UHdg", "UH"]
    assert dag_signals == ["UHdg", "UH", "UHdg"]
    assert dag_angles == [-angle for angle in reversed(raw_angles)]


def test_full_hamsim_qasm_snippet_uses_interleaved_phase_multiplexing():
    cos_record = {
        "component": "cos",
        "pyqsp_phases": [0.0, 0.0, 0.0],
        "qsvt_projector_phases": [0.1, 0.2, 0.3],
        "qasm_rz_angles": [1.0, 2.0, 3.0],
    }
    sin_record = {
        "component": "sin",
        "pyqsp_phases": [0.0, 0.0, 0.0, 0.0],
        "qsvt_projector_phases": [0.4, 0.5, 0.6, 0.7],
        "qasm_rz_angles": [4.0, 5.0, 6.0, 7.0],
    }

    snippet = full_hamsim_qasm_snippet(
        cos_record,
        sin_record,
        selector_qubit="q[0]",
        component_selector_qubit="q[1]",
        phase_qubit="q[2]",
        block_ancillas=["q[3]"],
        system_qubits=["q[4]"],
        signal_gate="UH",
        signal_gate_dagger="UHdg",
        controlled_signal_gate="cUH",
        controlled_signal_gate_dagger="cUHdg",
    )

    assert "// common phase 0" in snippet
    assert "// extra signed difference phase on sin branch only" in snippet
    assert "// common U" in snippet
    assert "// common U^\\dagger" in snippet
    assert "// sin-only final U" in snippet
    assert "cUH q[0], q[3], q[4];" in snippet
    assert "Cdag branch" not in snippet
    assert "// extract pyqsp imaginary response on the component selector" in snippet
    assert "h q[1];\ns q[1];\nx q[1];" in snippet
    assert "// combine 1/2 * (E_cos - i E_sin)" in snippet
    assert snippet.endswith("sdg q[0];\nh q[0];")


def test_full_hamsim_qasm_snippet_uses_degree_dependent_component_extraction():
    cos_record = {
        "component": "cos",
        "pyqsp_phases": [0.0] * 5,
        "qsvt_projector_phases": [0.1] * 5,
        "qasm_rz_angles": [1.0] * 5,
    }
    sin_record = {
        "component": "sin",
        "pyqsp_phases": [0.0] * 6,
        "qsvt_projector_phases": [0.2] * 6,
        "qasm_rz_angles": [2.0] * 6,
    }

    snippet = full_hamsim_qasm_snippet(
        cos_record,
        sin_record,
        selector_qubit="q[0]",
        component_selector_qubit="q[1]",
        phase_qubit="q[2]",
        block_ancillas=["q[3]"],
        system_qubits=["q[4]"],
        signal_gate="UH",
        signal_gate_dagger="UHdg",
        controlled_signal_gate="cUH",
        controlled_signal_gate_dagger="cUHdg",
    )

    assert "cUH q[0], q[3], q[4];" in snippet
    assert "h q[1];\nsdg q[1];\nx q[1];" in snippet


@pytest.mark.parametrize(
    ("example_dir", "block_ancillas", "selector_qubit", "system_qubits"),
    [
        (QSP_EXAMPLE_DIR, ["q[1]"], "q[2]", ["q[3]"]),
        (QSP_M2_EXAMPLE_DIR, ["q[1]", "q[2]"], "q[3]", ["q[4]"]),
    ],
)
def test_qsp_cos_sin_regression_passes_for_block_ancillas(example_dir, block_ancillas, selector_qubit, system_qubits):
    metadata = json.loads((example_dir / "expected_polynomials.json").read_text(encoding="utf-8"))

    assert metadata["phase_qubit"] == "q[0]"
    assert metadata["block_ancillas"] == block_ancillas
    assert metadata["selector_qubit"] == selector_qubit
    assert metadata["system_qubits"] == system_qubits

    _verify_qsp_metadata(metadata)


def test_hamsim_qsp_y_base_regression_uses_adjoint_not_complex_conjugate():
    y_base = "Y"

    raw_cos = _verify_qsp_file_on_base(
        QSP_EXAMPLE_DIR / "qsp_hamsim_cos_t05_eps1e-4_deg4.qasm",
        ancillas=(0, 1),
        systems=(2,),
        base=y_base,
        expected_polynomial=RAW_COS_COMPONENT_Y_RESPONSE,
    )
    assert pauli_expr_close(raw_cos.qsp_actual, parse_operator_expression(RAW_COS_COMPONENT_Y_RESPONSE))

    cos = _verify_qasm_dense_on_base(
        QSP_EXAMPLE_DIR / "qsp_hamsim_cos_selector_t05_eps1e-4_deg4.qasm",
        ancillas=(0, 1, 2),
        systems=(3,),
        base=_pauli_y_matrix(),
        expected_polynomial=COS_POLYNOMIAL,
    )

    sin = _verify_qasm_dense_on_base(
        QSP_EXAMPLE_DIR / "qsp_hamsim_sin_selector_t05_eps1e-4_deg5.qasm",
        ancillas=(0, 1, 2),
        systems=(3,),
        base=_pauli_y_matrix(),
        expected_polynomial=SIN_POLYNOMIAL,
    )

    final = _verify_qasm_dense_on_base(
        QSP_FULL_EXAMPLE_DIR / "qsp_hamsim_full_t05_eps01_deg3.qasm",
        ancillas=(0, 1, 2, 3),
        systems=(4,),
        base=_pauli_y_matrix(),
        expected_polynomial=FULL_HAMSIM_DEG3_POLYNOMIAL,
    )

    wrong_sin = _matrix_polynomial(SIN_POLYNOMIAL, _matrix_scale(_pauli_y_matrix(), -1))
    wrong_final = _matrix_polynomial(FULL_HAMSIM_DEG3_POLYNOMIAL, _matrix_scale(_pauli_y_matrix(), -1))
    assert not _matrix_close(sin, wrong_sin)
    assert not _matrix_close(final, wrong_final)

    mixed_base = _matrix_scale(_matrix_add(_matrix_add(_pauli_x_matrix(), _pauli_y_matrix()), _pauli_z_matrix()), 1 / 3)
    _verify_qasm_dense_on_base(
        QSP_EXAMPLE_DIR / "qsp_hamsim_cos_selector_t05_eps1e-4_deg4.qasm",
        ancillas=(0, 1, 2),
        systems=(3,),
        base=mixed_base,
        expected_polynomial=COS_POLYNOMIAL,
    )
    _verify_qasm_dense_on_base(
        QSP_EXAMPLE_DIR / "qsp_hamsim_sin_selector_t05_eps1e-4_deg5.qasm",
        ancillas=(0, 1, 2),
        systems=(3,),
        base=mixed_base,
        expected_polynomial=SIN_POLYNOMIAL,
    )
    _verify_qasm_dense_on_base(
        QSP_FULL_EXAMPLE_DIR / "qsp_hamsim_full_t05_eps01_deg3.qasm",
        ancillas=(0, 1, 2, 3),
        systems=(4,),
        base=mixed_base,
        expected_polynomial=FULL_HAMSIM_DEG3_POLYNOMIAL,
    )


def test_full_hamsim_xyz_example_verifies_mixed_hermitian_base():
    metadata = json.loads((QSP_FULL_XYZ_EXAMPLE_DIR / "expected_polynomial.json").read_text(encoding="utf-8"))
    ancillas = (
        metadata["selector_qubit"],
        metadata["component_selector_qubit"],
        metadata["phase_qubit"],
        *metadata["block_ancillas"],
    )

    result = verify_qasm_file(
        ROOT / metadata["qasm"],
        ancillas=tuple(_parse_qreg(qubit) for qubit in ancillas),
        systems=tuple(_parse_qreg(qubit) for qubit in metadata["system_qubits"]),
        base=metadata["base"],
        expected_polynomial=metadata["polynomial"],
        hermitian_base=True,
    )

    assert result.status == PASS
    assert result.qsp_polynomial_only is False
    assert pauli_expr_close(result.qsp_actual, parse_operator_expression(metadata["expected_operator"]))

    mixed_base = _matrix_scale(_matrix_add(_matrix_add(_pauli_x_matrix(), _pauli_y_matrix()), _pauli_z_matrix()), 1 / 3)
    _verify_qasm_dense_on_base(
        ROOT / metadata["qasm"],
        ancillas=(0, 1, 2, 3),
        systems=(4,),
        base=mixed_base,
        expected_polynomial=metadata["polynomial"],
    )


@pytest.mark.parametrize(
    ("qasm_path", "ancillas", "systems"),
    [
        (ROOT / "examples" / "qsp_t3_mcx_m1.qasm", (0, 1), (2,)),
        (ROOT / "examples" / "qsp_t3_mcx_m2.qasm", (0, 1, 2), (3,)),
    ],
)
def test_qsp_mcx_t3_passes_for_block_ancillas(qasm_path, ancillas, systems):
    result = verify_qasm_file(
        qasm_path,
        ancillas=ancillas,
        systems=systems,
        expected_polynomial="4*x^3 - 3*x",
        hermitian_base=True,
        compare_polynomial_only=True,
    )

    assert result.status == PASS
    assert polynomial_close(result.qsp_polynomial, "4*x^3 - 3*x")


def _parse_qreg(qubit: str) -> int:
    prefix = "q["
    assert qubit.startswith(prefix) and qubit.endswith("]")
    return int(qubit[len(prefix) : -1])


def _close_list(actual: list[float], expected: list[float]) -> bool:
    return len(actual) == len(expected) and all(math.isclose(a, b) for a, b in zip(actual, expected))


def _verify_qsp_metadata(metadata: dict) -> None:
    ancillas = (
        _parse_qreg(metadata["phase_qubit"]),
        *(_parse_qreg(qubit) for qubit in metadata["block_ancillas"]),
        _parse_qreg(metadata["selector_qubit"]),
    )
    systems = tuple(_parse_qreg(qubit) for qubit in metadata["system_qubits"])

    seen_components = set()
    for file_record in metadata["files"]:
        seen_components.add(file_record["component"])
        _verify_qasm_dense_on_base(
            ROOT / file_record["qasm"],
            ancillas=ancillas,
            systems=systems,
            expected_polynomial=file_record["polynomial"],
            base=_pauli_y_matrix(),
        )

    assert seen_components == {"cos", "sin"}


def _verify_qsp_file_on_base(
    qasm_path: Path,
    *,
    ancillas: tuple[int, ...],
    systems: tuple[int, ...],
    base: str,
    expected_polynomial: str,
):
    result = verify_qasm_file(
        qasm_path,
        ancillas=ancillas,
        systems=systems,
        base=base,
        expected_polynomial=expected_polynomial,
        hermitian_base=True,
    )

    assert result.status == PASS
    assert result.qsp_polynomial_only is False
    assert result.qsp_actual is not None
    assert result.qsp_expected is not None
    return result


Matrix = tuple[tuple[complex, ...], ...]


def _verify_qasm_dense_on_base(
    qasm_path: Path,
    *,
    ancillas: tuple[int, ...],
    systems: tuple[int, ...],
    base: Matrix,
    expected_polynomial: str,
) -> Matrix:
    assert len(systems) == 1
    gates = parse_qasm_file(qasm_path)
    num_qubits = max(qubit for gate in gates for qubit in gate.qubits) + 1
    unitary = _identity_matrix(2**num_qubits)

    for gate in gates:
        unitary = _apply_dense_gate(unitary, gate.name.lower(), gate.qubits, gate.parameter, num_qubits, systems[0], base)

    actual = _top_left_dense_block(unitary, num_qubits, ancillas, systems[0])
    expected = _matrix_polynomial(expected_polynomial, base)
    assert _matrix_close(actual, expected)
    return actual


def _apply_dense_gate(
    unitary: Matrix,
    name: str,
    qubits: tuple[int, ...],
    parameter: object,
    num_qubits: int,
    system_qubit: int,
    base: Matrix,
) -> Matrix:
    if name in {"x", "z", "h", "rz", "s", "sdg"}:
        matrix = {
            "x": ((0, 1), (1, 0)),
            "z": ((1, 0), (0, -1)),
            "h": ((1 / math.sqrt(2), 1 / math.sqrt(2)), (1 / math.sqrt(2), -1 / math.sqrt(2))),
            "s": ((1, 0), (0, 1j)),
            "sdg": ((1, 0), (0, -1j)),
        }.get(name)
        if name == "rz":
            angle = float(parameter)
            matrix = ((complex(math.cos(angle / 2), -math.sin(angle / 2)), 0), (0, complex(math.cos(angle / 2), math.sin(angle / 2))))
        return _apply_single_qubit_matrix(unitary, qubits[0], matrix, num_qubits)

    if name == "cx":
        return _apply_controlled_x(unitary, qubits[:-1], qubits[-1], num_qubits)
    if name == "mcx":
        return _apply_controlled_x(unitary, qubits[:-1], qubits[-1], num_qubits)

    control_count = _dense_controlled_uh_count(name)
    if name in {"uh", "uhdg"} or control_count is not None:
        controls = () if control_count is None else qubits[:control_count]
        operands = qubits if control_count is None else qubits[control_count:]
        block_qubits = tuple(qubit for qubit in operands if qubit != system_qubit)
        return _apply_dense_uh(unitary, controls, block_qubits, system_qubit, num_qubits, base)

    raise AssertionError(f"Unsupported dense gate in QSP test: {name}")


def _apply_single_qubit_matrix(unitary: Matrix, target: int, matrix: Matrix, num_qubits: int) -> Matrix:
    dim = len(unitary)
    out = [list(row) for row in unitary]
    bit = 1 << target
    for column in range(dim):
        for zero in range(dim):
            if zero & bit:
                continue
            one = zero | bit
            a0 = unitary[zero][column]
            a1 = unitary[one][column]
            out[zero][column] = matrix[0][0] * a0 + matrix[0][1] * a1
            out[one][column] = matrix[1][0] * a0 + matrix[1][1] * a1
    return tuple(tuple(row) for row in out)


def _apply_controlled_x(unitary: Matrix, controls: tuple[int, ...], target: int, num_qubits: int) -> Matrix:
    dim = len(unitary)
    out = [list(row) for row in unitary]
    target_bit = 1 << target
    control_mask = sum(1 << control for control in controls)
    for column in range(dim):
        for zero in range(dim):
            if zero & target_bit or (zero & control_mask) != control_mask:
                continue
            one = zero | target_bit
            out[zero][column] = unitary[one][column]
            out[one][column] = unitary[zero][column]
    return tuple(tuple(row) for row in out)


def _apply_dense_uh(
    unitary: Matrix,
    controls: tuple[int, ...],
    block_qubits: tuple[int, ...],
    system_qubit: int,
    num_qubits: int,
    base: Matrix,
) -> Matrix:
    if len(block_qubits) == 1:
        return _apply_dense_single_block_uh(unitary, controls, block_qubits[0], system_qubit, num_qubits, base)

    if _matrix_close(_matrix_mul(base, base), _identity_matrix(2)):
        return _apply_dense_unitary_base_on_zero_block(unitary, controls, block_qubits, system_qubit, num_qubits, base)

    raise AssertionError("Dense multi-block UH test helper only supports unitary bases")


def _apply_dense_single_block_uh(
    unitary: Matrix,
    controls: tuple[int, ...],
    block_qubit: int,
    system_qubit: int,
    num_qubits: int,
    base: Matrix,
) -> Matrix:
    dim = len(unitary)
    out = [list(row) for row in unitary]
    control_mask = sum(1 << control for control in controls)
    block_bit = 1 << block_qubit
    system_bit = 1 << system_qubit
    h2 = _matrix_mul(base, base)
    complement_scale = complex(math.sqrt(max(0.0, 1.0 - h2[0][0].real)))

    for column in range(dim):
        for base_index in range(dim):
            if base_index & block_bit or base_index & system_bit or (base_index & control_mask) != control_mask:
                continue
            indices = (base_index, base_index | system_bit, base_index | block_bit, base_index | block_bit | system_bit)
            values = [unitary[index][column] for index in indices]
            transformed = [
                base[0][0] * values[0] + base[0][1] * values[1] + complement_scale * values[2],
                base[1][0] * values[0] + base[1][1] * values[1] + complement_scale * values[3],
                complement_scale * values[0] - base[0][0] * values[2] - base[0][1] * values[3],
                complement_scale * values[1] - base[1][0] * values[2] - base[1][1] * values[3],
            ]
            for index, value in zip(indices, transformed):
                out[index][column] = value
    return tuple(tuple(row) for row in out)


def _apply_dense_unitary_base_on_zero_block(
    unitary: Matrix,
    controls: tuple[int, ...],
    block_qubits: tuple[int, ...],
    system_qubit: int,
    num_qubits: int,
    base: Matrix,
) -> Matrix:
    dim = len(unitary)
    out = [list(row) for row in unitary]
    control_mask = sum(1 << control for control in controls)
    block_mask = sum(1 << block for block in block_qubits)
    system_bit = 1 << system_qubit
    for column in range(dim):
        for zero in range(dim):
            if zero & system_bit or zero & block_mask or (zero & control_mask) != control_mask:
                continue
            one = zero | system_bit
            a0 = unitary[zero][column]
            a1 = unitary[one][column]
            out[zero][column] = base[0][0] * a0 + base[0][1] * a1
            out[one][column] = base[1][0] * a0 + base[1][1] * a1
    return tuple(tuple(row) for row in out)


def _top_left_dense_block(unitary: Matrix, num_qubits: int, ancillas: tuple[int, ...], system_qubit: int) -> Matrix:
    assert len(ancillas) + 1 == num_qubits
    zero_ancilla_mask = 0
    for bit in ancillas:
        zero_ancilla_mask |= 1 << bit
    assert zero_ancilla_mask | (1 << system_qubit) == (1 << num_qubits) - 1
    zero = 0
    one = 1 << system_qubit
    return ((unitary[zero][zero], unitary[zero][one]), (unitary[one][zero], unitary[one][one]))


def _matrix_polynomial(polynomial: str, base: Matrix) -> Matrix:
    x = sp.Symbol("x")
    parsed = parse_polynomial(polynomial)
    poly = sp.Poly(parsed, x)
    powers = {0: _identity_matrix(2)}
    result = _zero_matrix(2)
    for (degree,), coeff in poly.terms():
        while degree not in powers:
            powers[len(powers)] = _matrix_mul(powers[len(powers) - 1], base)
        result = _matrix_add(result, _matrix_scale(powers[degree], complex(sp.N(coeff, 30))))
    return result


def _dense_controlled_uh_count(name: str) -> int | None:
    if not name.endswith("uh") and not name.endswith("uhdg"):
        return None
    prefix = name.removesuffix("uhdg") if name.endswith("uhdg") else name.removesuffix("uh")
    if prefix and set(prefix) == {"c"}:
        return len(prefix)
    return None


def _pauli_x_matrix() -> Matrix:
    return ((0, 1), (1, 0))


def _pauli_y_matrix() -> Matrix:
    return ((0, -1j), (1j, 0))


def _pauli_z_matrix() -> Matrix:
    return ((1, 0), (0, -1))


def _identity_matrix(size: int) -> Matrix:
    return tuple(tuple(1 if row == column else 0 for column in range(size)) for row in range(size))


def _zero_matrix(size: int) -> Matrix:
    return tuple(tuple(0 for _ in range(size)) for _ in range(size))


def _matrix_add(left: Matrix, right: Matrix) -> Matrix:
    return tuple(tuple(a + b for a, b in zip(left_row, right_row)) for left_row, right_row in zip(left, right))


def _matrix_scale(matrix: Matrix, scalar: complex) -> Matrix:
    return tuple(tuple(scalar * value for value in row) for row in matrix)


def _matrix_mul(left: Matrix, right: Matrix) -> Matrix:
    width = len(right[0])
    return tuple(
        tuple(sum(left[row][inner] * right[inner][column] for inner in range(len(right))) for column in range(width))
        for row in range(len(left))
    )


def _matrix_close(left: Matrix, right: Matrix, tol: float = 1e-8) -> bool:
    return all(abs(a - b) <= tol for left_row, right_row in zip(left, right) for a, b in zip(left_row, right_row))
