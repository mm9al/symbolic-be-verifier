import json
import math
from pathlib import Path

import pytest

from symbolic.qsp import polynomial_expr, qasm_phase_data
from symbolic.verify import PASS, polynomial_close, verify_qasm_file


ROOT = Path(__file__).parents[1]
QSP_EXAMPLE_DIR = ROOT / "examples" / "qsp_hamsim_t05_eps1e-4"
QSP_M2_EXAMPLE_DIR = ROOT / "examples" / "qsp_hamsim_t05_eps1e-4_m2"


def test_phase_conversion_uses_qsvt_projector_convention():
    qsvt_phases, rz_angles = qasm_phase_data([0.0, 0.0, 0.0])

    assert _close_list(qsvt_phases, [math.pi / 4, math.pi / 2, math.pi / 4])
    assert _close_list(rz_angles, [-math.pi / 2, -math.pi, -math.pi / 2])


def test_polynomial_expr_formats_sparse_qsp_coefficients():
    assert polynomial_expr([0.5, 0.0, -0.125, 0.0, 0.01]) == "0.5 - 0.125*x^2 + 0.01*x^4"


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
        result = verify_qasm_file(
            ROOT / file_record["qasm"],
            ancillas=ancillas,
            systems=systems,
            expected_polynomial=file_record["polynomial"],
            hermitian_base=True,
            compare_polynomial_only=True,
        )

        assert result.status == PASS
        assert result.qsp_polynomial_only is True
        assert polynomial_close(result.qsp_polynomial, file_record["polynomial"])

    assert seen_components == {"cos", "sin"}
