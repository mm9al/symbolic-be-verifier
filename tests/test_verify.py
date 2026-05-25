from pathlib import Path

import pytest
import sympy as sp

from symbolic.branch_state import BranchState, Gate, UnsupportedGateError
from symbolic.expr import pauli
from symbolic.qasm_parser import parse_qasm_file, parse_qasm_text
from symbolic.verify import format_result, verify_qasm_file


def test_h_on_ancilla_update():
    state = BranchState.initial().apply(Gate("h", (0,)))
    expected = pauli("I").scale(sp.sqrt(2) / 2)

    assert state.b0.equals(expected)
    assert state.b1.equals(expected)


def test_lcu_x_plus_z_example_verifies():
    path = Path(__file__).parents[1] / "examples" / "lcu_x_plus_z.qasm"
    result = verify_qasm_file(path, expected="(X + Z)/2", keep_trace=True)

    assert result.success is True
    assert result.final_state.b0.equals((pauli("X") + pauli("Z")).scale(sp.Rational(1, 2)))
    assert len(result.trace) == len(parse_qasm_file(path)) + 1


def test_trace_output_can_show_steps():
    path = Path(__file__).parents[1] / "examples" / "lcu_x_plus_z.qasm"
    result = verify_qasm_file(path, expected="(X + Z)/2", keep_trace=True)
    output = format_result(result, show_trace=True)

    assert "Initial" in output
    assert "h q[0]" in output
    assert "Final B0 = (X + Z)/2" in output
    assert "PASS" in output


def test_rotation_is_parsed_but_not_supported():
    gates = parse_qasm_text(
        """
        OPENQASM 2.0;
        qreg q[2];
        rx(pi/2) q[0];
        """
    )

    with pytest.raises(UnsupportedGateError, match="Parameterized gate"):
        BranchState.initial().apply(gates[0])
