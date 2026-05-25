from pathlib import Path

import sympy as sp

from symbolic.branch_state import BranchState, Gate
from symbolic.expr import parse_operator_expression, pauli
from symbolic.qasm_parser import parse_qasm_file, parse_qasm_text
from symbolic.verify import FAIL, PASS, PASS_UP_TO_GLOBAL_PHASE, VerificationResult, format_result, proportional_phase, verify_qasm_file


def test_h_on_ancilla_update():
    state = BranchState.initial().apply(Gate("h", (0,)))
    expected = pauli("I").scale(sp.sqrt(2) / 2)

    assert state.b0.equals(expected)
    assert state.b1.equals(expected)


def test_lcu_x_plus_z_example_verifies():
    path = Path(__file__).parents[1] / "examples" / "lcu_x_plus_z.qasm"
    result = verify_qasm_file(path, expected="(X + Z)/2", keep_trace=True)

    assert result.success is True
    assert result.status == PASS
    assert result.final_state.b0.equals((pauli("X") + pauli("Z")).scale(sp.Rational(1, 2)))
    assert len(result.trace) == len(parse_qasm_file(path)) + 1


def test_lcu_xx_plus_zz_example_verifies():
    path = Path(__file__).parents[1] / "examples" / "lcu_xx_plus_zz.qasm"
    result = verify_qasm_file(path, expected="(X0 X1 + Z0 Z1)/2", systems=(1, 2), keep_trace=True)
    expected = (
        pauli("X", index=0, num_qubits=2) * pauli("X", index=1, num_qubits=2)
        + pauli("Z", index=0, num_qubits=2) * pauli("Z", index=1, num_qubits=2)
    ).scale(sp.Rational(1, 2))

    assert result.success is True
    assert result.status == PASS
    assert result.final_state.b0.equals(expected)
    assert str(result.final_state.b0) == "(X_0 X_1 + Z_0 Z_1)/2"
    assert len(result.trace) == len(parse_qasm_file(path)) + 1


def test_trace_output_can_show_steps():
    path = Path(__file__).parents[1] / "examples" / "lcu_x_plus_z.qasm"
    result = verify_qasm_file(path, expected="(X + Z)/2", keep_trace=True)
    output = format_result(result, show_trace=True)

    assert "Initial" in output
    assert "h q[0]" in output
    assert "Final B0 = (X + Z)/2" in output
    assert "PASS" in output


def test_system_rx_pi_is_supported():
    path = Path(__file__).parents[1] / "examples" / "rx_system_pi.qasm"
    result = verify_qasm_file(path, expected="-i X")

    assert result.success is True
    assert result.final_state.b0.equals(pauli("X").scale(-sp.I))
    assert result.final_state.b1.is_zero()


def test_ancilla_rz_theta_sandwich_is_supported():
    path = Path(__file__).parents[1] / "examples" / "rz_ancilla_theta_sandwich.qasm"
    theta = sp.Symbol("theta")
    result = verify_qasm_file(path, expected="cos(theta/2) I", keep_trace=True)

    assert result.success is True
    assert result.final_state.b0.equals(parse_operator_expression("cos(theta/2) I"))
    assert result.final_state.b1.equals(pauli("I").scale(-sp.I * sp.sin(theta / 2)))


def test_system_rotation_update_rules_all_axes():
    theta = sp.Symbol("theta")
    c = sp.cos(theta / 2)
    minus_i_s = -sp.I * sp.sin(theta / 2)

    rx_state = BranchState.initial(2).apply(Gate("rx", (2,), parameter=theta), systems=(1, 2))
    ry_state = BranchState.initial(2).apply(Gate("ry", (2,), parameter=theta), systems=(1, 2))
    rz_state = BranchState.initial(2).apply(Gate("rz", (2,), parameter=theta), systems=(1, 2))

    assert rx_state.b0.equals(pauli("I", num_qubits=2).scale(c) + pauli("X", index=1, num_qubits=2).scale(minus_i_s))
    assert ry_state.b0.equals(pauli("I", num_qubits=2).scale(c) + pauli("Y", index=1, num_qubits=2).scale(minus_i_s))
    assert rz_state.b0.equals(pauli("I", num_qubits=2).scale(c) + pauli("Z", index=1, num_qubits=2).scale(minus_i_s))


def test_ancilla_rz_uses_exponential_branch_phases():
    theta = sp.Symbol("theta")
    state = BranchState.initial().apply(Gate("rz", (0,), parameter=theta))

    assert state.b0.equals(pauli("I").scale(sp.exp(-sp.I * theta / 2)))
    assert state.b1.is_zero()


def test_ancilla_s_and_sdg_update_rules():
    prepared = BranchState.initial().apply(Gate("h", (0,)))

    s_state = prepared.apply(Gate("s", (0,)))
    sdg_state = prepared.apply(Gate("sdg", (0,)))

    assert s_state.b0.equals(pauli("I").scale(sp.sqrt(2) / 2))
    assert s_state.b1.equals(pauli("I").scale(sp.I * sp.sqrt(2) / 2))
    assert sdg_state.b0.equals(pauli("I").scale(sp.sqrt(2) / 2))
    assert sdg_state.b1.equals(pauli("I").scale(-sp.I * sp.sqrt(2) / 2))


def test_system_s_and_sdg_update_rules():
    s_state = BranchState.initial().apply(Gate("s", (1,)))
    sdg_state = BranchState.initial().apply(Gate("sdg", (1,)))

    expected_s = pauli("I").scale((1 + sp.I) / 2) + pauli("Z").scale((1 - sp.I) / 2)
    expected_sdg = pauli("I").scale((1 - sp.I) / 2) + pauli("Z").scale((1 + sp.I) / 2)

    assert s_state.b0.equals(expected_s)
    assert s_state.b1.is_zero()
    assert sdg_state.b0.equals(expected_sdg)
    assert sdg_state.b1.is_zero()


def test_system_s_uses_multi_system_index_mapping():
    state = BranchState.initial(2).apply(Gate("s", (2,)), systems=(1, 2))
    expected = pauli("I", num_qubits=2).scale((1 + sp.I) / 2) + pauli("Z", index=1, num_qubits=2).scale((1 - sp.I) / 2)

    assert state.b0.equals(expected)
    assert str(state.b0) == "(1/2 + i/2)*I + (1/2 - i/2)*Z_1"


def test_proportional_phase_detects_unit_global_phase():
    expected = (pauli("X") + pauli("Z")).scale(sp.Rational(1, 2))
    phase = (-1 - sp.I) / sp.sqrt(2)
    actual = expected.scale(phase)

    ok, found_phase = proportional_phase(actual, expected)

    assert ok is True
    assert sp.simplify(found_phase - phase) == 0


def test_verification_status_pass_up_to_global_phase():
    expected = (pauli("X") + pauli("Z")).scale(sp.Rational(1, 2))
    phase = (-1 - sp.I) / sp.sqrt(2)
    actual = expected.scale(phase)
    result = VerificationResult(final_state=BranchState(actual, pauli("I").scale(0)), trace=[], expected=expected)
    output = format_result(result)

    assert result.success is True
    assert result.status == PASS_UP_TO_GLOBAL_PHASE
    assert sp.simplify(result.global_phase - phase) == 0
    assert "PASS_UP_TO_GLOBAL_PHASE" in output
    assert "phase = sqrt(2)*(-1 - i)/2" in output


def test_proportional_but_non_unit_scale_fails():
    expected = pauli("X") + pauli("Z")
    actual = expected.scale(2)
    result = VerificationResult(final_state=BranchState(actual, pauli("I").scale(0)), trace=[], expected=expected)

    ok, phase = proportional_phase(actual, expected)

    assert ok is False
    assert phase == 2
    assert result.success is False
    assert result.status == FAIL


def test_ancilla_rx_and_ry_update_rules():
    gates = parse_qasm_text(
        """
        OPENQASM 2.0;
        qreg q[2];
        rx(theta) q[0];
        ry(-theta) q[0];
        """
    )
    theta = sp.Symbol("theta")

    state = BranchState.initial().apply(gates[0]).apply(gates[1])

    assert state.b0.equals(pauli("I").scale(sp.cos(theta / 2) ** 2 - sp.I * sp.sin(theta / 2) ** 2))
    assert state.b1.equals(pauli("I").scale(-sp.sin(theta / 2) * sp.cos(theta / 2) - sp.I * sp.sin(theta / 2) * sp.cos(theta / 2)))
