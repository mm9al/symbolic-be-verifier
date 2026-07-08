from pathlib import Path

import pytest
import sympy as sp

from symbolic.branch_state import BranchState, Gate, UnsupportedGateError
from symbolic.expr import parse_operator_expression, pauli
from symbolic.qasm_parser import parse_qasm_file, parse_qasm_text
from symbolic.verify import (
    FAIL,
    FAIL_GARBAGE,
    PASS,
    PASS_UP_TO_SCALE,
    VerificationResult,
    format_gate_profiles_csv,
    format_result,
    pauli_expr_close,
    polynomial_close,
    proportional_scale,
    scalar_close,
    verify_polynomial_approximates_exp,
    verify_qasm_file,
)
from symbolic.word_expr import WordExpr
from symbolic.word_expr import atom as word_atom


def test_h_on_ancilla_update():
    state = BranchState.initial().apply(Gate("h", (0,)))
    expected = pauli("I").scale(sp.sqrt(2) / 2)

    assert state.b0.equals(expected)
    assert state.b1.equals(expected)


def test_lcu_x_minus_z_example_verifies():
    path = Path(__file__).parents[1] / "examples" / "lcu_x_minus_z.qasm"
    result = verify_qasm_file(path, expected="(X - Z)/2", keep_trace=True)

    assert result.success is True
    assert result.status == PASS
    assert result.final_state.b0.equals((pauli("X") - pauli("Z")).scale(sp.Rational(1, 2)))
    assert len(result.trace) == len(parse_qasm_file(path)) + 1


def test_gate_profile_records_branch_and_term_counters():
    result = verify_qasm_file(
        Path(__file__).parents[1] / "examples" / "lcu_x_minus_z.qasm",
        expected="(X - Z)/2",
        profile_gates=True,
    )
    output = format_gate_profiles_csv(result)

    assert result.gate_profiles
    assert "gate_id,gate_name,gate,num_nonzero_branches,total_operator_terms,max_terms_per_branch" in output
    assert "time_this_gate,time_simplify,time_combine" in output
    assert result.gate_profiles[0].num_nonzero_branches >= 1


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
    path = Path(__file__).parents[1] / "examples" / "lcu_x_minus_z.qasm"
    result = verify_qasm_file(path, expected="(X - Z)/2", keep_trace=True)
    output = format_result(result, show_trace=True)

    assert "Initial" in output
    assert "h q[0]" in output
    assert "Final B0 = X/2 - Z/2" in output
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


def test_proportional_scale_detects_unit_global_phase():
    expected = (pauli("X") + pauli("Z")).scale(sp.Rational(1, 2))
    scale = (-1 - sp.I) / sp.sqrt(2)
    actual = expected.scale(scale)

    ok, found_scale = proportional_scale(actual, expected)

    assert ok is True
    assert sp.simplify(found_scale - scale) == 0


def test_verification_status_pass_up_to_scale_for_unit_global_phase():
    expected = (pauli("X") + pauli("Z")).scale(sp.Rational(1, 2))
    scale = (-1 - sp.I) / sp.sqrt(2)
    actual = expected.scale(scale)
    result = VerificationResult(final_state=BranchState(1, 1, {(0,): actual}), trace=[], expected=expected)
    output = format_result(result)

    assert result.success is True
    assert result.status == PASS_UP_TO_SCALE
    assert sp.simplify(result.scale - scale) == 0
    assert "PASS_UP_TO_SCALE" in output
    assert "scale = sqrt(2)*(-1 - i)/2" in output


def test_proportional_non_unit_scale_passes():
    expected = pauli("X") + pauli("Z")
    actual = expected.scale(2)
    result = VerificationResult(final_state=BranchState(1, 1, {(0,): actual}), trace=[], expected=expected)

    ok, scale = proportional_scale(actual, expected)

    assert ok is True
    assert scale == 2
    assert result.success is True
    assert result.status == PASS_UP_TO_SCALE


def test_numeric_tolerance_accepts_exponential_and_decimal_coefficients():
    actual_coeff = sp.exp(sp.I * sp.Rational(1160258681, 10**12))
    expected_coeff = sp.Float("0.9999993269", 20) + sp.I * sp.Float("0.001160258420", 20)
    actual = pauli("I").scale(actual_coeff)
    expected = pauli("I").scale(expected_coeff)

    assert actual.equals(expected) is False
    assert scalar_close(actual_coeff, expected_coeff)
    assert pauli_expr_close(actual, expected)

    result = VerificationResult(final_state=BranchState(1, 1, {(0,): actual}), trace=[], expected=expected)

    assert result.success is True
    assert result.status == PASS


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


def test_multi_ancilla_h_gates_create_sparse_bitstring_branches():
    state = BranchState.initial(num_system_qubits=1, num_ancillas=2)
    state = state.apply(Gate("h", (0,)), ancillas=(0, 1), systems=(2,))
    state = state.apply(Gate("h", (1,)), ancillas=(0, 1), systems=(2,))

    expected = pauli("I").scale(sp.Rational(1, 2))

    assert set(state.branches) == {(0, 0), (0, 1), (1, 0), (1, 1)}
    assert all(branch.equals(expected) for branch in state.branches.values())
    assert state.top_left().equals(expected)


def test_multi_ancilla_trace_output_lists_branch_keys():
    state = BranchState.initial(num_system_qubits=1, num_ancillas=2)
    result = VerificationResult(final_state=state, trace=[], expected=pauli("I"), ancilla_qubits=(0, 1), system_qubits=(2,))
    output = format_result(result, show_trace=True)

    assert "Ancilla qubits:" in output
    assert "q[1] -> ancilla[1]" in output
    assert "System qubits:" in output
    assert "q[2] -> system[0]" in output
    assert "Final B[00] = I" in output


def test_cx_ancilla_to_ancilla_permutes_branch_keys():
    state = BranchState(2, 1, {(1, 0): pauli("I")})
    state = state.apply(Gate("cx", (0, 1)), ancillas=(0, 1), systems=(2,))

    assert set(state.branches) == {(1, 1)}
    assert state.branch((1, 1)).equals(pauli("I"))


def test_cz_ancilla_to_ancilla_adds_minus_phase_on_11():
    state = BranchState(2, 1, {(1, 1): pauli("I"), (1, 0): pauli("X")})
    state = state.apply(Gate("cz", (0, 1)), ancillas=(0, 1), systems=(2,))

    assert state.branch((1, 1)).equals(pauli("I").scale(-1))
    assert state.branch((1, 0)).equals(pauli("X"))


def test_cx_ancilla_to_system_conditionally_left_multiplies_x():
    state = BranchState(2, 1, {(0, 0): pauli("I"), (1, 0): pauli("Z")})
    state = state.apply(Gate("cx", (0, 2)), ancillas=(0, 1), systems=(2,))

    assert state.branch((0, 0)).equals(pauli("I"))
    assert state.branch((1, 0)).equals(pauli("X") * pauli("Z"))


def test_cz_system_ancilla_is_symmetric_conditionally_left_multiplies_z():
    state = BranchState(2, 1, {(0, 1): pauli("X"), (0, 0): pauli("I")})
    state = state.apply(Gate("cz", (2, 1)), ancillas=(0, 1), systems=(2,))

    assert state.branch((0, 0)).equals(pauli("I"))
    assert state.branch((0, 1)).equals(pauli("Z") * pauli("X"))


def test_cx_system_to_ancilla_uses_projector_rule():
    state = BranchState(1, 1, {(0,): pauli("I")})
    state = state.apply(Gate("cx", (1, 0)), ancillas=(0,), systems=(1,))

    expected_stay = (pauli("I") + pauli("Z")).scale(sp.Rational(1, 2))
    expected_flip = (pauli("I") - pauli("Z")).scale(sp.Rational(1, 2))

    assert state.branch((0,)).equals(expected_stay)
    assert state.branch((1,)).equals(expected_flip)


def test_system_system_controlled_gates_are_unsupported_in_v04():
    state = BranchState.initial(num_system_qubits=2, num_ancillas=1)

    with pytest.raises(UnsupportedGateError, match="system-system cx"):
        state.apply(Gate("cx", (1, 2)), ancillas=(0,), systems=(1, 2))

    with pytest.raises(UnsupportedGateError, match="system-system cz"):
        state.apply(Gate("cz", (1, 2)), ancillas=(0,), systems=(1, 2))


def test_uhdg_uh_opaque_reduces_to_identity():
    path = Path(__file__).parents[1] / "examples" / "uhdg_uh_opaque.qasm"
    result = verify_qasm_file(
        path,
        ancillas=(0, 1),
        systems=(2,),
        base="(X - Z)/2",
        expected_polynomial="1",
    )

    assert result.success is True
    assert result.status == PASS
    assert result.qsp_normalized == WordExpr.identity()
    assert result.qsp_polynomial == 1


def test_qsp_polynomial_check_requires_hermitian_base_to_rewrite_hdg(tmp_path):
    path = tmp_path / "uhdg_only.qasm"
    path.write_text(
        """
        OPENQASM 2.0;
        opaque UHdg a, s;
        qreg q[2];
        UHdg q[0], q[1];
        """,
        encoding="utf-8",
    )

    result = verify_qasm_file(
        path,
        ancillas=(0,),
        systems=(1,),
        expected_polynomial="x",
        compare_polynomial_only=True,
    )
    hermitian_result = verify_qasm_file(
        path,
        ancillas=(0,),
        systems=(1,),
        expected_polynomial="x",
        hermitian_base=True,
        compare_polynomial_only=True,
    )

    assert result.status == FAIL_GARBAGE
    assert result.qsp_normalized == word_atom("Hd")
    assert hermitian_result.status == PASS
    assert hermitian_result.qsp_normalized == word_atom("H")


def test_controlled_uh_only_updates_selected_selector_branch():
    state = BranchState(
        2,
        1,
        {
            (0, 0): WordExpr.identity(),
            (0, 1): word_atom("H"),
            (1, 0): WordExpr.identity(),
            (1, 1): word_atom("H"),
        },
        expression_kind="word",
    )

    state = state.apply(Gate("cUH", (0, 1, 2)), ancillas=(0, 1), systems=(2,))

    assert state.branch((0, 0)).equals(WordExpr.identity())
    assert state.branch((0, 1)).equals(word_atom("H"))
    assert state.branch((1, 0)).equals(word_atom("H") + word_atom("G") * word_atom("H"))
    assert state.branch((1, 1)).equals(word_atom("A") + word_atom("C") * word_atom("H"))


def test_controlled_uhdg_uses_dagger_blocks_on_selected_selector_branch():
    state = BranchState(
        2,
        1,
        {
            (0, 0): WordExpr.identity(),
            (1, 0): WordExpr.identity(),
            (1, 1): word_atom("H"),
        },
        expression_kind="word",
    )

    state = state.apply(Gate("cUHdg", (0, 1, 2)), ancillas=(0, 1), systems=(2,))

    assert state.branch((0, 0)).equals(WordExpr.identity())
    assert state.branch((1, 0)).equals(word_atom("Hd") + word_atom("Ad") * word_atom("H"))
    assert state.branch((1, 1)).equals(word_atom("Gd") + word_atom("Cd") * word_atom("H"))


def test_parser_accepts_controlled_opaque_uh_declarations():
    gates = parse_qasm_text(
        """
        OPENQASM 2.0;
        opaque cUH c, a, s;
        opaque cUHdg c, a, s;
        qreg q[3];
        cUH q[0], q[1], q[2];
        cUHdg q[0], q[1], q[2];
        """
    )

    assert [gate.name for gate in gates] == ["cuh", "cuhdg"]
    assert gates[0].qubits == (0, 1, 2)


def test_qsp_t3_opaque_verifies_chebyshev_polynomial_on_hermitian_base():
    path = Path(__file__).parents[1] / "examples" / "qsp_t3_opaque.qasm"
    result = verify_qasm_file(
        path,
        ancillas=(0, 1),
        systems=(2,),
        base="(X - Z)/2",
        expected_polynomial="4*x^3 - 3*x",
        hermitian_base=True,
        keep_trace=True,
    )

    expected_operator = (pauli("X") - pauli("Z")).scale(sp.Rational(-1, 2))
    expected_polynomial = 4 * sp.Symbol("x") ** 3 - 3 * sp.Symbol("x")

    assert result.success is True
    assert result.status == PASS
    assert result.qsp_normalized == WordExpr({("H", "H", "H"): 4, ("H",): -3})
    assert sp.expand(result.qsp_polynomial - expected_polynomial) == 0
    assert result.qsp_actual.equals(expected_operator)
    assert "Normalized B[00] = -3*H + 4*H H H" in format_result(result, show_trace=False)


def test_qsp_polynomial_only_check_does_not_require_base():
    path = Path(__file__).parents[1] / "examples" / "qsp_t3_opaque.qasm"
    result = verify_qasm_file(
        path,
        ancillas=(0, 1),
        systems=(2,),
        expected_polynomial="4*x^3 - 3*x",
        hermitian_base=True,
        compare_polynomial_only=True,
    )

    assert result.success is True
    assert result.status == PASS
    assert result.qsp_polynomial_only is True
    assert result.qsp_actual is None
    assert result.qsp_expected is None


def test_qsp_polynomial_can_be_extracted_without_verification_route():
    path = Path(__file__).parents[1] / "examples" / "qsp_t3_opaque.qasm"
    result = verify_qasm_file(
        path,
        ancillas=(0, 1),
        systems=(2,),
        hermitian_base=True,
        extract_qsp_polynomial=True,
    )

    assert result.status is None
    assert result.success is None
    assert sp.expand(result.qsp_polynomial - (4 * sp.Symbol("x") ** 3 - 3 * sp.Symbol("x"))) == 0


def test_qsp_target_exp_check_can_drive_status_without_expected_polynomial():
    path = Path(__file__).parents[1] / "examples" / "uhdg_uh_opaque.qasm"
    result = verify_qasm_file(
        path,
        ancillas=(0, 1),
        systems=(2,),
        target_exp_tau=0.0,
        target_exp_epsilon=1e-8,
        hermitian_base=True,
    )

    assert result.status == PASS
    assert result.qsp_polynomial == 1
    assert result.qsp_approximation is not None
    assert result.qsp_approximation.success is True


def test_qsp_target_exp_grid_uses_markov_degree_bound():
    check = verify_polynomial_approximates_exp(
        "1 - i*x + 0.5*x^3",
        tau=0.5,
        epsilon=0.1,
        scale="0.25",
    )

    assert check.polynomial_degree == 3
    assert check.polynomial_derivative_bound == 3
    assert check.target_lipschitz == pytest.approx(0.125)
    assert check.num_grid_points == 100
    assert check.spacing == pytest.approx(float(sp.pi / 99))


def test_qsp_target_exp_grid_guard_fails_fast():
    with pytest.raises(ValueError, match="Approximation grid requires 100 Chebyshev points"):
        verify_polynomial_approximates_exp(
            "1 - i*x + 0.5*x^3",
            tau=0.5,
            epsilon=0.1,
            scale="0.25",
            max_grid_points=99,
        )


def test_qsp_target_exp_check_fails_when_grid_error_exceeds_half_epsilon():
    path = Path(__file__).parents[1] / "examples" / "qsp_t3_opaque.qasm"
    result = verify_qasm_file(
        path,
        ancillas=(0, 1),
        systems=(2,),
        target_exp_tau=0.0,
        target_exp_epsilon=0.1,
        hermitian_base=True,
    )

    assert result.status == FAIL
    assert result.qsp_approximation is not None
    assert result.qsp_approximation.max_grid_error > 0.05


def test_verification_routes_are_mutually_exclusive():
    path = Path(__file__).parents[1] / "examples" / "qsp_t3_opaque.qasm"

    with pytest.raises(ValueError, match="Choose exactly one verification route"):
        verify_qasm_file(
            path,
            ancillas=(0, 1),
            systems=(2,),
            expected_polynomial="4*x^3 - 3*x",
            target_exp_tau=0.0,
            target_exp_epsilon=0.1,
            hermitian_base=True,
        )
