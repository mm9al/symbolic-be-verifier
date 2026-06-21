import sympy as sp

from symbolic.qasm_parser import parse_qasm_text


def test_parse_qasm2_gate_subset_and_spacing():
    gates = parse_qasm_text(
        """
        OPENQASM 2.0;
        include "qelib1.inc";
        qreg q[2];

        h q[0];
        cx q[0],q[1];
        cz q[0], q[1];
        """
    )

    assert [gate.name for gate in gates] == ["h", "cx", "cz"]
    assert [gate.qubits for gate in gates] == [(0,), (0, 1), (0, 1)]


def test_parse_qasm3_register_and_parameterized_gate():
    gates = parse_qasm_text(
        """
        OPENQASM 3;
        include "stdgates.inc";
        qubit[2] qb;
        rx(pi*0.5) qb[0];
        """
    )

    assert len(gates) == 1
    assert gates[0].name == "rx"
    assert gates[0].qubits == (0,)
    assert sp.simplify(gates[0].parameter - sp.pi / 2) == 0


def test_parse_symbolic_rotation_parameter():
    gates = parse_qasm_text(
        """
        OPENQASM 2.0;
        qreg q[2];
        rz(-theta) q[0];
        """
    )

    assert gates[0].name == "rz"
    assert gates[0].parameter == -sp.Symbol("theta")


def test_parse_opaque_declarations_are_ignored_and_uh_gates_are_kept():
    gates = parse_qasm_text(
        """
        OPENQASM 2.0;
        qreg q[3];
        opaque UH a, s;
        opaque UHdg a, s;
        UH q[1], q[2];
        UHdg q[1], q[2];
        """
    )

    assert [gate.name for gate in gates] == ["uh", "uhdg"]
    assert [gate.qubits for gate in gates] == [(1, 2), (1, 2)]


def test_parse_abstract_mcx_with_multiple_controls():
    gates = parse_qasm_text(
        """
        OPENQASM 2.0;
        qreg q[4];
        mcx q[1], q[2], q[0];
        """
    )

    assert len(gates) == 1
    assert gates[0].name == "mcx"
    assert gates[0].qubits == (1, 2, 0)
    assert gates[0].parameter is None
