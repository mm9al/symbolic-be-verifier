import sympy as sp

from symbolic.expr import parse_operator_expression, pauli


def test_pauli_multiplication_rules():
    assert (pauli("X") * pauli("X")).equals(pauli("I"))
    assert (pauli("Y") * pauli("Y")).equals(pauli("I"))
    assert (pauli("Z") * pauli("Z")).equals(pauli("I"))

    assert (pauli("X") * pauli("Y")).equals(pauli("Z").scale(sp.I))
    assert (pauli("Y") * pauli("Z")).equals(pauli("X").scale(sp.I))
    assert (pauli("Z") * pauli("X")).equals(pauli("Y").scale(sp.I))

    assert (pauli("Y") * pauli("X")).equals(pauli("Z").scale(-sp.I))
    assert (pauli("Z") * pauli("Y")).equals(pauli("X").scale(-sp.I))
    assert (pauli("X") * pauli("Z")).equals(pauli("Y").scale(-sp.I))


def test_parse_operator_expression():
    parsed = parse_operator_expression("(X + Z)/2")
    expected = (pauli("X") + pauli("Z")).scale(sp.Rational(1, 2))

    assert parsed.equals(expected)


def test_parse_multiplied_paulis_uses_noncommutative_rules():
    assert parse_operator_expression("ZX").equals(pauli("Y").scale(sp.I))
    assert parse_operator_expression("XZ").equals(pauli("Y").scale(-sp.I))


def test_multi_qubit_pauli_string_multiplication():
    left = pauli("X", index=0, num_qubits=2) * pauli("Z", index=1, num_qubits=2)
    right = pauli("Z", index=0, num_qubits=2) * pauli("X", index=1, num_qubits=2)

    assert (left * right).equals((pauli("Y", 0, 2) * pauli("Y", 1, 2)).scale(1))


def test_parse_multi_qubit_operator_expression():
    parsed = parse_operator_expression("(X0 X1 + Z_0 Z_1)/2", num_qubits=2)
    expected = (
        pauli("X", index=0, num_qubits=2) * pauli("X", index=1, num_qubits=2)
        + pauli("Z", index=0, num_qubits=2) * pauli("Z", index=1, num_qubits=2)
    ).scale(sp.Rational(1, 2))

    assert parsed.equals(expected)
    assert str(parsed) == "(X_0 X_1 + Z_0 Z_1)/2"
