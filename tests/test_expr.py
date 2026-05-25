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
