import sympy as sp

from symbolic.verify import normalize_qsp_expression, verify_qasm_file
from symbolic.word_expr import WordExpr


def test_hermitian_replacement_adds_colliding_words():
    expr = WordExpr({("H",): 1 + 2 * sp.I, ("Hd",): 1 - 2 * sp.I})
    normalized = normalize_qsp_expression(expr, hermitian_base=True)
    assert normalized.terms == {("H",): 2}


def test_hermitian_replacement_cancels_colliding_words():
    expr = WordExpr({("H", "Hd"): 3, ("Hd", "H"): -3})
    assert normalize_qsp_expression(expr, hermitian_base=True).is_zero()


def test_real_extraction_preserves_full_linear_polynomial(tmp_path):
    path = tmp_path / "real_extraction.qasm"
    path.write_text("""OPENQASM 2.0;
include "qelib1.inc";
opaque cUH c,a,s;
opaque cUHdg c,a,s;
qreg q[3];
h q[0];
cUH q[0],q[1],q[2];
x q[0];
cUHdg q[0],q[1],q[2];
x q[0];
h q[0];
""")
    result = verify_qasm_file(path, ancillas=[0, 1], systems=[2],
                              expected_polynomial="x", compare_polynomial_only=True,
                              hermitian_base=True)
    assert result.status == "PASS"
