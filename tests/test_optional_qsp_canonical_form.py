from pathlib import Path

import pytest

from symbolic.branch_state import UnsupportedGateError
from symbolic.verify import main, verify_qasm_file


@pytest.fixture
def noncanonical_qasm(tmp_path):
    path = tmp_path / "noncanonical.qasm"
    # Z after UH does not change the selected zero block H, but addresses one
    # individual block ancilla and is outside the canonical input grammar.
    path.write_text("""OPENQASM 2.0;
qreg q[3];
UH q[0], q[1], q[2];
z q[1];
""")
    return path


def test_canonical_check_remains_enabled_by_default(noncanonical_qasm):
    with pytest.raises(UnsupportedGateError, match="distinguishes individual nonzero"):
        verify_qasm_file(noncanonical_qasm, ancillas=(0, 1), systems=(2,))


def test_api_can_skip_only_the_canonical_check(noncanonical_qasm):
    result = verify_qasm_file(
        noncanonical_qasm, ancillas=(0, 1), systems=(2,),
        expected_polynomial="x", hermitian_base=True, compare_polynomial_only=True,
        check_qsp_canonical_form=False,
    )
    assert result.status == "PASS"


def test_cli_can_skip_the_canonical_check(noncanonical_qasm, capsys):
    status = main([
        str(noncanonical_qasm), "--ancillas", "0", "1", "--systems", "2",
        "--expected-polynomial=x", "--hermitian-base", "--compare-polynomial-only",
        "--skip-qsp-canonical-check", "--result-only",
    ])
    assert status == 0
    assert capsys.readouterr().out == "PASS\n"


def test_skipping_canonical_check_still_rejects_unsupported_gates(noncanonical_qasm):
    with noncanonical_qasm.open("a") as output:
        output.write("unsupported_gate q[0];\n")
    with pytest.raises(UnsupportedGateError, match="Unsupported ancilla gate"):
        verify_qasm_file(
            noncanonical_qasm, ancillas=(0, 1), systems=(2,),
            check_qsp_canonical_form=False,
        )


def test_canonical_polynomial_is_identical_with_check_on_or_off():
    path = Path(__file__).parents[1] / "examples" / "qsp_t3_mcx_m2.qasm"
    options = dict(ancillas=(0, 1, 2), systems=(3,),
                   expected_polynomial="4*x^3-3*x", hermitian_base=True,
                   compare_polynomial_only=True)
    checked = verify_qasm_file(path, **options)
    unchecked = verify_qasm_file(path, check_qsp_canonical_form=False, **options)
    assert checked.status == unchecked.status == "PASS"
    assert checked.qsp_polynomial == unchecked.qsp_polynomial
