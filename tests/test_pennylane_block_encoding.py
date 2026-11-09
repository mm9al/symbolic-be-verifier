"""Native inputs and rigorous exact/tolerant Hermiticity decisions."""
import json
from pathlib import Path
import shutil

import pytest

from evaluation.pennylane import common, generate_benchmarks, run_block_encoding


def row_for(n):
    # Users can regenerate the sweep manifest with a different size subset.
    # Small saved artifacts remain independent regression inputs.
    directory = common.BENCHMARKS/"ising"/f"ising_periodic_n{n}"
    metadata_path = directory/"metadata.json"
    metadata = json.loads(metadata_path.read_text())
    keys = ("benchmark_id", "model", "graph", "n_system", "n_terms", "target_lcu_terms",
            "n_ancilla", "nominal_alpha", "gate_count")
    return dict(**{key: str(metadata[key]) for key in keys},
                qasm_path=str((directory/"block_encoding.qasm").relative_to(common.BENCHMARKS)),
                metadata_path=str(metadata_path.relative_to(common.BENCHMARKS)),
                metadata_sha256=common.sha256(metadata_path))


def require_verifier():
    pytest.importorskip("sympy")
    pytest.importorskip("flint")


def test_native_two_qubit_task1_accepts(tmp_path):
    require_verifier()
    result = run_block_encoding.evaluate(common.BENCHMARKS, row_for(2), tmp_path, "1e-8", 100)
    assert result["status"] == "ACCEPT", result
    assert result["proportionality_pass"]
    assert result["hermiticity_status"] == "PASS_STRUCTURE"
    assert abs(float(result["inferred_alpha"])-4) < 1e-12


@pytest.mark.parametrize("delta,status,hermiticity", [
    ("1e-12", "ACCEPT", "PASS_WITHIN_DELTA"),
    ("0", "REJECT", "FAIL_NON_HERMITIAN"),
])
def test_native_three_qubit_tolerant_and_exact_modes(tmp_path, delta, status, hermiticity):
    require_verifier()
    result = run_block_encoding.evaluate(common.BENCHMARKS, row_for(3), tmp_path, "1e-8", 100, delta)
    assert result["status"] == status, result
    assert result["proportionality_pass"]
    assert result["hermiticity_status"] == hermiticity
    assert result["exact_hermiticity_status"] == "FAIL_NON_HERMITIAN"
    assert result["hermiticity_delta"] == delta
    report = json.loads((tmp_path/result["report_path"]).read_text())
    assert report["hermiticity"]["exact"]["violations"]
    assert report["hermiticity"]["approximate_acceptance"] == (delta != "0")


def test_exact_hermiticity_does_not_override_failed_proportionality(tmp_path):
    require_verifier()
    result = run_block_encoding.evaluate(common.BENCHMARKS, row_for(2), tmp_path, "1e-18", 100)
    assert result["status"] == "REJECT", result
    assert not result["proportionality_pass"]
    assert result["hermiticity_status"] == "PASS_STRUCTURE"


def test_zero_containing_imaginary_ball_is_not_a_hermiticity_proof():
    require_verifier()
    from symbolic.expr import OpExpr
    from symbolic.interval import acb, arb
    actual = OpExpr({("X",): acb(1, arb("0 +/- 1e-40"))})
    result = run_block_encoding.classify_hermiticity(actual, dict(status="UNPROVEN"), "0")
    assert result["status"] == "UNPROVEN"


@pytest.mark.parametrize("imaginary,expected", [
    ("1/2048", "PASS_WITHIN_DELTA"),
    ("1/1024", "FAIL_DELTA"),
    ("1/512", "FAIL_DELTA"),
])
def test_delta_comparison_is_strict(imaginary, expected):
    require_verifier()
    from symbolic.expr import OpExpr
    from symbolic.interval import acb, arb
    actual = OpExpr({("X",): acb(1, arb(imaginary))})
    result = run_block_encoding.classify_hermiticity(actual, dict(status="UNPROVEN"), "1/1024")
    assert result["status"] == expected
    assert result["exact_status"] == "FAIL_NON_HERMITIAN"


def test_delta_uses_interval_upper_bound_instead_of_midpoint():
    require_verifier()
    from symbolic.expr import OpExpr
    from symbolic.interval import acb, arb
    actual = OpExpr({("X",): acb(1, arb("0 +/- 2e-12"))})
    result = run_block_encoding.classify_hermiticity(actual, dict(status="UNPROVEN"), "1e-12")
    assert result["status"] == "FAIL_DELTA"
    assert result["exact_status"] == "UNPROVEN"


def test_small_zero_containing_interval_can_pass_delta():
    require_verifier()
    from symbolic.expr import OpExpr
    from symbolic.interval import acb, arb
    actual = OpExpr({("X",): acb(1, arb("0 +/- 1e-14"))})
    result = run_block_encoding.classify_hermiticity(actual, dict(status="UNPROVEN"), "1e-12")
    assert result["status"] == "PASS_WITHIN_DELTA"
    assert result["exact_status"] == "UNPROVEN"


@pytest.mark.parametrize("delta", ["-1e-12", "nan", "inf"])
def test_invalid_delta_is_rejected_by_cli(monkeypatch, delta):
    monkeypatch.setattr("sys.argv", ["run_block_encoding.py", "--hermiticity-delta", delta])
    with pytest.raises(SystemExit) as exc:
        run_block_encoding.main()
    assert exc.value.code == 2


def test_changed_qasm_is_rejected_before_propagation(tmp_path):
    row = row_for(2)
    source = common.BENCHMARKS/Path(row["metadata_path"]).parent
    destination = tmp_path/Path(row["metadata_path"]).parent
    shutil.copytree(source, destination)
    qasm = tmp_path/row["qasm_path"]
    qasm.write_text(qasm.read_text()+"x q[2];\n")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        run_block_encoding.validate_inputs(tmp_path, row)


def test_inconsistent_structure_is_rejected_even_with_updated_hashes(tmp_path):
    row = row_for(2).copy()
    source = common.BENCHMARKS/Path(row["metadata_path"]).parent
    destination = tmp_path/Path(row["metadata_path"]).parent
    shutil.copytree(source, destination)
    structure_path = destination/"structure.json"
    structure = json.loads(structure_path.read_text())
    structure["select"] = []
    common.write_json(structure_path, structure)
    metadata_path = destination/"metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["files_sha256"]["structure.json"] = common.sha256(structure_path)
    common.write_json(metadata_path, metadata)
    row["metadata_sha256"] = common.sha256(metadata_path)
    with pytest.raises(ValueError, match="actual QASM"):
        run_block_encoding.validate_inputs(tmp_path, row)


def test_timeout_terminates_the_worker(tmp_path):
    result = run_block_encoding.run_case(common.BENCHMARKS, row_for(2), tmp_path, "1e-8", 100, 0.001)
    assert result["status"] == "TIMEOUT"


@pytest.mark.parametrize("n,source_terms", [(2, 3), (3, 6)])
def test_pennylane_generation_preserves_native_inputs(tmp_path, n, source_terms):
    pytest.importorskip("pennylane")
    row = generate_benchmarks.generate_case(tmp_path, n)
    metadata, structure, _ = run_block_encoding.validate_inputs(tmp_path, row)
    assert metadata["n_terms"] == source_terms
    assert metadata["target_lcu_terms"] == 2*n
    assert metadata["numerical_crosscheck"]["status"] == "PASS"
    assert not metadata["select_reconstructed"]
    assert not metadata["pi_snapping"]
    assert common.signatures(structure["unprepare"]) == common.signatures(common.inverse(structure["prepare"]))
