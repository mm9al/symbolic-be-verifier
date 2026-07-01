from pathlib import Path

from symbolic.verify import PASS, verify_qasm_file
from tools.benchmarks import make_hamiltonian, write_block_encoding_qasm


def test_ising_block_encoding_generator_verifies(tmp_path: Path):
    benchmark = make_hamiltonian("ising", 3)
    qasm_path = tmp_path / "ising_periodic_n3.qasm"
    write_block_encoding_qasm(benchmark, qasm_path)

    assert benchmark.graph == "periodic"
    assert len(benchmark.terms) == 6
    assert "ry(" in qasm_path.read_text(encoding="utf-8")

    result = verify_qasm_file(
        qasm_path,
        expected=str(benchmark.expected_operator()),
        ancillas=benchmark.ancillas,
        systems=benchmark.systems,
    )

    assert result.status == PASS


def test_maxcut_cycle_block_encoding_generator_verifies_negative_terms(tmp_path: Path):
    benchmark = make_hamiltonian("maxcut", 3)
    qasm_path = tmp_path / "maxcut_cycle_n3.qasm"
    write_block_encoding_qasm(benchmark, qasm_path)

    assert len(benchmark.terms) == 6
    assert "ry(" in qasm_path.read_text(encoding="utf-8")

    result = verify_qasm_file(
        qasm_path,
        expected=str(benchmark.expected_operator()),
        ancillas=benchmark.ancillas,
        systems=benchmark.systems,
    )

    assert result.status == PASS


def test_heisenberg_block_encoding_generator_verifies_y_terms(tmp_path: Path):
    benchmark = make_hamiltonian("heisenberg", 3)
    qasm_path = tmp_path / "heisenberg_n3.qasm"
    write_block_encoding_qasm(benchmark, qasm_path)

    result = verify_qasm_file(
        qasm_path,
        expected=str(benchmark.expected_operator()),
        ancillas=benchmark.ancillas,
        systems=benchmark.systems,
    )

    assert result.status == PASS
