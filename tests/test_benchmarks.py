from pathlib import Path
import re

import mpmath as mp

from symbolic.verify import PASS_EXACT, verify_qasm_file
from tools.benchmarks import hamsim_suite, make_hamiltonian, write_block_encoding_qasm


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

    assert result.status == PASS_EXACT


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

    assert result.status == PASS_EXACT


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

    assert result.status == PASS_EXACT


def test_block_encoding_generator_can_emit_high_precision_angles(tmp_path: Path):
    benchmark = make_hamiltonian("ising", 3)
    qasm_path = tmp_path / "ising_periodic_n3.qasm"
    write_block_encoding_qasm(benchmark, qasm_path, angle_precision=60)

    match = re.search(r"ry\((-?[0-9]+\.[0-9]+)", qasm_path.read_text(encoding="utf-8"))
    angles = re.findall(r"ry\((-?[0-9]+\.[0-9]+)", qasm_path.read_text(encoding="utf-8"))

    assert match is not None
    assert len(match.group(1).replace("-", "").replace(".", "")) >= 55
    with mp.workdps(80):
        assert abs(mp.mpf(angles[1]) - mp.pi / 4) < mp.mpf("1e-55")


def test_hamsim_suite_pins_the_intended_axis_variables():
    benchmarks = hamsim_suite()
    by_axis = {
        axis: [benchmark for benchmark in benchmarks if benchmark.axis == axis]
        for axis in ("vary_t", "vary_epsilon")
    }

    assert [benchmark.tau for benchmark in by_axis["vary_t"]] == [
        0.1,
        0.5,
        1.0,
        2.0,
        4.0,
        5.0,
        6.0,
        8.0,
        9.0,
        10.0,
        11.0,
        13.0,
        14.0,
        15.0,
        17.0,
        18.0,
        20.0,
        22.0,
        23.0,
        24.0,
        26.0,
        27.0,
        28.0,
        30.0,
        32.0,
        33.0,
        34.0,
        36.0,
        38.0,
        39.0,
        40.0,
        42.0,
        43.0,
        44.0,
        46.0,
        48.0,
        49.0,
        50.0,
    ]
    assert [benchmark.degree for benchmark in by_axis["vary_t"]] == [
        3,
        5,
        7,
        9,
        11,
        13,
        15,
        17,
        19,
        21,
        23,
        25,
        27,
        29,
        31,
        33,
        35,
        37,
        39,
        41,
        43,
        45,
        47,
        49,
        51,
        53,
        55,
        57,
        59,
        61,
        63,
        65,
        67,
        69,
        71,
        73,
        75,
        77,
    ]
    assert {benchmark.epsilon for benchmark in by_axis["vary_t"]} == {1e-4}
    assert {benchmark.block_ancilla_count for benchmark in by_axis["vary_t"]} == {1}

    assert {benchmark.tau for benchmark in by_axis["vary_epsilon"]} == {0.5}
    assert {benchmark.block_ancilla_count for benchmark in by_axis["vary_epsilon"]} == {1}
    assert [benchmark.epsilon for benchmark in by_axis["vary_epsilon"]] == [
        1e-1,
        1e-2,
        1e-3,
        1e-4,
        1e-5,
        1e-6,
        1e-7,
        1e-8,
        1e-9,
        1e-10,
        1e-11,
        1e-12,
    ]
    assert [benchmark.degree for benchmark in by_axis["vary_epsilon"]] == [
        3,
        3,
        3,
        5,
        5,
        7,
        7,
        7,
        9,
        9,
        9,
        11,
    ]
