from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_WORK_DIR = ROOT / "evaluation" / "bug_detection" / "generated"
DEFAULT_RESULTS_DIR = ROOT / "evaluation" / "bug_detection" / "results"
DEFAULT_BLOCK_ACCURACY_TARGETS = (1e-1, 1e-3, 1e-5, 1e-8)
QSP_GATE_ANGLE_PRECISION = 12

sys.path.insert(0, str(ROOT))

from evaluation import run_block_encoding, run_hamsim  # noqa: E402
from tools.benchmarks import (  # noqa: E402
    DEFAULT_BLOCK_ENCODING_ANGLE_PRECISION,
    HamsimBenchmark,
    make_hamiltonian,
    write_block_encoding_qasm,
    write_hamsim_qasm,
)
from symbolic.qsp import generate_full_hamsim_record  # noqa: E402

DEFAULT_HAMSIM_BUG_T_VALUES = (0.1, 0.5, 1.0, 2.0)
DEFAULT_HAMSIM_BUG_EPSILON = 1e-4


@dataclass(frozen=True)
class BugCase:
    benchmark_id: str
    intended_bug: bool
    qasm_path: Path
    row: dict[str, str]


def main() -> int:
    args = _parse_args()
    if args.work_dir.exists():
        shutil.rmtree(args.work_dir)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.results_dir.mkdir(parents=True, exist_ok=True)

    stale_combined = args.results_dir / "bug_detection_results.csv"
    if stale_combined.exists():
        stale_combined.unlink()

    block_rows = _run_block_cases(_block_cases(args.work_dir))
    hamsim_rows = _run_hamsim_cases(_hamsim_cases(args.work_dir))

    block_output = args.results_dir / "block_encoding_bug_detection_results.csv"
    hamsim_output = args.results_dir / "hamsim_bug_detection_results.csv"
    _write_rows(block_output, run_block_encoding._fieldnames() + ["bug_detection_result"], block_rows)
    _write_rows(hamsim_output, run_hamsim._fieldnames() + ["bug_detection_result"], hamsim_rows)

    print(f"Wrote {len(block_rows)} RQ1 bug-detection rows to {_format_path(block_output)}")
    print(f"Wrote {len(hamsim_rows)} RQ2 bug-detection rows to {_format_path(hamsim_output)}")
    return 0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run bug-detection cases through the same RQ1/RQ2 evaluation methods."
    )
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    return parser.parse_args()


def _block_cases(work_dir: Path) -> list[BugCase]:
    cases: list[BugCase] = []
    specs = (
        ("ising_n3", make_hamiltonian("ising", 3)),
        ("maxcut_n3", make_hamiltonian("maxcut", 3)),
        ("heisenberg_n3", make_hamiltonian("heisenberg", 3)),
    )
    mutations: tuple[tuple[str, Callable[[str, object], str], bool], ...] = (
        ("clean", lambda text, _benchmark: text, False),
        ("prepare_angle_delta_1e-07", lambda text, _benchmark: _mutate_prepare_angle(text, "1e-7"), True),
        ("prepare_angle_delta_1e-05", lambda text, _benchmark: _mutate_prepare_angle(text, "1e-5"), True),
        ("prepare_angle_delta_1e-03", lambda text, _benchmark: _mutate_prepare_angle(text, "1e-3"), True),
        ("select_cz_to_cx", lambda text, _benchmark: _mutate_first_select_cz_to_cx(text), True),
        ("normalization_full_selector", _mutate_to_full_selector_superposition, True),
    )

    for label, benchmark in specs:
        clean_path = work_dir / "block_encoding" / label / "_clean_source.qasm"
        write_block_encoding_qasm(benchmark, clean_path)
        clean_text = clean_path.read_text(encoding="utf-8")
        clean_path.unlink()
        manifest_row = benchmark.manifest_row(clean_path)

        for mutation_name, mutate, intended_bug in mutations:
            benchmark_id = f"{label}__{mutation_name}"
            qasm_path = work_dir / "block_encoding" / label / f"{benchmark_id}.qasm"
            qasm_path.write_text(mutate(clean_text, benchmark), encoding="utf-8")
            row = dict(manifest_row)
            row["benchmark_id"] = benchmark_id
            row["qasm_path"] = _format_path(qasm_path)
            cases.append(BugCase(benchmark_id, intended_bug, qasm_path, row))
    return cases


def _hamsim_cases(work_dir: Path) -> list[BugCase]:
    cases: list[BugCase] = []
    for tau in DEFAULT_HAMSIM_BUG_T_VALUES:
        epsilon = DEFAULT_HAMSIM_BUG_EPSILON
        label = f"hamsim_t{_float_token(tau)}_eps{_float_token(epsilon)}"
        record = generate_full_hamsim_record(
            tau=tau,
            epsilon=epsilon,
            ensure_bounded=True,
            method="sym_qsp",
            signal_operator="Wx",
        )
        benchmark = HamsimBenchmark(
            benchmark_id=label,
            axis="bug_detection",
            tau=tau,
            epsilon=epsilon,
            block_ancilla_count=1,
            record=record,
        )
        clean_path = work_dir / "hamsim" / label / "_clean_source.qasm"
        write_hamsim_qasm(benchmark, clean_path)
        clean_text = clean_path.read_text(encoding="utf-8")
        clean_path.unlink()
        manifest_row = benchmark.manifest_row(clean_path)

        for mutation_name, mutate, intended_bug in _hamsim_mutations(epsilon):
            benchmark_id = f"{label}__{mutation_name}"
            qasm_path = work_dir / "hamsim" / label / f"{benchmark_id}.qasm"
            qasm_path.write_text(mutate(clean_text), encoding="utf-8")
            row = dict(manifest_row)
            row["benchmark_id"] = benchmark_id
            row["qasm_path"] = _format_path(qasm_path)
            cases.append(BugCase(benchmark_id, intended_bug, qasm_path, row))
    return cases


def _hamsim_mutations(epsilon: float) -> tuple[tuple[str, Callable[[str], str], bool], ...]:
    small = f"{epsilon / 100:.1e}"
    near = f"{epsilon / 2:.1e}"
    large = f"{2 * epsilon:.1e}"
    return (
        ("clean", lambda text: text, False),
        (f"qsp_phase_delta_{small}", lambda text: _mutate_first_qsp_rz(text, small, precision=QSP_GATE_ANGLE_PRECISION), True),
        (f"qsp_phase_delta_{near}", lambda text: _mutate_first_qsp_rz(text, near, precision=QSP_GATE_ANGLE_PRECISION), True),
        (f"qsp_phase_delta_{large}", lambda text: _mutate_first_qsp_rz(text, large, precision=QSP_GATE_ANGLE_PRECISION), True),
        ("oracle_order_moved", _mutate_oracle_order, True),
        ("normalization_drop_final_h", _mutate_qsp_final_normalization, True),
    )


def _run_block_cases(cases: list[BugCase]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for case in cases:
        for epsilon in DEFAULT_BLOCK_ACCURACY_TARGETS:
            result = run_block_encoding._run_benchmark(
                case.row,
                timeout_sec=run_block_encoding.DEFAULT_TIMEOUT_SEC,
                profile_dir=None,
                check_mode="paper",
                block_encoding_epsilon=epsilon,
                block_encoding_residual_tolerance=0.0,
            )
            result["bug_detection_result"] = _block_bug_detection_result(case, result)
            rows.append(result)
    return rows


def _run_hamsim_cases(cases: list[BugCase]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for case in cases:
        result = run_hamsim._run_benchmark(
            case.row,
            timeout_sec=run_hamsim.DEFAULT_TIMEOUT_SEC,
            profile_dir=None,
            check_mode="approximation",
            max_approx_grid_points=1_000_000,
        )
        result["bug_detection_result"] = _hamsim_bug_detection_result(case, result)
        rows.append(result)
    return rows


def _block_bug_detection_result(case: BugCase, result: dict[str, str]) -> str:
    return _bug_detection_result(case, result)


def _hamsim_bug_detection_result(case: BugCase, result: dict[str, str]) -> str:
    return _bug_detection_result(case, result)


def _bug_detection_result(case: BugCase, result: dict[str, str]) -> str:
    detected_bug = result["status"] == "FAIL"
    if case.intended_bug:
        return "True Positive" if detected_bug else "False Negative"
    return "False Positive" if detected_bug else "True Negative"


def _write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _mutate_prepare_angle(text: str, delta: str) -> str:
    lines = text.splitlines()
    prep_start = _find_line(lines, "// PREP: uniform state over Hamiltonian terms") + 1
    select_start = _find_line(lines, "// SELECT: branch-conditioned signed Pauli terms")
    dagger_start = _find_line(lines, "// PREP dagger") + 1
    first = _find_matching_gate(lines, r"ry\(([^)]+)\) q\[0\];", start=prep_start, stop=select_start)
    last = _find_matching_gate(lines, r"ry\(([^)]+)\) q\[0\];", start=dagger_start, stop=len(lines), reverse=True)
    theta = _gate_angle(lines[first], precision=DEFAULT_BLOCK_ENCODING_ANGLE_PRECISION)
    mutated_theta = theta + _angle_delta(delta, precision=DEFAULT_BLOCK_ENCODING_ANGLE_PRECISION)
    lines[first] = _replace_gate_angle(
        lines[first],
        mutated_theta,
        precision=DEFAULT_BLOCK_ENCODING_ANGLE_PRECISION,
    )
    lines[last] = _replace_gate_angle(
        lines[last],
        -mutated_theta,
        precision=DEFAULT_BLOCK_ENCODING_ANGLE_PRECISION,
    )
    return "\n".join(lines) + "\n"


def _mutate_first_select_cz_to_cx(text: str) -> str:
    lines = text.splitlines()
    select_start = _find_line(lines, "// SELECT: branch-conditioned signed Pauli terms") + 1
    dagger_start = _find_line(lines, "// PREP dagger")
    index = _find_matching_gate(lines, r"cz q\[(\d+)\], q\[(\d+)\];", start=select_start, stop=dagger_start)
    lines[index] = lines[index].replace("cz ", "cx ", 1)
    return "\n".join(lines) + "\n"


def _mutate_to_full_selector_superposition(text: str, benchmark: object) -> str:
    selector_width = int(getattr(benchmark, "selector_ancillas"))
    lines = text.splitlines()
    prep_comment = _find_line(lines, "// PREP: uniform state over Hamiltonian terms")
    select_start = _find_line(lines, "// SELECT: branch-conditioned signed Pauli terms")
    dagger_comment = _find_line(lines, "// PREP dagger")
    prep = [f"h q[{index}];" for index in range(selector_width)]
    dagger = [f"h q[{index}];" for index in reversed(range(selector_width))]
    mutated = lines[: prep_comment + 1] + prep + lines[select_start: dagger_comment + 1] + dagger
    return "\n".join(mutated) + "\n"


def _mutate_first_qsp_rz(text: str, delta: str, *, precision: int) -> str:
    lines = text.splitlines()
    start = _find_line(lines, "// common phase 0")
    index = _find_matching_gate(lines, r"rz\(([^)]+)\) q\[2\];", start=start, stop=len(lines))
    angle = _gate_angle(lines[index], precision=precision)
    lines[index] = _replace_gate_angle(
        lines[index],
        angle + _angle_delta(delta, precision=precision),
        precision=precision,
    )
    return "\n".join(lines) + "\n"


def _mutate_oracle_order(text: str) -> str:
    lines = text.splitlines()
    oracle_indices = [
        index
        for index, line in enumerate(lines)
        if re.match(r"^(UH|UHdg) q\[\d+\], q\[\d+\];$", line.strip())
    ]
    if len(oracle_indices) < 2:
        raise ValueError("Expected at least two common QSP oracle calls")
    first, second = oracle_indices[:2]
    first_line = lines.pop(first)
    lines.insert(second - 1, first_line)
    return "\n".join(lines) + "\n"


def _mutate_qsp_final_normalization(text: str) -> str:
    old = "sdg q[0];\nh q[0];"
    new = "sdg q[0];\nrz(0.0) q[0];"
    if old not in text:
        raise ValueError("Could not find final selector normalization block")
    return text.replace(old, new, 1)


def _find_line(lines: list[str], needle: str) -> int:
    for index, line in enumerate(lines):
        if line.strip() == needle:
            return index
    raise ValueError(f"Could not find line: {needle}")


def _find_matching_gate(lines: list[str], pattern: str, *, start: int, stop: int, reverse: bool = False) -> int:
    indexes = range(stop - 1, start - 1, -1) if reverse else range(start, stop)
    regex = re.compile(pattern)
    for index in indexes:
        if regex.match(lines[index].strip()):
            return index
    raise ValueError(f"Could not find gate matching {pattern!r}")


def _gate_angle(line: str, *, precision: int) -> sp.Expr:
    match = re.search(r"\(([^)]+)\)", line)
    if match is None:
        raise ValueError(f"Gate has no angle: {line}")
    return sp.N(sp.sympify(match.group(1), locals={"pi": sp.pi}), precision)


def _angle_delta(value: str, *, precision: int) -> sp.Expr:
    return sp.Float(value, precision)


def _replace_gate_angle(line: str, value: sp.Expr, *, precision: int) -> str:
    if precision >= DEFAULT_BLOCK_ENCODING_ANGLE_PRECISION:
        angle_text = str(sp.N(value, precision))
    else:
        angle_text = f"{float(sp.N(value, 30)):.{precision}g}"
    return re.sub(r"\([^)]*\)", f"({angle_text})", line, count=1)


def _float_token(value: float) -> str:
    return f"{value:.0e}".replace("+", "").replace("-", "m")


def _format_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


if __name__ == "__main__":
    raise SystemExit(main())
