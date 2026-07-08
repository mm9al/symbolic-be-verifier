from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
import resource
import sys
import time
import tracemalloc
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "benchmarks" / "hamsim" / "manifest.csv"
RESULTS_DIR = ROOT / "evaluation" / "results"
RESULTS_PATH = RESULTS_DIR / "hamsim_results.csv"
DEFAULT_TIMEOUT_SEC = 1800.0

sys.path.insert(0, str(ROOT))

from symbolic.qasm_parser import parse_qasm_file  # noqa: E402
from symbolic.verify import (  # noqa: E402
    DEFAULT_MAX_APPROXIMATION_GRID_POINTS,
    PASS,
    gate_profile_fieldnames,
    gate_profile_rows,
    rescale_polynomial_for_target,
    verify_polynomial_approximates_exp,
    verify_qasm_file,
)


def main() -> int:
    args = _parse_args()
    manifest_path = args.manifest
    if not manifest_path.exists():
        raise SystemExit(f"Missing {manifest_path.relative_to(ROOT)}. Run evaluation/generate_hamsim_benchmarks.py first.")

    rows = _read_manifest(manifest_path)
    if args.axes:
        selected = set(args.axes)
        rows = [row for row in rows if row["axis"] in selected]

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    max_approx_grid_points = None if args.max_approx_grid_points == 0 else args.max_approx_grid_points
    timeout_axes: set[str] = set()
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_fieldnames())
        writer.writeheader()
        for row in rows:
            if row["axis"] in timeout_axes:
                result = _skipped_after_timeout(row, args.check_mode)
            else:
                result = _run_benchmark(
                    row,
                    timeout_sec=args.timeout_sec,
                    profile_dir=args.profile_dir,
                    check_mode=args.check_mode,
                    max_approx_grid_points=max_approx_grid_points,
                )
                if result["status"] == "TIMEOUT":
                    timeout_axes.add(row["axis"])
            writer.writerow(result)
            handle.flush()

    print(f"Wrote {_format_path(args.output)}")
    return 0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run RQ2 Hamiltonian-simulation QSP symbolic verification benchmarks.")
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--output", type=Path, default=RESULTS_PATH)
    parser.add_argument("--axes", nargs="+", choices=("vary_t", "vary_epsilon"))
    parser.add_argument(
        "--check-mode",
        choices=("polynomial", "approximation", "both"),
        default="polynomial",
        help=(
            "Use exact generated-polynomial comparison, numerical exp(-iHt) approximation checking, "
            "or both checks from one symbolic polynomial."
        ),
    )
    parser.add_argument(
        "--timeout-sec",
        type=float,
        default=DEFAULT_TIMEOUT_SEC,
        help=f"Per-benchmark timeout in seconds. Use 0 to disable. Default: {DEFAULT_TIMEOUT_SEC:g}.",
    )
    parser.add_argument(
        "--profile-dir",
        type=Path,
        help="Optional directory for per-benchmark gate-level profile CSV files.",
    )
    parser.add_argument(
        "--max-approx-grid-points",
        type=int,
        default=DEFAULT_MAX_APPROXIMATION_GRID_POINTS,
        help=(
            "Fail fast if an approximation grid would exceed this many points. "
            "Use 0 to disable."
        ),
    )
    return parser.parse_args()


def _read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _run_benchmark(
    row: dict[str, str],
    *,
    timeout_sec: float,
    profile_dir: Path | None,
    check_mode: str,
    max_approx_grid_points: int | None,
) -> dict[str, str]:
    if timeout_sec <= 0:
        return _run_benchmark_inner(
            row,
            profile_dir=profile_dir,
            check_mode=check_mode,
            max_approx_grid_points=max_approx_grid_points,
        )

    queue: mp.Queue = mp.Queue(maxsize=1)
    process = mp.Process(target=_worker, args=(row, queue, profile_dir, check_mode, max_approx_grid_points))
    started = time.perf_counter()
    process.start()
    process.join(timeout_sec)

    if process.is_alive():
        process.terminate()
        process.join()
        return {
            **_metadata(row, check_mode),
            "runtime_sec": f"{time.perf_counter() - started:.9f}",
            "symbolic_runtime_sec": "",
            "approx_runtime_sec": "",
            "tracemalloc_peak_mb": "",
            "process_maxrss_mb": "",
            "status": "TIMEOUT",
            "success": "False",
            "approx_max_grid_error": "",
            "approx_worst_x": "",
            "approx_grid_points": "",
            "error": f"timeout after {timeout_sec:g}s",
        }

    if not queue.empty():
        return queue.get()

    return {
        **_metadata(row, check_mode),
        "runtime_sec": f"{time.perf_counter() - started:.9f}",
        "symbolic_runtime_sec": "",
        "approx_runtime_sec": "",
        "tracemalloc_peak_mb": "",
        "process_maxrss_mb": "",
        "status": "ERROR",
        "success": "False",
        "approx_max_grid_error": "",
        "approx_worst_x": "",
        "approx_grid_points": "",
        "error": f"worker exited with code {process.exitcode}",
    }


def _worker(
    row: dict[str, str],
    queue: "mp.Queue",
    profile_dir: Path | None,
    check_mode: str,
    max_approx_grid_points: int | None,
) -> None:
    queue.put(
        _run_benchmark_inner(
            row,
            profile_dir=profile_dir,
            check_mode=check_mode,
            max_approx_grid_points=max_approx_grid_points,
        )
    )


def _run_benchmark_inner(
    row: dict[str, str],
    *,
    profile_dir: Path | None,
    check_mode: str,
    max_approx_grid_points: int | None,
) -> dict[str, str]:
    error = ""
    status = "ERROR"
    success = "False"
    approx_max_grid_error = ""
    approx_worst_x = ""
    approx_grid_points = ""
    symbolic_runtime_sec = ""
    approx_runtime_sec = ""
    qasm_path = ROOT / row["qasm_path"]
    ancillas = _parse_int_list(row["ancillas"])
    systems = _parse_int_list(row["systems"])

    tracemalloc.start()
    started = time.perf_counter()
    try:
        kwargs = {
            "ancillas": ancillas,
            "systems": systems,
            "keep_trace": False,
            "profile_gates": profile_dir is not None,
            "hermitian_base": True,
        }
        if check_mode in {"polynomial", "both"}:
            kwargs.update(
                {
                    "expected_polynomial": row["expected_polynomial"],
                    "compare_polynomial_only": True,
                }
            )
        else:
            kwargs["extract_qsp_polynomial"] = True

        symbolic_started = time.perf_counter()
        result = verify_qasm_file(qasm_path, **kwargs)
        symbolic_runtime_sec = f"{time.perf_counter() - symbolic_started:.9f}"
        if profile_dir is not None:
            _write_gate_profile(row, result, profile_dir)
        status = result.status or "NO_EXPECTED"
        success = str(result.success)
        if check_mode in {"approximation", "both"}:
            if result.qsp_polynomial is None:
                raise ValueError(f"{check_mode} check mode expected a generated QSP polynomial")
            approx_started = time.perf_counter()
            approximation = verify_polynomial_approximates_exp(
                rescale_polynomial_for_target(result.qsp_polynomial, row["target_scale"]),
                tau=float(row["tau"]),
                epsilon=float(row["epsilon"]),
                max_grid_points=max_approx_grid_points,
            )
            approx_runtime_sec = f"{time.perf_counter() - approx_started:.9f}"
            approx_max_grid_error = f"{approximation.max_grid_error:.12g}"
            approx_worst_x = f"{approximation.worst_x:.12g}"
            approx_grid_points = str(approximation.num_grid_points)
            if check_mode == "both":
                status = PASS if result.status == PASS and approximation.success else "FAIL"
            else:
                status = PASS if approximation.success else "FAIL"
            success = str(status == PASS)
        if result.qsp_approximation is not None:
            approx_max_grid_error = f"{result.qsp_approximation.max_grid_error:.12g}"
            approx_worst_x = f"{result.qsp_approximation.worst_x:.12g}"
            approx_grid_points = str(result.qsp_approximation.num_grid_points)
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}"
    runtime_sec = time.perf_counter() - started
    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return {
        **_metadata(row, check_mode),
        "runtime_sec": f"{runtime_sec:.9f}",
        "symbolic_runtime_sec": symbolic_runtime_sec,
        "approx_runtime_sec": approx_runtime_sec,
        "tracemalloc_peak_mb": f"{_bytes_to_mib(peak_bytes):.6f}",
        "process_maxrss_mb": f"{_process_maxrss_mb():.6f}",
        "status": status,
        "success": success,
        "approx_max_grid_error": approx_max_grid_error,
        "approx_worst_x": approx_worst_x,
        "approx_grid_points": approx_grid_points,
        "error": error,
    }


def _metadata(row: dict[str, str], check_mode: str) -> dict[str, str]:
    return {
        "benchmark_id": row["benchmark_id"],
        "rq": row["rq"],
        "stage": row["stage"],
        "axis": row["axis"],
        "check_mode": check_mode,
        "tau": row["tau"],
        "epsilon": row["epsilon"],
        "degree": row["degree"],
        "cos_degree": row["cos_degree"],
        "sin_degree": row["sin_degree"],
        "uh_ancillas": row["uh_ancillas"],
        "n_ancilla": row["n_ancilla"],
        "n_system": row["n_system"],
        "target_scale": row["target_scale"],
        "gate_count": str(_gate_count(row)),
        "qasm_path": row["qasm_path"],
    }


def _write_gate_profile(row: dict[str, str], result, profile_dir: Path) -> None:
    profile_dir.mkdir(parents=True, exist_ok=True)
    path = profile_dir / f"{row['benchmark_id']}_gate_profile.csv"
    fieldnames = ["benchmark_id", *gate_profile_fieldnames()]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for profile_row in gate_profile_rows(result):
            writer.writerow({"benchmark_id": row["benchmark_id"], **profile_row})


def _skipped_after_timeout(row: dict[str, str], check_mode: str) -> dict[str, str]:
    return {
        **_metadata(row, check_mode),
        "runtime_sec": "",
        "symbolic_runtime_sec": "",
        "approx_runtime_sec": "",
        "tracemalloc_peak_mb": "",
        "process_maxrss_mb": "",
        "status": "SKIPPED_AFTER_TIMEOUT",
        "success": "",
        "approx_max_grid_error": "",
        "approx_worst_x": "",
        "approx_grid_points": "",
        "error": "previous benchmark in this axis timed out",
    }


def _fieldnames() -> list[str]:
    return [
        "benchmark_id",
        "rq",
        "stage",
        "axis",
        "check_mode",
        "tau",
        "epsilon",
        "degree",
        "cos_degree",
        "sin_degree",
        "uh_ancillas",
        "n_ancilla",
        "n_system",
        "target_scale",
        "gate_count",
        "qasm_path",
        "runtime_sec",
        "symbolic_runtime_sec",
        "approx_runtime_sec",
        "tracemalloc_peak_mb",
        "process_maxrss_mb",
        "status",
        "success",
        "approx_max_grid_error",
        "approx_worst_x",
        "approx_grid_points",
        "error",
    ]


def _parse_int_list(text: str) -> tuple[int, ...]:
    return tuple(int(piece) for piece in text.split())


def _gate_count(row: dict[str, str]) -> int:
    return len(parse_qasm_file(ROOT / row["qasm_path"]))


def _bytes_to_mib(value: int) -> float:
    return value / (1024 * 1024)


def _process_maxrss_mb() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return rss / (1024 * 1024)
    return rss / 1024


def _format_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


if __name__ == "__main__":
    raise SystemExit(main())
