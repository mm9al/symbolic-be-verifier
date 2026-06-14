from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
import resource
import sys
import time
import tracemalloc
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "benchmarks" / "manifest.csv"
RESULTS_DIR = ROOT / "evaluation" / "results"
RESULTS_PATH = RESULTS_DIR / "symbolic_results.csv"
DEFAULT_TIMEOUT_SEC = 300.0

sys.path.insert(0, str(ROOT))

from symbolic.expr import OpExpr, identity, pauli  # noqa: E402
from symbolic.verify import verify_qasm_file  # noqa: E402


def main() -> int:
    args = _parse_args()
    if not MANIFEST_PATH.exists():
        raise SystemExit(f"Missing {MANIFEST_PATH.relative_to(ROOT)}. Run evaluation/generate_benchmarks.py first.")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    timeout_families: set[str] = set()
    with RESULTS_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_fieldnames())
        writer.writeheader()
        for row in _read_manifest():
            if row["family"] in timeout_families:
                result = _skipped_after_timeout(row)
            else:
                result = _run_benchmark(row, timeout_sec=args.timeout_sec)
                if result["status"] == "TIMEOUT":
                    timeout_families.add(row["family"])
            writer.writerow(result)
            handle.flush()
    print(f"Wrote {RESULTS_PATH.relative_to(ROOT)}")
    return 0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run symbolic verification benchmarks with a per-benchmark timeout.")
    parser.add_argument(
        "--timeout-sec",
        type=float,
        default=DEFAULT_TIMEOUT_SEC,
        help=f"Per-benchmark timeout in seconds. Use 0 to disable. Default: {DEFAULT_TIMEOUT_SEC:g}.",
    )
    return parser.parse_args()


def _read_manifest() -> list[dict[str, str]]:
    with MANIFEST_PATH.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _run_benchmark(row: dict[str, str], *, timeout_sec: float) -> dict[str, str]:
    if timeout_sec <= 0:
        return _run_benchmark_inner(row)

    queue: mp.Queue = mp.Queue(maxsize=1)
    process = mp.Process(target=_run_benchmark_worker, args=(row, queue))
    started = time.perf_counter()
    process.start()
    process.join(timeout_sec)

    if process.is_alive():
        process.terminate()
        process.join()
        runtime_sec = time.perf_counter() - started
        return {
            **_metadata(row),
            "runtime_sec": f"{runtime_sec:.9f}",
            "tracemalloc_peak_mb": "",
            "process_maxrss_mb": "",
            "status": "TIMEOUT",
            "success": "False",
            "error": f"timeout after {timeout_sec:g}s",
        }

    if not queue.empty():
        return queue.get()

    runtime_sec = time.perf_counter() - started
    return {
        **_metadata(row),
        "runtime_sec": f"{runtime_sec:.9f}",
        "tracemalloc_peak_mb": "",
        "process_maxrss_mb": "",
        "status": "ERROR",
        "success": "False",
        "error": f"worker exited with code {process.exitcode}",
    }


def _run_benchmark_worker(row: dict[str, str], queue: "mp.Queue") -> None:
    queue.put(_run_benchmark_inner(row))


def _run_benchmark_inner(row: dict[str, str]) -> dict[str, str]:
    qasm_path = ROOT / row["qasm_path"]
    ancillas = _parse_int_list(row["ancillas"])
    systems = _parse_int_list(row["systems"])

    error = ""
    status = "ERROR"
    success = "False"

    tracemalloc.start()
    started = time.perf_counter()
    try:
        result = verify_qasm_file(
            qasm_path,
            expected=_expected(row),
            ancillas=ancillas,
            systems=systems,
            keep_trace=False,
        )
        status = result.status or "NO_EXPECTED"
        success = str(result.success)
    except Exception as exc:  # keep the evaluation run going
        error = f"{type(exc).__name__}: {exc}"
    runtime_sec = time.perf_counter() - started
    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return {
        **_metadata(row),
        "runtime_sec": f"{runtime_sec:.9f}",
        "tracemalloc_peak_mb": f"{_bytes_to_mib(peak_bytes):.6f}",
        "process_maxrss_mb": f"{_process_maxrss_mb():.6f}",
        "status": status,
        "success": success,
        "error": error,
    }


def _metadata(row: dict[str, str]) -> dict[str, str]:
    return {
        "benchmark_id": row["benchmark_id"],
        "kind": row["kind"],
        "family": row["family"],
        "n_system": row["n_system"],
        "n_ancilla": row["n_ancilla"],
        "qasm_path": row["qasm_path"],
    }


def _skipped_after_timeout(row: dict[str, str]) -> dict[str, str]:
    return {
        **_metadata(row),
        "runtime_sec": "",
        "tracemalloc_peak_mb": "",
        "process_maxrss_mb": "",
        "status": "SKIPPED_AFTER_TIMEOUT",
        "success": "",
        "error": "previous benchmark in this family timed out",
    }


def _expected(row: dict[str, str]) -> str | OpExpr:
    if row["family"] == "branch_term_product":
        return _branch_term_expected(int(row["n_system"]))
    return row["expected"]


def _branch_term_expected(n_system: int) -> OpExpr:
    expected = identity(n_system)
    identity_op = identity(n_system)
    for index in range(n_system):
        factor = (identity_op + pauli("X", index=index, num_qubits=n_system)).scale(sp.Rational(1, 2))
        expected = factor * expected
    return expected


def _parse_int_list(text: str) -> tuple[int, ...]:
    return tuple(int(piece) for piece in text.split())


def _bytes_to_mib(value: int) -> float:
    return value / (1024 * 1024)


def _process_maxrss_mb() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return rss / (1024 * 1024)
    return rss / 1024


def _fieldnames() -> list[str]:
    return [
        "benchmark_id",
        "kind",
        "family",
        "n_system",
        "n_ancilla",
        "qasm_path",
        "runtime_sec",
        "tracemalloc_peak_mb",
        "process_maxrss_mb",
        "status",
        "success",
        "error",
    ]


if __name__ == "__main__":
    raise SystemExit(main())
