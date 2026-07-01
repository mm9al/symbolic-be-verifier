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
MANIFEST_PATH = ROOT / "benchmarks" / "block_encoding" / "manifest.csv"
RESULTS_DIR = ROOT / "evaluation" / "results"
RESULTS_PATH = RESULTS_DIR / "block_encoding_results.csv"
DEFAULT_TIMEOUT_SEC = 300.0

sys.path.insert(0, str(ROOT))

from symbolic.verify import gate_profile_fieldnames, gate_profile_rows, verify_qasm_file  # noqa: E402


def main() -> int:
    args = _parse_args()
    manifest_path = args.manifest
    if not manifest_path.exists():
        raise SystemExit(f"Missing {manifest_path.relative_to(ROOT)}. Run evaluation/generate_benchmarks.py first.")

    rows = _read_manifest(manifest_path)
    if args.models:
        selected = set(args.models)
        rows = [row for row in rows if row["model"] in selected]

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    timeout_models: set[str] = set()
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_fieldnames())
        writer.writeheader()
        for row in rows:
            if row["model"] in timeout_models:
                result = _skipped_after_timeout(row)
            else:
                result = _run_benchmark(row, timeout_sec=args.timeout_sec, profile_dir=args.profile_dir)
                if result["status"] == "TIMEOUT":
                    timeout_models.add(row["model"])
            writer.writerow(result)
            handle.flush()

    print(f"Wrote {_format_path(args.output)}")
    return 0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run RQ1 block-encoding symbolic verification benchmarks.")
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--output", type=Path, default=RESULTS_PATH)
    parser.add_argument("--models", nargs="+", choices=("ising", "maxcut", "heisenberg"))
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
    return parser.parse_args()


def _read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _run_benchmark(row: dict[str, str], *, timeout_sec: float, profile_dir: Path | None) -> dict[str, str]:
    if timeout_sec <= 0:
        return _run_benchmark_inner(row, profile_dir=profile_dir)

    queue: mp.Queue = mp.Queue(maxsize=1)
    process = mp.Process(target=_worker, args=(row, queue, profile_dir))
    started = time.perf_counter()
    process.start()
    process.join(timeout_sec)

    if process.is_alive():
        process.terminate()
        process.join()
        return {
            **_metadata(row),
            "runtime_sec": f"{time.perf_counter() - started:.9f}",
            "tracemalloc_peak_mb": "",
            "process_maxrss_mb": "",
            "status": "TIMEOUT",
            "success": "False",
            "error": f"timeout after {timeout_sec:g}s",
        }

    if not queue.empty():
        return queue.get()

    return {
        **_metadata(row),
        "runtime_sec": f"{time.perf_counter() - started:.9f}",
        "tracemalloc_peak_mb": "",
        "process_maxrss_mb": "",
        "status": "ERROR",
        "success": "False",
        "error": f"worker exited with code {process.exitcode}",
    }


def _worker(row: dict[str, str], queue: "mp.Queue", profile_dir: Path | None) -> None:
    queue.put(_run_benchmark_inner(row, profile_dir=profile_dir))


def _run_benchmark_inner(row: dict[str, str], *, profile_dir: Path | None) -> dict[str, str]:
    error = ""
    status = "ERROR"
    success = "False"
    qasm_path = ROOT / row["qasm_path"]
    ancillas = _parse_int_list(row["ancillas"])
    systems = _parse_int_list(row["systems"])

    tracemalloc.start()
    started = time.perf_counter()
    try:
        result = verify_qasm_file(
            qasm_path,
            expected=row["expected"],
            ancillas=ancillas,
            systems=systems,
            keep_trace=False,
            profile_gates=profile_dir is not None,
        )
        if profile_dir is not None:
            _write_gate_profile(row, result, profile_dir)
        status = result.status or "NO_EXPECTED"
        success = str(result.success)
    except Exception as exc:
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
        "rq": row["rq"],
        "stage": row["stage"],
        "model": row["model"],
        "graph": row["graph"],
        "n_system": row["n_system"],
        "n_terms": row["n_terms"],
        "alpha": row["alpha"],
        "selector_ancillas": row["selector_ancillas"],
        "n_ancilla": row["n_ancilla"],
        "max_locality": row["max_locality"],
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


def _skipped_after_timeout(row: dict[str, str]) -> dict[str, str]:
    return {
        **_metadata(row),
        "runtime_sec": "",
        "tracemalloc_peak_mb": "",
        "process_maxrss_mb": "",
        "status": "SKIPPED_AFTER_TIMEOUT",
        "success": "",
        "error": "previous benchmark in this model timed out",
    }


def _fieldnames() -> list[str]:
    return [
        "benchmark_id",
        "rq",
        "stage",
        "model",
        "graph",
        "n_system",
        "n_terms",
        "alpha",
        "selector_ancillas",
        "n_ancilla",
        "max_locality",
        "qasm_path",
        "runtime_sec",
        "tracemalloc_peak_mb",
        "process_maxrss_mb",
        "status",
        "success",
        "error",
    ]


def _parse_int_list(text: str) -> tuple[int, ...]:
    return tuple(int(piece) for piece in text.split())


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
