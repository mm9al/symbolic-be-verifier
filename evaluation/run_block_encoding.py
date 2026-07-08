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
DEFAULT_TIMEOUT_SEC = 1800.0
DEFAULT_BLOCK_ENCODING_EPSILON = 1e-8

sys.path.insert(0, str(ROOT))

from symbolic.qasm_parser import parse_qasm_file  # noqa: E402
from symbolic.verify import (  # noqa: E402
    DEFAULT_BLOCK_ENCODING_RESIDUAL_TOLERANCE,
    gate_profile_fieldnames,
    gate_profile_rows,
    parse_scalar_expression,
    verify_qasm_file,
)


def main() -> int:
    args = _parse_args()
    if args.check_mode == "paper" and args.block_encoding_epsilon <= 0:
        raise SystemExit("--block-encoding-epsilon must be positive")
    if args.check_mode == "paper" and args.block_encoding_residual_tolerance < 0:
        raise SystemExit("--block-encoding-residual-tolerance must be nonnegative")
    block_encoding_epsilon = args.block_encoding_epsilon if args.check_mode == "paper" else None
    block_encoding_residual_tolerance = (
        args.block_encoding_residual_tolerance if args.check_mode == "paper" else None
    )
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
                result = _skipped_after_timeout(
                    row,
                    check_mode=args.check_mode,
                    block_encoding_epsilon=block_encoding_epsilon,
                    block_encoding_residual_tolerance=block_encoding_residual_tolerance,
                )
            else:
                result = _run_benchmark(
                    row,
                    timeout_sec=args.timeout_sec,
                    profile_dir=args.profile_dir,
                    check_mode=args.check_mode,
                    block_encoding_epsilon=block_encoding_epsilon,
                    block_encoding_residual_tolerance=block_encoding_residual_tolerance,
                )
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
        "--check-mode",
        choices=("paper", "exact"),
        default="paper",
        help=(
            "Verification mode. 'paper' checks approximate proportionality against H; "
            "'exact' compares against the normalized top-left block from the manifest. Default: paper."
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
        "--block-encoding-epsilon",
        type=float,
        default=DEFAULT_BLOCK_ENCODING_EPSILON,
        help=(
            "Tolerance for --check-mode paper. The manifest's normalized expected block "
            f"is multiplied by alpha to recover the target H. Default: {DEFAULT_BLOCK_ENCODING_EPSILON:g}."
        ),
    )
    parser.add_argument(
        "--block-encoding-residual-tolerance",
        type=float,
        default=DEFAULT_BLOCK_ENCODING_RESIDUAL_TOLERANCE,
        help=(
            "Numerical floor for paper-mode coefficient residual checks with decimal QASM angles. "
            f"Default: {DEFAULT_BLOCK_ENCODING_RESIDUAL_TOLERANCE:g}."
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
    block_encoding_epsilon: float | None,
    block_encoding_residual_tolerance: float | None,
) -> dict[str, str]:
    if timeout_sec <= 0:
        return _run_benchmark_inner(
            row,
            profile_dir=profile_dir,
            check_mode=check_mode,
            block_encoding_epsilon=block_encoding_epsilon,
            block_encoding_residual_tolerance=block_encoding_residual_tolerance,
        )

    queue: mp.Queue = mp.Queue(maxsize=1)
    process = mp.Process(
        target=_worker,
        args=(row, queue, profile_dir, check_mode, block_encoding_epsilon, block_encoding_residual_tolerance),
    )
    started = time.perf_counter()
    process.start()
    process.join(timeout_sec)

    if process.is_alive():
        process.terminate()
        process.join()
        return {
            **_metadata(row),
            **_check_metadata(check_mode, block_encoding_epsilon, block_encoding_residual_tolerance, None),
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
        **_check_metadata(check_mode, block_encoding_epsilon, block_encoding_residual_tolerance, None),
        "runtime_sec": f"{time.perf_counter() - started:.9f}",
        "tracemalloc_peak_mb": "",
        "process_maxrss_mb": "",
        "status": "ERROR",
        "success": "False",
        "error": f"worker exited with code {process.exitcode}",
    }


def _worker(
    row: dict[str, str],
    queue: "mp.Queue",
    profile_dir: Path | None,
    check_mode: str,
    block_encoding_epsilon: float | None,
    block_encoding_residual_tolerance: float | None,
) -> None:
    queue.put(
        _run_benchmark_inner(
            row,
            profile_dir=profile_dir,
            check_mode=check_mode,
            block_encoding_epsilon=block_encoding_epsilon,
            block_encoding_residual_tolerance=block_encoding_residual_tolerance,
        )
    )


def _run_benchmark_inner(
    row: dict[str, str],
    *,
    profile_dir: Path | None,
    check_mode: str,
    block_encoding_epsilon: float | None,
    block_encoding_residual_tolerance: float | None,
) -> dict[str, str]:
    error = ""
    status = "ERROR"
    success = "False"
    qasm_path = ROOT / row["qasm_path"]
    ancillas = _parse_int_list(row["ancillas"])
    systems = _parse_int_list(row["systems"])
    expected = row["expected"]
    if check_mode == "paper":
        alpha = parse_scalar_expression(row["alpha"])
        expected = f"({alpha})*({expected})"

    tracemalloc.start()
    started = time.perf_counter()
    try:
        result = verify_qasm_file(
            qasm_path,
            expected=expected,
            ancillas=ancillas,
            systems=systems,
            keep_trace=False,
            profile_gates=profile_dir is not None,
            block_encoding_epsilon=block_encoding_epsilon,
            block_encoding_residual_tolerance=(
                block_encoding_residual_tolerance
                if block_encoding_residual_tolerance is not None
                else DEFAULT_BLOCK_ENCODING_RESIDUAL_TOLERANCE
            ),
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
        **_check_metadata(
            check_mode,
            block_encoding_epsilon,
            block_encoding_residual_tolerance,
            result.block_encoding_check if "result" in locals() else None,
        ),
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


def _skipped_after_timeout(
    row: dict[str, str],
    *,
    check_mode: str,
    block_encoding_epsilon: float | None,
    block_encoding_residual_tolerance: float | None,
) -> dict[str, str]:
    return {
        **_metadata(row),
        **_check_metadata(check_mode, block_encoding_epsilon, block_encoding_residual_tolerance, None),
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
        "gate_count",
        "qasm_path",
        "check_mode",
        "check_epsilon",
        "inferred_alpha",
        "residual_norm",
        "residual_threshold",
        "residual_numerical_tolerance",
        "residual_acceptance_threshold",
        "runtime_sec",
        "tracemalloc_peak_mb",
        "process_maxrss_mb",
        "status",
        "success",
        "error",
    ]


def _check_metadata(
    check_mode: str,
    block_encoding_epsilon: float | None,
    block_encoding_residual_tolerance: float | None,
    check,
) -> dict[str, str]:
    if check_mode != "paper":
        return {
            "check_mode": check_mode,
            "check_epsilon": "",
            "inferred_alpha": "",
            "residual_norm": "",
            "residual_threshold": "",
            "residual_numerical_tolerance": "",
            "residual_acceptance_threshold": "",
        }
    if check is None:
        return {
            "check_mode": check_mode,
            "check_epsilon": f"{block_encoding_epsilon:.12g}",
            "inferred_alpha": "",
            "residual_norm": "",
            "residual_threshold": "",
            "residual_numerical_tolerance": f"{block_encoding_residual_tolerance:.12g}",
            "residual_acceptance_threshold": "",
        }
    return {
        "check_mode": check_mode,
        "check_epsilon": f"{block_encoding_epsilon:.12g}",
        "inferred_alpha": f"{check.alpha:.12g}",
        "residual_norm": f"{check.residual_norm:.12g}",
        "residual_threshold": f"{check.threshold:.12g}",
        "residual_numerical_tolerance": f"{check.numerical_tolerance:.12g}",
        "residual_acceptance_threshold": f"{check.acceptance_threshold:.12g}",
    }


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
