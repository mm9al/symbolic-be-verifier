from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "benchmarks" / "manifest.csv"
RESULTS_PATH = ROOT / "results" / "rq2_branch_term_growth_rerun.csv"

sys.path.insert(0, str(ROOT))

from evaluation import generate_benchmarks  # noqa: E402
from evaluation.run_symbolic import _expected  # noqa: E402
from symbolic.verify import verify_qasm_file  # noqa: E402


def main() -> int:
    args = _parse_args()
    if not MANIFEST_PATH.exists():
        generate_benchmarks.main()

    rows = [
        row
        for row in _read_manifest()
        if row["family"] == "branch_term_product" and int(row["n_ancilla"]) <= args.max_m
    ]
    if not rows:
        raise SystemExit("No RQ2 benchmarks selected.")

    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    with RESULTS_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_fieldnames())
        writer.writeheader()
        for row in rows:
            writer.writerow(_run(row))
            handle.flush()

    print(f"Wrote {RESULTS_PATH.relative_to(ROOT)}")
    return 0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a small RQ2 branch/term-growth scalability rerun.")
    parser.add_argument(
        "--max-m",
        type=int,
        default=4,
        help="Largest branch/term-growth ancilla count to rerun. Default: 4.",
    )
    return parser.parse_args()


def _read_manifest() -> list[dict[str, str]]:
    with MANIFEST_PATH.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _run(row: dict[str, str]) -> dict[str, str]:
    started = time.perf_counter()
    error = ""
    try:
        result = verify_qasm_file(
            ROOT / row["qasm_path"],
            expected=_expected(row),
            ancillas=_parse_int_list(row["ancillas"]),
            systems=_parse_int_list(row["systems"]),
        )
        status = result.status or "NO_EXPECTED"
        success = str(result.success)
    except Exception as exc:
        status = "ERROR"
        success = "False"
        error = f"{type(exc).__name__}: {exc}"
    runtime_sec = time.perf_counter() - started
    return {
        "benchmark_id": row["benchmark_id"],
        "family": row["family"],
        "n_system": row["n_system"],
        "n_ancilla": row["n_ancilla"],
        "expected_terms": row["expected_terms"],
        "qasm_path": row["qasm_path"],
        "runtime_sec": f"{runtime_sec:.9f}",
        "status": status,
        "success": success,
        "error": error,
    }


def _parse_int_list(text: str) -> tuple[int, ...]:
    return tuple(int(piece) for piece in text.split())


def _fieldnames() -> list[str]:
    return [
        "benchmark_id",
        "family",
        "n_system",
        "n_ancilla",
        "expected_terms",
        "qasm_path",
        "runtime_sec",
        "status",
        "success",
        "error",
    ]


if __name__ == "__main__":
    raise SystemExit(main())
