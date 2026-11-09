"""Rebuild saved comparisons and README tables without rerunning circuits."""
from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from evaluation.pennylane.common import sha256
from symbolic.interval import acb, arb, ctx


def load_results(directory: Path) -> tuple[dict, dict]:
    """Reject mixed CSV/JSON generations before using their timings."""
    rows, reports = {}, {}
    with (directory/"block_encoding_results.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            n = int(row["n_system"])
            if n in rows:
                raise ValueError(f"Duplicate size {n}: {directory}")
            rows[n] = row
            if not row.get("report_path"):
                continue
            report = json.loads((directory/row["report_path"]).read_text())
            if any(value != str(report["summary"].get(key, ""))
                   for key, value in row.items() if value != ""):
                raise ValueError(f"CSV/JSON mismatch: {directory}, n={n}")
            reports[n] = report
    return rows, reports


def parse_coefficient(text: str):
    """Read Acb's recorded real + imaginary display without binary64 casts."""
    if not text.endswith("j"):
        return acb(arb(text))
    value = text[:-1]
    if " + " in value:
        real, imaginary = value.rsplit(" + ", 1)
        return acb(arb(real), arb(imaginary))
    if " - " in value:
        real, imaginary = value.rsplit(" - ", 1)
        return acb(arb(real), -arb(imaginary))
    return acb(0, arb(value))


def check_pair(left: dict, right: dict) -> None:
    for key in ("qasm_sha256", "metadata_sha256", "target_H", "ancillas", "systems", "global_phase_gamma"):
        if left[key] != right[key]:
            raise ValueError(f"Comparison input mismatch: {key}")
    for key in ("precision", "check_epsilon", "hermiticity_delta"):
        if left["summary"][key] != right["summary"][key]:
            raise ValueError(f"Comparison setting mismatch: {key}")
    with ctx.workdps(int(right["summary"]["precision"])):
        error = arb(right.get("bounded_pruning", {}).get("operator_error_upper", "0")).abs_upper()
        a, b = left["actual_block_coefficients"], right["actual_block_coefficients"]
        for word in a.keys() | b.keys():
            difference = parse_coefficient(a.get(word, "0"))-parse_coefficient(b.get(word, "0"))
            compatible = difference.contains(0) if error.is_zero() else difference.abs_lower() <= error
            if not compatible:
                raise ValueError(f"Coefficient enclosures incompatible with recorded error: {word}")


def snapshot_checks(directory: Path) -> dict:
    metadata = json.loads((directory/"run_metadata.json").read_text())
    snapshots = metadata["source_snapshots"]
    if not snapshots:
        raise ValueError("No recorded source snapshots")
    checks = {name: sha256(directory/path) == metadata["sources_sha256"][name]
              for name, path in snapshots.items()}
    if not all(checks.values()):
        raise ValueError("Source snapshot hash differs from measurement metadata")
    return checks


def csv_text(rows: list[dict]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def build_outputs(base: Path = HERE) -> dict[Path, str]:
    raw_dir = base/"results_local_summaries"/"baseline"
    local_dir, sparse_dir = base/"results_local_summaries", base/"results_bounded_pruning"
    raw_rows, raw = load_results(raw_dir)
    local_rows, local = load_results(local_dir)
    sparse_rows, sparse = load_results(sparse_dir)
    comparisons, pruning = [], []
    paired = {"local_summaries": [], "bounded_pruning": []}
    for n in sorted(raw.keys() & local.keys()):
        a, b = raw_rows[n], local_rows[n]
        check_pair(raw[n], local[n])
        paired["local_summaries"].append(n)
        comparisons.append(dict(n_system=n, source_gate_count=b["gate_count"], executed_steps=b["executed_steps"],
            local_summary_count=b["local_summary_count"], baseline_runtime_sec=a["runtime_sec"],
            optimized_runtime_sec=b["runtime_sec"], speedup=float(a["runtime_sec"])/float(b["runtime_sec"]),
            baseline_status=a["status"], optimized_status=b["status"], qasm_sha256=local[n]["qasm_sha256"]))
    for n in sorted(sparse):
        a, b = local_rows.get(n, {}), sparse_rows[n]
        comparable = n in local
        if comparable:
            check_pair(local[n], sparse[n])
            paired["bounded_pruning"].append(n)
        pruning.append(dict(n=n, local_summary_sec=a["runtime_sec"] if comparable else "",
            bounded_sec=b["runtime_sec"], speedup=float(a["runtime_sec"])/float(b["runtime_sec"]) if comparable else "",
            local_summary_final_terms=len(local[n]["actual_block_coefficients"]) if comparable else "",
            bounded_final_terms=b["final_pauli_terms"], status=b["status"], local_summary_status=a.get("status", ""),
            propagation_error_upper=b["truncation_error_upper"], same_input_and_settings=True if comparable else "",
            coefficient_enclosures_compatible_with_error_bound=True if comparable else ""))
    outputs = {local_dir/"comparison.csv": csv_text(comparisons), sparse_dir/"comparison.csv": csv_text(pruning)}
    checks = dict(generator="evaluation/pennylane/refresh_comparisons.py", reran_circuits=False,
        timing_scope="Saved per-case runtime_sec; one run per mode, no repetitions; baseline and optimized runs are separate batches",
        matched_sizes=paired, csv_json_summaries_match=True, paired_inputs_and_settings_match=True,
        paired_coefficient_enclosures_compatible=True,
        coefficient_check_scope="Regression compatibility of serialized coefficient balls with eta; not a new equivalence proof",
        source_snapshot_hashes_match_measurement_metadata=snapshot_checks(sparse_dir),
        input_artifacts_sha256={})
    for directory in (raw_dir, local_dir, sparse_dir):
        for path in [directory/"block_encoding_results.csv", directory/"run_metadata.json",
                     *(directory/row["report_path"] for row in load_results(directory)[0].values() if row.get("report_path"))]:
            checks["input_artifacts_sha256"][str(path.relative_to(base))] = sha256(path)
    validation_path = sparse_dir/"validation_checks.json"
    validation = json.loads(validation_path.read_text())
    validation["source_snapshot_scope"] = "The two saved mode sources match the hashes of this directory's recorded measurement run"
    validation["source_snapshot_hashes_match_measurement_metadata"] = True
    validation["source_snapshot_checks"] = checks["source_snapshot_hashes_match_measurement_metadata"]
    validation["comparison_checks"] = "../comparison_checks.json"
    outputs[validation_path] = json.dumps(validation, indent=2, allow_nan=False)+"\n"
    outputs[base/"comparison_checks.json"] = json.dumps(checks, indent=2, allow_nan=False)+"\n"
    local_table = ["| n | 原本逐閘 (s) | 局部摘要 (s) | 加速 | 原始 gates → 執行 steps |",
                   "|---|---:|---:|---:|---:|"]
    for row in comparisons:
        local_table.append(f"| {row['n_system']} | {float(row['baseline_runtime_sec']):.3f} | {float(row['optimized_runtime_sec']):.3f} | {row['speedup']:.1f}x | {row['source_gate_count']} → {row['executed_steps']} |")
    sparse_table = ["| n | 既有局部摘要 (s) | 新模式 (s) | 最後 Pauli 項數：既有 → 新模式 |", "|---|---:|---:|---:|"]
    for row in pruning:
        if row["n"] < 5:
            continue
        baseline = f"{float(row['local_summary_sec']):.3f}" if row["local_summary_sec"] else row["local_summary_status"]
        sparse_table.append(f"| {row['n']} | {baseline} | {float(row['bounded_sec']):.3f} | {row['local_summary_final_terms'] or '—'} → {row['bounded_final_terms']} |")
    readme_path = base/"README.md"
    readme = readme_path.read_text()
    for name, table in (("local-summary-comparison", local_table), ("bounded-pruning-comparison", sparse_table)):
        begin, end = f"<!-- BEGIN {name} -->", f"<!-- END {name} -->"
        before, rest = readme.split(begin, 1)
        _, after = rest.split(end, 1)
        readme = before+begin+"\n"+"\n".join(table)+"\n"+end+after
    outputs[readme_path] = readme
    return outputs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Fail if checked-in comparisons or README tables are stale")
    args = parser.parse_args()
    outputs = build_outputs()
    stale = [path for path, text in outputs.items() if not path.exists() or path.read_bytes() != text.encode()]
    if args.check:
        if stale:
            print("Stale outputs: "+", ".join(str(path.relative_to(HERE)) for path in stale))
            return 1
    else:
        for path, text in outputs.items():
            path.write_text(text)
    print("Comparisons, source snapshots, CSV/JSON summaries and README tables are consistent")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
