"""TASK 1 evaluation of native PennyLane periodic-Ising block encodings."""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import importlib.metadata
import json
import math
import multiprocessing as mp
from pathlib import Path
import platform
import re
import resource
import sys
import time
import tracemalloc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from evaluation.pennylane.common import (BENCHMARKS, RESULTS, ROOT, gate, inverse, ising_terms, parse_sizes,
                                        qasm_body, sha256, signatures, target_expression, write_json)


DEFAULT_HERMITICITY_DELTA = "1e-12"

FIELDS = ["benchmark_id", "model", "n_system", "n_terms", "target_lcu_terms", "n_ancilla",
          "nominal_alpha", "gate_count", "precision", "check_epsilon", "inferred_alpha",
          "residual_norm_upper", "residual_threshold", "operator_error_upper",
          "core_proportionality_pass", "proportionality_pass", "hermiticity_delta",
          "max_imaginary_coefficient_upper", "max_imaginary_coefficient_pauli",
          "exact_hermiticity_status", "hermiticity_status",
          "runtime_sec", "symbolic_runtime_sec", "hermiticity_runtime_sec",
          "tracemalloc_peak_mb", "process_maxrss_mb", "status", "success", "error", "report_path"]
SUMMARY_FIELDS = ["execution_mode", "local_summary_count", "summarized_source_gates", "executed_steps",
                  "removed_identity_gates", "summary_compilation_sec"]
PRUNING_FIELDS = ["pruning_threshold", "pruning_budget", "truncation_error_upper",
                  "scaled_truncation_error_upper", "pruned_terms", "final_pauli_terms",
                  "peak_total_branch_terms", "peak_terms_per_branch", "roundoff_error_upper"]


def validate_inputs(base: Path, row: dict) -> tuple[dict, dict, Path]:
    metadata_path = base / row["metadata_path"]
    if sha256(metadata_path) != row["metadata_sha256"]:
        raise ValueError("Metadata SHA-256 mismatch")
    metadata = json.loads(metadata_path.read_text())
    n = int(row["n_system"])
    expected = dict(benchmark_id=f"ising_periodic_n{n}", model="ising", graph="periodic",
                    n_system=n, target_lcu_terms=2*n, target_H=target_expression(n),
                    independent_terms=ising_terms(n), nominal_alpha=str(2*n))
    if any(metadata.get(key) != value for key, value in expected.items()):
        raise ValueError("Inputs differ from the independent periodic-Ising specification")
    for name, digest in metadata["files_sha256"].items():
        if sha256(metadata_path.parent / name) != digest:
            raise ValueError(f"Input SHA-256 mismatch: {name}")
    m = metadata["n_ancilla"]
    if (metadata["ancillas"] != list(range(m)) or metadata["systems"] != list(range(m, m+n))
            or m != (metadata["n_terms"] - 1).bit_length()):
        raise ValueError("Unexpected ancilla/system registers")
    qasm = base / row["qasm_path"]
    if qasm.resolve() != (metadata_path.parent / "block_encoding.qasm").resolve():
        raise ValueError("QASM path differs from the hashed artifact")
    structure = json.loads((metadata_path.parent / "structure.json").read_text())
    if qasm_body(qasm) != structure["prepare"] + structure["select"] + structure["unprepare"]:
        raise ValueError("Structure segments do not match the actual QASM")
    gamma = sum(map(Fraction, structure["global_phases"]), Fraction())
    if (gamma != Fraction(structure["gamma"]) or gamma != Fraction(metadata["global_phase_gamma"])
            or gamma != sum((sum(map(Fraction, group), Fraction())
                             for group in structure["segment_global_phases"]), Fraction())):
        raise ValueError("Global-phase metadata is inconsistent")
    return metadata, structure, qasm


def structural_hermiticity(structure: dict, ancillas: list[int], systems: list[int]) -> dict:
    """Certify the actual decimal P†SP using exact SELECT branch execution."""
    import sympy as sp
    from symbolic.qasm_parser import parse_qasm_text
    from symbolic.verify import run_circuit

    prepare, select, unprepare = (structure[key] for key in ("prepare", "select", "unprepare"))
    if signatures(unprepare) != signatures(inverse(prepare)):
        return dict(status="UNPROVEN", reason="Native unprepare is not the exact inverse of PREPARE")
    for line in prepare:
        wires = set(map(int, re.findall(r"q\[(\d+)\]", line)))
        if not wires <= set(ancillas) or not line.startswith(("ry(", "cx ")):
            return dict(status="UNPROVEN", reason="PREPARE is outside the ancilla-only RY/CX proof")
    phases = [sum(map(Fraction, group), Fraction()) for group in structure["segment_global_phases"]]
    if len(phases) != 3 or phases[0] + phases[2] or phases[1]:
        return dict(status="UNPROVEN", reason="Global phases do not certify a Hermitian sandwich")
    gates = parse_qasm_text("\n".join(select))
    if any(g.parameter is not None and g.parameter != 0 for g in gates):
        return dict(status="UNPROVEN", reason="Native SELECT contains nonzero decimal rotations")
    zero = (0,) * len(ancillas)
    branches = []
    for index in range(2 ** len(ancillas)):
        bits = format(index, f"0{len(ancillas)}b")
        flips = [gate("x", [wire]) for bit, wire in zip(bits, ancillas) if bit == "1"]
        state = run_circuit(parse_qasm_text("\n".join(flips+select+list(reversed(flips)))),
                            ancillas=ancillas, systems=systems).final_state
        if any(key != zero and not value.is_zero() for key, value in state.branches.items()):
            return dict(status="UNPROVEN", reason="SELECT does not preserve every ancilla basis branch")
        actual = state.top_left()
        if any(sp.simplify(sp.im(coefficient)) != 0 for coefficient in actual.terms.values()):
            return dict(status="UNPROVEN", reason="SELECT branch is not exactly Hermitian")
        branches.append(dict(index=index, operator=str(actual)))
    return dict(status="PASS_STRUCTURE", proof="Exact native P†SP with Hermitian diagonal SELECT blocks",
                select_branches=branches)


def classify_exact_hermiticity(actual, structure_proof: dict) -> dict:
    from symbolic.interval import to_acb
    violations = [{"pauli": "".join(word), "imaginary_coefficient": str(to_acb(c).imag)}
                  for word, c in actual.terms.items() if not to_acb(c).imag.contains(0)]
    if violations:
        return dict(status="FAIL_NON_HERMITIAN", proof="A Pauli coefficient's imaginary ball excludes zero",
                    violations=violations, structure=structure_proof)
    if structure_proof["status"] == "PASS_STRUCTURE":
        return structure_proof
    if all(to_acb(c).imag.is_zero() for c in actual.terms.values()):
        return dict(status="PASS_EXACT_COEFFICIENTS", proof="Every imaginary coefficient is exactly zero",
                    structure=structure_proof)
    return dict(status="UNPROVEN", reason="Zero-containing imaginary intervals do not prove exact Hermiticity",
                structure=structure_proof)


def classify_hermiticity(actual, structure_proof: dict,
                         delta: str = DEFAULT_HERMITICITY_DELTA) -> dict:
    """Accept exact Hermiticity or a rigorous max-Pauli-coefficient bound < delta."""
    from symbolic.interval import arb, to_acb

    delta_fraction = Fraction(delta)
    if delta_fraction < 0:
        raise ValueError("Hermiticity delta must be nonnegative")
    delta_ball = arb(str(delta_fraction))
    bounds = {"".join(word): to_acb(c).imag.abs_upper() for word, c in actual.terms.items()}
    worst_pauli = max(bounds, key=bounds.get) if bounds else ""
    maximum = bounds.get(worst_pauli, arb(0))
    exact = classify_exact_hermiticity(actual, structure_proof)
    within_delta = bool(delta_fraction > 0 and maximum < delta_ball.lower())
    if exact["status"].startswith("PASS_"):
        status = exact["status"]
    elif within_delta:
        status = "PASS_WITHIN_DELTA"
    elif delta_fraction > 0:
        status = "FAIL_DELTA"
    else:
        status = exact["status"]
    return dict(status=status, exact_status=exact["status"], exact=exact,
                delta=delta, delta_ball=str(delta_ball), within_delta=within_delta,
                max_imaginary_coefficient_upper=str(maximum),
                max_imaginary_coefficient_pauli=worst_pauli,
                criterion="Exact Hermiticity OR max_P |Im(c_P)| < delta, using rigorous upper bounds on B",
                approximate_acceptance=status == "PASS_WITHIN_DELTA")


def evaluate(base: Path, row: dict, output: Path, epsilon: str, precision: int,
             hermiticity_delta: str = DEFAULT_HERMITICITY_DELTA,
             local_summaries: bool = False, bounded_pruning: bool = False,
             pruning_threshold: str = "1e-24", pruning_budget: str = "1e-14") -> dict:
    from symbolic.expr import parse_operator_expression
    from symbolic.interval import acb, arb, ctx, to_acb
    from symbolic.verify import block_encoding_proportionality_check, verify_qasm_file

    started = time.perf_counter()
    tracemalloc.start()
    summary = {key: row[key] for key in FIELDS if key in row}
    summary.update(status="ERROR", success=False, check_epsilon=epsilon,
                   hermiticity_delta=hermiticity_delta)
    report = {}
    try:
        metadata, structure, qasm = validate_inputs(base, row)
        n = metadata["n_system"]
        precision = max(100, precision, 2*n)
        summary["precision"] = precision
        with ctx.workdps(precision):
            target = parse_operator_expression(target_expression(n), num_qubits=n)
            tick = time.perf_counter()
            summary_report = None
            pruning_report = None
            if local_summaries or bounded_pruning:
                if bounded_pruning:
                    from evaluation.pennylane.bounded_pruning import verify_with_pruning
                    result, summary_report, pruning_report = verify_with_pruning(
                        qasm, ancillas=metadata["ancillas"], systems=metadata["systems"], precision=precision,
                        threshold=pruning_threshold, budget=pruning_budget)
                else:
                    from evaluation.pennylane.local_summaries import verify_with_summaries
                    result, summary_report = verify_with_summaries(
                        qasm, ancillas=metadata["ancillas"], systems=metadata["systems"], precision=precision)
                summary.update(execution_mode="bounded_pruning" if bounded_pruning else "exact_local_summaries",
                    local_summary_count=summary_report["summary_count"],
                    summarized_source_gates=summary_report["summarized_source_gates"],
                    executed_steps=summary_report["executed_steps"],
                    removed_identity_gates=summary_report["removed_identity_gates"],
                    summary_compilation_sec=summary_report["summary_compilation_sec"])
            else:
                result = verify_qasm_file(qasm, ancillas=metadata["ancillas"], systems=metadata["systems"],
                                          keep_trace=False, qsp_calculation_precision=precision)
            gamma = arb(str(Fraction(structure["gamma"])))
            actual = result.final_state.top_left().scale((-acb(0, 1)*gamma).exp())
            summary["symbolic_runtime_sec"] = time.perf_counter()-tick
            fit = block_encoding_proportionality_check(actual, target, epsilon=float(Fraction(epsilon)))
            # Recheck a concrete decimal alpha, rather than accepting its midpoint
            # or using the generator's nominal normalization as a certificate.
            alpha_text = fit.alpha_ball.mid().str(precision, radius=False)
            alpha = arb(alpha_text)
            difference = actual.scale(alpha)-target
            residual = sum((to_acb(c).abs_upper()**2 for c in difference.terms.values()), arb(0)).sqrt().abs_upper()
            dimension_factor = arb(2)**n
            dimension_factor = dimension_factor.sqrt()
            threshold = arb(str(Fraction(epsilon)))/dimension_factor
            proportionality_pass = bool(alpha.lower() > 0 and residual <= threshold.lower())
            operator_error = (dimension_factor*residual).abs_upper()
            truncation_error = arb(0)
            if pruning_report is not None:
                from evaluation.pennylane.bounded_pruning import residual_certificate
                truncation_error = arb(pruning_report["operator_error_upper"]).abs_upper()
                certificate = residual_certificate(residual, alpha, dimension_factor, epsilon, truncation_error)
                proportionality_pass = certificate["passed"]
                operator_error = certificate["operator_error_upper"]
                summary.update(pruning_threshold=pruning_threshold, pruning_budget=pruning_budget,
                    truncation_error_upper=str(truncation_error),
                    scaled_truncation_error_upper=str(certificate["scaled_truncation_error_upper"]),
                    pruned_terms=pruning_report["dropped_terms"],
                    roundoff_error_upper=pruning_report["roundoff_error_upper"],
                    final_pauli_terms=pruning_report["final_pauli_terms"],
                    peak_total_branch_terms=pruning_report["peak_total_branch_terms"],
                    peak_terms_per_branch=pruning_report["peak_terms_per_branch"])
            tick = time.perf_counter()
            proof = structural_hermiticity(structure, metadata["ancillas"], metadata["systems"])
            if pruning_report is not None:
                from evaluation.pennylane.bounded_pruning import classify_with_pruning
                hermiticity = classify_with_pruning(actual, proof, hermiticity_delta, truncation_error)
            else:
                hermiticity = classify_hermiticity(actual, proof, hermiticity_delta)
            summary["hermiticity_runtime_sec"] = time.perf_counter()-tick
            success = proportionality_pass and hermiticity["status"].startswith("PASS_")
            summary.update(status="ACCEPT" if success else "REJECT", success=success,
                inferred_alpha=alpha_text, residual_norm_upper=str(residual), residual_threshold=str(threshold),
                operator_error_upper=str(operator_error),
                core_proportionality_pass=fit.success, proportionality_pass=proportionality_pass,
                max_imaginary_coefficient_upper=hermiticity["max_imaginary_coefficient_upper"],
                max_imaginary_coefficient_pauli=hermiticity["max_imaginary_coefficient_pauli"],
                exact_hermiticity_status=hermiticity["exact_status"],
                hermiticity_status=hermiticity["status"])
            if not success:
                summary["error"] = (f"Hermiticity policy failed (delta={hermiticity_delta})" if not hermiticity["status"].startswith("PASS_")
                                    else "Coefficient proportionality threshold exceeded")
            report = dict(target_H=target_expression(n), ancillas=metadata["ancillas"], systems=metadata["systems"],
                hermiticity=hermiticity, alpha_fit_ball=str(fit.alpha_ball),
                alpha_policy="Infer alpha, choose a concrete decimal candidate, recompute the residual",
                acceptance_rule="(Exact Hermiticity OR max_P |Im(c_P)| < delta) AND alpha>0 AND ||alpha*c-h||_2 <= epsilon_BE/2^(n/2)",
                actual_block_coefficients={"".join(word): str(to_acb(c)) for word, c in actual.terms.items()},
                global_phase_gamma=structure["gamma"], metadata_sha256=row["metadata_sha256"],
                qasm_sha256=sha256(qasm))
            if summary_report is not None:
                report["local_summaries"] = summary_report
            if pruning_report is not None:
                report["bounded_pruning"] = pruning_report
                report["coefficient_scope"] = "Sparse approximate block B_tilde; ||B_original-B_tilde|| <= eta (pruning plus charged recentering)"
                report["acceptance_rule"] = ("(Exact structure OR max_P |Im(c_tilde_P)|+eta < delta) AND alpha>0 "
                    "AND 2^(n/2)*||alpha*c_tilde-h||_2 + alpha*eta <= epsilon_BE")
                report["core_proportionality_scope"] = "Fit/check on B_tilde only; final proportionality_pass also includes alpha*eta"
    except Exception as exc:
        summary["error"] = f"{type(exc).__name__}: {exc}"
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    summary.update(runtime_sec=time.perf_counter()-started, tracemalloc_peak_mb=peak/(1024**2),
                   process_maxrss_mb=rss/(1024**2 if sys.platform == "darwin" else 1024))
    path = output / "cases" / f"{row['benchmark_id']}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    summary["report_path"] = str(path.relative_to(output))
    write_json(path, dict(**report, summary=summary))
    return summary


def worker(connection, base, row, output, epsilon, precision,
           hermiticity_delta=DEFAULT_HERMITICITY_DELTA, local_summaries=False,
           bounded_pruning=False, pruning_threshold="1e-24", pruning_budget="1e-14"):
    try:
        connection.send(evaluate(base, row, output, epsilon, precision, hermiticity_delta, local_summaries,
                                 bounded_pruning, pruning_threshold, pruning_budget))
    finally:
        connection.close()


def run_case(base, row, output, epsilon, precision, timeout,
             hermiticity_delta=DEFAULT_HERMITICITY_DELTA, local_summaries=False,
             bounded_pruning=False, pruning_threshold="1e-24", pruning_budget="1e-14"):
    if timeout == 0:
        return evaluate(base, row, output, epsilon, precision, hermiticity_delta, local_summaries,
                        bounded_pruning, pruning_threshold, pruning_budget)
    context = mp.get_context("spawn")
    reader, sender = context.Pipe(duplex=False)
    process = context.Process(target=worker, args=(sender, base, row, output, epsilon, precision,
                                                 hermiticity_delta, local_summaries, bounded_pruning,
                                                 pruning_threshold, pruning_budget))
    started = time.perf_counter()
    process.start()
    sender.close()
    try:
        if reader.poll(timeout):
            try:
                return reader.recv()
            except EOFError:
                return dict(benchmark_id=row["benchmark_id"], status="ERROR", success=False,
                            error="Worker exited before returning a result")
        return dict(**{key: row[key] for key in FIELDS if key in row}, status="TIMEOUT", success=False,
                    check_epsilon=epsilon, hermiticity_delta=hermiticity_delta,
                    runtime_sec=time.perf_counter()-started,
                    error=f"timeout after {timeout:g}s")
    finally:
        process.join(timeout=0.1)
        if process.is_alive():
            process.terminate()
        process.join()
        reader.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=BENCHMARKS/"manifest.csv")
    parser.add_argument("--sizes", help="Optional subset of generated sizes")
    parser.add_argument("--output-dir", type=Path,
                        help="Default: results/; summaries: results_local_summaries/; pruning: results_bounded_pruning/")
    parser.add_argument("--local-summaries", action="store_true",
                        help="Opt-in exact CCX/CCZ subcircuit summaries; original QASM and core verifier are unchanged")
    parser.add_argument("--bounded-pruning", action="store_true",
                        help="Opt-in certified sparse propagation; includes exact local summaries; separate output folder")
    parser.add_argument("--pruning-threshold", default="1e-24",
                        help="Only discard coefficients with rigorous magnitude upper bound <= this value")
    parser.add_argument("--pruning-budget", default="1e-14",
                        help="Shared operator norm budget on B for pruning and charged ball recentering; 0 disables both")
    parser.add_argument("--block-encoding-epsilon", default="1e-8", help="TASK 1 epsilon_BE (no TASK 2 budget split)")
    parser.add_argument("--hermiticity-delta", default=DEFAULT_HERMITICITY_DELTA,
                        help="Accept max_P |Im(c_P)| < delta on B; default 1e-12; 0 restores exact Hermiticity")
    parser.add_argument("--precision", type=int, default=100, help="Decimal digits; at least max(100,2*n) is used")
    parser.add_argument("--timeout-sec", type=float, default=1800, help="Per-case timeout; 0 disables isolation/timeouts")
    parser.add_argument("--continue-after-timeout", action="store_true")
    args = parser.parse_args()
    try:
        budget = Fraction(args.block_encoding_epsilon)
        if budget <= 0 or not math.isfinite(float(budget)) or float(budget) <= 0:
            raise ValueError("--block-encoding-epsilon must be finite and positive")
        if Fraction(args.hermiticity_delta) < 0:
            raise ValueError("--hermiticity-delta must be nonnegative")
        if Fraction(args.pruning_threshold) <= 0 or Fraction(args.pruning_budget) < 0:
            raise ValueError("Require positive --pruning-threshold and nonnegative --pruning-budget")
        sizes = None if args.sizes is None else parse_sizes(args.sizes)
    except (ValueError, ZeroDivisionError, OverflowError) as exc:
        parser.error(str(exc))
    if args.precision < 1 or args.timeout_sec < 0 or not math.isfinite(args.timeout_sec):
        parser.error("Require positive precision and a finite, nonnegative timeout")
    if not args.manifest.exists():
        parser.error(f"Missing {args.manifest}; run generate_benchmarks.py first")
    with args.manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if sizes is not None:
        missing = set(sizes)-{int(row["n_system"]) for row in rows}
        if missing:
            parser.error(f"Sizes absent from manifest: {sorted(missing)}")
        rows = [row for row in rows if int(row["n_system"]) in sizes]
    if not rows:
        parser.error("Manifest contains no selected inputs")
    rows.sort(key=lambda row: int(row["n_system"]))
    output = (args.output_dir or (HERE/"results_bounded_pruning" if args.bounded_pruning else
                                 HERE/"results_local_summaries" if args.local_summaries else RESULTS)).resolve()
    output.mkdir(parents=True, exist_ok=True)
    source_paths = [HERE/"run_block_encoding.py", HERE/"common.py", *sorted((ROOT/"symbolic").glob("*.py"))]
    if args.local_summaries or args.bounded_pruning:
        source_paths.append(HERE/"local_summaries.py")
    if args.bounded_pruning:
        source_paths.append(HERE/"bounded_pruning.py")
    write_json(output/"run_metadata.json", dict(manifest_sha256=sha256(args.manifest),
        command=sys.argv, python=platform.python_version(), platform=platform.platform(),
        block_encoding_epsilon=args.block_encoding_epsilon,
        execution_mode="bounded_pruning" if args.bounded_pruning else "exact_local_summaries" if args.local_summaries else "raw_gates",
        bounded_pruning=(dict(threshold=args.pruning_threshold, budget=args.pruning_budget,
            charged_ball_recentering=True, acceptance_includes_propagation_error=True) if args.bounded_pruning else None),
        hermiticity_policy=dict(delta=args.hermiticity_delta,
            criterion="Exact Hermiticity OR max_P |Im(c_P)| < delta on the final unscaled block B",
            rigorous_upper_bounds=True, comparison="strictly less than",
            exact_status_recorded_separately=True),
        packages={name: importlib.metadata.version(name) for name in ("sympy", "python-flint")},
        sources_sha256={str(p.relative_to(ROOT)): sha256(p) for p in source_paths},
        timing_scope="Per-case TASK 1 wall time, including input checks, symbolic execution, alpha/residual and Hermiticity; excludes generation and process startup",
        certificate_scope="Exact serialized decimal QASM with the recorded unconditional global phase"))
    counts, timed_out = {}, False
    with (output/"block_encoding_results.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS+
                                (SUMMARY_FIELDS if args.local_summaries or args.bounded_pruning else [])+
                                (PRUNING_FIELDS if args.bounded_pruning else []),
                                lineterminator="\n")
        writer.writeheader()
        for row in rows:
            if timed_out and not args.continue_after_timeout:
                result = dict(**{key: row[key] for key in FIELDS if key in row},
                              status="SKIPPED_AFTER_TIMEOUT", success=False,
                              hermiticity_delta=args.hermiticity_delta, error="Earlier Ising size timed out")
            else:
                result = run_case(args.manifest.resolve().parent, row, output, args.block_encoding_epsilon,
                                  args.precision, args.timeout_sec, args.hermiticity_delta, args.local_summaries,
                                  args.bounded_pruning, args.pruning_threshold, args.pruning_budget)
            timed_out |= result["status"] == "TIMEOUT"
            counts[result["status"]] = counts.get(result["status"], 0)+1
            writer.writerow(result)
            handle.flush()
            print(f"{row['benchmark_id']}: {result['status']} {result.get('hermiticity_status', '')}", flush=True)
    print(f"Results: {output/'block_encoding_results.csv'}; {counts}")
    return 1 if any(counts.get(key) for key in ("ERROR", "TIMEOUT")) else 0


if __name__ == "__main__":
    raise SystemExit(main())
