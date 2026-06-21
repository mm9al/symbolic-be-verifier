#!/usr/bin/env python3
"""Generate pyqsp Hamiltonian-simulation polynomials and QSP phases.

Examples:
    python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --component cos
    python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --component both --format json
    python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --component exp --no-angles
    python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --component cos-selector --qasm-snippet
    python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --write-examples
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from symbolic.qsp import (  # noqa: E402
    generate_component,
    generate_exp_component,
    generate_selector_component,
    polynomial_expr,
    qasm_snippet,
    selector_qasm_snippet,
    write_selector_examples,
)


def _print_text(
    records: list[dict[str, Any]],
    *,
    variable: str,
    precision: int,
    include_logs: bool,
    include_qasm_snippet: bool,
    qasm_args: argparse.Namespace,
) -> None:
    for index, record in enumerate(records):
        if index:
            print()

        print(f"[{record['component']}]")
        print(f"tau = {record['tau']}")
        print(f"epsilon = {record['epsilon']}")
        print(f"ensure_bounded = {record['ensure_bounded']}")
        print(f"scale = {record['scale']}")
        print(f"degree = {record['degree']}")
        if record["component"] == "exp":
            _print_exp_record(record, variable=variable, precision=precision, include_logs=include_logs)
            if include_qasm_snippet:
                print("QASM snippet skipped: pyqsp hamsim provides separate cos/sin phase sequences, not one exp sequence.")
            continue

        if "parity" in record:
            print(f"parity = {record['parity']}")

        cheb = record["chebyshev_coefficients"]
        mono = record["monomial_coefficients"]
        print(f"Chebyshev coefficients = {cheb}")
        print(f"monomial coefficients  = {mono}")
        print("monomial polynomial   = " + polynomial_expr(mono, variable, precision))

        if "pyqsp_phases" in record:
            print(f"phase count = {record['phase_count']}")
            print(f"QSP degree = {record['qsp_degree']}")
            print(f"pyqsp phases = {record['pyqsp_phases']}")
            print(f"QSVT projector phases = {record['qsvt_projector_phases']}")
            print(f"pyqsp response part = {record['pyqsp_response_part']}")
            if "selector_output" in record:
                print(f"selector output = {record['selector_output']}")
            print(f"qasm rz angles = {record['qasm_rz_angles']}")

        _print_pyqsp_log(record, include_logs=include_logs)

        if include_qasm_snippet:
            _print_qasm_snippet(record, qasm_args)


def _print_exp_record(record: dict[str, Any], *, variable: str, precision: int, include_logs: bool) -> None:
    cos_mono = record["cos"]["monomial_coefficients"]
    sin_mono = record["sin"]["monomial_coefficients"]
    cos_expr = polynomial_expr(cos_mono, variable, precision)
    sin_expr = polynomial_expr(sin_mono, variable, precision)
    print("target = scale * exp(-i * tau * x)")
    print(f"cos monomial coefficients = {cos_mono}")
    print(f"sin monomial coefficients = {sin_mono}")
    print(f"complex monomial coefficients = {record['monomial_coefficients_complex']}")
    print(f"exp polynomial = ({cos_expr}) - I*({sin_expr})")
    _print_pyqsp_log(record, include_logs=include_logs)


def _print_pyqsp_log(record: dict[str, Any], *, include_logs: bool) -> None:
    if include_logs and record["pyqsp_log"]:
        print("pyqsp log:")
        for line in record["pyqsp_log"]:
            print(f"  {line}")


def _print_qasm_snippet(record: dict[str, Any], qasm_args: argparse.Namespace) -> None:
    if "pyqsp_phases" not in record:
        print("QASM snippet skipped because --no-angles was used.")
        return

    block_ancillas = _block_ancillas_from_args(qasm_args)
    system_qubits = _system_qubits_from_args(qasm_args)
    print("QASM snippet:")
    if "selector_component" in record:
        print(
            selector_qasm_snippet(
                record["pyqsp_phases"],
                record["qsvt_projector_phases"],
                record["qasm_rz_angles"],
                selector_qubit=qasm_args.selector_qubit,
                phase_qubit=qasm_args.phase_qubit,
                block_ancillas=block_ancillas,
                system_qubits=system_qubits,
                signal_gate=qasm_args.signal_gate,
                signal_gate_dagger=qasm_args.signal_gate_dagger,
            )
        )
        return

    print(
        qasm_snippet(
            record["pyqsp_phases"],
            record["qsvt_projector_phases"],
            record["qasm_rz_angles"],
            phase_qubit=qasm_args.phase_qubit,
            block_ancillas=block_ancillas,
            system_qubits=system_qubits,
            signal_gate=qasm_args.signal_gate,
            signal_gate_dagger=qasm_args.signal_gate_dagger,
        )
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate Hamiltonian-simulation polynomials and QSP phase angles using pyqsp."
    )
    parser.add_argument("--tau", type=float, required=True, help="Simulation time.")
    parser.add_argument(
        "--epsilon",
        type=float,
        default=1e-4,
        help="Hamiltonian-simulation polynomial truncation target.",
    )
    parser.add_argument(
        "--component",
        choices=["cos", "sin", "both", "exp", "cos-selector", "sin-selector"],
        default="cos",
        help="Which Hamiltonian-simulation component to generate.",
    )
    parser.add_argument(
        "--unbounded",
        action="store_true",
        help="Disable pyqsp's default 0.5 rescaling of the target polynomial.",
    )
    parser.add_argument(
        "--method",
        choices=["sym_qsp", "laurent"],
        default="sym_qsp",
        help="QSP angle-finding method.",
    )
    parser.add_argument(
        "--signal-operator",
        choices=["Wx", "Wz"],
        default="Wx",
        help="QSP signal convention for angle generation.",
    )
    parser.add_argument("--no-angles", action="store_true", help="Only generate polynomial coefficients.")
    parser.add_argument("--format", choices=["text", "json"], default="text", help="Output format.")
    parser.add_argument(
        "--precision",
        type=int,
        default=17,
        help="Significant digits for the printed polynomial expression.",
    )
    parser.add_argument("--variable", default="x", help="Variable name used in the printed monomial polynomial.")
    parser.add_argument(
        "--show-pyqsp-log",
        action="store_true",
        help="Show pyqsp truncation and iteration diagnostics in text output.",
    )
    parser.add_argument(
        "--qasm-snippet",
        action="store_true",
        help="Print a QASM phase/U sequence snippet after the data.",
    )
    parser.add_argument(
        "--write-examples",
        action="store_true",
        help="Write selector-wrapped cos/sin QASM examples and expected polynomial metadata under --examples-dir.",
    )
    parser.add_argument(
        "--examples-dir",
        type=Path,
        default=Path("examples"),
        help="Directory where --write-examples creates qsp_hamsim_* folders.",
    )
    parser.add_argument("--phase-qubit", default="q[0]")
    parser.add_argument("--block-ancilla", default="q[1]")
    parser.add_argument("--block-ancillas", nargs="+", help="Block-encoding ancilla qubits, e.g. q[1] q[2].")
    parser.add_argument("--system-qubit", default="q[2]")
    parser.add_argument("--system-qubits", nargs="+", help="System qubits, e.g. q[4] q[5].")
    parser.add_argument("--selector-qubit", default="q[3]")
    parser.add_argument("--signal-gate", default="UH")
    parser.add_argument("--signal-gate-dagger", default="UHdg")

    args = parser.parse_args()
    if args.epsilon <= 0:
        parser.error("--epsilon must be positive")
    if not math.isfinite(args.tau):
        parser.error("--tau must be finite")
    if args.precision <= 0:
        parser.error("--precision must be positive")
    if args.qasm_snippet and args.component in {"both", "exp"}:
        parser.error("--qasm-snippet requires --component cos, sin, cos-selector, or sin-selector")
    if args.write_examples and args.no_angles:
        parser.error("--write-examples requires angle generation")
    return args


def main() -> int:
    args = parse_args()
    if args.write_examples:
        block_ancillas = _block_ancillas_from_args(args)
        system_qubits = _system_qubits_from_args(args)
        metadata = write_selector_examples(
            tau=args.tau,
            epsilon=args.epsilon,
            ensure_bounded=not args.unbounded,
            method=args.method,
            signal_operator=args.signal_operator,
            precision=args.precision,
            examples_dir=args.examples_dir,
            selector_qubit=args.selector_qubit,
            phase_qubit=args.phase_qubit,
            block_ancillas=block_ancillas,
            system_qubits=system_qubits,
            signal_gate=args.signal_gate,
            signal_gate_dagger=args.signal_gate_dagger,
        )
        if args.format == "json":
            print(json.dumps(metadata, indent=2, sort_keys=True))
        else:
            print(f"Wrote metadata: {metadata['metadata']}")
            for file_record in metadata["files"]:
                print(f"Wrote {file_record['component']}: {file_record['qasm']}")
                print(f"  polynomial = {file_record['polynomial']}")
        return 0

    records = _generate_records(args)
    if args.format == "json":
        output: Any = records[0] if len(records) == 1 else records
        print(json.dumps(output, indent=2, sort_keys=True))
    else:
        _print_text(
            records,
            variable=args.variable,
            precision=args.precision,
            include_logs=args.show_pyqsp_log,
            include_qasm_snippet=args.qasm_snippet,
            qasm_args=args,
        )
    return 0


def _generate_records(args: argparse.Namespace) -> list[dict[str, Any]]:
    if args.component in {"cos-selector", "sin-selector"}:
        return [
            generate_selector_component(
                component=args.component.removesuffix("-selector"),
                tau=args.tau,
                epsilon=args.epsilon,
                ensure_bounded=not args.unbounded,
                compute_angles=not args.no_angles,
                method=args.method,
                signal_operator=args.signal_operator,
            )
        ]

    if args.component == "exp":
        return [
            generate_exp_component(
                tau=args.tau,
                epsilon=args.epsilon,
                ensure_bounded=not args.unbounded,
            )
        ]

    components = ["cos", "sin"] if args.component == "both" else [args.component]
    return [
        generate_component(
            component=component,
            tau=args.tau,
            epsilon=args.epsilon,
            ensure_bounded=not args.unbounded,
            compute_angles=not args.no_angles,
            method=args.method,
            signal_operator=args.signal_operator,
        )
        for component in components
    ]


def _block_ancillas_from_args(args: argparse.Namespace) -> list[str]:
    return args.block_ancillas if args.block_ancillas is not None else [args.block_ancilla]


def _system_qubits_from_args(args: argparse.Namespace) -> list[str]:
    return args.system_qubits if args.system_qubits is not None else [args.system_qubit]


if __name__ == "__main__":
    raise SystemExit(main())
