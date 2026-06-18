#!/usr/bin/env python3
"""Generate pyqsp Hamiltonian-simulation polynomials and QSP phases.

Examples:
    python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --component cos
    python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --component both --format json
    python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --component exp --no-angles
    python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --component cos --qasm-snippet

The polynomial coefficients are printed in both Chebyshev and monomial bases.
Use monomial coefficients for symbolic verifier expectations such as
    c0*I + c2*H^2 + c4*H^4
or, in one-variable checks,
    c0 + c2*x^2 + c4*x^4.

QASM snippets keep the pyqsp phases as data, but write physical OpenQASM
Rz angles after converting pyqsp's Wx phases to QSVT projector phases:
    psi_0 = phi_0 + pi/4
    psi_j = phi_j + pi/2, 1 <= j <= d - 1
    psi_d = phi_d + pi/4 for even d, and phi_d - pi/4 for odd d
and then theta_rz = -2 * psi.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import math
import os
import sys
from pathlib import Path
from typing import Any

try:
    import numpy as np
except ModuleNotFoundError as exc:
    raise SystemExit(
        "Missing dependency 'numpy'. Run this tool with an environment that has "
        "pyqsp's dependencies installed, for example from the repository root: "
        "python3 symbolic-be-verifier/tools/hamsim_qsp.py --tau 0.5"
    ) from exc


def _configure_import_path() -> None:
    """Make the sibling pyqsp checkout importable when it is not installed."""

    # Avoid pyqsp/matplotlib trying to build caches under the user's home dir.
    tmp_mpl = Path(os.environ.get("TMPDIR", "/tmp")) / "symbolic_be_verifier_mpl"
    tmp_mpl.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(tmp_mpl))

    repo_parent = Path(__file__).resolve().parents[2]
    local_pyqsp = repo_parent / "pyqsp"
    if local_pyqsp.exists():
        sys.path.insert(0, str(local_pyqsp))


_configure_import_path()

try:
    with contextlib.redirect_stderr(io.StringIO()):
        from pyqsp.angle_sequence import QuantumSignalProcessingPhases  # noqa: E402
        from pyqsp.poly import PolyCosineTX, PolySineTX  # noqa: E402
except ModuleNotFoundError as exc:
    raise SystemExit(
        "Could not import pyqsp or one of its dependencies. Install pyqsp's "
        "requirements, or keep the local pyqsp checkout next to "
        "symbolic-be-verifier."
    ) from exc


def _clean_float(value: Any, *, zero_tol: float = 1e-14) -> float:
    value = float(value)
    if abs(value) < zero_tol:
        return 0.0
    return value


def _clean_list(values: Any) -> list[float]:
    return [_clean_float(v) for v in values]


def _pad(values: list[float], length: int) -> list[float]:
    return values + [0.0] * (length - len(values))


def _format_number(value: float, precision: int) -> str:
    if value == 0.0:
        return "0"
    return f"{value:.{precision}g}"


def _polynomial_expr(coeffs: list[float], variable: str, precision: int) -> str:
    terms: list[str] = []
    for degree, coeff in enumerate(coeffs):
        if coeff == 0.0:
            continue
        abs_coeff = abs(coeff)
        if degree == 0:
            body = _format_number(abs_coeff, precision)
        elif degree == 1:
            body = f"{_format_number(abs_coeff, precision)}*{variable}"
        else:
            body = f"{_format_number(abs_coeff, precision)}*{variable}^{degree}"

        if not terms:
            terms.append(body if coeff > 0 else f"-{body}")
        else:
            terms.append(("+ " if coeff > 0 else "- ") + body)
    return " ".join(terms) if terms else "0"


def _complex_coefficients(real_coeffs: list[float], imag_coeffs: list[float]) -> list[dict[str, float]]:
    length = max(len(real_coeffs), len(imag_coeffs))
    real_padded = _pad(real_coeffs, length)
    imag_padded = _pad(imag_coeffs, length)
    return [
        {"real": _clean_float(real), "imag": _clean_float(imag)}
        for real, imag in zip(real_padded, imag_padded)
    ]


def pyqsp_phase_to_qsvt_projector_phase(phi: float, *, index: int, degree: int) -> float:
    if index == 0:
        shift = np.pi / 4
    elif index == degree:
        shift = np.pi / 4 if degree % 2 == 0 else -np.pi / 4
    else:
        shift = np.pi / 2
    return phi + shift


def qsvt_projector_phase_to_qasm_rz(psi: float) -> float:
    return -2.0 * psi


def qasm_phase_data(phases: list[float]) -> tuple[list[float], list[float]]:
    degree = len(phases) - 1
    qsvt_phases = [
        pyqsp_phase_to_qsvt_projector_phase(phi, index=index, degree=degree)
        for index, phi in enumerate(phases)
    ]
    rz_angles = [qsvt_projector_phase_to_qasm_rz(psi) for psi in qsvt_phases]
    return _clean_list(qsvt_phases), _clean_list(rz_angles)


def _phase_block(angle: float, phase_qubit: str, block_ancilla: str) -> list[str]:
    angle_text = f"{angle:.12g}"
    return [
        f"h {phase_qubit};",
        f"cz {phase_qubit}, {block_ancilla};",
        f"h {phase_qubit};",
        f"rz({angle_text}) {phase_qubit};",
        f"h {phase_qubit};",
        f"cz {phase_qubit}, {block_ancilla};",
        f"h {phase_qubit};",
    ]


def _qasm_snippet(
    pyqsp_phases: list[float],
    qsvt_projector_phases: list[float],
    qasm_rz_angles: list[float],
    *,
    phase_qubit: str,
    block_ancilla: str,
    system_qubit: str,
    signal_gate: str,
    signal_gate_dagger: str,
) -> str:
    lines: list[str] = []
    for idx, (phase, psi, theta) in enumerate(zip(pyqsp_phases, qsvt_projector_phases, qasm_rz_angles)):
        lines.append(f"// phi_{idx} = {phase:.12g}")
        lines.append(f"// psi_{idx} = {psi:.12g}")
        lines.append(f"// theta_rz_{idx} = -2 * psi_{idx} = {theta:.12g}")
        lines.extend(_phase_block(theta, phase_qubit, block_ancilla))
        if idx == len(pyqsp_phases) - 1:
            continue

        use_dagger = idx % 2 == 1
        gate = signal_gate_dagger if use_dagger else signal_gate
        label = "U^\\dagger" if use_dagger else "U"
        lines.append("")
        lines.append(f"// {label}")
        lines.append(f"{gate} {block_ancilla}, {system_qubit};")
        lines.append("")
    return "\n".join(lines)


def generate_component(
    *,
    component: str,
    tau: float,
    epsilon: float,
    ensure_bounded: bool,
    compute_angles: bool,
    method: str,
    signal_operator: str,
) -> dict[str, Any]:
    generator_cls = {
        "cos": PolyCosineTX,
        "sin": PolySineTX,
    }[component]

    # pyqsp prints truncation diagnostics from the polynomial generator and
    # iteration logs from the phase solver. Capture them so this tool can emit
    # structured output cleanly.
    logs = io.StringIO()
    with contextlib.redirect_stdout(logs):
        if ensure_bounded:
            cheb_coeffs, scale = generator_cls().generate(
                tau,
                epsilon,
                return_coef=True,
                ensure_bounded=True,
                return_scale=True,
                chebyshev_basis=True,
            )
        else:
            cheb_coeffs = generator_cls().generate(
                tau,
                epsilon,
                return_coef=True,
                ensure_bounded=False,
                chebyshev_basis=True,
            )
            scale = 1.0

        phases = None
        reduced_phases = None
        parity = None
        if compute_angles:
            result = QuantumSignalProcessingPhases(
                cheb_coeffs,
                method=method,
                signal_operator=signal_operator,
                chebyshev_basis=(method == "sym_qsp"),
            )
            if method == "sym_qsp":
                phases, reduced_phases, parity = result
            else:
                phases = result

    cheb = _clean_list(cheb_coeffs)
    mono = _clean_list(np.polynomial.chebyshev.cheb2poly(cheb_coeffs))

    data: dict[str, Any] = {
        "component": component,
        "tau": tau,
        "epsilon": epsilon,
        "ensure_bounded": ensure_bounded,
        "scale": _clean_float(scale),
        "degree": len(cheb) - 1,
        "chebyshev_coefficients": cheb,
        "monomial_coefficients": mono,
        "pyqsp_log": [line for line in logs.getvalue().splitlines() if line],
    }

    if phases is not None:
        phase_list = _clean_list(phases)
        qsvt_projector_phases, qasm_rz_angles = qasm_phase_data(phase_list)
        data["pyqsp_phases"] = phase_list
        data["qsvt_projector_phases"] = qsvt_projector_phases
        data["qasm_rz_angles"] = qasm_rz_angles
        data["pyqsp_response_part"] = "imag" if method == "sym_qsp" else None
        data["phase_conversion"] = "psi_0=phi_0+pi/4; psi_j=phi_j+pi/2; psi_d=phi_d+pi/4; theta_rz=-2*psi"
        data["phase_count"] = len(phase_list)
        data["qsp_degree"] = len(phase_list) - 1
    if reduced_phases is not None:
        data["reduced_phases"] = _clean_list(reduced_phases)
    if parity is not None:
        data["parity"] = int(parity)

    return data


def generate_exp_component(
    *,
    tau: float,
    epsilon: float,
    ensure_bounded: bool,
) -> dict[str, Any]:
    cos_record = generate_component(
        component="cos",
        tau=tau,
        epsilon=epsilon,
        ensure_bounded=ensure_bounded,
        compute_angles=False,
        method="sym_qsp",
        signal_operator="Wx",
    )
    sin_record = generate_component(
        component="sin",
        tau=tau,
        epsilon=epsilon,
        ensure_bounded=ensure_bounded,
        compute_angles=False,
        method="sym_qsp",
        signal_operator="Wx",
    )

    cos_mono = cos_record["monomial_coefficients"]
    sin_mono = sin_record["monomial_coefficients"]
    cos_cheb = cos_record["chebyshev_coefficients"]
    sin_cheb = sin_record["chebyshev_coefficients"]

    return {
        "component": "exp",
        "tau": tau,
        "epsilon": epsilon,
        "ensure_bounded": ensure_bounded,
        "scale": cos_record["scale"],
        "degree": max(cos_record["degree"], sin_record["degree"]),
        "target": "scale * exp(-i * tau * x)",
        "chebyshev_coefficients_complex": _complex_coefficients(cos_cheb, [-v for v in sin_cheb]),
        "monomial_coefficients_complex": _complex_coefficients(cos_mono, [-v for v in sin_mono]),
        "cos": cos_record,
        "sin": sin_record,
        "pyqsp_log": cos_record["pyqsp_log"] + sin_record["pyqsp_log"],
    }


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
            cos_mono = record["cos"]["monomial_coefficients"]
            sin_mono = record["sin"]["monomial_coefficients"]
            cos_expr = _polynomial_expr(cos_mono, variable, precision)
            sin_expr = _polynomial_expr(sin_mono, variable, precision)
            print("target = scale * exp(-i * tau * x)")
            print(f"cos monomial coefficients = {cos_mono}")
            print(f"sin monomial coefficients = {sin_mono}")
            print(f"complex monomial coefficients = {record['monomial_coefficients_complex']}")
            print(f"exp polynomial = ({cos_expr}) - I*({sin_expr})")
            if include_logs and record["pyqsp_log"]:
                print("pyqsp log:")
                for line in record["pyqsp_log"]:
                    print(f"  {line}")
            if include_qasm_snippet:
                print("QASM snippet skipped: pyqsp hamsim provides separate cos/sin phase sequences, not one exp sequence.")
            continue
        if "parity" in record:
            print(f"parity = {record['parity']}")

        cheb = record["chebyshev_coefficients"]
        mono = record["monomial_coefficients"]
        print(f"Chebyshev coefficients = {cheb}")
        print(f"monomial coefficients  = {mono}")
        print(
            "monomial polynomial   = "
            + _polynomial_expr(mono, variable, precision)
        )

        if "pyqsp_phases" in record:
            print(f"phase count = {record['phase_count']}")
            print(f"QSP degree = {record['qsp_degree']}")
            print(f"pyqsp phases = {record['pyqsp_phases']}")
            print(f"QSVT projector phases = {record['qsvt_projector_phases']}")
            print(f"pyqsp response part = {record['pyqsp_response_part']}")
            print(f"qasm rz angles = {record['qasm_rz_angles']}")

        if include_logs and record["pyqsp_log"]:
            print("pyqsp log:")
            for line in record["pyqsp_log"]:
                print(f"  {line}")

        if include_qasm_snippet:
            if "pyqsp_phases" not in record:
                print("QASM snippet skipped because --no-angles was used.")
            else:
                print("QASM snippet:")
                print(
                    _qasm_snippet(
                        record["pyqsp_phases"],
                        record["qsvt_projector_phases"],
                        record["qasm_rz_angles"],
                        phase_qubit=qasm_args.phase_qubit,
                        block_ancilla=qasm_args.block_ancilla,
                        system_qubit=qasm_args.system_qubit,
                        signal_gate=qasm_args.signal_gate,
                        signal_gate_dagger=qasm_args.signal_gate_dagger,
                    )
                )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate Hamiltonian-simulation polynomials and QSP phase angles "
            "using pyqsp."
        )
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
        choices=["cos", "sin", "both", "exp"],
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
    parser.add_argument(
        "--no-angles",
        action="store_true",
        help="Only generate polynomial coefficients.",
    )
    parser.add_argument(
        "--format",
        choices=["text", "json"],
        default="text",
        help="Output format.",
    )
    parser.add_argument(
        "--precision",
        type=int,
        default=12,
        help="Significant digits for the printed polynomial expression.",
    )
    parser.add_argument(
        "--variable",
        default="x",
        help="Variable name used in the printed monomial polynomial.",
    )
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
    parser.add_argument("--phase-qubit", default="q[0]")
    parser.add_argument("--block-ancilla", default="q[1]")
    parser.add_argument("--system-qubit", default="q[2]")
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
        parser.error("--qasm-snippet requires --component cos or --component sin")
    return args


def main() -> int:
    args = parse_args()
    if args.component == "exp":
        records = [
            generate_exp_component(
                tau=args.tau,
                epsilon=args.epsilon,
                ensure_bounded=not args.unbounded,
            )
        ]
    else:
        components = ["cos", "sin"] if args.component == "both" else [args.component]
        records = [
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


if __name__ == "__main__":
    raise SystemExit(main())
