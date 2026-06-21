from __future__ import annotations

import contextlib
import io
import json
import os
import re
import sys
from pathlib import Path
from typing import Any


def clean_float(value: Any, *, zero_tol: float = 1e-14) -> float:
    value = float(value)
    if abs(value) < zero_tol:
        return 0.0
    return value


def clean_list(values: Any) -> list[float]:
    return [clean_float(v) for v in values]


def polynomial_expr(coeffs: list[float], variable: str = "x", precision: int = 17) -> str:
    terms: list[str] = []
    for degree, coeff in enumerate(coeffs):
        if coeff == 0.0:
            continue
        abs_coeff = abs(coeff)
        if degree == 0:
            body = format_number(abs_coeff, precision)
        elif degree == 1:
            body = f"{format_number(abs_coeff, precision)}*{variable}"
        else:
            body = f"{format_number(abs_coeff, precision)}*{variable}^{degree}"

        if not terms:
            terms.append(body if coeff > 0 else f"-{body}")
        else:
            terms.append(("+ " if coeff > 0 else "- ") + body)
    return " ".join(terms) if terms else "0"


def format_number(value: float, precision: int) -> str:
    if value == 0.0:
        return "0"
    return f"{value:.{precision}g}"


def example_tag(tau: float, epsilon: float) -> str:
    return f"t{_format_float_token(tau)}_eps{_format_float_token(epsilon)}"


def pyqsp_phase_to_qsvt_projector_phase(phi: float, *, index: int, degree: int) -> float:
    import math

    if index == 0:
        shift = math.pi / 4
    elif index == degree:
        shift = math.pi / 4 if degree % 2 == 0 else -math.pi / 4
    else:
        shift = math.pi / 2
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
    return clean_list(qsvt_phases), clean_list(rz_angles)


def qasm_snippet(
    pyqsp_phases: list[float],
    qsvt_projector_phases: list[float],
    qasm_rz_angles: list[float],
    *,
    phase_qubit: str,
    block_ancillas: list[str],
    system_qubits: list[str],
    signal_gate: str,
    signal_gate_dagger: str,
) -> str:
    block_ancillas = _normalize_qubit_list(block_ancillas, label="block ancilla")
    system_qubits = _normalize_qubit_list(system_qubits, label="system qubit")
    lines: list[str] = []
    for idx, (phase, psi, theta) in enumerate(zip(pyqsp_phases, qsvt_projector_phases, qasm_rz_angles)):
        lines.append(f"// phi_{idx} = {phase:.12g}")
        lines.append(f"// psi_{idx} = {psi:.12g}")
        lines.append(f"// theta_rz_{idx} = -2 * psi_{idx} = {theta:.12g}")
        lines.extend(_phase_block(theta, phase_qubit, block_ancillas))
        if idx == len(pyqsp_phases) - 1:
            continue

        use_dagger = idx % 2 == 1
        gate = signal_gate_dagger if use_dagger else signal_gate
        label = "U^\\dagger" if use_dagger else "U"
        lines.append("")
        lines.append(f"// {label}")
        lines.append(f"{gate} {_format_gate_operands([*block_ancillas, *system_qubits])};")
        lines.append("")
    return "\n".join(lines)


def selector_qasm_snippet(
    pyqsp_phases: list[float],
    qsvt_projector_phases: list[float],
    qasm_rz_angles: list[float],
    *,
    selector_qubit: str,
    phase_qubit: str,
    block_ancillas: list[str],
    system_qubits: list[str],
    signal_gate: str,
    signal_gate_dagger: str,
) -> str:
    block_ancillas = _normalize_qubit_list(block_ancillas, label="block ancilla")
    system_qubits = _normalize_qubit_list(system_qubits, label="system qubit")
    lines: list[str] = []
    lines.append(f"h {selector_qubit};")
    lines.append("")
    for idx, (phase, psi, theta) in enumerate(zip(pyqsp_phases, qsvt_projector_phases, qasm_rz_angles)):
        lines.append(f"// phi_{idx} = {phase:.12g}")
        lines.append(f"// psi_{idx} = {psi:.12g}")
        lines.append(f"// solid-control theta_{idx} = {theta:.12g}")
        lines.append(f"// open-control SELECT_RZ theta_{idx} = {-theta:.12g}")
        lines.extend(_selector_phase_block(theta, selector_qubit, phase_qubit, block_ancillas))
        if idx == len(pyqsp_phases) - 1:
            continue

        use_dagger = idx % 2 == 1
        gate = signal_gate_dagger if use_dagger else signal_gate
        label = "U^\\dagger" if use_dagger else "U"
        lines.append("")
        lines.append(f"// {label}")
        lines.append(f"{gate} {_format_gate_operands([*block_ancillas, *system_qubits])};")
        lines.append("")

    lines.append("")
    lines.append(f"h {selector_qubit};")
    lines.append(f"sdg {selector_qubit};")
    lines.append(f"x {selector_qubit};")
    return "\n".join(lines)


def selector_qasm_file(
    record: dict[str, Any],
    *,
    precision: int,
    selector_qubit: str,
    phase_qubit: str,
    block_ancillas: list[str],
    system_qubits: list[str],
    signal_gate: str,
    signal_gate_dagger: str,
) -> str:
    block_ancillas = _normalize_qubit_list(block_ancillas, label="block ancilla")
    system_qubits = _normalize_qubit_list(system_qubits, label="system qubit")
    component = record["selector_component"]
    qreg_size = _qreg_size(selector_qubit, phase_qubit, *block_ancillas, *system_qubits)
    polynomial = polynomial_expr(record["monomial_coefficients"], "x", precision)
    snippet = selector_qasm_snippet(
        record["pyqsp_phases"],
        record["qsvt_projector_phases"],
        record["qasm_rz_angles"],
        selector_qubit=selector_qubit,
        phase_qubit=phase_qubit,
        block_ancillas=block_ancillas,
        system_qubits=system_qubits,
        signal_gate=signal_gate,
        signal_gate_dagger=signal_gate_dagger,
    )
    phase_lines = "\n".join(f"//   phi_{idx} = {phase:.17g}" for idx, phase in enumerate(record["pyqsp_phases"]))
    theta_lines = "\n".join(f"//   theta_{idx} = {theta:.17g}" for idx, theta in enumerate(record["qasm_rz_angles"]))
    block_comment = "\n".join(f"// {qubit} = block-encoding ancilla[{index}]" for index, qubit in enumerate(block_ancillas))
    system_comment = "\n".join(f"// {qubit} = system qubit[{index}]" for index, qubit in enumerate(system_qubits))
    signal_signature = _opaque_signature([f"a{index}" for index in range(len(block_ancillas))], [f"s{index}" for index in range(len(system_qubits))])
    mcx_signature = _opaque_signature([f"c{index}" for index in range(len(block_ancillas))], ["t"])

    return f"""OPENQASM 2.0;
include "qelib1.inc";

opaque {signal_gate} {signal_signature};
opaque {signal_gate_dagger} {signal_signature};
opaque mcx {mcx_signature};

qreg q[{qreg_size}];

// {phase_qubit} = phase ancilla
{block_comment}
// {selector_qubit} = selector ancilla
{system_comment}

// Selector-wrapped Hamiltonian simulation {component} component.
// Generated by tools/hamsim_qsp.py with:
//   tau = {record["tau"]}
//   epsilon = {record["epsilon"]}
//   method = sym_qsp
//   signal_operator = Wx
//   ensure_bounded = {record["ensure_bounded"]}
//   scale = {record["scale"]}
//
// The base QSP component has Im B[00] = P_{component}(H). The selector wrapper
// coherently runs theta and -theta branches, takes their difference, multiplies
// by -i, and moves the result to the all-zero selector branch. Therefore the
// final all-zero branch is directly P_{component}(H), so the verifier compares
// the full polynomial.
//
// Open-control multi-controlled X gates on the block ancillas are expanded
// with x gates around abstract mcx. This matches the paper convention that
// phases act on the |0...0> block-ancilla branch.
//
// P_{component}(x) = {polynomial}
//
// QSP phases:
{phase_lines}
//
// Physical OpenQASM rz angles use theta_rz = -2 * psi:
{theta_lines}

{snippet}
"""


def verification_command(
    record: dict[str, Any],
    qasm_path: Path,
    *,
    precision: int,
    selector_qubit: str,
    phase_qubit: str,
    block_ancillas: list[str],
    system_qubits: list[str],
) -> list[str]:
    block_ancillas = _normalize_qubit_list(block_ancillas, label="block ancilla")
    system_qubits = _normalize_qubit_list(system_qubits, label="system qubit")
    polynomial = polynomial_expr(record["monomial_coefficients"], "x", precision)
    return [
        ".venv/bin/python",
        "-m",
        "symbolic.verify",
        str(qasm_path),
        "--ancillas",
        phase_qubit,
        *block_ancillas,
        selector_qubit,
        "--systems",
        *system_qubits,
        "--expected-polynomial",
        polynomial,
        "--hermitian-base",
        "--compare-polynomial-only",
    ]


def write_selector_examples(
    *,
    tau: float,
    epsilon: float,
    ensure_bounded: bool,
    method: str,
    signal_operator: str,
    precision: int,
    examples_dir: Path,
    selector_qubit: str,
    phase_qubit: str,
    block_ancillas: list[str],
    system_qubits: list[str],
    signal_gate: str,
    signal_gate_dagger: str,
) -> dict[str, Any]:
    if method != "sym_qsp" or signal_operator != "Wx":
        raise ValueError("--write-examples currently expects --method sym_qsp --signal-operator Wx")

    block_ancillas = _normalize_qubit_list(block_ancillas, label="block ancilla")
    system_qubits = _normalize_qubit_list(system_qubits, label="system qubit")
    tag = example_tag(tau, epsilon)
    suffix = "" if len(block_ancillas) == 1 else f"_m{len(block_ancillas)}"
    output_dir = examples_dir / f"qsp_hamsim_{tag}{suffix}"
    output_dir.mkdir(parents=True, exist_ok=True)

    records = [
        generate_selector_component(
            component=component,
            tau=tau,
            epsilon=epsilon,
            ensure_bounded=ensure_bounded,
            compute_angles=True,
            method=method,
            signal_operator=signal_operator,
        )
        for component in ("cos", "sin")
    ]

    files: list[dict[str, Any]] = []
    for record in records:
        component = record["selector_component"]
        filename = f"qsp_hamsim_{component}_selector_{tag}_deg{record['degree']}.qasm"
        qasm_path = output_dir / filename
        qasm_text = selector_qasm_file(
            record,
            precision=precision,
            selector_qubit=selector_qubit,
            phase_qubit=phase_qubit,
            block_ancillas=block_ancillas,
            system_qubits=system_qubits,
            signal_gate=signal_gate,
            signal_gate_dagger=signal_gate_dagger,
        )
        qasm_path.write_text(qasm_text, encoding="utf-8")
        polynomial = polynomial_expr(record["monomial_coefficients"], "x", precision)
        files.append(
            {
                "component": component,
                "qasm": str(qasm_path),
                "degree": record["degree"],
                "polynomial": polynomial,
                "monomial_coefficients": record["monomial_coefficients"],
                "verification_command": verification_command(
                    record,
                    qasm_path,
                    precision=precision,
                    selector_qubit=selector_qubit,
                    phase_qubit=phase_qubit,
                    block_ancillas=block_ancillas,
                    system_qubits=system_qubits,
                ),
            }
        )

    metadata = {
        "tau": tau,
        "epsilon": epsilon,
        "ensure_bounded": ensure_bounded,
        "scale": records[0]["scale"],
        "selector_qubit": selector_qubit,
        "phase_qubit": phase_qubit,
        "block_ancillas": block_ancillas,
        "system_qubits": system_qubits,
        "comparison": "full polynomial",
        "files": files,
    }
    metadata_path = output_dir / "expected_polynomials.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    metadata["metadata"] = str(metadata_path)
    return metadata


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
    np, quantum_signal_processing_phases, poly_cosine_tx, poly_sine_tx = _load_pyqsp()
    generator_cls = {
        "cos": poly_cosine_tx,
        "sin": poly_sine_tx,
    }[component]

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
            result = quantum_signal_processing_phases(
                cheb_coeffs,
                method=method,
                signal_operator=signal_operator,
                chebyshev_basis=(method == "sym_qsp"),
            )
            if method == "sym_qsp":
                phases, reduced_phases, parity = result
            else:
                phases = result

    cheb = clean_list(cheb_coeffs)
    mono = clean_list(np.polynomial.chebyshev.cheb2poly(cheb_coeffs))

    data: dict[str, Any] = {
        "component": component,
        "tau": tau,
        "epsilon": epsilon,
        "ensure_bounded": ensure_bounded,
        "scale": clean_float(scale),
        "degree": len(cheb) - 1,
        "chebyshev_coefficients": cheb,
        "monomial_coefficients": mono,
        "pyqsp_log": [line for line in logs.getvalue().splitlines() if line],
    }

    if phases is not None:
        phase_list = clean_list(phases)
        qsvt_projector_phases, qasm_rz_angles = qasm_phase_data(phase_list)
        data["pyqsp_phases"] = phase_list
        data["qsvt_projector_phases"] = qsvt_projector_phases
        data["qasm_rz_angles"] = qasm_rz_angles
        data["pyqsp_response_part"] = "imag" if method == "sym_qsp" else None
        data["phase_conversion"] = "psi_0=phi_0+pi/4; psi_j=phi_j+pi/2; psi_d=phi_d+pi/4; theta_rz=-2*psi"
        data["phase_count"] = len(phase_list)
        data["qsp_degree"] = len(phase_list) - 1
    if reduced_phases is not None:
        data["reduced_phases"] = clean_list(reduced_phases)
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


def generate_selector_component(
    *,
    component: str,
    tau: float,
    epsilon: float,
    ensure_bounded: bool,
    compute_angles: bool,
    method: str,
    signal_operator: str,
) -> dict[str, Any]:
    base = generate_component(
        component=component,
        tau=tau,
        epsilon=epsilon,
        ensure_bounded=ensure_bounded,
        compute_angles=compute_angles,
        method=method,
        signal_operator=signal_operator,
    )
    data = dict(base)
    data["component"] = f"{component}-selector"
    data["selector_component"] = component
    data["selector_output"] = "all-zero branch equals the former imaginary polynomial part"
    phase_conversion = data.get("phase_conversion")
    if phase_conversion:
        data["phase_conversion"] = phase_conversion + "; selector wrapper extracts imaginary part"
    else:
        data["phase_conversion"] = "selector wrapper extracts imaginary part"
    return data


def _configure_import_path() -> None:
    tmp_mpl = Path(os.environ.get("TMPDIR", "/tmp")) / "symbolic_be_verifier_mpl"
    tmp_mpl.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(tmp_mpl))

    repo_parent = Path(__file__).resolve().parents[2]
    local_pyqsp = repo_parent / "pyqsp"
    if local_pyqsp.exists() and str(local_pyqsp) not in sys.path:
        sys.path.insert(0, str(local_pyqsp))


def _load_pyqsp() -> tuple[Any, Any, Any, Any]:
    _configure_import_path()
    try:
        import numpy as np
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "Missing dependency 'numpy'. Run this tool with an environment that has "
            "pyqsp's dependencies installed."
        ) from exc

    try:
        with contextlib.redirect_stderr(io.StringIO()):
            from pyqsp.angle_sequence import QuantumSignalProcessingPhases
            from pyqsp.poly import PolyCosineTX, PolySineTX
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "Could not import pyqsp or one of its dependencies. Install pyqsp's "
            "requirements, or keep the local pyqsp checkout next to symbolic-be-verifier."
        ) from exc

    return np, QuantumSignalProcessingPhases, PolyCosineTX, PolySineTX


def _phase_block(angle: float, phase_qubit: str, block_ancillas: list[str]) -> list[str]:
    angle_text = f"{angle:.12g}"
    return [
        *_open_mcx_block(block_ancillas, phase_qubit),
        f"rz({angle_text}) {phase_qubit};",
        *_open_mcx_block(block_ancillas, phase_qubit),
    ]


def _select_rz(angle: float, selector_qubit: str, phase_qubit: str) -> list[str]:
    half = angle / 2.0
    half_text = f"{half:.12g}"
    neg_half_text = f"{-half:.12g}"
    return [
        f"x {selector_qubit};",
        f"rz({half_text}) {phase_qubit};",
        f"cx {selector_qubit}, {phase_qubit};",
        f"rz({neg_half_text}) {phase_qubit};",
        f"cx {selector_qubit}, {phase_qubit};",
        f"x {selector_qubit};",
        f"rz({neg_half_text}) {phase_qubit};",
        f"cx {selector_qubit}, {phase_qubit};",
        f"rz({half_text}) {phase_qubit};",
        f"cx {selector_qubit}, {phase_qubit};",
    ]


def _open_mcx_block(block_ancillas: list[str], phase_qubit: str) -> list[str]:
    controls = ", ".join([*block_ancillas, phase_qubit])
    return [
        *(f"x {block_ancilla};" for block_ancilla in block_ancillas),
        f"mcx {controls};",
        *(f"x {block_ancilla};" for block_ancilla in reversed(block_ancillas)),
    ]


def _selector_phase_block(angle: float, selector_qubit: str, phase_qubit: str, block_ancillas: list[str]) -> list[str]:
    open_control_angle = -angle
    return [
        *_open_mcx_block(block_ancillas, phase_qubit),
        *_select_rz(open_control_angle, selector_qubit, phase_qubit),
        *_open_mcx_block(block_ancillas, phase_qubit),
    ]


def _normalize_qubit_list(values: list[str], *, label: str) -> list[str]:
    normalized = [value.strip() for value in values]
    if not normalized:
        raise ValueError(f"At least one {label} is required")
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"Duplicate {label}s are not allowed: {normalized}")
    return normalized


def _format_gate_operands(qubits: list[str]) -> str:
    return ", ".join(qubits)


def _opaque_signature(ancilla_args: list[str], system_args: list[str]) -> str:
    return ", ".join([*ancilla_args, *system_args])


def _qubit_index(qubit: str) -> int:
    match = re.fullmatch(r"q\[(\d+)\]", qubit.strip())
    if not match:
        raise ValueError(f"Expected a q[index] qubit name, got {qubit!r}")
    return int(match.group(1))


def _qreg_size(*qubits: str) -> int:
    return max(_qubit_index(qubit) for qubit in qubits) + 1


def _format_float_token(value: float) -> str:
    if value != 0.0 and abs(value) < 1e-3:
        text = f"{value:.0e}"
        text = text.replace("e-0", "e-").replace("e+0", "e").replace("e+", "e")
        return text
    text = f"{value:g}".replace("-", "m").replace(".", "p")
    if text.startswith("0p"):
        return "0" + text[2:]
    return text


def _pad(values: list[float], length: int) -> list[float]:
    return values + [0.0] * (length - len(values))


def _complex_coefficients(real_coeffs: list[float], imag_coeffs: list[float]) -> list[dict[str, float]]:
    length = max(len(real_coeffs), len(imag_coeffs))
    real_padded = _pad(real_coeffs, length)
    imag_padded = _pad(imag_coeffs, length)
    return [
        {"real": clean_float(real), "imag": clean_float(imag)}
        for real, imag in zip(real_padded, imag_padded)
    ]
