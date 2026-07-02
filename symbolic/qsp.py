from __future__ import annotations

import contextlib
from dataclasses import dataclass
import io
import json
import os
import re
import sys
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class QasmGate:
    name: str
    qubits: tuple[str, ...]
    parameter: float | None = None


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


def hamsim_exp_polynomial_expr(
    cos_coeffs: list[float],
    sin_coeffs: list[float],
    *,
    precision: int = 17,
    scale: float = 1.0,
    variable: str = "x",
) -> str:
    terms: list[str] = []
    length = max(len(cos_coeffs), len(sin_coeffs))
    for degree in range(length):
        real = clean_float(scale * (cos_coeffs[degree] if degree < len(cos_coeffs) else 0.0))
        imag = clean_float(-scale * (sin_coeffs[degree] if degree < len(sin_coeffs) else 0.0))
        coeff = _format_complex_coeff(real, imag, precision)
        if coeff is None:
            continue

        if degree == 0:
            body = coeff
        elif degree == 1:
            body = f"{_coefficient_factor(coeff)}*{variable}"
        else:
            body = f"{_coefficient_factor(coeff)}*{variable}^{degree}"

        if not terms:
            terms.append(body)
        elif body.startswith("-"):
            terms.append("- " + body[1:])
        else:
            terms.append("+ " + body)
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


def build_raw_qsp_gate_list(
    pyqsp_phases: list[float],
    qsvt_projector_phases: list[float],
    qasm_rz_angles: list[float],
    *,
    phase_qubit: str,
    block_ancillas: list[str],
    system_qubits: list[str],
    signal_gate: str,
    signal_gate_dagger: str,
) -> list[QasmGate]:
    block_ancillas = _normalize_qubit_list(block_ancillas, label="block ancilla")
    system_qubits = _normalize_qubit_list(system_qubits, label="system qubit")
    gates: list[QasmGate] = []
    for idx, theta in enumerate(qasm_rz_angles):
        gates.extend(_phase_block_gates(-theta, phase_qubit, block_ancillas))
        if idx == len(pyqsp_phases) - 1:
            continue

        use_dagger = idx % 2 == 1
        gate = signal_gate_dagger if use_dagger else signal_gate
        gates.append(QasmGate(gate, tuple([*block_ancillas, *system_qubits])))
    return gates


def adjoint_gate(
    gate: QasmGate,
    *,
    signal_gate: str = "UH",
    signal_gate_dagger: str = "UHdg",
) -> QasmGate:
    name = gate.name.lower()
    if name in {"rz", "rx", "ry"}:
        if gate.parameter is None:
            raise ValueError(f"{gate.name} requires an angle")
        return QasmGate(gate.name, gate.qubits, -gate.parameter)

    if gate.name == signal_gate:
        return QasmGate(signal_gate_dagger, gate.qubits)
    if gate.name == signal_gate_dagger:
        return QasmGate(signal_gate, gate.qubits)

    if name == "s":
        return QasmGate("sdg", gate.qubits)
    if name == "sdg":
        return QasmGate("s", gate.qubits)

    if name in {"h", "x", "z", "cx", "cz", "mcx"}:
        return gate

    raise NotImplementedError(gate.name)


def daggerize_gate_list(
    gates: list[QasmGate],
    *,
    signal_gate: str = "UH",
    signal_gate_dagger: str = "UHdg",
) -> list[QasmGate]:
    return [
        adjoint_gate(gate, signal_gate=signal_gate, signal_gate_dagger=signal_gate_dagger)
        for gate in reversed(gates)
    ]


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
    controlled_signal_gate: str = "cUH",
    controlled_signal_gate_dagger: str = "cUHdg",
) -> str:
    block_ancillas = _normalize_qubit_list(block_ancillas, label="block ancilla")
    system_qubits = _normalize_qubit_list(system_qubits, label="system qubit")
    raw_gates = build_raw_qsp_gate_list(
        pyqsp_phases,
        qsvt_projector_phases,
        qasm_rz_angles,
        phase_qubit=phase_qubit,
        block_ancillas=block_ancillas,
        system_qubits=system_qubits,
        signal_gate=signal_gate,
        signal_gate_dagger=signal_gate_dagger,
    )
    dagger_gates = daggerize_gate_list(raw_gates, signal_gate=signal_gate, signal_gate_dagger=signal_gate_dagger)
    lines: list[str] = []
    lines.append(f"h {selector_qubit};")
    lines.append("")

    lines.append("// C branch: selector = 0")
    lines.extend(
        _controlled_gate_list_lines(
            raw_gates,
            controls=[(selector_qubit, 0)],
            signal_gate=signal_gate,
            signal_gate_dagger=signal_gate_dagger,
            controlled_signal_gate=controlled_signal_gate,
            controlled_signal_gate_dagger=controlled_signal_gate_dagger,
        )
    )
    lines.append("")
    lines.append("// Cdag branch: selector = 1")
    lines.extend(
        _controlled_gate_list_lines(
            dagger_gates,
            controls=[(selector_qubit, 1)],
            signal_gate=signal_gate,
            signal_gate_dagger=signal_gate_dagger,
            controlled_signal_gate=controlled_signal_gate,
            controlled_signal_gate_dagger=controlled_signal_gate_dagger,
        )
    )
    lines.append("")

    lines.append(f"z {selector_qubit};")
    lines.append(f"h {selector_qubit};")
    lines.append(f"rz(3.14159265359) {selector_qubit};")
    return "\n".join(lines)


def full_hamsim_qasm_snippet(
    cos_record: dict[str, Any],
    sin_record: dict[str, Any],
    *,
    selector_qubit: str,
    component_selector_qubit: str,
    phase_qubit: str,
    block_ancillas: list[str],
    system_qubits: list[str],
    signal_gate: str,
    signal_gate_dagger: str,
    controlled_signal_gate: str,
    controlled_signal_gate_dagger: str,
) -> str:
    block_ancillas = _normalize_qubit_list(block_ancillas, label="block ancilla")
    system_qubits = _normalize_qubit_list(system_qubits, label="system qubit")
    if len(system_qubits) != 1:
        raise ValueError("full Hamiltonian simulation currently expects one system qubit")
    if len(sin_record["qsvt_projector_phases"]) != len(cos_record["qsvt_projector_phases"]) + 1:
        raise ValueError("full Hamiltonian simulation expects one more sine phase than cosine phase")

    lines: list[str] = []
    lines.append(f"h {selector_qubit};")
    lines.append(f"h {component_selector_qubit};")
    lines.append("")

    for idx, cos_psi in enumerate(cos_record["qsvt_projector_phases"]):
        sin_psi = sin_record["qsvt_projector_phases"][idx]
        delta = clean_float(sin_psi - cos_psi)
        lines.append(f"// common phase {idx}")
        lines.append(f"// cos phi_{idx} = {cos_record['pyqsp_phases'][idx]:.12g}")
        lines.append(f"// cos psi_{idx} = {cos_psi:.12g}")
        lines.append(f"// sin phi_{idx} = {sin_record['pyqsp_phases'][idx]:.12g}")
        lines.append(f"// sin psi_{idx} = {sin_psi:.12g}")
        lines.append(f"// Delta_{idx} = sin psi_{idx} - cos psi_{idx} = {delta:.12g}")
        lines.extend(
            _base_signed_phase(
                cos_psi,
                component_selector_qubit=component_selector_qubit,
                phase_qubit=phase_qubit,
                block_ancillas=block_ancillas,
            )
        )
        lines.append("")
        lines.append("// extra signed difference phase on sin branch only")
        lines.extend(
            _controlled_signed_phase(
                delta,
                selector_qubit=selector_qubit,
                component_selector_qubit=component_selector_qubit,
                phase_qubit=phase_qubit,
                block_ancillas=block_ancillas,
            )
        )
        lines.append("")

        if idx == len(cos_record["qsvt_projector_phases"]) - 1:
            continue

        use_dagger = idx % 2 == 1
        gate = signal_gate_dagger if use_dagger else signal_gate
        label = "common U^\\dagger" if use_dagger else "common U"
        lines.append(f"// {label}")
        lines.append(f"{gate} {_format_gate_operands([*block_ancillas, *system_qubits])};")
        lines.append("")

    final_idx = len(cos_record["qsvt_projector_phases"])
    final_sin_psi = sin_record["qsvt_projector_phases"][final_idx]
    lines.append("// sin-only final U")
    lines.append(f"{controlled_signal_gate} {_format_gate_operands([selector_qubit, *block_ancillas, *system_qubits])};")
    lines.append("")
    lines.append(f"// final sin-only phase {final_idx}")
    lines.append(f"// sin phi_{final_idx} = {sin_record['pyqsp_phases'][final_idx]:.12g}")
    lines.append(f"// sin psi_{final_idx} = {final_sin_psi:.12g}")
    lines.extend(
        _controlled_signed_phase(
            final_sin_psi,
            selector_qubit=selector_qubit,
            component_selector_qubit=component_selector_qubit,
            phase_qubit=phase_qubit,
            block_ancillas=block_ancillas,
        )
    )
    lines.append("")

    lines.append("// extract pyqsp imaginary response on the component selector")
    lines.append(f"h {component_selector_qubit};")
    lines.append(f"{_component_extraction_phase_gate(int(sin_record.get('degree', len(sin_record['pyqsp_phases']) - 1)))} {component_selector_qubit};")
    lines.append(f"x {component_selector_qubit};")
    lines.append("")
    lines.append("// combine 1/2 * (E_cos - i E_sin)")
    lines.append(f"sdg {selector_qubit};")
    lines.append(f"h {selector_qubit};")
    return "\n".join(lines)


def full_hamsim_qasm_file(
    record: dict[str, Any],
    *,
    precision: int,
    selector_qubit: str,
    component_selector_qubit: str,
    phase_qubit: str,
    block_ancillas: list[str],
    system_qubits: list[str],
    signal_gate: str,
    signal_gate_dagger: str,
    controlled_signal_gate: str,
    controlled_signal_gate_dagger: str,
) -> str:
    block_ancillas = _normalize_qubit_list(block_ancillas, label="block ancilla")
    system_qubits = _normalize_qubit_list(system_qubits, label="system qubit")
    qreg_size = _qreg_size(selector_qubit, component_selector_qubit, phase_qubit, *block_ancillas, *system_qubits)
    snippet = full_hamsim_qasm_snippet(
        record["cos"],
        record["sin"],
        selector_qubit=selector_qubit,
        component_selector_qubit=component_selector_qubit,
        phase_qubit=phase_qubit,
        block_ancillas=block_ancillas,
        system_qubits=system_qubits,
        signal_gate=signal_gate,
        signal_gate_dagger=signal_gate_dagger,
        controlled_signal_gate=controlled_signal_gate,
        controlled_signal_gate_dagger=controlled_signal_gate_dagger,
    )
    cos_polynomial = polynomial_expr(record["cos"]["monomial_coefficients"], "x", precision)
    sin_polynomial = polynomial_expr(record["sin"]["monomial_coefficients"], "x", precision)
    exp_polynomial = hamsim_exp_polynomial_expr(
        record["cos"]["monomial_coefficients"],
        record["sin"]["monomial_coefficients"],
        precision=precision,
        scale=0.5,
    )
    cos_phase_lines = "\n".join(f"//   cos phi_{idx} = {phase:.17g}" for idx, phase in enumerate(record["cos"]["pyqsp_phases"]))
    sin_phase_lines = "\n".join(f"//   sin phi_{idx} = {phase:.17g}" for idx, phase in enumerate(record["sin"]["pyqsp_phases"]))
    block_comment = "\n".join(f"// {qubit} = block-encoding ancilla[{index}]" for index, qubit in enumerate(block_ancillas))
    system_comment = "\n".join(f"// {qubit} = system qubit[{index}]" for index, qubit in enumerate(system_qubits))
    controlled_signature = _opaque_signature(["c", *[f"a{index}" for index in range(len(block_ancillas))]], [f"s{index}" for index in range(len(system_qubits))])
    mcx_signature = _opaque_signature(["selector", *[f"c{index}" for index in range(len(block_ancillas))]], ["t"])

    return f"""OPENQASM 2.0;
include "qelib1.inc";

opaque {signal_gate} {controlled_signature.removeprefix("c, ")};
opaque {signal_gate_dagger} {controlled_signature.removeprefix("c, ")};
opaque {controlled_signal_gate} {controlled_signature};
opaque {controlled_signal_gate_dagger} {controlled_signature};
opaque mcx {mcx_signature};

qreg q[{qreg_size}];

// {selector_qubit} = Hamiltonian-simulation selector ancilla
// {component_selector_qubit} = QSP imaginary-part extraction selector ancilla
// {phase_qubit} = QSP phase ancilla
{block_comment}
{system_comment}

// Full Hamiltonian simulation QSP block generated by tools/hamsim_qsp.py.
// The outer selector multiplexes the QSP phases: each common layer applies the
// cosine signed phase, then selector |1> adds the sine-cosine phase difference.
// The common U_H/U_Hdg skeleton is shared; only the sine branch has the final
// extra U_H and phase. The component selector extracts pyqsp's imaginary
// response before the outer sdg+h combines P_cos(H) - i P_sin(H).
//
// P_cos(x) = {cos_polynomial}
// P_sin(x) = {sin_polynomial}
// Expected all-zero selector branch = {exp_polynomial}
//
// QSP phases:
{cos_phase_lines}
{sin_phase_lines}

{snippet}
"""


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
    controlled_signal_gate: str = "cUH",
    controlled_signal_gate_dagger: str = "cUHdg",
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
        controlled_signal_gate=controlled_signal_gate,
        controlled_signal_gate_dagger=controlled_signal_gate_dagger,
    )
    phase_lines = "\n".join(f"//   phi_{idx} = {phase:.17g}" for idx, phase in enumerate(record["pyqsp_phases"]))
    theta_lines = "\n".join(f"//   theta_{idx} = {theta:.17g}" for idx, theta in enumerate(record["qasm_rz_angles"]))
    block_comment = "\n".join(f"// {qubit} = block-encoding ancilla[{index}]" for index, qubit in enumerate(block_ancillas))
    system_comment = "\n".join(f"// {qubit} = system qubit[{index}]" for index, qubit in enumerate(system_qubits))
    signal_signature = _opaque_signature([f"a{index}" for index in range(len(block_ancillas))], [f"s{index}" for index in range(len(system_qubits))])
    controlled_signature = _opaque_signature(["c", *[f"a{index}" for index in range(len(block_ancillas))]], [f"s{index}" for index in range(len(system_qubits))])
    mcx_signature = _opaque_signature([f"c{index}" for index in range(len(block_ancillas))], ["t"])

    return f"""OPENQASM 2.0;
include "qelib1.inc";

opaque {signal_gate} {signal_signature};
opaque {signal_gate_dagger} {signal_signature};
opaque {controlled_signal_gate} {controlled_signature};
opaque {controlled_signal_gate_dagger} {controlled_signature};
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
// runs C on selector |0> and the true adjoint C^dagger on selector |1>, then
// extracts (C - C^dagger)/(2i). Therefore the final all-zero branch is directly
// P_{component}(H), so the verifier compares the full polynomial.
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


def full_hamsim_verification_command(
    record: dict[str, Any],
    qasm_path: Path,
    *,
    precision: int,
    selector_qubit: str,
    component_selector_qubit: str,
    phase_qubit: str,
    block_ancillas: list[str],
    system_qubits: list[str],
) -> list[str]:
    block_ancillas = _normalize_qubit_list(block_ancillas, label="block ancilla")
    system_qubits = _normalize_qubit_list(system_qubits, label="system qubit")
    target_scale = clean_float(0.5 * record["scale"])
    return [
        ".venv/bin/python",
        "-m",
        "symbolic.verify",
        str(qasm_path),
        "--ancillas",
        selector_qubit,
        component_selector_qubit,
        phase_qubit,
        *block_ancillas,
        "--systems",
        *system_qubits,
        "--hamsim-tau",
        format_number(record["tau"], precision),
        "--hamsim-epsilon",
        format_number(record["epsilon"], precision),
        "--hamsim-scale",
        format_number(target_scale, precision),
        "--hermitian-base",
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
    controlled_signal_gate: str = "cUH",
    controlled_signal_gate_dagger: str = "cUHdg",
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
        filename = f"qsp_hamsim_{component}_selector_{tag}{suffix}_deg{record['degree']}.qasm"
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
            controlled_signal_gate=controlled_signal_gate,
            controlled_signal_gate_dagger=controlled_signal_gate_dagger,
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
        "assumption": "Hermitian base: selector wrappers use C^dagger from reverse-order adjoint gates, and Hd is rewritten to H for polynomial evaluation",
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


def write_full_hamsim_example(
    *,
    tau: float,
    epsilon: float,
    ensure_bounded: bool,
    method: str,
    signal_operator: str,
    precision: int,
    examples_dir: Path,
    selector_qubit: str,
    component_selector_qubit: str,
    phase_qubit: str,
    block_ancillas: list[str],
    system_qubits: list[str],
    signal_gate: str,
    signal_gate_dagger: str,
    controlled_signal_gate: str,
    controlled_signal_gate_dagger: str,
) -> dict[str, Any]:
    if method != "sym_qsp" or signal_operator != "Wx":
        raise ValueError("--write-examples --component full currently expects --method sym_qsp --signal-operator Wx")

    block_ancillas = _normalize_qubit_list(block_ancillas, label="block ancilla")
    system_qubits = _normalize_qubit_list(system_qubits, label="system qubit")
    record = generate_full_hamsim_record(
        tau=tau,
        epsilon=epsilon,
        ensure_bounded=ensure_bounded,
        method=method,
        signal_operator=signal_operator,
    )
    tag = example_tag(tau, epsilon)
    suffix = "" if len(block_ancillas) == 1 else f"_m{len(block_ancillas)}"
    output_dir = examples_dir / f"qsp_hamsim_full_{tag}{suffix}"
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"qsp_hamsim_full_{tag}{suffix}_deg{record['degree']}.qasm"
    qasm_path = output_dir / filename
    qasm_text = full_hamsim_qasm_file(
        record,
        precision=precision,
        selector_qubit=selector_qubit,
        component_selector_qubit=component_selector_qubit,
        phase_qubit=phase_qubit,
        block_ancillas=block_ancillas,
        system_qubits=system_qubits,
        signal_gate=signal_gate,
        signal_gate_dagger=signal_gate_dagger,
        controlled_signal_gate=controlled_signal_gate,
        controlled_signal_gate_dagger=controlled_signal_gate_dagger,
    )
    qasm_path.write_text(qasm_text, encoding="utf-8")
    polynomial = hamsim_exp_polynomial_expr(
        record["cos"]["monomial_coefficients"],
        record["sin"]["monomial_coefficients"],
        precision=precision,
        scale=0.5,
    )

    metadata = {
        "tau": tau,
        "epsilon": epsilon,
        "ensure_bounded": ensure_bounded,
        "assumption": "Hermitian base: full block shares the cosine/sine QSP signal skeleton and multiplexes signed phase gadgets; UHdg is interpreted as an adjoint block and Hd is rewritten to H for polynomial evaluation",
        "scale": record["scale"],
        "selector_qubit": selector_qubit,
        "component_selector_qubit": component_selector_qubit,
        "phase_qubit": phase_qubit,
        "block_ancillas": block_ancillas,
        "system_qubits": system_qubits,
        "component": "full",
        "qasm": str(qasm_path),
        "degree": record["degree"],
        "polynomial": polynomial,
        "verification_command": full_hamsim_verification_command(
            record,
            qasm_path,
            precision=precision,
            selector_qubit=selector_qubit,
            component_selector_qubit=component_selector_qubit,
            phase_qubit=phase_qubit,
            block_ancillas=block_ancillas,
            system_qubits=system_qubits,
        ),
    }
    metadata_path = output_dir / "expected_polynomial.json"
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


def generate_full_hamsim_record(
    *,
    tau: float,
    epsilon: float,
    ensure_bounded: bool,
    method: str,
    signal_operator: str,
) -> dict[str, Any]:
    cos_record = generate_component(
        component="cos",
        tau=tau,
        epsilon=epsilon,
        ensure_bounded=ensure_bounded,
        compute_angles=True,
        method=method,
        signal_operator=signal_operator,
    )
    sin_record = generate_component(
        component="sin",
        tau=tau,
        epsilon=epsilon,
        ensure_bounded=ensure_bounded,
        compute_angles=True,
        method=method,
        signal_operator=signal_operator,
    )
    return {
        "component": "full",
        "tau": tau,
        "epsilon": epsilon,
        "ensure_bounded": ensure_bounded,
        "scale": cos_record["scale"],
        "degree": max(cos_record["degree"], sin_record["degree"]),
        "selector_output": "all-zero selector branch is 1/2 * (P_cos(H) - i P_sin(H))",
        "cos": cos_record,
        "sin": sin_record,
        "pyqsp_log": cos_record["pyqsp_log"] + sin_record["pyqsp_log"],
    }


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


def _phase_block_gates(angle: float, phase_qubit: str, block_ancillas: list[str]) -> list[QasmGate]:
    return [
        *_open_mcx_block_gates(block_ancillas, phase_qubit),
        QasmGate("rz", (phase_qubit,), angle),
        *_open_mcx_block_gates(block_ancillas, phase_qubit),
    ]


def _controlled_gate_list_lines(
    gates: list[QasmGate],
    *,
    controls: list[tuple[str, int]],
    signal_gate: str,
    signal_gate_dagger: str,
    controlled_signal_gate: str,
    controlled_signal_gate_dagger: str,
) -> list[str]:
    lines: list[str] = []
    for gate in gates:
        lines.extend(
            _controlled_gate_lines(
                gate,
                controls=controls,
                signal_gate=signal_gate,
                signal_gate_dagger=signal_gate_dagger,
                controlled_signal_gate=controlled_signal_gate,
                controlled_signal_gate_dagger=controlled_signal_gate_dagger,
            )
        )
    return lines


def _controlled_gate_lines(
    gate: QasmGate,
    *,
    controls: list[tuple[str, int]],
    signal_gate: str,
    signal_gate_dagger: str,
    controlled_signal_gate: str,
    controlled_signal_gate_dagger: str,
) -> list[str]:
    if not controls:
        return [_format_qasm_gate(gate)]

    name = gate.name.lower()
    open_controls = [qubit for qubit, value in controls if value == 0]
    solid_controls = [qubit for qubit, value in controls if value == 1]

    if name in {"rz", "rx", "ry"}:
        if name != "rz":
            raise NotImplementedError(f"Controlled {gate.name} emission is not supported")
        if gate.parameter is None:
            raise ValueError("rz requires an angle")
        return _multi_controlled_rz(gate.parameter, open_controls=open_controls, solid_controls=solid_controls, target=gate.qubits[0])

    if name == "x":
        target = gate.qubits[0]
        return _controlled_x_lines(open_controls=open_controls, solid_controls=solid_controls, target=target)

    if name == "cx":
        target = gate.qubits[1]
        return _controlled_x_lines(
            open_controls=open_controls,
            solid_controls=[*solid_controls, gate.qubits[0]],
            target=target,
        )

    if name == "mcx":
        target = gate.qubits[-1]
        return _controlled_x_lines(
            open_controls=open_controls,
            solid_controls=[*solid_controls, *gate.qubits[:-1]],
            target=target,
        )

    if gate.name in {signal_gate, signal_gate_dagger}:
        dagger = gate.name == signal_gate_dagger
        controlled_name = _controlled_signal_gate_name(
            len(controls),
            dagger=dagger,
            signal_gate=signal_gate,
            signal_gate_dagger=signal_gate_dagger,
            controlled_signal_gate=controlled_signal_gate,
            controlled_signal_gate_dagger=controlled_signal_gate_dagger,
        )
        line = f"{controlled_name} {_format_gate_operands([*solid_controls, *open_controls, *gate.qubits])};"
        return [
            *(f"x {control};" for control in open_controls),
            line,
            *(f"x {control};" for control in reversed(open_controls)),
        ]

    raise NotImplementedError(f"Controlled {gate.name} emission is not supported")


def _controlled_x_lines(*, open_controls: list[str], solid_controls: list[str], target: str) -> list[str]:
    controls = [*solid_controls, *open_controls]
    if not controls:
        return [f"x {target};"]
    gate = "cx" if len(controls) == 1 else "mcx"
    return [
        *(f"x {control};" for control in open_controls),
        f"{gate} {_format_gate_operands([*controls, target])};",
        *(f"x {control};" for control in reversed(open_controls)),
    ]


def _controlled_signal_gate_name(
    control_count: int,
    *,
    dagger: bool,
    signal_gate: str,
    signal_gate_dagger: str,
    controlled_signal_gate: str,
    controlled_signal_gate_dagger: str,
) -> str:
    if control_count == 1:
        return controlled_signal_gate_dagger if dagger else controlled_signal_gate
    base = signal_gate_dagger if dagger else signal_gate
    return f"{'c' * control_count}{base}"


def _format_qasm_gate(gate: QasmGate) -> str:
    if gate.parameter is None:
        return f"{gate.name} {_format_gate_operands(list(gate.qubits))};"
    return f"{gate.name}({gate.parameter:.12g}) {_format_gate_operands(list(gate.qubits))};"


def _base_signed_phase(
    psi: float,
    *,
    component_selector_qubit: str,
    phase_qubit: str,
    block_ancillas: list[str],
) -> list[str]:
    signed_angle = -2.0 * psi
    return [
        "// base cos signed phase",
        *_controlled_x_lines(
            open_controls=[component_selector_qubit],
            solid_controls=[],
            target=phase_qubit,
        ),
        "",
        *_phase_block(signed_angle, phase_qubit, block_ancillas),
        "",
        *_controlled_x_lines(
            open_controls=[component_selector_qubit],
            solid_controls=[],
            target=phase_qubit,
        ),
    ]


def _controlled_signed_phase(
    psi: float,
    *,
    selector_qubit: str,
    component_selector_qubit: str,
    phase_qubit: str,
    block_ancillas: list[str],
) -> list[str]:
    signed_angle = -2.0 * psi
    return [
        *_controlled_x_lines(
            open_controls=[component_selector_qubit],
            solid_controls=[selector_qubit],
            target=phase_qubit,
        ),
        "",
        *_controlled_phase_block(
            signed_angle,
            selector_qubit=selector_qubit,
            phase_qubit=phase_qubit,
            block_ancillas=block_ancillas,
        ),
        "",
        *_controlled_x_lines(
            open_controls=[component_selector_qubit],
            solid_controls=[selector_qubit],
            target=phase_qubit,
        ),
    ]


def _controlled_phase_block(
    angle: float,
    *,
    selector_qubit: str,
    phase_qubit: str,
    block_ancillas: list[str],
) -> list[str]:
    return [
        *_controlled_x_lines(
            open_controls=block_ancillas,
            solid_controls=[selector_qubit],
            target=phase_qubit,
        ),
        *_multi_controlled_rz(
            angle,
            open_controls=[],
            solid_controls=[selector_qubit],
            target=phase_qubit,
        ),
        *_controlled_x_lines(
            open_controls=block_ancillas,
            solid_controls=[selector_qubit],
            target=phase_qubit,
        ),
    ]


def _controlled_selector_qsp_execution(
    record: dict[str, Any],
    *,
    selector_qubit: str,
    component_selector_qubit: str,
    selector_value: int,
    phase_qubit: str,
    block_ancillas: list[str],
    system_qubits: list[str],
    controlled_signal_gate: str,
    controlled_signal_gate_dagger: str,
) -> list[str]:
    component = record["component"]
    control_label = f"selector = {selector_value}"
    lines = [f"// {component} branch ({control_label})"]
    for idx, (phase, psi, theta) in enumerate(
        zip(record["pyqsp_phases"], record["qsvt_projector_phases"], record["qasm_rz_angles"])
    ):
        lines.append(f"// {component} phi_{idx} = {phase:.12g}")
        lines.append(f"// {component} psi_{idx} = {psi:.12g}")
        lines.append(f"// {component} theta_rz_{idx} = {theta:.12g}")
        lines.extend(
            _controlled_selector_phase_block(
                theta,
                selector_qubit=selector_qubit,
                component_selector_qubit=component_selector_qubit,
                selector_value=selector_value,
                phase_qubit=phase_qubit,
                block_ancillas=block_ancillas,
            )
        )
        if idx == len(record["pyqsp_phases"]) - 1:
            continue

        use_dagger = idx % 2 == 1
        gate = controlled_signal_gate_dagger if use_dagger else controlled_signal_gate
        label = "controlled U^\\dagger" if use_dagger else "controlled U"
        lines.append("")
        lines.append(f"// {component} {label}")
        lines.extend(
            _controlled_signal_gate_lines(
                gate,
                selector_qubit=selector_qubit,
                selector_value=selector_value,
                block_ancillas=block_ancillas,
                system_qubits=system_qubits,
            )
        )
        lines.append("")
    return lines


def _controlled_signal_gate_lines(
    gate: str,
    *,
    selector_qubit: str,
    selector_value: int,
    block_ancillas: list[str],
    system_qubits: list[str],
) -> list[str]:
    gate_line = f"{gate} {_format_gate_operands([selector_qubit, *block_ancillas, *system_qubits])};"
    if selector_value == 1:
        return [gate_line]
    return [
        f"x {selector_qubit};",
        gate_line,
        f"x {selector_qubit};",
    ]


def _component_extraction_phase_gate(sin_degree: int) -> str:
    if sin_degree % 2 != 1:
        raise ValueError("full Hamiltonian simulation expects the sine component to have odd degree")
    half_degree = (sin_degree - 1) // 2
    return "s" if half_degree % 2 == 1 else "sdg"


def _full_hamsim_needs_global_sign_flip(sin_degree: int) -> bool:
    if sin_degree % 2 != 1:
        raise ValueError("full Hamiltonian simulation expects the sine component to have odd degree")
    return ((sin_degree - 1) // 2) % 2 == 1


def _controlled_selector_phase_block(
    angle: float,
    *,
    selector_qubit: str,
    component_selector_qubit: str,
    selector_value: int,
    phase_qubit: str,
    block_ancillas: list[str],
) -> list[str]:
    open_control_angle = -angle
    open_controls = [*block_ancillas]
    solid_controls: list[str] = []
    if selector_value == 0:
        open_controls.append(selector_qubit)
    else:
        solid_controls.append(selector_qubit)
    return [
        *_open_mcx_controls(open_controls, solid_controls, phase_qubit),
        *_controlled_select_rz(
            open_control_angle,
            selector_qubit=selector_qubit,
            component_selector_qubit=component_selector_qubit,
            selector_value=selector_value,
            phase_qubit=phase_qubit,
        ),
        *_open_mcx_controls(open_controls, solid_controls, phase_qubit),
    ]


def _controlled_select_rz(
    angle: float,
    *,
    selector_qubit: str,
    component_selector_qubit: str,
    selector_value: int,
    phase_qubit: str,
) -> list[str]:
    outer_open: list[str] = []
    outer_solid: list[str] = []
    if selector_value == 0:
        outer_open.append(selector_qubit)
    else:
        outer_solid.append(selector_qubit)

    return [
        *_multi_controlled_rz(
            angle,
            open_controls=[*outer_open, component_selector_qubit],
            solid_controls=outer_solid,
            target=phase_qubit,
        ),
        *_multi_controlled_rz(
            -angle,
            open_controls=outer_open,
            solid_controls=[*outer_solid, component_selector_qubit],
            target=phase_qubit,
        ),
    ]


def _multi_controlled_rz(angle: float, *, open_controls: list[str], solid_controls: list[str], target: str) -> list[str]:
    half_text = f"{angle / 2.0:.12g}"
    neg_half_text = f"{-angle / 2.0:.12g}"
    return [
        *_open_mcx_controls(open_controls, solid_controls, target),
        f"rz({neg_half_text}) {target};",
        *_open_mcx_controls(open_controls, solid_controls, target),
        f"rz({half_text}) {target};",
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
    return _open_mcx_controls(block_ancillas, [], phase_qubit)


def _open_mcx_block_gates(block_ancillas: list[str], phase_qubit: str) -> list[QasmGate]:
    return _open_mcx_controls_gates(block_ancillas, [], phase_qubit)


def _open_mcx_controls(open_controls: list[str], solid_controls: list[str], target: str) -> list[str]:
    controls = ", ".join([*open_controls, *solid_controls, target])
    return [
        *(f"x {control};" for control in open_controls),
        f"mcx {controls};",
        *(f"x {control};" for control in reversed(open_controls)),
    ]


def _open_mcx_controls_gates(open_controls: list[str], solid_controls: list[str], target: str) -> list[QasmGate]:
    return [
        *(QasmGate("x", (control,)) for control in open_controls),
        QasmGate("mcx", tuple([*open_controls, *solid_controls, target])),
        *(QasmGate("x", (control,)) for control in reversed(open_controls)),
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


def _format_complex_coeff(real: float, imag: float, precision: int) -> str | None:
    if real == 0.0 and imag == 0.0:
        return None
    if imag == 0.0:
        return format_number(real, precision)
    if real == 0.0:
        return _format_imaginary(imag, precision)

    sign = "+" if imag > 0 else "-"
    return f"{format_number(real, precision)} {sign} {_format_imaginary(abs(imag), precision)}"


def _format_imaginary(value: float, precision: int) -> str:
    if value == 1.0:
        return "i"
    if value == -1.0:
        return "-i"
    if value < 0:
        return f"-{format_number(abs(value), precision)}*i"
    return f"{format_number(value, precision)}*i"


def _coefficient_factor(coeff: str) -> str:
    if " + " in coeff or " - " in coeff:
        return f"({coeff})"
    return coeff
