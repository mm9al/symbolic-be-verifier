"""Export actual PennyLane PrepSelPrep circuits for the periodic ZZ+X chain."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from fractions import Fraction
import importlib.metadata
import math
from pathlib import Path
import platform
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from evaluation.pennylane.common import (BENCHMARKS, DEFAULT_SIZES, gate, ising_terms, parse_sizes,
                                        sha256, target_expression, write_json)


# Fully decompose controlled gates, including Toffoli: the current verifier's
# mcx primitive only permits ancilla targets. No custom SELECT is substituted.
NAMES = {"Hadamard": "h", "PauliX": "x", "PauliY": "y", "PauliZ": "z",
         "S": "s", "Adjoint(S)": "sdg", "T": "t", "Adjoint(T)": "tdg",
         "RX": "rx", "RY": "ry", "RZ": "rz", "CNOT": "cx", "CY": "cy", "CZ": "cz"}


def lower(ops, stop):
    import pennylane as qml
    return qml.transforms.decompose(qml.tape.QuantumScript(ops),
                                    stopping_condition=stop)[0][0].operations


def export(ops):
    body, phases = [], []
    for op in lower(ops, lambda o: o.name in NAMES or o.name == "GlobalPhase"):
        if op.name == "GlobalPhase":
            phases.append(format(float(op.data[0]), ".17g"))
        elif op.name in NAMES:
            angle = format(float(op.data[0]), ".17g") if op.data else None
            if angle is not None and not math.isfinite(float(angle)):
                raise ValueError("Nonfinite source angle")
            body.append(gate(NAMES[op.name], op.wires, angle))
        else:
            raise ValueError(f"PennyLane could not lower {op.name} into the export basis")
    return body, phases


def numerical_crosscheck(be, h, body, gamma, blocks, systems):
    """Small export regression only; this is not the TASK 1 certificate."""
    import numpy as np
    import pennylane as qml
    import re
    constructors = {"h": qml.Hadamard, "x": qml.X, "y": qml.Y, "z": qml.Z,
                    "s": qml.S, "sdg": qml.adjoint(qml.S), "t": qml.T,
                    "tdg": qml.adjoint(qml.T), "rx": qml.RX, "ry": qml.RY,
                    "rz": qml.RZ, "cx": qml.CNOT, "cy": qml.CY, "cz": qml.CZ}
    ops = []
    for line in body:
        name = line.split("(")[0].split()[0]
        wires = list(map(int, re.findall(r"q\[(\d+)\]", line)))
        parameter = re.search(r"\(([^)]+)\)", line)
        args = [] if parameter is None else [float(parameter.group(1))]
        ops.append(constructors[name](*args, wires=wires))
    wires = blocks + systems
    source = qml.matrix(be, wire_order=wires)
    exported = np.exp(-1j * float(gamma)) * qml.matrix(qml.tape.QuantumScript(ops), wire_order=wires)
    target = qml.matrix(h, wire_order=systems) / sum(abs(c) for c in h.terms()[0])
    dimension = 2 ** len(systems)
    errors = dict(export_max_abs=float(np.max(np.abs(exported - source))),
                  source_block_max_abs=float(np.max(np.abs(source[:dimension, :dimension] - target))))
    if any(not math.isfinite(error) or error > 1e-10 for error in errors.values()):
        raise ValueError(f"Numerical source/export cross-check failed: {errors}")
    return dict(status="PASS", threshold=1e-10, **errors)


def generate_case(output: Path, n: int, matrix_max_qubits: int = 8) -> dict:
    import pennylane as qml
    from pennylane.ops.op_math import Adjoint
    started = time.perf_counter()
    # qml.dot may merge duplicate terms (n=2 has 2*ZZ + X0 + X1).
    # Determine the register size from the actual PennyLane LCU.
    ops = [qml.Z(i) @ qml.Z((i + 1) % n) for i in range(n)] + [qml.X(i) for i in range(n)]
    h = qml.dot([1.0] * (2 * n), ops)
    coeffs, _ = h.terms()
    m = (len(coeffs) - 1).bit_length()
    blocks, systems = list(range(m)), list(range(m, m + n))
    h = qml.map_wires(h, dict(enumerate(systems)))
    be = qml.PrepSelPrep(h, control=blocks)

    def boundary(op):
        base = op.base if isinstance(op, Adjoint) else op
        return isinstance(base, (qml.AmplitudeEmbedding, qml.Select))

    parts = lower([be], boundary)
    if (len(parts) != 3 or not isinstance(parts[0], qml.AmplitudeEmbedding)
            or not isinstance(parts[1], qml.Select)
            or not isinstance(parts[2], Adjoint)
            or not qml.equal(parts[2].base, parts[0])):
        raise ValueError("Unexpected native PrepSelPrep boundaries")
    segments = [export([part]) for part in parts]
    body, phases = export([be])
    if body != sum((segment[0] for segment in segments), []):
        raise ValueError("Segment exports differ from the full native decomposition")
    gamma = sum(map(Fraction, phases), Fraction())
    segment_gamma = sum((sum(map(Fraction, item[1]), Fraction()) for item in segments), Fraction())
    if gamma != segment_gamma:
        raise ValueError("Segment global phases differ from the full export")

    benchmark_id = f"ising_periodic_n{n}"
    directory = output / "ising" / benchmark_id
    directory.mkdir(parents=True, exist_ok=True)
    qasm = directory / "block_encoding.qasm"
    qasm.write_text("\n".join(['OPENQASM 2.0;', 'include "qelib1.inc";',
        f"qreg q[{m+n}];", f"// PennyLane PrepSelPrep: {benchmark_id}",
        f"// target H = {target_expression(n)}", "// PREPARE",
        *segments[0][0], "// SELECT", *segments[1][0], "// PREPARE dagger", *segments[2][0]]) + "\n")
    structure = directory / "structure.json"
    write_json(structure, dict(prepare=segments[0][0], select=segments[1][0], unprepare=segments[2][0],
        global_phases=phases, segment_global_phases=[item[1] for item in segments],
        gamma=str(gamma), source_select=repr(parts[1]),
        origin="Unmodified native PennyLane decomposition of the same PrepSelPrep object"))
    crosscheck = (numerical_crosscheck(be, h, body, gamma, blocks, systems)
                  if m + n <= matrix_max_qubits else
                  dict(status="NOT_RUN", reason=f"More than {matrix_max_qubits} total qubits"))
    metadata = dict(benchmark_id=benchmark_id, model="ising", graph="periodic", n_system=n,
        target_lcu_terms=2*n, n_terms=len(coeffs), n_ancilla=m, ancillas=blocks, systems=systems,
        target_H=target_expression(n), independent_terms=ising_terms(n), nominal_alpha=str(2*n),
        source_coefficients=[format(float(c), ".17g") for c in coeffs],
        source_operators=[repr(o) for o in h.terms()[1]], gate_count=len(body),
        gate_counts=dict(Counter(line.split("(")[0].split()[0] for line in body)),
        global_phase_gamma=str(gamma), decimal_angle_digits=17,
        added_workspace_wires=0, select_reconstructed=False, pi_snapping=False,
        numerical_crosscheck=crosscheck, generation_sec=time.perf_counter()-started,
        python=platform.python_version(),
        packages={name: importlib.metadata.version(name) for name in ("pennylane", "numpy")},
        generator_sources_sha256={name: sha256(HERE/name) for name in ("common.py", "generate_benchmarks.py")},
        files_sha256={p.name: sha256(p) for p in (qasm, structure)})
    write_json(directory / "metadata.json", metadata)
    print(f"Generated {benchmark_id}: {len(coeffs)} source terms, {m} ancillas, {len(body)} gates", flush=True)
    return dict(benchmark_id=benchmark_id, model="ising", graph="periodic", n_system=n,
        n_terms=len(coeffs), target_lcu_terms=2*n, n_ancilla=m, nominal_alpha=2*n, gate_count=len(body),
        qasm_path=str(qasm.relative_to(output)), metadata_path=str((directory/"metadata.json").relative_to(output)),
        metadata_sha256=sha256(directory/"metadata.json"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", default=DEFAULT_SIZES, help="Comma-separated sizes or inclusive ranges, e.g. 2..8")
    parser.add_argument("--out-dir", type=Path, default=BENCHMARKS)
    parser.add_argument("--matrix-max-qubits", type=int, default=8, help="Limit small numerical export checks; 0 disables them")
    args = parser.parse_args()
    try:
        sizes = parse_sizes(args.sizes)
    except ValueError as exc:
        parser.error(str(exc))
    if args.matrix_max_qubits < 0:
        parser.error("--matrix-max-qubits must be nonnegative")
    output = args.out_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    with (output/"manifest.csv").open("w", newline="") as handle:
        writer = None
        for n in sizes:
            row = generate_case(output, n, args.matrix_max_qubits)
            if writer is None:
                writer = csv.DictWriter(handle, fieldnames=list(row), lineterminator="\n")
                writer.writeheader()
            writer.writerow(row)
            handle.flush()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
