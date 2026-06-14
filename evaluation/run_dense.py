from __future__ import annotations

import argparse
import csv
import math
import multiprocessing as mp
import re
import resource
import sys
import time
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "benchmarks" / "manifest.csv"
RESULTS_DIR = ROOT / "evaluation" / "results"
RESULTS_PATH = RESULTS_DIR / "dense_results.csv"
DEFAULT_MAX_UNITARY_MEMORY_MB = 4096.0
DEFAULT_TIMEOUT_SEC = 300.0
TOLERANCE = 1e-9

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - exercised only on missing deps
    raise SystemExit("run_dense.py requires numpy. Install it with: pip install numpy") from exc


@dataclass(frozen=True)
class Gate:
    name: str
    qubits: tuple[int, ...]
    parameter: object = None


_GATE_RE = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)\s*(?:\(([^)]*)\))?\s+(.+);$")
_QUBIT_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*\[(\d+)\]$")


def main() -> int:
    args = _parse_args()
    if not MANIFEST_PATH.exists():
        raise SystemExit(f"Missing {MANIFEST_PATH.relative_to(ROOT)}. Run evaluation/generate_benchmarks.py first.")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rows = [
        _run_benchmark(
            row,
            max_unitary_memory_mb=args.max_unitary_memory_mb,
            max_n_system=args.max_n_system,
            timeout_sec=args.timeout_sec,
        )
        for row in _read_manifest()
    ]
    _write_results(rows)
    print(f"Wrote {RESULTS_PATH.relative_to(ROOT)}")
    return 0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run dense-unitary baselines for benchmarks that fit a memory budget.")
    parser.add_argument(
        "--max-unitary-memory-mb",
        type=float,
        default=DEFAULT_MAX_UNITARY_MEMORY_MB,
        help=(
            "Maximum estimated full-unitary matrix size to run, in MiB. "
            f"Default: {DEFAULT_MAX_UNITARY_MEMORY_MB:g}. "
            "This is a lower-bound estimate; temporary arrays can need more memory."
        ),
    )
    parser.add_argument(
        "--max-n-system",
        type=int,
        default=None,
        help="Optional hard cap on system qubits, in addition to the memory budget.",
    )
    parser.add_argument(
        "--timeout-sec",
        type=float,
        default=DEFAULT_TIMEOUT_SEC,
        help=f"Per-benchmark timeout in seconds. Use 0 to disable. Default: {DEFAULT_TIMEOUT_SEC:g}.",
    )
    return parser.parse_args()


def _read_manifest() -> list[dict[str, str]]:
    with MANIFEST_PATH.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _run_benchmark(
    row: dict[str, str],
    *,
    max_unitary_memory_mb: float,
    max_n_system: int | None,
    timeout_sec: float,
) -> dict[str, str]:
    return _run_benchmark_with_timeout(
        row,
        max_unitary_memory_mb=max_unitary_memory_mb,
        max_n_system=max_n_system,
        timeout_sec=timeout_sec,
    )


def _run_benchmark_with_timeout(
    row: dict[str, str],
    *,
    max_unitary_memory_mb: float,
    max_n_system: int | None,
    timeout_sec: float,
) -> dict[str, str]:
    base = _dense_metadata(row)
    skip_reason = _skip_reason(
        n_system=int(row["n_system"]),
        memory_estimate_mb=float(base["dense_memory_estimate_mb"]),
        max_unitary_memory_mb=max_unitary_memory_mb,
        max_n_system=max_n_system,
    )
    if skip_reason is not None:
        return {
            **base,
            "runtime_sec": "",
            "process_maxrss_mb": "",
            "max_abs_error": "",
            "status": "SKIPPED",
            "success": "",
            "error": skip_reason,
        }

    if timeout_sec <= 0:
        return _run_benchmark_inner(row)

    queue: mp.Queue = mp.Queue(maxsize=1)
    process = mp.Process(target=_run_benchmark_worker, args=(row, queue))
    started = time.perf_counter()
    process.start()
    process.join(timeout_sec)
    if process.is_alive():
        process.terminate()
        process.join()
        runtime_sec = time.perf_counter() - started
        return {
            **base,
            "runtime_sec": f"{runtime_sec:.9f}",
            "process_maxrss_mb": "",
            "max_abs_error": "",
            "status": "TIMEOUT",
            "success": "False",
            "error": f"timeout after {timeout_sec:g}s",
        }
    if not queue.empty():
        return queue.get()
    runtime_sec = time.perf_counter() - started
    return {
        **base,
        "runtime_sec": f"{runtime_sec:.9f}",
        "process_maxrss_mb": "",
        "max_abs_error": "",
        "status": "ERROR",
        "success": "False",
        "error": f"worker exited with code {process.exitcode}",
    }


def _run_benchmark_worker(row: dict[str, str], queue: "mp.Queue") -> None:
    queue.put(_run_benchmark_inner(row))


def _dense_metadata(row: dict[str, str]) -> dict[str, str]:
    n_system = int(row["n_system"])
    n_ancilla = int(row["n_ancilla"])
    num_qubits = n_system + n_ancilla
    dim = 1 << num_qubits
    memory_estimate_mb = _dense_unitary_estimate_mb(num_qubits)
    return {
        "benchmark_id": row["benchmark_id"],
        "kind": row["kind"],
        "family": row["family"],
        "n_system": row["n_system"],
        "n_ancilla": row["n_ancilla"],
        "qasm_path": row["qasm_path"],
        "matrix_dim": str(dim),
        "dense_memory_estimate_mb": f"{memory_estimate_mb:.6f}",
    }


def _run_benchmark_inner(row: dict[str, str]) -> dict[str, str]:
    n_system = int(row["n_system"])
    num_qubits = n_system + int(row["n_ancilla"])
    base = _dense_metadata(row)
    started = time.perf_counter()
    try:
        unitary = _simulate_unitary(ROOT / row["qasm_path"], num_qubits)
        actual = _top_left_block(unitary, num_qubits, _parse_int_list(row["ancillas"]))
        expected = _expected_matrix(row["family"], n_system)
        max_abs_error = float(np.max(np.abs(actual - expected)))
        success = max_abs_error <= TOLERANCE
        status = "PASS" if success else "FAIL"
        error = ""
    except Exception as exc:  # keep the evaluation run going
        max_abs_error = math.nan
        success = False
        status = "ERROR"
        error = f"{type(exc).__name__}: {exc}"
    runtime_sec = time.perf_counter() - started

    return {
        **base,
        "runtime_sec": f"{runtime_sec:.9f}",
        "process_maxrss_mb": f"{_process_maxrss_mb():.6f}",
        "max_abs_error": "" if math.isnan(max_abs_error) else f"{max_abs_error:.9e}",
        "status": status,
        "success": str(success),
        "error": error,
    }


def _skip_reason(
    *,
    n_system: int,
    memory_estimate_mb: float,
    max_unitary_memory_mb: float,
    max_n_system: int | None,
) -> str | None:
    if max_n_system is not None and n_system > max_n_system:
        return f"n_system > {max_n_system}"
    if memory_estimate_mb > max_unitary_memory_mb:
        return f"dense unitary estimate {memory_estimate_mb:.1f} MiB > budget {max_unitary_memory_mb:.1f} MiB"
    return None


def _simulate_unitary(qasm_path: Path, num_qubits: int) -> "np.ndarray":
    dim = 1 << num_qubits
    unitary = np.eye(dim, dtype=np.complex128)
    for gate in _parse_qasm_file(qasm_path):
        name = gate.name.lower()
        if len(gate.qubits) == 1:
            unitary = _apply_single_qubit_gate(unitary, _single_qubit_matrix(name, gate.parameter), gate.qubits[0], num_qubits)
        elif len(gate.qubits) == 2:
            unitary = _apply_two_qubit_gate(unitary, name, gate.qubits[0], gate.qubits[1], num_qubits)
        else:
            raise ValueError(f"Unsupported gate arity for dense baseline: {gate}")
    return unitary


def _parse_qasm_file(path: Path) -> list[Gate]:
    gates = []
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw_line.split("//", 1)[0].strip()
        if not line:
            continue
        if (
            line.startswith("OPENQASM")
            or line.startswith("include")
            or line.startswith("qreg ")
            or line.startswith("qubit[")
            or line.startswith("creg ")
            or line.startswith("barrier ")
        ):
            continue
        match = _GATE_RE.match(line)
        if match is None:
            raise ValueError(f"Unsupported QASM statement on line {line_number}: {raw_line.strip()}")
        name, parameter_text, operands_text = match.groups()
        qubits = tuple(_parse_qubit(operand.strip(), line_number) for operand in operands_text.split(","))
        parameter = _parse_parameter(parameter_text) if parameter_text is not None else None
        gates.append(Gate(name=name.lower(), qubits=qubits, parameter=parameter))
    return gates


def _parse_qubit(text: str, line_number: int) -> int:
    match = _QUBIT_RE.match(text)
    if match is None:
        raise ValueError(f"Expected qubit operand on line {line_number}: {text!r}")
    return int(match.group(1))


def _parse_parameter(text: str) -> float:
    cleaned = text.strip().replace("^", "**")
    allowed = {"pi": math.pi}
    return float(eval(cleaned, {"__builtins__": {}}, allowed))


def _single_qubit_matrix(name: str, theta: object) -> "np.ndarray":
    if name == "h":
        return np.array([[1, 1], [1, -1]], dtype=np.complex128) / np.sqrt(2)
    if name == "x":
        return np.array([[0, 1], [1, 0]], dtype=np.complex128)
    if name == "z":
        return np.array([[1, 0], [0, -1]], dtype=np.complex128)
    if name == "s":
        return np.array([[1, 0], [0, 1j]], dtype=np.complex128)
    if name == "sdg":
        return np.array([[1, 0], [0, -1j]], dtype=np.complex128)
    if name in {"rx", "ry", "rz"}:
        angle = complex(theta)
        c = np.cos(angle / 2)
        s = np.sin(angle / 2)
        if name == "rx":
            return np.array([[c, -1j * s], [-1j * s, c]], dtype=np.complex128)
        if name == "ry":
            return np.array([[c, -s], [s, c]], dtype=np.complex128)
        return np.array([[np.exp(-0.5j * angle), 0], [0, np.exp(0.5j * angle)]], dtype=np.complex128)
    raise ValueError(f"Unsupported dense single-qubit gate: {name}")


def _apply_single_qubit_gate(unitary: "np.ndarray", matrix: "np.ndarray", target: int, num_qubits: int) -> "np.ndarray":
    dim = unitary.shape[0]
    step = 1 << (num_qubits - 1 - target)
    period = step * 2
    for base in range(0, dim, period):
        for offset in range(step):
            row0 = base + offset
            row1 = row0 + step
            old0 = unitary[row0, :].copy()
            old1 = unitary[row1, :].copy()
            unitary[row0, :] = matrix[0, 0] * old0 + matrix[0, 1] * old1
            unitary[row1, :] = matrix[1, 0] * old0 + matrix[1, 1] * old1
    return unitary


def _apply_two_qubit_gate(unitary: "np.ndarray", name: str, left: int, right: int, num_qubits: int) -> "np.ndarray":
    if name == "cx":
        for row0 in range(unitary.shape[0]):
            if _bit(row0, left, num_qubits) and not _bit(row0, right, num_qubits):
                row1 = _flip_bit(row0, right, num_qubits)
                tmp = unitary[row0, :].copy()
                unitary[row0, :] = unitary[row1, :]
                unitary[row1, :] = tmp
        return unitary
    if name == "cz":
        for row in range(unitary.shape[0]):
            if _bit(row, left, num_qubits) and _bit(row, right, num_qubits):
                unitary[row, :] *= -1
        return unitary
    raise ValueError(f"Unsupported dense two-qubit gate: {name}")


def _top_left_block(unitary: "np.ndarray", num_qubits: int, ancillas: tuple[int, ...]) -> "np.ndarray":
    indices = [index for index in range(unitary.shape[0]) if all(_bit(index, ancilla, num_qubits) == 0 for ancilla in ancillas)]
    return unitary[np.ix_(indices, indices)]


def _expected_matrix(family: str, n_system: int) -> "np.ndarray":
    if family == "x_plus_z":
        return (_pauli_product("X") + _pauli_product("Z")) / 2
    if family == "xx_plus_zz":
        return (_pauli_product("XX") + _pauli_product("ZZ")) / 2
    if family == "scaling_xstring":
        return (np.eye(1 << n_system, dtype=np.complex128) + _pauli_product("X" * n_system)) / 2
    if family == "branch_term_product":
        single = (np.eye(2, dtype=np.complex128) + _pauli_product("X")) / 2
        out = np.array([[1]], dtype=np.complex128)
        for _ in range(n_system):
            out = np.kron(out, single)
        return out
    raise ValueError(f"Unsupported dense expected family: {family}")


def _pauli_product(ops: str) -> "np.ndarray":
    matrices = {
        "I": np.eye(2, dtype=np.complex128),
        "X": np.array([[0, 1], [1, 0]], dtype=np.complex128),
        "Y": np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
        "Z": np.array([[1, 0], [0, -1]], dtype=np.complex128),
    }
    out = np.array([[1]], dtype=np.complex128)
    for op in ops:
        out = np.kron(out, matrices[op])
    return out


def _bit(index: int, qubit: int, num_qubits: int) -> int:
    return (index >> (num_qubits - 1 - qubit)) & 1


def _flip_bit(index: int, qubit: int, num_qubits: int) -> int:
    return index ^ (1 << (num_qubits - 1 - qubit))


def _parse_int_list(text: str) -> tuple[int, ...]:
    return tuple(int(piece) for piece in text.split())


def _dense_unitary_estimate_mb(num_qubits: int) -> float:
    dim = 1 << num_qubits
    return dim * dim * np.dtype(np.complex128).itemsize / (1024 * 1024)


def _process_maxrss_mb() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return rss / (1024 * 1024)
    return rss / 1024


def _write_results(rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "benchmark_id",
        "kind",
        "family",
        "n_system",
        "n_ancilla",
        "qasm_path",
        "matrix_dim",
        "dense_memory_estimate_mb",
        "runtime_sec",
        "process_maxrss_mb",
        "max_abs_error",
        "status",
        "success",
        "error",
    ]
    with RESULTS_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
