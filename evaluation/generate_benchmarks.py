from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS_DIR = ROOT / "benchmarks"
HANDWRITTEN_DIR = BENCHMARKS_DIR / "handwritten"
GENERATED_DIR = BENCHMARKS_DIR / "generated"
MANIFEST_PATH = BENCHMARKS_DIR / "manifest.csv"

SCALING_N_SYSTEM = (1, 2, 4, 6, 8, 10, 12, 16, 20, 32)
BRANCH_TERM_M = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
DEFAULT_DENSE_MAX_N = 12


def main() -> int:
    HANDWRITTEN_DIR.mkdir(parents=True, exist_ok=True)
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    _clear_generated_benchmarks()

    rows = []
    rows.append(_write_handwritten_x_plus_z())
    rows.append(_write_handwritten_xx_plus_zz())

    for n_system in SCALING_N_SYSTEM:
        rows.append(_write_scaling_xstring(n_system))

    for m in BRANCH_TERM_M:
        rows.append(_write_branch_term_product(m))

    _write_manifest(rows)
    print(f"Wrote {len(rows)} benchmarks to {BENCHMARKS_DIR.relative_to(ROOT)}")
    print(f"Wrote manifest to {MANIFEST_PATH.relative_to(ROOT)}")
    return 0


def _write_handwritten_x_plus_z() -> dict[str, str]:
    path = HANDWRITTEN_DIR / "x_plus_z.qasm"
    path.write_text(
        """OPENQASM 2.0;
include "qelib1.inc";

qreg q[2];

// q[0] = ancilla
// q[1] = system

h q[0];

// SELECT X on |0>
x q[0];
cx q[0], q[1];
x q[0];

// SELECT Z on |1>
cz q[0], q[1];

h q[0];
""",
        encoding="utf-8",
    )
    return _manifest_row(
        benchmark_id="x_plus_z",
        kind="handwritten",
        family="x_plus_z",
        qasm_path=path,
        n_system=1,
        systems=(1,),
        expected="(X + Z)/2",
    )


def _write_handwritten_xx_plus_zz() -> dict[str, str]:
    path = HANDWRITTEN_DIR / "xx_plus_zz.qasm"
    path.write_text(
        """OPENQASM 2.0;
include "qelib1.inc";

qreg q[3];

// q[0] = ancilla
// q[1], q[2] = system qubits

h q[0];

// SELECT XX on |0>
x q[0];
cx q[0], q[1];
cx q[0], q[2];
x q[0];

// SELECT ZZ on |1>
cz q[0], q[1];
cz q[0], q[2];

h q[0];
""",
        encoding="utf-8",
    )
    return _manifest_row(
        benchmark_id="xx_plus_zz",
        kind="handwritten",
        family="xx_plus_zz",
        qasm_path=path,
        n_system=2,
        systems=(1, 2),
        expected="(X0 X1 + Z0 Z1)/2",
    )


def _clear_generated_benchmarks() -> None:
    for pattern in ("scaling_xstring_n*.qasm", "branch_term_product_m*.qasm"):
        for path in GENERATED_DIR.glob(pattern):
            path.unlink()


def _write_scaling_xstring(n_system: int) -> dict[str, str]:
    path = GENERATED_DIR / f"scaling_xstring_n{n_system}.qasm"
    system_lines = "\n".join(f"cx q[0], q[{index}];" for index in range(1, n_system + 1))
    path.write_text(
        f"""OPENQASM 2.0;
include "qelib1.inc";

qreg q[{n_system + 1}];

// q[0] = ancilla
// q[1]..q[{n_system}] = system qubits
// B0 = (I + X0 X1 ... X{n_system - 1}) / 2

h q[0];
{system_lines}
h q[0];
""",
        encoding="utf-8",
    )
    x_string = " ".join(f"X{index}" for index in range(n_system))
    return _manifest_row(
        benchmark_id=f"scaling_xstring_n{n_system}",
        kind="generated",
        family="scaling_xstring",
        qasm_path=path,
        n_system=n_system,
        systems=tuple(range(1, n_system + 1)),
        expected=f"(I + {x_string})/2",
    )


def _write_branch_term_product(m: int) -> dict[str, str]:
    path = GENERATED_DIR / f"branch_term_product_m{m}.qasm"
    lines = []
    for index in range(m):
        ancilla = index
        system = m + index
        lines.extend(
            [
                f"h q[{ancilla}];",
                f"cx q[{ancilla}], q[{system}];",
                f"h q[{ancilla}];",
                "",
            ]
        )
    path.write_text(
        f"""OPENQASM 2.0;
include "qelib1.inc";

qreg q[{2 * m}];

// q[0]..q[{m - 1}] = ancilla qubits
// q[{m}]..q[{2 * m - 1}] = system qubits
// B[0...0] = product_j (I + X_j) / 2

{chr(10).join(lines).rstrip()}
""",
        encoding="utf-8",
    )
    return _manifest_row(
        benchmark_id=f"branch_term_product_m{m}",
        kind="generated",
        family="branch_term_product",
        qasm_path=path,
        n_system=m,
        n_ancilla=m,
        ancillas=tuple(range(m)),
        systems=tuple(range(m, 2 * m)),
        expected=_branch_term_expected(m),
        expected_terms=2**m,
    )


def _manifest_row(
    *,
    benchmark_id: str,
    kind: str,
    family: str,
    qasm_path: Path,
    n_system: int,
    systems: tuple[int, ...],
    expected: str,
    n_ancilla: int = 1,
    ancillas: tuple[int, ...] = (0,),
    expected_terms: int | None = None,
) -> dict[str, str]:
    return {
        "benchmark_id": benchmark_id,
        "kind": kind,
        "family": family,
        "qasm_path": str(qasm_path.relative_to(ROOT)),
        "n_system": str(n_system),
        "n_ancilla": str(n_ancilla),
        "ancillas": " ".join(str(ancilla) for ancilla in ancillas),
        "systems": " ".join(str(system) for system in systems),
        "expected": expected,
        "expected_terms": str(expected_terms or _expected_terms(family, n_system)),
        "dense_max_n": str(DEFAULT_DENSE_MAX_N),
    }


def _branch_term_expected(m: int) -> str:
    numerator = " * ".join(f"(I + X{index})" for index in range(m))
    return f"({numerator})/{2**m}"


def _expected_terms(family: str, n_system: int) -> int:
    if family == "scaling_xstring":
        return 2
    if family in {"x_plus_z", "xx_plus_zz"}:
        return 2
    if family == "branch_term_product":
        return 2**n_system
    return 0


def _write_manifest(rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "benchmark_id",
        "kind",
        "family",
        "qasm_path",
        "n_system",
        "n_ancilla",
        "ancillas",
        "systems",
        "expected",
        "expected_terms",
        "dense_max_n",
    ]
    with MANIFEST_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
