from __future__ import annotations

import csv
import os
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "evaluation" / "results"
SYMBOLIC_RESULTS_PATH = RESULTS_DIR / "symbolic_results.csv"
DENSE_RESULTS_PATH = RESULTS_DIR / "dense_results.csv"
SYSTEM_RUNTIME_PNG_PATH = RESULTS_DIR / "system_runtime_vs_n.png"
SYSTEM_MEMORY_PNG_PATH = RESULTS_DIR / "system_memory_vs_n.png"
BRANCH_RUNTIME_PNG_PATH = RESULTS_DIR / "branch_term_runtime_vs_m.png"
BRANCH_MEMORY_PNG_PATH = RESULTS_DIR / "branch_term_memory_vs_m.png"
LEGACY_OUTPUT_PATHS = (
    RESULTS_DIR / "runtime_vs_n.png",
    RESULTS_DIR / "dense_memory_vs_n.png",
)


def main() -> int:
    if not SYMBOLIC_RESULTS_PATH.exists():
        raise SystemExit(f"Missing {SYMBOLIC_RESULTS_PATH.relative_to(ROOT)}. Run evaluation/run_symbolic.py first.")
    if not DENSE_RESULTS_PATH.exists():
        raise SystemExit(f"Missing {DENSE_RESULTS_PATH.relative_to(ROOT)}. Run evaluation/run_dense.py first.")

    plt, mticker = _load_matplotlib()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    symbolic = _read_csv(SYMBOLIC_RESULTS_PATH)
    dense = _read_csv(DENSE_RESULTS_PATH)

    written = [
        _plot_runtime(
            plt,
            mticker,
            symbolic,
            dense,
            family="scaling_xstring",
            x_field="n_system",
            output_path=SYSTEM_RUNTIME_PNG_PATH,
            title="System-Size Scalability: Runtime",
            x_label="system qubits",
        ),
        _plot_memory(
            plt,
            mticker,
            symbolic,
            dense,
            family="scaling_xstring",
            x_field="n_system",
            output_path=SYSTEM_MEMORY_PNG_PATH,
            title="System-Size Scalability: Memory",
            x_label="system qubits",
        ),
        _plot_runtime(
            plt,
            mticker,
            symbolic,
            dense,
            family="branch_term_product",
            x_field="n_ancilla",
            output_path=BRANCH_RUNTIME_PNG_PATH,
            title="Branch/Term-Growth Scalability: Runtime",
            x_label="ancilla qubits m",
        ),
        _plot_memory(
            plt,
            mticker,
            symbolic,
            dense,
            family="branch_term_product",
            x_field="n_ancilla",
            output_path=BRANCH_MEMORY_PNG_PATH,
            title="Branch/Term-Growth Scalability: Memory",
            x_label="ancilla qubits m",
        ),
    ]
    _remove_legacy_outputs()

    for path in written:
        print(f"Wrote {path.relative_to(ROOT)}")
    return 0


def _load_matplotlib():
    cache_root = Path(tempfile.gettempdir()) / "symbolic-be-verifier-plot-cache"
    matplotlib_cache = cache_root / "matplotlib"
    xdg_cache = cache_root / "xdg"
    matplotlib_cache.mkdir(parents=True, exist_ok=True)
    xdg_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))
    os.environ.setdefault("XDG_CACHE_HOME", str(xdg_cache))

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mticker
    except ImportError as exc:
        raise SystemExit("plot_results.py requires matplotlib. Install it with: pip install matplotlib") from exc
    return plt, mticker


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _points(
    rows: list[dict[str, str]],
    *,
    family: str,
    x_field: str,
    y_field: str,
    statuses: set[str] | None = None,
) -> list[tuple[int, float]]:
    points = []
    for row in rows:
        if row["family"] != family:
            continue
        if statuses is not None and row.get("status") not in statuses:
            continue
        value = row.get(y_field, "")
        if not value:
            continue
        try:
            y = float(value)
        except ValueError:
            continue
        if y > 0:
            points.append((int(row[x_field]), y))
    return sorted(points)


def _plot_runtime(
    plt,
    mticker,
    symbolic: list[dict[str, str]],
    dense: list[dict[str, str]],
    *,
    family: str,
    x_field: str,
    output_path: Path,
    title: str,
    x_label: str,
) -> Path:
    symbolic_runtime = _points(symbolic, family=family, x_field=x_field, y_field="runtime_sec", statuses={"PASS"})
    max_symbolic_x = _max_x(symbolic_runtime)
    dense_runtime = _filter_x_at_most(
        _points(dense, family=family, x_field=x_field, y_field="runtime_sec", statuses={"PASS"}),
        max_symbolic_x,
    )
    x_ticks = _all_x_ticks(symbolic_runtime, dense_runtime)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    _plot_series(ax, symbolic_runtime, marker="o", label="symbolic verifier", color="#1f77b4")
    _plot_series(ax, dense_runtime, marker="s", label="dense baseline", color="#d62728")

    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel("runtime (seconds)")
    ax.set_yscale("log")
    ax.set_xticks(x_ticks)
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%d"))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(_format_seconds))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.grid(True, which="major", alpha=0.28)
    ax.grid(True, which="minor", axis="y", alpha=0.12)
    ax.legend()

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def _plot_memory(
    plt,
    mticker,
    symbolic: list[dict[str, str]],
    dense: list[dict[str, str]],
    *,
    family: str,
    x_field: str,
    output_path: Path,
    title: str,
    x_label: str,
) -> Path:
    symbolic_memory = _mib_points_to_bytes(
        _points(symbolic, family=family, x_field=x_field, y_field="tracemalloc_peak_mb", statuses={"PASS"})
    )
    dense_measured_memory = _mib_points_to_bytes(
        _points(dense, family=family, x_field=x_field, y_field="process_maxrss_mb", statuses={"PASS"})
    )
    dense_skipped_estimate = _mib_points_to_bytes(
        _points(dense, family=family, x_field=x_field, y_field="dense_memory_estimate_mb", statuses={"SKIPPED", "TIMEOUT"})
    )
    max_symbolic_x = _max_x(symbolic_memory)
    dense_skipped_estimate = _filter_x_at_most(dense_skipped_estimate, max_symbolic_x)
    dense_estimate_continuation = _estimate_continuation(dense_measured_memory, dense_skipped_estimate)
    x_ticks = _all_x_ticks(symbolic_memory, dense_measured_memory, dense_skipped_estimate)
    y_values = [y for points in (symbolic_memory, dense_measured_memory, dense_estimate_continuation) for _, y in points]

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    _plot_series(ax, symbolic_memory, marker="o", label="symbolic peak memory", color="#1f77b4")
    _plot_series(ax, dense_measured_memory, marker="s", label="dense measured memory", color="#d62728")
    _plot_series(
        ax,
        dense_estimate_continuation,
        marker="s",
        label="dense estimate (not run)",
        color="#d62728",
        linestyle="--",
    )

    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel("memory")
    ax.set_yscale("log")
    ax.set_xticks(x_ticks)
    ax.set_yticks(_binary_memory_ticks(y_values))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%d"))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(_format_memory))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.grid(True, which="major", alpha=0.28)
    ax.grid(True, which="minor", axis="y", alpha=0.12)
    ax.legend()

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def _plot_series(
    ax,
    points: list[tuple[int, float]],
    *,
    marker: str,
    label: str,
    color: str,
    linestyle: str = "-",
) -> None:
    if not points:
        return
    x_values, y_values = zip(*points)
    ax.plot(x_values, y_values, marker=marker, linewidth=2, markersize=5, label=label, color=color, linestyle=linestyle)


def _all_x_ticks(*series: list[tuple[int, float]]) -> list[int]:
    return sorted({x for points in series for x, _ in points})


def _max_x(points: list[tuple[int, float]]) -> int | None:
    if not points:
        return None
    return max(x for x, _ in points)


def _filter_x_at_most(points: list[tuple[int, float]], max_x: int | None) -> list[tuple[int, float]]:
    if max_x is None:
        return points
    return [(x, y) for x, y in points if x <= max_x]


def _estimate_continuation(
    measured: list[tuple[int, float]],
    estimated: list[tuple[int, float]],
) -> list[tuple[int, float]]:
    if not measured or not estimated:
        return estimated
    return [measured[-1], *estimated]


def _mib_points_to_bytes(points: list[tuple[int, float]]) -> list[tuple[int, float]]:
    return [(x, y * 1024 * 1024) for x, y in points]


def _binary_memory_ticks(values: list[float]) -> list[float]:
    if not values:
        return []
    landmarks = set(_lower_memory_ticks(min(values)))
    landmarks.update(1024**power for power in range(0, 9))
    low = min(values) / 3
    high = max(values) * 3
    return [float(tick) for tick in sorted(landmarks) if low <= tick <= high]


def _lower_memory_ticks(min_value: float) -> list[float]:
    units = [1024**power for power in range(0, 9)]
    lower_units = [unit for unit in units if unit <= min_value]
    if not lower_units:
        return [1.0]
    base = lower_units[-1]
    return [base, 4 * base, 16 * base]


def _format_seconds(value: float, _position: int) -> str:
    if value >= 1:
        return f"{value:g}s"
    if value >= 1e-3:
        return f"{value * 1e3:g}ms"
    return f"{value * 1e6:g}us"


def _format_memory(value_bytes: float, _position: int) -> str:
    units = ("B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB")
    value = value_bytes
    unit = units[0]
    for unit in units:
        if abs(value) < 1024 or unit == units[-1]:
            break
        value /= 1024
    if value >= 100:
        return f"{value:.0f} {unit}"
    if value >= 10:
        return f"{value:.1f} {unit}"
    return f"{value:.2g} {unit}"


def _remove_legacy_outputs() -> None:
    for path in LEGACY_OUTPUT_PATHS:
        if path.exists():
            path.unlink()


if __name__ == "__main__":
    raise SystemExit(main())
