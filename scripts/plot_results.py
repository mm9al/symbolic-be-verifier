from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data" / "raw"
RESULTS_DIR = ROOT / "results"

sys.path.insert(0, str(ROOT))

from evaluation import plot_results as evaluation_plot_results  # noqa: E402


def main() -> int:
    symbolic_csv = DATA_RAW / "symbolic_results.csv"
    dense_csv = DATA_RAW / "dense_results.csv"
    if not symbolic_csv.exists() or not dense_csv.exists():
        raise SystemExit("Missing raw CSV data. Expected data/raw/symbolic_results.csv and data/raw/dense_results.csv.")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    evaluation_plot_results.RESULTS_DIR = RESULTS_DIR
    evaluation_plot_results.SYMBOLIC_RESULTS_PATH = symbolic_csv
    evaluation_plot_results.DENSE_RESULTS_PATH = dense_csv
    evaluation_plot_results.SYSTEM_RUNTIME_PNG_PATH = RESULTS_DIR / "system_runtime_vs_n.png"
    evaluation_plot_results.SYSTEM_MEMORY_PNG_PATH = RESULTS_DIR / "system_memory_vs_n.png"
    evaluation_plot_results.BRANCH_RUNTIME_PNG_PATH = RESULTS_DIR / "branch_term_runtime_vs_m.png"
    evaluation_plot_results.BRANCH_MEMORY_PNG_PATH = RESULTS_DIR / "branch_term_memory_vs_m.png"
    evaluation_plot_results.LEGACY_OUTPUT_PATHS = (
        RESULTS_DIR / "runtime_vs_n.png",
        RESULTS_DIR / "dense_memory_vs_n.png",
    )

    return evaluation_plot_results.main()


if __name__ == "__main__":
    raise SystemExit(main())
