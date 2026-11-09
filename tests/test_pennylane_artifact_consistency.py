"""Reject the stale or mixed artifacts that can invalidate timing comparisons."""
import csv
import json
from pathlib import Path
import shutil

import pytest

from evaluation.pennylane import refresh_comparisons as reports


def test_saved_comparisons_checks_and_readme_match_the_recorded_runs():
    for path, expected in reports.build_outputs().items():
        assert path.read_bytes() == expected.encode(), str(path)


def test_mixed_csv_and_case_json_cannot_supply_a_timing(tmp_path):
    source = reports.HERE/"results_bounded_pruning"
    with (source/"block_encoding_results.csv").open(newline="") as handle:
        row = next(csv.DictReader(handle))
    destination = tmp_path/row["report_path"]
    destination.parent.mkdir(parents=True)
    shutil.copyfile(source/row["report_path"], destination)
    row["runtime_sec"] = "999"
    with (tmp_path/"block_encoding_results.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    with pytest.raises(ValueError, match="CSV/JSON mismatch"):
        reports.load_results(tmp_path)


def test_changed_source_snapshot_cannot_be_marked_as_matching(tmp_path):
    source = reports.HERE/"results_bounded_pruning"
    metadata = json.loads((source/"run_metadata.json").read_text())
    name, relative = next(iter(metadata["source_snapshots"].items()))
    snapshot = tmp_path/relative
    snapshot.parent.mkdir(parents=True)
    snapshot.write_text("# Different source version\n")
    metadata["source_snapshots"] = {name: relative}
    (tmp_path/"run_metadata.json").write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="Source snapshot hash differs"):
        reports.snapshot_checks(tmp_path)


@pytest.mark.parametrize("change", ["input", "setting", "coefficient"])
def test_incompatible_reports_cannot_produce_a_speedup(change):
    path = Path("cases/ising_periodic_n2.json")
    left = json.loads((reports.HERE/"results_local_summaries"/path).read_text())
    right = json.loads((reports.HERE/"results_bounded_pruning"/path).read_text())
    if change == "input":
        right["qasm_sha256"] = "different-input"
    elif change == "setting":
        right["summary"]["hermiticity_delta"] = "1e-8"
    else:
        right["actual_block_coefficients"]["XX"] = "1"
    with pytest.raises(ValueError, match="mismatch|incompatible"):
        reports.check_pair(left, right)
