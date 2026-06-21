# Artifact Guide: Symbolic Block-Encoding Verifier

Artifact DOI: `10.5281/zenodo.20781657`

This artifact accompanies the final report **"Symbolic Verification of Block-Encoding Circuits"**.

The artifact contains the prototype symbolic verifier, OpenQASM examples, raw evaluation data, expected plots, and scripts for reproducing the main experimental plots and running smaller reviewer-friendly reruns.

## 1. How to Obtain the Artifact

Open the Zenodo DOI record:

```text
10.5281/zenodo.20781657
```

Download the archived artifact package from Zenodo and extract it.

After extraction, enter the extracted directory that contains this `ARTIFACT.md` file and `requirements.txt`.

For example:

```bash
cd <extracted-artifact-directory>
```

No Git checkout or branch selection is needed. The Zenodo archive is the fixed artifact version for evaluation.

## 2. Artifact Overview

The verifier checks the top-left block of block-encoding circuits by symbolically tracking ancilla branches:

```text
|Phi> = sum_b |b>_a B_b |psi>
```

The main checked property is whether the final all-zero branch `B[00...0]` equals the expected target operator.

The artifact supports these claims from the report:

* The verifier can prove correctness of small LCU-style block-encoding circuits.
* The verifier supports multiple ancilla and system qubits.
* The verifier avoids constructing the full dense unitary matrix.
* The evaluation data and scripts reproduce the reported scalability trends for system size and branch/term growth.

## 3. Requirements

Tested environment:

* Python 3.10 or later
* Linux or macOS
* Python packages listed in `requirements.txt`

No GPU is required.

## 4. Installation

From the extracted artifact directory, create a fresh Python environment and install the dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

All commands below assume that the virtual environment is activated and that the current working directory is the extracted artifact directory.

## 5. Smoke Test

Run:

```bash
bash scripts/run_smoke_tests.sh
```

Expected result:

```text
All smoke tests passed.
```

The smoke test covers one-ancilla LCU verification, two-ancilla branch tracking, two-system-qubit operators, and a symbolic rotation-angle example.

## 6. Reproducing the Evaluation

The full raw CSV data used to generate the report figures is included under `data/raw/`.

Recreate the plots from those CSVs with:

```bash
python scripts/plot_results.py
```

This writes:

```text
results/system_runtime_vs_n.png
results/system_memory_vs_n.png
results/branch_term_runtime_vs_m.png
results/branch_term_memory_vs_m.png
```

The expected versions of the plots are also archived under:

```text
results/expected/
```

Reviewer-friendly reruns are available separately. They intentionally use smaller defaults than the full historical run.

### RQ1: System-size scalability

Run:

```bash
python scripts/reproduce_rq1.py
```

Output:

```text
results/rq1_system_size_rerun.csv
```

### RQ2: Branch / term-growth scalability

Run:

```bash
python scripts/reproduce_rq2.py
```

Output:

```text
results/rq2_branch_term_growth_rerun.csv
```

The default RQ2 rerun stops at `m=4` because larger branch/term-growth cases can take minutes to hours on a laptop. Use `--max-m` or `--max-n-system` to change the rerun sizes.

## 7. Running the Unit Tests

Run:

```bash
python -m pytest
```

Expected result:

```text
All tests should pass.
```

In the archived artifact version, the test suite contains 34 tests.

## 8. Expected Runtime

The smoke tests should finish within a few seconds.

Plot recreation from `data/raw/` should also finish within a few seconds.

The reviewer-friendly RQ1 and RQ2 reruns should be suitable for a laptop. The full branch/term-growth benchmark is much longer; the included raw CSV preserves the full run used for the report.

Runtime numbers may vary across machines, but the qualitative trends should remain the same.

## 9. Directory Structure

```text
symbolic/            verifier implementation
examples/            small OpenQASM examples for smoke tests
benchmarks/          generated and handwritten benchmark circuits
evaluation/          original benchmark generation, dense baseline, and plotting code
scripts/             artifact-facing smoke, rerun, and plotting entry points
data/raw/            raw CSV evaluation data used in the report
results/             reviewer-generated CSVs and plots
results/expected/    expected plots archived with the artifact
tests/               unit tests
```

## 10. Verification Results

The verifier reports:

* `PASS` if the final all-zero branch exactly matches the expected operator.
* `PASS_UP_TO_GLOBAL_PHASE` if the result differs only by a unit global phase.
* `FAIL` otherwise.

## 11. Suggested Evaluation Checklist

A reviewer can evaluate the artifact by running the following commands from the extracted artifact directory:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

bash scripts/run_smoke_tests.sh
python scripts/plot_results.py
python scripts/reproduce_rq1.py
python scripts/reproduce_rq2.py
python -m pytest
```

The expected outcome is that all smoke tests pass, all plots are regenerated, the reviewer-friendly rerun CSVs are produced, and the unit tests pass.
