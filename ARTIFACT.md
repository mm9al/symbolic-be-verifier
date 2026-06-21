# Artifact Guide: Symbolic Block-Encoding Verifier

This artifact accompanies the final report "Symbolic Verification of Block-Encoding Circuits".
It contains the prototype verifier, OpenQASM examples, raw evaluation data, and scripts for
reproducing the main experimental plots and running smaller reviewer-friendly reruns.

## 1. Artifact Overview

The verifier checks the top-left block of block-encoding circuits by symbolically tracking
ancilla branches:

```text
|Phi> = sum_b |b>_a B_b |psi>
```

The main checked property is whether the final all-zero branch `B[00...0]` equals the
expected target operator.

The artifact supports these claims from the report:

- The verifier can prove correctness of small LCU-style block-encoding circuits.
- The verifier supports multiple ancilla and system qubits.
- The verifier avoids constructing the full dense unitary matrix.
- The evaluation data and scripts reproduce the reported scalability trends for system size
  and branch/term growth.

## 2. Requirements

Tested environment:

- Python 3.10 or later
- Linux or macOS
- Python packages listed in `requirements.txt`

No GPU is required.

## 3. Installation

From a fresh clone:

```bash
git clone https://github.com/mm9al/symbolic-be-verifier.git
cd symbolic-be-verifier

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## 4. Smoke Test

Run:

```bash
scripts/run_smoke_tests.sh
```

Expected result:

```text
All smoke tests passed.
```

The smoke test covers one-ancilla LCU verification, two-ancilla branch tracking,
two-system-qubit operators, and a symbolic rotation-angle example.

## 5. Reproducing the Evaluation

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

Reviewer-friendly reruns are available separately. They intentionally use smaller defaults
than the full historical run.

RQ1, system-size scalability:

```bash
python scripts/reproduce_rq1.py
```

Output:

```text
results/rq1_system_size_rerun.csv
```

RQ2, branch/term-growth scalability:

```bash
python scripts/reproduce_rq2.py
```

Output:

```text
results/rq2_branch_term_growth_rerun.csv
```

The default RQ2 rerun stops at `m=4` because larger branch/term-growth cases can take
minutes to hours on a laptop. Use `--max-m` or `--max-n-system` to change the rerun sizes.

## 6. Expected Runtime

The smoke tests should finish within a few seconds.
Plot recreation from `data/raw/` should also finish within a few seconds.
The full branch/term-growth benchmark is much longer; the included raw CSV preserves the
full run used for the report.

Runtime numbers may vary across machines, but the qualitative trends should remain the same.

## 7. Directory Structure

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

## 8. Verification Results

The verifier reports:

- `PASS` if the final all-zero branch exactly matches the expected operator.
- `PASS_UP_TO_GLOBAL_PHASE` if the result differs only by a unit global phase.
- `FAIL` otherwise.

## 9. DOI Release Checklist

Before submitting the artifact:

1. Push this repository to GitHub and make it public.
2. Connect the repository in Zenodo's GitHub integration.
3. Create a GitHub release, for example `v1.0.0-artifact`.
4. Zenodo will archive that release and assign a DOI.
5. Replace the DOI placeholder in `README.md`, add the DOI to `CITATION.cff`, and update
   the report artifact section with the Zenodo DOI.
