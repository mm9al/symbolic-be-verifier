# Bug-Detection Evaluation

This directory contains an independent bug-detection suite for the symbolic
verifier. It generates fresh block-encoding and full Hamiltonian-simulation QSP
circuits under `generated/`, injects controlled defects, and runs them through
the same RQ1 and RQ2 evaluation paths used by the main evaluation scripts.

Run:

```bash
.venv/bin/python evaluation/bug_detection/run_bug_detection.py
```

Outputs:

```text
evaluation/bug_detection/results/block_encoding_bug_detection_results.csv
evaluation/bug_detection/results/hamsim_bug_detection_results.csv
```

The RQ1 CSV uses the same fields as `evaluation/results/block_encoding_results.csv`.
The RQ2 CSV uses the same fields as `evaluation/results/hamsim_results.csv`.
Both append only one final field, `bug_detection_result`. Its value is
`success` only when the case has an injected bug and the verifier returns
`FAIL`; all other cases are marked `unsuccess`.

Covered mutation families:

- `prepare_coefficients`: perturb PREPARE state-preparation angles.
- `select_pauli`: replace a SELECT-controlled Pauli gate.
- `qsp_phase`: perturb QSP `rz` phase angles.
- `oracle_order`: move a QSP signal oracle across a phase layer.
- `normalization`: use a wrong selector support/normalization or final selector
  extraction.

RQ1 uses paper-mode block-encoding verification. RQ2 uses approximation-mode
Hamiltonian-simulation QSP verification.
RQ1 PREPARE angles use the same 60-digit generation precision as the main RQ1
benchmarks, including the mutated PREPARE-angle cases. RQ2 uses the same QSP
QASM angle formatting and 17-digit polynomial metadata precision as the main
RQ2 generator.
The default RQ2 bug-detection suite follows the paper-style fixed-epsilon
degree sweep: `epsilon = 1e-4`, with `t = 0.1, 0.5, 1.0, 2.0`, giving QSP
degrees `3, 5, 7, 9`.
