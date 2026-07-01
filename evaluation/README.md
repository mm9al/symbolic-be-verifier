# Evaluation

RQ1 evaluates block-encoding verification for three Hamiltonian families:

- transverse-field Ising chain: `(n - 1)` nearest-neighbor `ZZ` terms plus `n` single-qubit `X` terms
- MaxCut on the cycle graph `C_n`: one `+I/2` and one `-Z_i Z_j/2` LCU term per edge
- isotropic Heisenberg chain: `XX + YY + ZZ` on each open-chain edge

Generate the QASM benchmark suite and manifest:

```bash
.venv/bin/python evaluation/generate_benchmarks.py
```

Run symbolic block-encoding verification:

```bash
.venv/bin/python evaluation/run_block_encoding.py
```

The result CSV is written to:

```text
evaluation/results/block_encoding_results.csv
```

The default sizes are `n = 4, 8, 16, 32`. Use `--sizes 4,8` on the generator
for a quick smoke run.
