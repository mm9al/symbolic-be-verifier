# Evaluation

RQ1 evaluates block-encoding verification for three Hamiltonian families:

- periodic transverse-field Ising chain: `n` nearest-neighbor `ZZ` terms plus `n` single-qubit `X` terms
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

Write per-gate profiling CSV files while running selected benchmarks:

```bash
.venv/bin/python evaluation/run_block_encoding.py --models heisenberg --profile-dir evaluation/results/gate_profiles
```

Each profile row records `gate_id`, `gate_name`, nonzero branch count, total
operator terms, max terms per branch, total gate time, scalar simplify time,
and Pauli-term combine time.

The result CSV is written to:

```text
evaluation/results/block_encoding_results.csv
```

The default sizes are the consecutive integers `n = 2..128`.
All three families use the same `ry`-based uniform PREP generator, including
cases where the term count happens to be a power of two. Use `--sizes 2..8` on
the generator for a quick smoke run, or pass a custom list such as
`--sizes 2,4,8,16`.
