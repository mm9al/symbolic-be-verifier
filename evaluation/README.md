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

RQ2 evaluates full Hamiltonian-simulation QSP verification along two axes:

- vary simulation time `t`, with `epsilon = 1e-4` and one `U_H` ancilla
- vary target `epsilon`, with `t = 0.5` and one `U_H` ancilla

Generate the QASM benchmark suite and manifest:

```bash
.venv/bin/python evaluation/generate_hamsim_benchmarks.py
```

Run symbolic Hamiltonian-simulation QSP verification:

```bash
.venv/bin/python evaluation/run_hamsim.py
```

The default check compares the generated all-zero branch against the exact
polynomial in the manifest. To check the polynomial numerically against
`target_scale * exp(-i*x*t)` instead on a fixed-epsilon axis, use:

```bash
.venv/bin/python evaluation/run_hamsim.py --axes vary_t --check-mode approximation
```

To keep the exact polynomial check as the status gate while also filling
`approx_max_grid_error`, `approx_worst_x`, and `approx_grid_points`, use:

```bash
.venv/bin/python evaluation/run_hamsim.py --axes vary_t --check-mode both
```

Avoid approximation mode for the default `vary_epsilon` axis unless you
intentionally want the very fine numerical grid implied by tiny epsilon values.
The approximation grid uses Chebyshev nodes `x_m = cos(m*pi/M)` with
`M = ceil(pi * (t + degree) / epsilon)` and evaluates `M + 1` points. The
approximation check first rescales the extracted branch by `1 / target_scale`,
so the reported error is for `B(x) / target_scale - exp(-i*x*t)`. By default
the runner fails fast above 10,000,000 points; pass `--max-approx-grid-points 0`
only when you intentionally want to disable that guard.
In default `polynomial` mode, the `approx_*` fields are blank because no grid
approximation check is run. The `error` field is blank on successful rows and
filled only for exceptions, timeouts, or skipped rows.
The `runtime_sec` field is total wall time. For RQ2 rows, `symbolic_runtime_sec`
measures symbolic execution plus polynomial extraction/comparison, and
`approx_runtime_sec` measures only the final numerical comparison against
`exp(-i*x*t)`.

The result CSV is written to:

```text
evaluation/results/hamsim_results.csv
```

The default RQ2 suite uses `t = 0.1, 0.5, 1.0, 2.0, 4.0`,
`epsilon = 1e-1, 1e-4, 1e-6, 1e-10, 1e-12`, and `m = 1`.
For `t = 0.5`, the epsilon axis gives full QSP degrees `3, 5, 7, 9, 11`.
The degree columns record the generated polynomial degree; on the `vary_t`
axis, epsilon and `m` stay fixed even though the polynomial degree changes with
the requested simulation time.
