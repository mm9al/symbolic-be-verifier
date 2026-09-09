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

RQ1 uses the paper's final proportionality check by default. For each row, the
runner multiplies the manifest's normalized expected block by the declared
`alpha` to recover the target Hamiltonian `H`, infers the best real positive
`alpha`, and accepts when the coefficient residual is at most
`epsilon / 2^(n/2)`. The default `--block-encoding-epsilon` is `1e-8`, and
the optional `--block-encoding-residual-tolerance` numerical floor defaults to
`0` so the default run follows the paper threshold exactly. Use
`--check-mode exact` only for the older exact comparison against the normalized
top-left block.

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

The RQ1 CSV includes `check_mode`, `check_epsilon`, `inferred_alpha`,
`residual_norm`, `residual_threshold`, `residual_numerical_tolerance`, and
`residual_acceptance_threshold` so paper-mode results record the final
proportionality check data explicitly.

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

By default, the runner uses the paper-style approximation check on the fixed
epsilon axis:

```bash
.venv/bin/python evaluation/run_hamsim.py --axes vary_t --check-mode approximation
```

The approximation grid uses
`D = ceil(e*|t|/2 + log(32/epsilon))` and `M = 4*max(D, degree)`, then
evaluates the Chebyshev roots `x_m = cos((2m - 1)*pi/(2M))` for `m = 1..M`.
The approximation check computes `beta = <B,e> / ||B||^2` on those grid values
and reports the error for `beta*B(x) - exp(-i*x*t)`. The `approx_beta` CSV
field records the fitted scale.

To compare the generated all-zero branch against the exact polynomial in the
manifest instead, use:

```bash
.venv/bin/python evaluation/run_hamsim.py --check-mode polynomial
```

To run the exact polynomial check and also fill `approx_max_grid_error`,
`approx_worst_x`, `approx_grid_points`, and `approx_beta`, use:

```bash
.venv/bin/python evaluation/run_hamsim.py --axes vary_t --check-mode both
```

Pass `--axes vary_epsilon` explicitly if you want the epsilon axis. Be careful
with approximation mode there: small epsilon values imply very fine grids. By
default the runner fails fast above 10,000,000 points; pass
`--max-approx-grid-points 0` only when you intentionally want to disable that
guard. In `polynomial` mode, the `approx_*` fields are blank because no grid
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

The generated RQ2 suite uses `t = 0.1, 0.5, 1.0, 2.0, 4.0`,
`epsilon = 1e-1, 1e-4, 1e-6, 1e-10, 1e-12`, and `m = 1`.
For `t = 0.5`, the epsilon axis gives full QSP degrees `3, 5, 7, 9, 11`.
The degree columns record the generated polynomial degree; on the `vary_t`
axis, epsilon and `m` stay fixed even though the polynomial degree changes with
the requested simulation time.
