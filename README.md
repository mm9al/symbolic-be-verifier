# Symbolic Block-Encoding Verifier

This project is a prototype symbolic verifier for block-encoding circuits.

Given an OpenQASM 2.0 circuit and specified ancilla qubits, the verifier tracks the intermediate state in the form

$$
|\Phi\rangle = \sum_b |b\rangle_a B_b |\psi\rangle.
$$

Here, each $B_b$ is a symbolic operator acting on the system qubits.
In v0.4, the symbolic operators are represented in the Pauli-string basis, for
example `X_0`, `Z_1`, and `X_0 Y_1 Z_2`.
Scalar coefficients are SymPy expressions, so symbolic angles such as `theta`
and exact angles such as `pi/2` are preserved and simplified separately from
Pauli-string multiplication.

The goal is to verify the top-left block of a block-encoding circuit by checking the final all-zero branch $B_{00\cdots0}$.

## Reproducible Environment

The full verifier, benchmark generators, evaluation runners, and test suite are
tested with Python 3.11. Create one environment containing both the core
verifier and the QSP/Hamiltonian-simulation dependencies:

```bash
python3.11 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
```

Verify the installation before reproducing evaluation results:

```bash
.venv/bin/python -m pytest -q
```

The dependency versions in `requirements.txt` are pinned to the tested
artifact environment. Python 3.9 is not supported because the required
`python-flint==0.9.0` wheels are unavailable for that interpreter.

## Current Goal

Support a small OpenQASM subset:

- h
- s
- sdg
- t
- tdg
- x
- z
- rx
- ry
- rz
- cx
- cz
- cy

Rotation gates support symbolic and exact angles such as `theta`, `-theta`,
`pi`, `pi/2`, `0.5*pi`, and `2*pi`.

`t` and `tdg` use exact algebraic coefficients, including their absolute phase,
on either ancillas or systems. `cy` supports distinct known control/target
qubits in either register, including both-system operands in concrete operator
mode. QSP word mode still rejects explicit system operations and retains its
canonical input checks.

## Single-Circuit Block-Encoding Verification

Use the paper-style proportionality check for a single block-encoding circuit.
Pass the unnormalized target Hamiltonian `H` to `--expected`; the verifier
infers the real positive subnormalization factor `alpha` and checks the
coefficient residual against `epsilon / 2^(n/2)`:

```bash
.venv/bin/python -m symbolic.verify examples/lcu_x_minus_z.qasm \
  --expected "X - Z" \
  --block-encoding-epsilon 1e-8
```

## RQ1 Block-Encoding Evaluation

The checked-in RQ1 benchmark suite targets block-encoding verification for
three Hamiltonian families:

- periodic transverse-field Ising, with `2n` Pauli terms
- MaxCut on the deterministic cycle graph `C_n`, with `2n` LCU terms
- open-chain isotropic Heisenberg, with `3(n - 1)` Pauli terms including `Y`

Generate the default consecutive-size benchmark suite, `n = 2..128`:

```bash
.venv/bin/python evaluation/generate_benchmarks.py
```

Run the symbolic block-encoding evaluation:

```bash
.venv/bin/python evaluation/run_block_encoding.py
```

The block-encoding evaluation uses the paper's final proportionality check by
default: it treats the manifest Hamiltonian as `H`, infers the real positive
subnormalization factor `alpha`, and checks the coefficient residual with
`epsilon = 1e-8` against the paper threshold `epsilon / 2^(n/2)`. This is the
only approximate acceptance threshold. Pass `--check-mode exact` only to
reproduce the older normalized top-left-block comparison.

Results are written to `evaluation/results/block_encoding_results.csv`. For a
quick smoke suite, pass a smaller size list such as:

```bash
.venv/bin/python evaluation/generate_benchmarks.py --sizes 2..8
```

For multiple ancilla and system qubits, pass QASM indices in the order they
should appear in branch keys and Pauli strings. Quote `q[...]` arguments in
shells such as zsh, or pass plain integers.

```bash
.venv/bin/python -m symbolic.verify examples/two_ancilla_hh.qasm \
  --ancillas 'q[0]' 'q[1]' \
  --systems 'q[2]' \
  --expected "I" \
  --block-encoding-epsilon 1e-8 \
  --trace
```

Bitstring convention:

```text
--ancillas q[0] q[1] q[2]
```

corresponds to branch key `(b0, b1, b2)`, where `b0` is the bit of `q[0]`,
`b1` is the bit of `q[1]`, and `b2` is the bit of `q[2]`. The top-left block is
always the all-zero branch, e.g. `B[000]`.

For multiple system qubits, `--systems q[2] q[3]` means `q[2] -> system[0]`
and `q[3] -> system[1]`, so `X_0 Z_1` means `X` on `q[2]` and `Z` on `q[3]`.

```bash
.venv/bin/python -m symbolic.verify examples/lcu_xx_plus_zz.qasm \
  --systems 1 2 \
  --expected "X0 X1 + Z0 Z1" \
  --block-encoding-epsilon 1e-8
```

The legacy exact-mode symbolic-angle command is retained below only as a
commented reference. It checks the normalized top-left block directly and is
not part of the paper-style RQ1 evaluation:

```bash
# Legacy exact-mode command (not the paper RQ1 check):
# .venv/bin/python -m symbolic.verify examples/rz_ancilla_theta_sandwich.qasm \
#   --expected "cos(theta/2) I" \
#   --trace
```

Print every intermediate update:

```bash
.venv/bin/python -m symbolic.verify examples/lcu_x_minus_z.qasm \
  --expected "X - Z" \
  --block-encoding-epsilon 1e-8 \
  --trace
```

For multi-ancilla circuits, trace output is sparse and branch-oriented:

```text
0 Initial
  B[00] = I
1 h q[0]
  B[00] = sqrt(2)/2*I
  B[10] = sqrt(2)/2*I
```

The verifier starts from

$$B_{00\cdots0} = I$$

In the paper-style block-encoding route, it treats `--expected` as the target
Hamiltonian `H` and verifies the final all-zero branch with the proportionality
and residual-threshold check described above. Direct comparison against an
expected normalized top-left block belongs only to the legacy exact mode.

Verification results:

- `PASS`: with `--block-encoding-epsilon`, the final all-zero branch is
  proportional to the target Hamiltonian within the requested tolerance, using
  the inferred real positive subnormalization factor `alpha`.
- `PASS_EXACT`: in the legacy exact mode, the final all-zero branch exactly
  equals the expected normalized operator.
- `PASS_UP_TO_SCALE`: in the legacy exact mode, the final all-zero branch is
  `scale * expected` for a scalar `scale`.
- `FAIL`: final all-zero branch is neither exactly equal nor equal up to a scalar.

Numeric coefficient comparisons use a default absolute tolerance of `1e-8`.
Override it with `--tolerance` when running `symbolic.verify`. Decimal QASM
parameters and every inexact coefficient derived from them are propagated with
python-flint Arb/Acb ball arithmetic (at least 100 decimal digits by default,
automatically raised to two digits per system qubit for large circuits). A numeric
comparison passes only when the rigorous upper bound on the residual is no
larger than the rigorous lower bound on the acceptance threshold. If the two
intervals overlap, the verifier conservatively returns `FAIL`; it never accepts
from the midpoint alone.

Ball arithmetic certifies rounding error in the implemented computation. It
does not remove modelling error, finite precision already present in generated
QASM angle literals, or the approximation/discretization assumptions of a
verification route.

## Verification Routes

The command-line verifier supports three separate routes:

- Block encoding: compare the final all-zero branch to an operator passed with
  `--expected`. With `--block-encoding-epsilon`, treat `--expected` as the
  target Hamiltonian `H`, infer the real positive `alpha` minimizing
  `||alpha*B - H||`, and accept when the coefficient residual is at most
  `epsilon / 2^(n/2)`.
- Hamiltonian simulation: compute the final symbolic polynomial and check it
  directly against `scale * exp(-i*x*tau)` on `[-1, 1]` with `--hamsim-tau`,
  `--hamsim-epsilon`, and `--hamsim-scale`.
- General QSP: compute the final symbolic polynomial and compare it to a target
  polynomial passed with `--expected-polynomial`.

## General QSP Verification

For QSP word-mode circuits, pass the QSP phase ancilla, block-encoding
ancillas, and system qubits explicitly. The examples use this layout:

```text
q[0]       = QSP phase ancilla
q[1..m]    = block-encoding ancillas of U_H
q[m+1]     = selector / realification ancilla
q[m+2..]   = system qubits
```

QSP word mode accepts QASM in the QSP canonical form used by the paper. Every
`UH`, `UHdg`, `cUH`, and `cUHdg` invocation must use the same block-ancilla and
system registers. With two or more block ancillas, a phase conditioned on the
all-zero block state must use the complete open-control sandwich: apply `x` to
every block ancilla, use all of them as controls, and restore every `x` in
reverse order. The verifier rejects direct gates that distinguish individual
nonzero block states and controls that cover only part of the block register.

Canonical-form validation is enabled by default. For experimental imports,
pass `--skip-qsp-canonical-check`, or call
`verify_qasm_file(..., check_qsp_canonical_form=False)` in Python. This only
skips input validation; gate support and the zero/nonzero block abstraction
are unchanged. For multiple block ancillas, a result on noncanonical input
requires a separate justification that the abstraction represents that circuit.

The QASM files still keep pyqsp phases separate from physical OpenQASM rotation
angles. They first convert pyqsp's Wx phases to QSVT projector phases, then
write the physical `rz(theta_rz)` angle:

```text
psi_0 = phi_0 + pi/4
psi_j = phi_j + pi/2, 1 <= j <= d - 1
psi_d = phi_d + pi/4 for even d, and phi_d - pi/4 for odd d
theta_rz = -2 * psi
```

Use `--compare-polynomial-only` to compare polynomial coefficients directly
instead of evaluating the polynomial on `--base`.

For QSP polynomial checks, the verifier eagerly applies the oriented
block-unitarity rewrite rules after every QSP signal step, preventing unreduced
word expressions from accumulating across the circuit. It normalizes the
all-zero branch once more at the end. Pass `--defer-qsp-rewrite` only when an
unreduced intermediate trace is needed. With `--hermitian-base`, the verifier
also rewrites `Hd` to `H`. The normalized expression must then be a polynomial
in `H` only, i.e. a linear combination of `I`, `H`, `H H`, `H H H`, and so on.
If any `A`, `G`, `C`, adjoint garbage block, or unreduced `Hd` remains, the
check fails before comparing coefficients with the target polynomial.

QSP verification uses Arb/Acb balls with a 100-decimal-digit minimum for parsing
physical QASM rotation angles, trigonometric coefficients, coefficient
propagation, and the final Hamiltonian-simulation approximation grid. This
retains a rounding-error enclosure through high-degree coefficient propagation
and grid evaluation. Override the minimum with
`--qsp-calculation-precision` when deliberately testing another precision.
Generated physical QSP rotations contain 17 significant digits, matching the
meaningful precision available from pyqsp's binary64 phase synthesis; the
ball verifier preserves those input values without inventing extra phase
digits.

Ordinary QSP verification example:

```bash
.venv/bin/python -m symbolic.verify examples/qsp_t3_opaque.qasm \
  --ancillas 'q[0]' 'q[1]' \
  --systems 'q[2]' \
  --base '(X - Z)/2' \
  --expected-polynomial '4*x^3 - 3*x' \
  --hermitian-base \
  --result-only
```

Polynomial-only check:

```bash
.venv/bin/python -m symbolic.verify examples/qsp_t3_opaque.qasm \
  --ancillas 'q[0]' 'q[1]' \
  --systems 'q[2]' \
  --expected-polynomial '4*x^3 - 3*x' \
  --hermitian-base \
  --compare-polynomial-only \
  --result-only
```

The verifier treats `UHdg` as an adjoint block, not as elementwise complex
conjugation. With `--hermitian-base`, word-mode normalization rewrites `Hd` to
`H`, so the final syntactic polynomial check only relies on `H^\dagger = H`; it
does not require a real block encoding with `H^* = H`.

## Hamiltonian Simulation Verification

QSP generation lives in `symbolic.qsp`; `tools/hamsim_qsp.py` is a command-line
wrapper around that module. Its pyqsp scientific Python dependencies (`numpy`,
`scipy`, and `matplotlib`) are included in the reproducible `.venv` setup
above.

Generate a full Hamiltonian-simulation QASM file and metadata:

```bash
python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --component full --write-examples \
  --selector-qubit 'q[0]' \
  --component-selector-qubit 'q[1]' \
  --phase-qubit 'q[2]' \
  --block-ancilla 'q[3]' \
  --system-qubit 'q[4]'
```

For a block encoding whose `U_H` has multiple block ancillas, pass the whole
block register with `--block-ancillas`:

```bash
python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --component full --write-examples \
  --selector-qubit 'q[0]' \
  --component-selector-qubit 'q[1]' \
  --phase-qubit 'q[2]' \
  --block-ancillas 'q[3]' 'q[4]' \
  --system-qubit 'q[5]'
```

This writes:

```text
examples/qsp_hamsim_full_t05_eps1e-4/
  qsp_hamsim_full_t05_eps1e-4_deg5.qasm
  expected_polynomial.json
```

The full block shares the cosine/sine QSP signal skeleton and multiplexes the
signed phase gadgets. Relative to the generated component polynomials, the
all-zero selector branch is:

```text
1/2 * (P_cos(H) - i P_sin(H))
```

For the checked-in pyqsp Hamiltonian-simulation examples, `P_cos` and `P_sin`
are already half-scaled. Therefore the final expected polynomial is
approximately:

```text
1/4 * exp(-i H tau)
```

The hamsim route checks the final symbolic polynomial numerically against
`exp(-i*x*tau)` on `[-1, 1]`. The verifier first chooses the scalar
`beta = <p,e> / ||p||^2` on the grid, where `e_m = exp(-i*x_m*tau)`, then uses
the paper grid
`D = ceil(e*|tau|/2 + log(32/epsilon))` and
`M = 4*max(D, degree)` with Chebyshev roots
`x_m = cos((2m - 1)*pi/(2M))` for `m = 1..M`. Arb certifies `D`, the roots,
the fitted `beta`, and each complex residual. It accepts only when the maximum
rigorous residual upper bound is at most the lower bound of `epsilon/2`.

The checked-in degree-3 full Hamiltonian simulation fixture uses:

```text
q[0] = Hamiltonian-simulation selector ancilla
q[1] = QSP imaginary-part extraction selector ancilla
q[2] = QSP phase ancilla
q[3] = block-encoding ancilla
q[4] = system qubit
```

Verify the degree-3 fixture with:

```bash
.venv/bin/python -m symbolic.verify \
  examples/qsp_hamsim_full_t05_eps01/qsp_hamsim_full_t05_eps01_deg3.qasm \
  --ancillas 'q[0]' 'q[1]' 'q[2]' 'q[3]' \
  --systems 'q[4]' \
  --hamsim-tau 0.5 \
  --hamsim-epsilon 0.1 \
  --hermitian-base \
  --result-only
```

Verify the degree-5 fixture with:

```bash
.venv/bin/python -m symbolic.verify \
  examples/qsp_hamsim_full_t05_eps1e-4/qsp_hamsim_full_t05_eps1e-4_deg5.qasm \
  --ancillas 'q[0]' 'q[1]' 'q[2]' 'q[3]' \
  --systems 'q[4]' \
  --hamsim-tau 0.5 \
  --hamsim-epsilon 1e-4 \
  --hermitian-base \
  --result-only
```

Use `--result-only` when you only want the verifier status instead of the full
symbolic branches.

```bash
time .venv/bin/python -m symbolic.verify \
  examples/qsp_hamsim_full_t05_eps01/qsp_hamsim_full_t05_eps01_deg3.qasm \
  --ancillas 'q[0]' 'q[1]' 'q[2]' 'q[3]' \
  --systems 'q[4]' \
  --hamsim-tau 0.5 \
  --hamsim-epsilon 0.1 \
  --hermitian-base \
  --result-only
```

## Abstract Multi-Control QSP

The verifier supports a verifier-only abstract multi-controlled X gate:

```qasm
mcx q[control_0], q[control_1], ..., q[target];
```

The final operand is the target ancilla. All preceding operands are control
ancillas. This gate is not decomposed into OpenQASM basis gates yet; it directly
permutes symbolic branches by flipping the target exactly when all controls are
`1`.

Multi-block-ancilla QSP keeps the same circuit shape as
`examples/qsp_t3_opaque.qasm`. The zero-control phase sandwich is written with
explicit `x` gates around `mcx`:

```qasm
x q[block_0];
x q[block_1];
mcx q[block_0], q[block_1], q[phase];
rz(theta) q[phase];
mcx q[block_0], q[block_1], q[phase];
x q[block_1];
x q[block_0];
```

This implements the QSP phase split:

```text
block register = |0...0>  -> Rz(-theta) on the phase ancilla
otherwise                 -> Rz(theta) on the phase ancilla
```

Multi-ancilla `UH`/`UHdg` gates are accepted as abstract block encodings:

```qasm
UH q[block_0], q[block_1], ..., q[system_0], ...;
UHdg q[block_0], q[block_1], ..., q[system_0], ...;
cUH q[selector], q[block_0], q[block_1], ..., q[system_0], ...;
cUHdg q[selector], q[block_0], q[block_1], ..., q[system_0], ...;
```

In word mode, the verifier tracks the all-zero block subspace and abstracts the
orthogonal block-ancilla subspace as one complement. This abstraction is sound
for the canonical form above because the circuit never distinguishes states
inside that complement. The verifier checks this restriction before symbolic
execution when `m >= 2`. For `m = 1`, the complement contains only the single
state `|1>`, so ordinary branch operations remain exact.

Run the abstract QSP regressions with:

```bash
.venv/bin/python -m pytest tests/test_qsp.py::test_qsp_mcx_t3_passes_for_block_ancillas
```
