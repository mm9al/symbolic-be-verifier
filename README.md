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

## Current Goal

Support a small OpenQASM subset:

- h
- s
- sdg
- x
- z
- rx
- ry
- rz
- cx
- cz

Rotation gates support symbolic and exact angles such as `theta`, `-theta`,
`pi`, `pi/2`, `0.5*pi`, and `2*pi`.

## Usage

Print only the final symbolic branches:

```bash
.venv/bin/python -m symbolic.verify examples/lcu_x_plus_z.qasm --expected "(X + Z)/2"
```

For multiple ancilla and system qubits, pass QASM indices in the order they
should appear in branch keys and Pauli strings. Quote `q[...]` arguments in
shells such as zsh, or pass plain integers.

```bash
.venv/bin/python -m symbolic.verify examples/two_ancilla_hh.qasm \
  --ancillas 'q[0]' 'q[1]' \
  --systems 'q[2]' \
  --expected "I/2" \
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
  --expected "(X0 X1 + Z0 Z1)/2"
```

Rotation-gate example with a symbolic angle:

```bash
.venv/bin/python -m symbolic.verify examples/rz_ancilla_theta_sandwich.qasm \
  --expected "cos(theta/2) I" \
  --trace
```

Print every intermediate update:

```bash
.venv/bin/python -m symbolic.verify examples/lcu_x_plus_z.qasm --expected "(X + Z)/2" --trace
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

and verifies the final all-zero branch against the expected top-left block.

Verification results:

- `PASS`: final all-zero branch exactly equals the expected operator.
- `PASS_UP_TO_GLOBAL_PHASE`: final all-zero branch is `phase * expected` and `|phase| = 1`.
- `FAIL`: final all-zero branch is neither exactly equal nor equal up to a global phase.

## Hamiltonian Simulation QSP Data

QSP generation lives in `symbolic.qsp`; `tools/hamsim_qsp.py` is a command-line
wrapper around that module. Use the tool to generate selector-wrapped
Hamiltonian-simulation cosine and sine QASM examples. Synthesis needs pyqsp's
scientific Python dependencies (`numpy`, `scipy`, `matplotlib`), so run it with
an environment that has pyqsp available, not the minimal verifier `.venv`.

```bash
python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --write-examples
```

This writes both selector examples into one folder, plus polynomial metadata:

```text
examples/qsp_hamsim_t05_eps1e-4/
  qsp_hamsim_cos_selector_t05_eps1e-4_deg4.qasm
  qsp_hamsim_sin_selector_t05_eps1e-4_deg5.qasm
  expected_polynomials.json
```

The original pyqsp component circuits place the target polynomial in the
imaginary part of the top-left block. The selector examples add `q[3]`, run the
`theta` and `-theta` branches coherently, extract their difference, multiply by
`-i`, and move the result to the all-zero branch. The final block is therefore
directly `P_cos(H)` or `P_sin(H)`, so the verifier compares the full polynomial.

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

```bash
.venv/bin/python -m symbolic.verify examples/qsp_hamsim_t05_eps1e-4/qsp_hamsim_cos_selector_t05_eps1e-4_deg4.qasm \
  --ancillas 'q[3]' 'q[0]' 'q[1]' \
  --systems 'q[2]' \
  --expected-polynomial '0.49999966355545339 - 0.062493938728279574*x^2 + 0.0012858918109143005*x^4' \
  --hermitian-base \
  --compare-polynomial-only
```

```bash
.venv/bin/python -m symbolic.verify examples/qsp_hamsim_t05_eps1e-4/qsp_hamsim_sin_selector_t05_eps1e-4_deg5.qasm \
  --ancillas 'q[3]' 'q[0]' 'q[1]' \
  --systems 'q[2]' \
  --expected-polynomial '0.24999991579484243*x - 0.010415992523176125*x^3 + 0.00012885803586171963*x^5' \
  --hermitian-base \
  --compare-polynomial-only
```

The checked-in `examples/qsp_hamsim_t05_eps1e-4` fixture is the regression case
for one block ancilla:

```text
phase ancilla     q[0]
block ancilla     q[1]
system qubit      q[2]
selector ancilla  q[3]
```

Run the regression with:

```bash
.venv/bin/python -m pytest tests/test_qsp.py
```

This verifies that both checked-in selector circuits still synthesize and
verify as:

```text
m = 1 block ancilla -> qsp cos/sin synthesis PASS
```

The checked-in `examples/qsp_hamsim_t05_eps1e-4_m2` fixture is the first
multi-block-ancilla cos/sin regression. It uses the fixed multi-ancilla QSP
layout:

```text
q[0]       = QSP phase ancilla
q[1..m]    = block-encoding ancillas of U_H
q[m+1]     = selector / realification ancilla
q[m+2..]   = system qubits
```

For `m = 2`, this is:

```text
phase ancilla     q[0]
block ancillas    q[1], q[2]
selector ancilla  q[3]
system qubit      q[4]
```

Run the multi-ancilla cos/sin regression with:

```bash
.venv/bin/python -m pytest tests/test_qsp.py::test_multi_block_ancilla_qsp_cos_sin_regression_passes
```

Generate fresh multi-ancilla selector examples with pyqsp available:

```bash
python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --write-examples \
  --block-ancillas 'q[1]' 'q[2]' \
  --selector-qubit 'q[3]' \
  --system-qubits 'q[4]'
```

The current multi-ancilla cos/sin verifier path assumes a real block encoding,
so `U_H = U_H^*`. Gates that require complex conjugation behavior, such as
system-level `Y`, `S`, or `Sdg` inside `U_H`, are intentionally out of scope for
this checkpoint.

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
```

In word mode, the verifier tracks the all-zero block subspace and abstracts the
orthogonal block-ancilla subspace as one complement. This is enough for the
first multi-ancilla QSP regression before adding a physical decomposition pass.

Run the abstract QSP regressions with:

```bash
.venv/bin/python -m pytest tests/test_qsp.py::test_qsp_mcx_m1_regression_matches_t3 \
  tests/test_qsp.py::test_qsp_mcx_m2_t3_passes
```

For debugging the raw component data without writing files:

```bash
python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --component both --format json
python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --component cos-selector --qasm-snippet
python3 tools/hamsim_qsp.py --tau 0.5 --epsilon 1e-4 --component sin-selector --qasm-snippet
```
