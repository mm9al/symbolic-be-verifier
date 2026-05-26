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
