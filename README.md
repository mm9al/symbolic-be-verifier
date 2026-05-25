# Symbolic Block-Encoding Verifier

This project is a prototype symbolic verifier for block-encoding circuits.

Given an OpenQASM 2.0 circuit and a specified ancilla qubit, the verifier tracks the intermediate state in the form

\[
|\Phi\rangle = |0\rangle_a B_0 |\psi\rangle + |1\rangle_a B_1 |\psi\rangle.
\]

Here, \(B_0\) and \(B_1\) are symbolic operators acting on the system qubits.
In v0.3, the symbolic operators are represented in the Pauli-string basis, for
example `X_0`, `Z_1`, and `X_0 Y_1 Z_2`.
Scalar coefficients are SymPy expressions, so symbolic angles such as `theta`
and exact angles such as `pi/2` are preserved and simplified separately from
Pauli-string multiplication.

The goal is to verify the top-left block of a block-encoding circuit by checking the final operator \(B_0\).

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

For multiple system qubits, pass the system QASM indices in the order they
should appear in Pauli strings. For example, `q[1] -> system 0` and
`q[2] -> system 1`:

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

The verifier starts from

\[
B_0 = I,\quad B_1 = 0
\]

and verifies the final `B0` branch against the expected top-left block.

Verification results:

- `PASS`: final `B0` exactly equals the expected operator.
- `PASS_UP_TO_GLOBAL_PHASE`: final `B0 = phase * expected` and `|phase| = 1`.
- `FAIL`: final `B0` is neither exactly equal nor equal up to a global phase.
