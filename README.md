# Symbolic Block-Encoding Verifier

This project is a prototype symbolic verifier for block-encoding circuits.

Given an OpenQASM 2.0 circuit and a specified ancilla qubit, the verifier tracks the intermediate state in the form

\[
|\Phi\rangle = |0\rangle_a B_0 |\psi\rangle + |1\rangle_a B_1 |\psi\rangle.
\]

Here, \(B_0\) and \(B_1\) are symbolic operators acting on the system qubits.

The goal is to verify the top-left block of a block-encoding circuit by checking the final operator \(B_0\).

## Current Goal

Support a small OpenQASM subset:

- h
- x
- z
- rx
- ry
- rz
- cx
- cz