OPENQASM 2.0;
include "qelib1.inc";

qreg q[2];

// q[0] = ancilla
// q[1] = system qubit

rx(pi) q[1];
