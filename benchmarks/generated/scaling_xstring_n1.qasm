OPENQASM 2.0;
include "qelib1.inc";

qreg q[2];

// q[0] = ancilla
// q[1]..q[1] = system qubits
// B0 = (I + X0 X1 ... X0) / 2

h q[0];
cx q[0], q[1];
h q[0];
