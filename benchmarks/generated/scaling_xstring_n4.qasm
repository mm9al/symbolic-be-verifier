OPENQASM 2.0;
include "qelib1.inc";

qreg q[5];

// q[0] = ancilla
// q[1]..q[4] = system qubits
// B0 = (I + X0 X1 ... X3) / 2

h q[0];
cx q[0], q[1];
cx q[0], q[2];
cx q[0], q[3];
cx q[0], q[4];
h q[0];
