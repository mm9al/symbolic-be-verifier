OPENQASM 2.0;
include "qelib1.inc";

qreg q[7];

// q[0] = ancilla
// q[1]..q[6] = system qubits
// B0 = (I + X0 X1 ... X5) / 2

h q[0];
cx q[0], q[1];
cx q[0], q[2];
cx q[0], q[3];
cx q[0], q[4];
cx q[0], q[5];
cx q[0], q[6];
h q[0];
