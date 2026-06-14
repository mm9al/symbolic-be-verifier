OPENQASM 2.0;
include "qelib1.inc";

qreg q[9];

// q[0] = ancilla
// q[1]..q[8] = system qubits
// B0 = (I + X0 X1 ... X7) / 2

h q[0];
cx q[0], q[1];
cx q[0], q[2];
cx q[0], q[3];
cx q[0], q[4];
cx q[0], q[5];
cx q[0], q[6];
cx q[0], q[7];
cx q[0], q[8];
h q[0];
