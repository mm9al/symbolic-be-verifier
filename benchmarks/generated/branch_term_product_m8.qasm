OPENQASM 2.0;
include "qelib1.inc";

qreg q[16];

// q[0]..q[7] = ancilla qubits
// q[8]..q[15] = system qubits
// B[0...0] = product_j (I + X_j) / 2

h q[0];
cx q[0], q[8];
h q[0];

h q[1];
cx q[1], q[9];
h q[1];

h q[2];
cx q[2], q[10];
h q[2];

h q[3];
cx q[3], q[11];
h q[3];

h q[4];
cx q[4], q[12];
h q[4];

h q[5];
cx q[5], q[13];
h q[5];

h q[6];
cx q[6], q[14];
h q[6];

h q[7];
cx q[7], q[15];
h q[7];
