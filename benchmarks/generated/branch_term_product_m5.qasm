OPENQASM 2.0;
include "qelib1.inc";

qreg q[10];

// q[0]..q[4] = ancilla qubits
// q[5]..q[9] = system qubits
// B[0...0] = product_j (I + X_j) / 2

h q[0];
cx q[0], q[5];
h q[0];

h q[1];
cx q[1], q[6];
h q[1];

h q[2];
cx q[2], q[7];
h q[2];

h q[3];
cx q[3], q[8];
h q[3];

h q[4];
cx q[4], q[9];
h q[4];
