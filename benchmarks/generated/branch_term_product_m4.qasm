OPENQASM 2.0;
include "qelib1.inc";

qreg q[8];

// q[0]..q[3] = ancilla qubits
// q[4]..q[7] = system qubits
// B[0...0] = product_j (I + X_j) / 2

h q[0];
cx q[0], q[4];
h q[0];

h q[1];
cx q[1], q[5];
h q[1];

h q[2];
cx q[2], q[6];
h q[2];

h q[3];
cx q[3], q[7];
h q[3];
