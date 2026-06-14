OPENQASM 2.0;
include "qelib1.inc";

qreg q[4];

// q[0]..q[1] = ancilla qubits
// q[2]..q[3] = system qubits
// B[0...0] = product_j (I + X_j) / 2

h q[0];
cx q[0], q[2];
h q[0];

h q[1];
cx q[1], q[3];
h q[1];
