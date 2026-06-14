OPENQASM 2.0;
include "qelib1.inc";

qreg q[2];

// q[0]..q[0] = ancilla qubits
// q[1]..q[1] = system qubits
// B[0...0] = product_j (I + X_j) / 2

h q[0];
cx q[0], q[1];
h q[0];
