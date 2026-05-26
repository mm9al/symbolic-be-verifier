OPENQASM 2.0;
include "qelib1.inc";

qreg q[3];

// q[0], q[1] = ancilla qubits
// q[2] = system qubit

h q[0];
h q[1];
