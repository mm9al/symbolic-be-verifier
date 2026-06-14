OPENQASM 2.0;
include "qelib1.inc";

qreg q[3];

// q[0] = ancilla
// q[1], q[2] = system qubits

h q[0];

// SELECT XX on |0>
x q[0];
cx q[0], q[1];
cx q[0], q[2];
x q[0];

// SELECT ZZ on |1>
cz q[0], q[1];
cz q[0], q[2];

h q[0];
