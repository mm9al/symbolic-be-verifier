OPENQASM 2.0;
include "qelib1.inc";

qreg q[2];

// q[0] = ancilla
// q[1] = system

h q[0];

// SELECT X on |0>
x q[0];
cx q[0], q[1];
x q[0];

// SELECT Z on |1>
cz q[0], q[1];

h q[0];
