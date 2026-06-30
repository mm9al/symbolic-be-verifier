OPENQASM 2.0;
include "qelib1.inc";

qreg q[2];

// PREPARE
x q[0];
h q[0];

// SELECT X
x q[0];
cx q[0], q[1];
x q[0];

// SELECT Z
cz q[0], q[1];

// PREPARE dagger
h q[0];
