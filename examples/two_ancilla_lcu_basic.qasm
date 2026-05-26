OPENQASM 2.0;
include "qelib1.inc";

qreg q[3];

// q[0], q[1] = ancillas
// q[2] = system[0]

// PREPARE: uniform superposition over 00, 01, 10, 11
h q[0];
h q[1];

// SELECT
// if q[0] = 1, apply X on system
cx q[0], q[2];

// if q[1] = 1, apply Z on system
cz q[1], q[2];

// PREPARE dagger
h q[1];
h q[0];