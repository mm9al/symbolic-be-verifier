OPENQASM 2.0;
include "qelib1.inc";

opaque UH a, s;
opaque UHdg a, s;

qreg q[3];

// q[0] = phase ancilla
// q[1] = block-encoding ancilla
// q[2] = system qubit

// This is a degree-4 QSP/QSVT-style test.
// The phase angles are written as numerical constants.
// Target polynomial:
// P(x) ~= cos(0.5 x)
//      ~= 0.9999993269 - 0.1249878775 x^2 + 0.0102871345 x^4

// phase phi_0
h q[0];
cz q[0], q[1];
h q[0];
rz(1.5707963268) q[0];
h q[0];
cz q[0], q[1];
h q[0];

// U
UH q[1], q[2];

// phase phi_1
h q[0];
cz q[0], q[1];
h q[0];
rz(-0.2526802551) q[0];
h q[0];
cz q[0], q[1];
h q[0];

// U^\dagger
UHdg q[1], q[2];

// phase phi_2
h q[0];
cz q[0], q[1];
h q[0];
rz(0.5053605103) q[0];
h q[0];
cz q[0], q[1];
h q[0];

// U
UH q[1], q[2];

// phase phi_3
h q[0];
cz q[0], q[1];
h q[0];
rz(-0.2526802551) q[0];
h q[0];
cz q[0], q[1];
h q[0];

// U^\dagger
UHdg q[1], q[2];

// phase phi_4
h q[0];
cz q[0], q[1];
h q[0];
rz(1.5707963268) q[0];
h q[0];
cz q[0], q[1];
h q[0];