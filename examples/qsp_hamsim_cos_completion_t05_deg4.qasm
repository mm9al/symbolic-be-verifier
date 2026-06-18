OPENQASM 2.0;
include "qelib1.inc";

opaque UH a, s;
opaque UHdg a, s;

qreg q[3];

// q[0] = phase ancilla
// q[1] = block-encoding ancilla
// q[2] = system qubit

// This is a valid degree-4 QSP-style completion.
// Re P(x) ~= cos(0.5 x).

// phi_0 = -2.971548155103
h q[0];
cz q[0], q[1];
h q[0];
rz(-2.971548155103) q[0];
h q[0];
cz q[0], q[1];
h q[0];

// U
UH q[1], q[2];

// phi_1 = -0.013239571518
h q[0];
cz q[0], q[1];
h q[0];
rz(-0.013239571518) q[0];
h q[0];
cz q[0], q[1];
h q[0];

// U^\dagger
UHdg q[1], q[2];

// phi_2 = 2.462143062415
h q[0];
cz q[0], q[1];
h q[0];
rz(2.462143062415) q[0];
h q[0];
cz q[0], q[1];
h q[0];

// U
UH q[1], q[2];

// phi_3 = 0.498060453743
h q[0];
cz q[0], q[1];
h q[0];
rz(0.498060453743) q[0];
h q[0];
cz q[0], q[1];
h q[0];

// U^\dagger
UHdg q[1], q[2];

// phi_4 = 0.991905457551
h q[0];
cz q[0], q[1];
h q[0];
rz(0.991905457551) q[0];
h q[0];
cz q[0], q[1];
h q[0];