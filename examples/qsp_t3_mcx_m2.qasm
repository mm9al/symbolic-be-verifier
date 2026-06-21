OPENQASM 2.0;
include "qelib1.inc";

opaque UH a0, a1, s;
opaque UHdg a0, a1, s;
opaque mcx c0, c1, t;

qreg q[4];

// q[0] = QSP phase ancilla
// q[1], q[2] = block-encoding ancillas
// q[3] = system qubit

UH q[1], q[2], q[3];

x q[1];
x q[2];
mcx q[1], q[2], q[0];
rz(pi) q[0];
mcx q[1], q[2], q[0];
x q[2];
x q[1];

UHdg q[1], q[2], q[3];

x q[1];
x q[2];
mcx q[1], q[2], q[0];
rz(pi) q[0];
mcx q[1], q[2], q[0];
x q[2];
x q[1];

UH q[1], q[2], q[3];

rz(2*pi) q[0];
