OPENQASM 2.0;
include "qelib1.inc";

opaque UH a, s;
opaque UHdg a, s;
opaque mcx c, t;

qreg q[3];

// q[0] = QSP phase ancilla
// q[1] = block-encoding ancilla
// q[2] = system qubit

UH q[1], q[2];

x q[1];
mcx q[1], q[0];
rz(pi) q[0];
mcx q[1], q[0];
x q[1];

UHdg q[1], q[2];

x q[1];
mcx q[1], q[0];
rz(pi) q[0];
mcx q[1], q[0];
x q[1];

UH q[1], q[2];

rz(2*pi) q[0];
