OPENQASM 2.0;
include "qelib1.inc";

opaque UH a, s;
opaque UHdg a, s;

qreg q[3];

UH q[1], q[2];

x q[1];
cx q[1], q[0];
rz(pi) q[0];
cx q[1], q[0];
x q[1];

UHdg q[1], q[2];

x q[1];
cx q[1], q[0];
rz(pi) q[0];
cx q[1], q[0];
x q[1];

UH q[1], q[2];

rz(2*pi) q[0];
