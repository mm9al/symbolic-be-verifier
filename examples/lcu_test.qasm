// Generated from Cirq v1.4.0

OPENQASM 2.0;
include "qelib1.inc";


// Qubits: [selection, target]
qreg q[2];


ry(pi*0.5) q[0];
rz(pi*1.5) q[1];
rx(pi*1.0) q[0];
ry(pi*-0.5) q[1];
x q[0];
s q[0];
cz q[0],q[1];
ry(pi*0.5) q[1];
rz(pi*1.0) q[1];
rx(pi*0.5) q[1];
s q[1];
ry(pi*-0.5) q[1];
cz q[0],q[1];
ry(pi*0.5) q[1];
x q[0];
ry(pi*-0.5) q[1];
s q[0];
z q[1];
rz(pi*1.5) q[1];
ry(pi*-0.5) q[1];
cz q[0],q[1];
ry(pi*0.5) q[1];
rz(pi*1.5) q[1];
ry(pi*-0.5) q[1];
cz q[0],q[1];
ry(pi*0.5) q[1];
ry(pi*0.25) q[0];
z q[1];
rx(pi*-1.0) q[0];
ry(pi*-0.25) q[0];
