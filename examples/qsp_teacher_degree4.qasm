OPENQASM 2.0;
include "qelib1.inc";

opaque UH a, s;
opaque UHdg a, s;

qreg q[3];

// q[0] = QSP phase ancilla
// q[1] = block-encoding ancilla of U_H
// q[2] = system qubit
//
// Hand-written degree-4 QSP example with projector phases
// phi = (2.30, -1.05, 1.07, -1.05, 0.73).
//
// Each phase is implemented by a zero-controlled X sandwich on q[0]:
//     x q[1]; cx q[1], q[0]; x q[1];
//     rz(2 * phi_j) q[0];
//     x q[1]; cx q[1], q[0]; x q[1];
// so the block-ancilla-zero branch gets exp(+i phi_j) and the
// block-ancilla-one branch gets exp(-i phi_j).

// phase phi_0 = 2.30
x q[1];
cx q[1], q[0];
x q[1];
rz(4.60) q[0];
x q[1];
cx q[1], q[0];
x q[1];

UH q[1], q[2];

// phase phi_1 = -1.05
x q[1];
cx q[1], q[0];
x q[1];
rz(-2.10) q[0];
x q[1];
cx q[1], q[0];
x q[1];

UHdg q[1], q[2];

// phase phi_2 = 1.07
x q[1];
cx q[1], q[0];
x q[1];
rz(2.14) q[0];
x q[1];
cx q[1], q[0];
x q[1];

UH q[1], q[2];

// phase phi_3 = -1.05
x q[1];
cx q[1], q[0];
x q[1];
rz(-2.10) q[0];
x q[1];
cx q[1], q[0];
x q[1];

UHdg q[1], q[2];

// phase phi_4 = 0.73
x q[1];
cx q[1], q[0];
x q[1];
rz(1.46) q[0];
x q[1];
cx q[1], q[0];
x q[1];
