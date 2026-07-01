OPENQASM 2.0;
include "qelib1.inc";

qreg q[5];

// RQ1 block-encoding benchmark: heisenberg_n2
// q[0]..q[1] = selector ancillas
// q[2] = work ancilla
// q[3]..q[4] = system qubits
// top-left block = H / alpha, alpha = 3

// PREP: uniform state over Hamiltonian terms
ry(1.2309594173407745) q[0];
x q[0];
ry(0.78539816339744839) q[1];
cx q[0], q[1];
ry(-0.78539816339744839) q[1];
cx q[0], q[1];
x q[0];

// SELECT: branch-conditioned signed Pauli terms
// term 0: +X0 X1
x q[0];
x q[1];
mcx q[0], q[1], q[2];
cx q[2], q[3];
cx q[2], q[4];
mcx q[0], q[1], q[2];
x q[0];
x q[1];
// term 1: +Y0 Y1
x q[0];
mcx q[0], q[1], q[2];
sdg q[2];
cx q[2], q[3];
cz q[2], q[3];
sdg q[2];
cx q[2], q[4];
cz q[2], q[4];
mcx q[0], q[1], q[2];
x q[0];
// term 2: +Z0 Z1
x q[1];
mcx q[0], q[1], q[2];
cz q[2], q[3];
cz q[2], q[4];
mcx q[0], q[1], q[2];
x q[1];

// PREP dagger
x q[0];
cx q[0], q[1];
ry(0.78539816339744839) q[1];
cx q[0], q[1];
ry(-0.78539816339744839) q[1];
x q[0];
ry(-1.2309594173407745) q[0];
