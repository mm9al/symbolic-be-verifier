OPENQASM 2.0;
include "qelib1.inc";

qreg q[7];

// RQ1 block-encoding benchmark: heisenberg_n3
// q[0]..q[2] = selector ancillas
// q[3] = work ancilla
// q[4]..q[6] = system qubits
// top-left block = H / alpha, alpha = 6

// PREP: uniform state over Hamiltonian terms
ry(1.2309594173407745) q[0];
x q[0];
ry(0.78539816339744839) q[1];
cx q[0], q[1];
ry(-0.78539816339744839) q[1];
cx q[0], q[1];
x q[0];
x q[0];
x q[1];
mcx q[0], q[1], q[3];
ry(0.78539816339744839) q[2];
cx q[3], q[2];
ry(-0.78539816339744839) q[2];
cx q[3], q[2];
mcx q[0], q[1], q[3];
x q[0];
x q[1];
x q[0];
mcx q[0], q[1], q[3];
ry(0.78539816339744839) q[2];
cx q[3], q[2];
ry(-0.78539816339744839) q[2];
cx q[3], q[2];
mcx q[0], q[1], q[3];
x q[0];
x q[1];
mcx q[0], q[1], q[3];
ry(0.78539816339744839) q[2];
cx q[3], q[2];
ry(-0.78539816339744839) q[2];
cx q[3], q[2];
mcx q[0], q[1], q[3];
x q[1];

// SELECT: branch-conditioned signed Pauli terms
// term 0: +X0 X1
x q[0];
x q[1];
x q[2];
mcx q[0], q[1], q[2], q[3];
cx q[3], q[4];
cx q[3], q[5];
mcx q[0], q[1], q[2], q[3];
x q[0];
x q[1];
x q[2];
// term 1: +Y0 Y1
x q[0];
x q[1];
mcx q[0], q[1], q[2], q[3];
sdg q[3];
cx q[3], q[4];
cz q[3], q[4];
sdg q[3];
cx q[3], q[5];
cz q[3], q[5];
mcx q[0], q[1], q[2], q[3];
x q[0];
x q[1];
// term 2: +Z0 Z1
x q[0];
x q[2];
mcx q[0], q[1], q[2], q[3];
cz q[3], q[4];
cz q[3], q[5];
mcx q[0], q[1], q[2], q[3];
x q[0];
x q[2];
// term 3: +X1 X2
x q[0];
mcx q[0], q[1], q[2], q[3];
cx q[3], q[5];
cx q[3], q[6];
mcx q[0], q[1], q[2], q[3];
x q[0];
// term 4: +Y1 Y2
x q[1];
x q[2];
mcx q[0], q[1], q[2], q[3];
sdg q[3];
cx q[3], q[5];
cz q[3], q[5];
sdg q[3];
cx q[3], q[6];
cz q[3], q[6];
mcx q[0], q[1], q[2], q[3];
x q[1];
x q[2];
// term 5: +Z1 Z2
x q[1];
mcx q[0], q[1], q[2], q[3];
cz q[3], q[5];
cz q[3], q[6];
mcx q[0], q[1], q[2], q[3];
x q[1];

// PREP dagger
x q[1];
mcx q[0], q[1], q[3];
cx q[3], q[2];
ry(0.78539816339744839) q[2];
cx q[3], q[2];
ry(-0.78539816339744839) q[2];
mcx q[0], q[1], q[3];
x q[1];
x q[0];
mcx q[0], q[1], q[3];
cx q[3], q[2];
ry(0.78539816339744839) q[2];
cx q[3], q[2];
ry(-0.78539816339744839) q[2];
mcx q[0], q[1], q[3];
x q[0];
x q[1];
x q[0];
mcx q[0], q[1], q[3];
cx q[3], q[2];
ry(0.78539816339744839) q[2];
cx q[3], q[2];
ry(-0.78539816339744839) q[2];
mcx q[0], q[1], q[3];
x q[1];
x q[0];
x q[0];
cx q[0], q[1];
ry(0.78539816339744839) q[1];
cx q[0], q[1];
ry(-0.78539816339744839) q[1];
x q[0];
ry(-1.2309594173407745) q[0];
