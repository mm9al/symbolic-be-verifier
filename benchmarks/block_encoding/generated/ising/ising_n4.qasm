OPENQASM 2.0;
include "qelib1.inc";

qreg q[8];

// RQ1 block-encoding benchmark: ising_n4
// q[0]..q[2] = selector ancillas
// q[3] = work ancilla
// q[4]..q[7] = system qubits
// top-left block = H / alpha, alpha = 7

// PREP: uniform state over Hamiltonian terms
ry(1.4274487578895312) q[0];
x q[0];
ry(0.78539816339744839) q[1];
cx q[0], q[1];
ry(-0.78539816339744839) q[1];
cx q[0], q[1];
x q[0];
ry(0.61547970867038726) q[1];
cx q[0], q[1];
ry(-0.61547970867038726) q[1];
cx q[0], q[1];
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
// term 0: +Z0 Z1
x q[0];
x q[1];
x q[2];
mcx q[0], q[1], q[2], q[3];
cz q[3], q[4];
cz q[3], q[5];
mcx q[0], q[1], q[2], q[3];
x q[0];
x q[1];
x q[2];
// term 1: +Z1 Z2
x q[0];
x q[1];
mcx q[0], q[1], q[2], q[3];
cz q[3], q[5];
cz q[3], q[6];
mcx q[0], q[1], q[2], q[3];
x q[0];
x q[1];
// term 2: +Z2 Z3
x q[0];
x q[2];
mcx q[0], q[1], q[2], q[3];
cz q[3], q[6];
cz q[3], q[7];
mcx q[0], q[1], q[2], q[3];
x q[0];
x q[2];
// term 3: +X0
x q[0];
mcx q[0], q[1], q[2], q[3];
cx q[3], q[4];
mcx q[0], q[1], q[2], q[3];
x q[0];
// term 4: +X1
x q[1];
x q[2];
mcx q[0], q[1], q[2], q[3];
cx q[3], q[5];
mcx q[0], q[1], q[2], q[3];
x q[1];
x q[2];
// term 5: +X2
x q[1];
mcx q[0], q[1], q[2], q[3];
cx q[3], q[6];
mcx q[0], q[1], q[2], q[3];
x q[1];
// term 6: +X3
x q[2];
mcx q[0], q[1], q[2], q[3];
cx q[3], q[7];
mcx q[0], q[1], q[2], q[3];
x q[2];

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
cx q[0], q[1];
ry(0.61547970867038726) q[1];
cx q[0], q[1];
ry(-0.61547970867038726) q[1];
x q[0];
cx q[0], q[1];
ry(0.78539816339744839) q[1];
cx q[0], q[1];
ry(-0.78539816339744839) q[1];
x q[0];
ry(-1.4274487578895312) q[0];
