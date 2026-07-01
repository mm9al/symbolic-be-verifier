OPENQASM 2.0;
include "qelib1.inc";

qreg q[5];

// RQ1 block-encoding benchmark: maxcut_cycle_n2
// q[0]..q[1] = selector ancillas
// q[2] = work ancilla
// q[3]..q[4] = system qubits
// top-left block = H / alpha, alpha = 2

// PREP: uniform state over Hamiltonian terms
ry(1.5707963267948968) q[0];
x q[0];
ry(0.78539816339744839) q[1];
cx q[0], q[1];
ry(-0.78539816339744839) q[1];
cx q[0], q[1];
x q[0];
ry(0.78539816339744839) q[1];
cx q[0], q[1];
ry(-0.78539816339744839) q[1];
cx q[0], q[1];

// SELECT: branch-conditioned signed Pauli terms
// term 0: +I
// term 1: -Z0 Z1
x q[0];
mcx q[0], q[1], q[2];
z q[2];
cz q[2], q[3];
cz q[2], q[4];
mcx q[0], q[1], q[2];
x q[0];
// term 2: +I
// term 3: -Z0 Z1
mcx q[0], q[1], q[2];
z q[2];
cz q[2], q[3];
cz q[2], q[4];
mcx q[0], q[1], q[2];

// PREP dagger
cx q[0], q[1];
ry(0.78539816339744839) q[1];
cx q[0], q[1];
ry(-0.78539816339744839) q[1];
x q[0];
cx q[0], q[1];
ry(0.78539816339744839) q[1];
cx q[0], q[1];
ry(-0.78539816339744839) q[1];
x q[0];
ry(-1.5707963267948968) q[0];
