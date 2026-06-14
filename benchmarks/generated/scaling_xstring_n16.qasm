OPENQASM 2.0;
include "qelib1.inc";

qreg q[17];

// q[0] = ancilla
// q[1]..q[16] = system qubits
// B0 = (I + X0 X1 ... X15) / 2

h q[0];
cx q[0], q[1];
cx q[0], q[2];
cx q[0], q[3];
cx q[0], q[4];
cx q[0], q[5];
cx q[0], q[6];
cx q[0], q[7];
cx q[0], q[8];
cx q[0], q[9];
cx q[0], q[10];
cx q[0], q[11];
cx q[0], q[12];
cx q[0], q[13];
cx q[0], q[14];
cx q[0], q[15];
cx q[0], q[16];
h q[0];
