OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
rx(-pi) q[1];
rx(-pi) q[3];
h q[4];
crz(pi/2) q[3],q[4];
h q[3];
crz(pi/4) q[2],q[4];
crz(pi/2) q[2],q[3];
h q[2];
crz(pi/8) q[1],q[4];
crz(pi/4) q[1],q[3];
crz(pi/2) q[1],q[2];
h q[1];
crz(pi/16) q[0],q[4];
crz(pi/8) q[0],q[3];
crz(pi/4) q[0],q[2];
crz(pi/2) q[0],q[1];
h q[0];
p(pi) q[0];
p(pi/2) q[1];
p(pi/4) q[2];
p(pi/8) q[3];
p(pi/16) q[4];
h q[0];
crz(-pi/2) q[0],q[1];
crz(-pi/4) q[0],q[2];
crz(-pi/8) q[0],q[3];
crz(-pi/16) q[0],q[4];
h q[1];
crz(-pi/2) q[1],q[2];
crz(-pi/4) q[1],q[3];
crz(-pi/8) q[1],q[4];
h q[2];
crz(-pi/2) q[2],q[3];
crz(-pi/4) q[2],q[4];
h q[3];
crz(-pi/2) q[3],q[4];
h q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
