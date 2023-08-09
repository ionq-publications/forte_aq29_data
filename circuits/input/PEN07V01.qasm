OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
ry(-pi/2) q[0];
rz(-pi) q[0];
ry(-pi/2) q[1];
rz(-pi) q[1];
ry(-pi/2) q[2];
rz(-pi) q[2];
ry(-pi/2) q[3];
rz(-pi) q[3];
ry(-pi/2) q[4];
rz(-pi) q[4];
ry(-pi/2) q[5];
rz(-pi) q[5];
rx(-pi) q[6];
p(11*pi) q[1];
p(11*pi/2) q[2];
p(11*pi/4) q[3];
p(11*pi/8) q[4];
p(11*pi/16) q[5];
cx q[5],q[6];
p(-11*pi/16) q[6];
cx q[5],q[6];
p(11*pi/16) q[6];
cx q[4],q[6];
p(-11*pi/8) q[6];
cx q[4],q[6];
p(11*pi/8) q[6];
cx q[3],q[6];
p(-11*pi/4) q[6];
cx q[3],q[6];
p(11*pi/4) q[6];
cx q[2],q[6];
p(-11*pi/2) q[6];
cx q[2],q[6];
p(11*pi/2) q[6];
cx q[1],q[6];
p(-11*pi) q[6];
cx q[1],q[6];
p(11*pi) q[6];
ry(-pi/2) q[0];
rz(-pi) q[0];
p(-pi/4) q[1];
cx q[0],q[1];
p(pi/4) q[1];
cx q[0],q[1];
p(-pi/8) q[2];
cx q[0],q[2];
p(pi/8) q[2];
cx q[0],q[2];
p(-pi/16) q[3];
cx q[0],q[3];
p(pi/16) q[3];
cx q[0],q[3];
p(-pi/32) q[4];
cx q[0],q[4];
p(pi/32) q[4];
cx q[0],q[4];
p(-pi/64) q[5];
cx q[0],q[5];
p(pi/64) q[5];
cx q[0],q[5];
ry(-pi/2) q[1];
rz(-pi) q[1];
p(-pi/4) q[2];
cx q[1],q[2];
p(pi/4) q[2];
cx q[1],q[2];
p(-pi/8) q[3];
cx q[1],q[3];
p(pi/8) q[3];
cx q[1],q[3];
p(-pi/16) q[4];
cx q[1],q[4];
p(pi/16) q[4];
cx q[1],q[4];
p(-pi/32) q[5];
cx q[1],q[5];
p(pi/32) q[5];
cx q[1],q[5];
ry(-pi/2) q[2];
rz(-pi) q[2];
p(-pi/4) q[3];
cx q[2],q[3];
p(pi/4) q[3];
cx q[2],q[3];
p(-pi/8) q[4];
cx q[2],q[4];
p(pi/8) q[4];
cx q[2],q[4];
p(-pi/16) q[5];
cx q[2],q[5];
p(pi/16) q[5];
cx q[2],q[5];
ry(-pi/2) q[3];
rz(-pi) q[3];
p(-pi/4) q[4];
cx q[3],q[4];
p(pi/4) q[4];
cx q[3],q[4];
p(-pi/8) q[5];
cx q[3],q[5];
p(pi/8) q[5];
cx q[3],q[5];
ry(-pi/2) q[4];
rz(-pi) q[4];
p(-pi/4) q[5];
cx q[4],q[5];
p(pi/4) q[5];
cx q[4],q[5];
ry(-pi/2) q[5];
rz(-pi) q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
