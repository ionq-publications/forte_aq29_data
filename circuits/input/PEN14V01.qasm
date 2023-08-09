OPENQASM 2.0;
include "qelib1.inc";
qreg q[14];
creg c[14];
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
ry(-pi/2) q[6];
rz(-pi) q[6];
ry(-pi/2) q[7];
rz(-pi) q[7];
ry(-pi/2) q[8];
rz(-pi) q[8];
ry(-pi/2) q[9];
rz(-pi) q[9];
ry(-pi/2) q[10];
rz(-pi) q[10];
ry(-pi/2) q[11];
rz(-pi) q[11];
ry(-pi/2) q[12];
rz(-pi) q[12];
rx(-pi) q[13];
p(1961*pi) q[0];
p(3080.3315968447923) q[1];
p(6.016272650087485) q[10];
p(3.0081363250437425) q[11];
p(1.5040681625218713) q[12];
cx q[12],q[13];
p(-1.5040681625218713) q[13];
cx q[12],q[13];
p(1.5040681625218713) q[13];
cx q[11],q[13];
p(-3.0081363250437425) q[13];
cx q[11],q[13];
p(3.0081363250437425) q[13];
cx q[10],q[13];
p(-6.016272650087485) q[13];
cx q[10],q[13];
p(6.016272650087485) q[13];
p(1540.1657984223962) q[2];
p(770.0828992111981) q[3];
p(385.04144960559904) q[4];
p(192.52072480279952) q[5];
p(96.26036240139976) q[6];
p(48.13018120069988) q[7];
p(24.06509060034994) q[8];
p(12.03254530017497) q[9];
cx q[9],q[13];
p(-12.03254530017497) q[13];
cx q[9],q[13];
p(12.03254530017497) q[13];
cx q[8],q[13];
p(-24.06509060034994) q[13];
cx q[8],q[13];
p(24.06509060034994) q[13];
cx q[7],q[13];
p(-48.13018120069988) q[13];
cx q[7],q[13];
p(48.13018120069988) q[13];
cx q[6],q[13];
p(-96.26036240139976) q[13];
cx q[6],q[13];
p(96.26036240139976) q[13];
cx q[5],q[13];
p(-192.52072480279952) q[13];
cx q[5],q[13];
p(192.52072480279952) q[13];
cx q[4],q[13];
p(-385.04144960559904) q[13];
cx q[4],q[13];
p(385.04144960559904) q[13];
cx q[3],q[13];
p(-770.0828992111981) q[13];
cx q[3],q[13];
p(770.0828992111981) q[13];
cx q[2],q[13];
p(-1540.1657984223962) q[13];
cx q[2],q[13];
p(1540.1657984223962) q[13];
cx q[1],q[13];
p(-3080.3315968447923) q[13];
cx q[1],q[13];
p(3080.3315968447923) q[13];
cx q[0],q[13];
p(-1961*pi) q[13];
cx q[0],q[13];
p(1961*pi) q[13];
ry(-pi/2) q[0];
rz(-pi) q[0];
p(-pi/4) q[1];
cx q[0],q[1];
p(pi/4) q[1];
cx q[0],q[1];
p(-pi/2048) q[10];
p(-pi/4096) q[11];
p(-pi/8192) q[12];
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
p(-pi/128) q[6];
cx q[0],q[6];
p(pi/128) q[6];
cx q[0],q[6];
p(-pi/256) q[7];
cx q[0],q[7];
p(pi/256) q[7];
cx q[0],q[7];
p(-pi/512) q[8];
cx q[0],q[8];
p(pi/512) q[8];
cx q[0],q[8];
p(-pi/1024) q[9];
cx q[0],q[9];
p(pi/1024) q[9];
cx q[0],q[9];
cx q[0],q[10];
p(pi/2048) q[10];
cx q[0],q[10];
cx q[0],q[11];
p(pi/4096) q[11];
cx q[0],q[11];
cx q[0],q[12];
p(pi/8192) q[12];
cx q[0],q[12];
ry(-pi/2) q[1];
rz(-pi) q[1];
p(-pi/1024) q[10];
p(-pi/2048) q[11];
p(-pi/4096) q[12];
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
p(-pi/64) q[6];
cx q[1],q[6];
p(pi/64) q[6];
cx q[1],q[6];
p(-pi/128) q[7];
cx q[1],q[7];
p(pi/128) q[7];
cx q[1],q[7];
p(-pi/256) q[8];
cx q[1],q[8];
p(pi/256) q[8];
cx q[1],q[8];
p(-pi/512) q[9];
cx q[1],q[9];
p(pi/512) q[9];
cx q[1],q[9];
cx q[1],q[10];
p(pi/1024) q[10];
cx q[1],q[10];
cx q[1],q[11];
p(pi/2048) q[11];
cx q[1],q[11];
cx q[1],q[12];
p(pi/4096) q[12];
cx q[1],q[12];
p(-pi/512) q[10];
p(-pi/1024) q[11];
p(-pi/2048) q[12];
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
p(-pi/32) q[6];
cx q[2],q[6];
p(pi/32) q[6];
cx q[2],q[6];
p(-pi/64) q[7];
cx q[2],q[7];
p(pi/64) q[7];
cx q[2],q[7];
p(-pi/128) q[8];
cx q[2],q[8];
p(pi/128) q[8];
cx q[2],q[8];
p(-pi/256) q[9];
cx q[2],q[9];
p(pi/256) q[9];
cx q[2],q[9];
cx q[2],q[10];
p(pi/512) q[10];
cx q[2],q[10];
cx q[2],q[11];
p(pi/1024) q[11];
cx q[2],q[11];
cx q[2],q[12];
p(pi/2048) q[12];
cx q[2],q[12];
p(-pi/256) q[10];
p(-pi/512) q[11];
p(-pi/1024) q[12];
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
p(-pi/16) q[6];
cx q[3],q[6];
p(pi/16) q[6];
cx q[3],q[6];
p(-pi/32) q[7];
cx q[3],q[7];
p(pi/32) q[7];
cx q[3],q[7];
p(-pi/64) q[8];
cx q[3],q[8];
p(pi/64) q[8];
cx q[3],q[8];
p(-pi/128) q[9];
cx q[3],q[9];
p(pi/128) q[9];
cx q[3],q[9];
cx q[3],q[10];
p(pi/256) q[10];
cx q[3],q[10];
cx q[3],q[11];
p(pi/512) q[11];
cx q[3],q[11];
cx q[3],q[12];
p(pi/1024) q[12];
cx q[3],q[12];
p(-pi/128) q[10];
p(-pi/256) q[11];
p(-pi/512) q[12];
ry(-pi/2) q[4];
rz(-pi) q[4];
p(-pi/4) q[5];
cx q[4],q[5];
p(pi/4) q[5];
cx q[4],q[5];
p(-pi/8) q[6];
cx q[4],q[6];
p(pi/8) q[6];
cx q[4],q[6];
p(-pi/16) q[7];
cx q[4],q[7];
p(pi/16) q[7];
cx q[4],q[7];
p(-pi/32) q[8];
cx q[4],q[8];
p(pi/32) q[8];
cx q[4],q[8];
p(-pi/64) q[9];
cx q[4],q[9];
p(pi/64) q[9];
cx q[4],q[9];
cx q[4],q[10];
p(pi/128) q[10];
cx q[4],q[10];
cx q[4],q[11];
p(pi/256) q[11];
cx q[4],q[11];
cx q[4],q[12];
p(pi/512) q[12];
cx q[4],q[12];
p(-pi/64) q[10];
p(-pi/128) q[11];
p(-pi/256) q[12];
ry(-pi/2) q[5];
rz(-pi) q[5];
p(-pi/4) q[6];
cx q[5],q[6];
p(pi/4) q[6];
cx q[5],q[6];
p(-pi/8) q[7];
cx q[5],q[7];
p(pi/8) q[7];
cx q[5],q[7];
p(-pi/16) q[8];
cx q[5],q[8];
p(pi/16) q[8];
cx q[5],q[8];
p(-pi/32) q[9];
cx q[5],q[9];
p(pi/32) q[9];
cx q[5],q[9];
cx q[5],q[10];
p(pi/64) q[10];
cx q[5],q[10];
cx q[5],q[11];
p(pi/128) q[11];
cx q[5],q[11];
cx q[5],q[12];
p(pi/256) q[12];
cx q[5],q[12];
p(-pi/32) q[10];
p(-pi/64) q[11];
p(-pi/128) q[12];
ry(-pi/2) q[6];
rz(-pi) q[6];
p(-pi/4) q[7];
cx q[6],q[7];
p(pi/4) q[7];
cx q[6],q[7];
p(-pi/8) q[8];
cx q[6],q[8];
p(pi/8) q[8];
cx q[6],q[8];
p(-pi/16) q[9];
cx q[6],q[9];
p(pi/16) q[9];
cx q[6],q[9];
cx q[6],q[10];
p(pi/32) q[10];
cx q[6],q[10];
cx q[6],q[11];
p(pi/64) q[11];
cx q[6],q[11];
cx q[6],q[12];
p(pi/128) q[12];
cx q[6],q[12];
p(-pi/16) q[10];
p(-pi/32) q[11];
p(-pi/64) q[12];
ry(-pi/2) q[7];
rz(-pi) q[7];
p(-pi/4) q[8];
cx q[7],q[8];
p(pi/4) q[8];
cx q[7],q[8];
p(-pi/8) q[9];
cx q[7],q[9];
p(pi/8) q[9];
cx q[7],q[9];
cx q[7],q[10];
p(pi/16) q[10];
cx q[7],q[10];
cx q[7],q[11];
p(pi/32) q[11];
cx q[7],q[11];
cx q[7],q[12];
p(pi/64) q[12];
cx q[7],q[12];
p(-pi/8) q[10];
p(-pi/16) q[11];
p(-pi/32) q[12];
ry(-pi/2) q[8];
rz(-pi) q[8];
p(-pi/4) q[9];
cx q[8],q[9];
p(pi/4) q[9];
cx q[8],q[9];
cx q[8],q[10];
p(pi/8) q[10];
cx q[8],q[10];
cx q[8],q[11];
p(pi/16) q[11];
cx q[8],q[11];
cx q[8],q[12];
p(pi/32) q[12];
cx q[8],q[12];
p(-pi/4) q[10];
p(-pi/8) q[11];
p(-pi/16) q[12];
ry(-pi/2) q[9];
rz(-pi) q[9];
cx q[9],q[10];
p(pi/4) q[10];
cx q[9],q[10];
cx q[9],q[11];
p(pi/8) q[11];
cx q[9],q[11];
cx q[9],q[12];
p(pi/16) q[12];
cx q[9],q[12];
ry(-pi/2) q[10];
rz(-pi) q[10];
p(-pi/4) q[11];
cx q[10],q[11];
p(pi/4) q[11];
cx q[10],q[11];
p(-pi/8) q[12];
cx q[10],q[12];
p(pi/8) q[12];
cx q[10],q[12];
ry(-pi/2) q[11];
rz(-pi) q[11];
p(-pi/4) q[12];
cx q[11],q[12];
p(pi/4) q[12];
cx q[11],q[12];
ry(-pi/2) q[12];
rz(-pi) q[12];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
measure q[8] -> c[8];
measure q[9] -> c[9];
measure q[10] -> c[10];
measure q[11] -> c[11];
measure q[12] -> c[12];
measure q[13] -> c[13];
