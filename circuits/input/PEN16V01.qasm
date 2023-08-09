OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
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
ry(-pi/2) q[13];
rz(-pi) q[13];
ry(-pi/2) q[14];
rz(-pi) q[14];
rx(-pi) q[15];
p(23.120158435012385) q[10];
p(11.560079217506193) q[11];
p(5.780039608753096) q[12];
p(2.890019804376548) q[13];
p(1.445009902188274) q[14];
cx q[14],q[15];
p(-1.445009902188274) q[15];
cx q[14],q[15];
p(1.445009902188274) q[15];
cx q[13],q[15];
p(-2.890019804376548) q[15];
cx q[13],q[15];
p(2.890019804376548) q[15];
cx q[12],q[15];
p(-5.780039608753096) q[15];
cx q[12],q[15];
p(5.780039608753096) q[15];
cx q[11],q[15];
p(-11.560079217506193) q[15];
cx q[11],q[15];
p(11.560079217506193) q[15];
cx q[10],q[15];
p(-23.120158435012385) q[15];
cx q[10],q[15];
p(23.120158435012385) q[15];
p(471*pi) q[4];
p(739.8450699203963) q[5];
p(369.92253496019816) q[6];
p(184.96126748009908) q[7];
p(92.48063374004954) q[8];
p(46.24031687002477) q[9];
cx q[9],q[15];
p(-46.24031687002477) q[15];
cx q[9],q[15];
p(46.24031687002477) q[15];
cx q[8],q[15];
p(-92.48063374004954) q[15];
cx q[8],q[15];
p(92.48063374004954) q[15];
cx q[7],q[15];
p(-184.96126748009908) q[15];
cx q[7],q[15];
p(184.96126748009908) q[15];
cx q[6],q[15];
p(-369.92253496019816) q[15];
cx q[6],q[15];
p(369.92253496019816) q[15];
cx q[5],q[15];
p(-739.8450699203963) q[15];
cx q[5],q[15];
p(739.8450699203963) q[15];
cx q[4],q[15];
p(-471*pi) q[15];
cx q[4],q[15];
p(471*pi) q[15];
ry(-pi/2) q[0];
rz(-pi) q[0];
p(-pi/4) q[1];
cx q[0],q[1];
p(pi/4) q[1];
cx q[0],q[1];
p(-pi/2048) q[10];
p(-pi/4096) q[11];
p(-pi/8192) q[12];
p(-pi/16384) q[13];
p(-pi/32768) q[14];
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
cx q[0],q[13];
p(pi/16384) q[13];
cx q[0],q[13];
cx q[0],q[14];
p(pi/32768) q[14];
cx q[0],q[14];
ry(-pi/2) q[1];
rz(-pi) q[1];
p(-pi/1024) q[10];
p(-pi/2048) q[11];
p(-pi/4096) q[12];
p(-pi/8192) q[13];
p(-pi/16384) q[14];
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
cx q[1],q[13];
p(pi/8192) q[13];
cx q[1],q[13];
cx q[1],q[14];
p(pi/16384) q[14];
cx q[1],q[14];
p(-pi/512) q[10];
p(-pi/1024) q[11];
p(-pi/2048) q[12];
p(-pi/4096) q[13];
p(-pi/8192) q[14];
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
cx q[2],q[13];
p(pi/4096) q[13];
cx q[2],q[13];
cx q[2],q[14];
p(pi/8192) q[14];
cx q[2],q[14];
p(-pi/256) q[10];
p(-pi/512) q[11];
p(-pi/1024) q[12];
p(-pi/2048) q[13];
p(-pi/4096) q[14];
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
cx q[3],q[13];
p(pi/2048) q[13];
cx q[3],q[13];
cx q[3],q[14];
p(pi/4096) q[14];
cx q[3],q[14];
p(-pi/128) q[10];
p(-pi/256) q[11];
p(-pi/512) q[12];
p(-pi/1024) q[13];
p(-pi/2048) q[14];
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
cx q[4],q[13];
p(pi/1024) q[13];
cx q[4],q[13];
cx q[4],q[14];
p(pi/2048) q[14];
cx q[4],q[14];
p(-pi/64) q[10];
p(-pi/128) q[11];
p(-pi/256) q[12];
p(-pi/512) q[13];
p(-pi/1024) q[14];
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
cx q[5],q[13];
p(pi/512) q[13];
cx q[5],q[13];
cx q[5],q[14];
p(pi/1024) q[14];
cx q[5],q[14];
p(-pi/32) q[10];
p(-pi/64) q[11];
p(-pi/128) q[12];
p(-pi/256) q[13];
p(-pi/512) q[14];
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
cx q[6],q[13];
p(pi/256) q[13];
cx q[6],q[13];
cx q[6],q[14];
p(pi/512) q[14];
cx q[6],q[14];
p(-pi/16) q[10];
p(-pi/32) q[11];
p(-pi/64) q[12];
p(-pi/128) q[13];
p(-pi/256) q[14];
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
cx q[7],q[13];
p(pi/128) q[13];
cx q[7],q[13];
cx q[7],q[14];
p(pi/256) q[14];
cx q[7],q[14];
p(-pi/8) q[10];
p(-pi/16) q[11];
p(-pi/32) q[12];
p(-pi/64) q[13];
p(-pi/128) q[14];
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
cx q[8],q[13];
p(pi/64) q[13];
cx q[8],q[13];
cx q[8],q[14];
p(pi/128) q[14];
cx q[8],q[14];
p(-pi/4) q[10];
p(-pi/8) q[11];
p(-pi/16) q[12];
p(-pi/32) q[13];
p(-pi/64) q[14];
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
cx q[9],q[13];
p(pi/32) q[13];
cx q[9],q[13];
cx q[9],q[14];
p(pi/64) q[14];
cx q[9],q[14];
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
p(-pi/16) q[13];
cx q[10],q[13];
p(pi/16) q[13];
cx q[10],q[13];
p(-pi/32) q[14];
cx q[10],q[14];
p(pi/32) q[14];
cx q[10],q[14];
ry(-pi/2) q[11];
rz(-pi) q[11];
p(-pi/4) q[12];
cx q[11],q[12];
p(pi/4) q[12];
cx q[11],q[12];
p(-pi/8) q[13];
cx q[11],q[13];
p(pi/8) q[13];
cx q[11],q[13];
p(-pi/16) q[14];
cx q[11],q[14];
p(pi/16) q[14];
cx q[11],q[14];
ry(-pi/2) q[12];
rz(-pi) q[12];
p(-pi/4) q[13];
cx q[12],q[13];
p(pi/4) q[13];
cx q[12],q[13];
p(-pi/8) q[14];
cx q[12],q[14];
p(pi/8) q[14];
cx q[12],q[14];
ry(-pi/2) q[13];
rz(-pi) q[13];
p(-pi/4) q[14];
cx q[13],q[14];
p(pi/4) q[14];
cx q[13],q[14];
ry(-pi/2) q[14];
rz(-pi) q[14];
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
measure q[14] -> c[14];
measure q[15] -> c[15];
