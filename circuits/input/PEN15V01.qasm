OPENQASM 2.0;
include "qelib1.inc";
qreg q[15];
creg c[15];
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
rx(-pi) q[14];
p(5.878214379177777) q[10];
p(2.9391071895888885) q[11];
p(1.4695535947944443) q[12];
p(0.7347767973972221) q[13];
cx q[13],q[14];
p(-0.7347767973972221) q[14];
cx q[13],q[14];
p(0.7347767973972221) q[14];
cx q[12],q[14];
p(-1.4695535947944443) q[14];
cx q[12],q[14];
p(1.4695535947944443) q[14];
cx q[11],q[14];
p(-2.9391071895888885) q[14];
cx q[11],q[14];
p(2.9391071895888885) q[14];
cx q[10],q[14];
p(-5.878214379177777) q[14];
cx q[10],q[14];
p(5.878214379177777) q[14];
p(479*pi) q[2];
p(752.4114405347555) q[3];
p(376.2057202673777) q[4];
p(188.10286013368886) q[5];
p(94.05143006684443) q[6];
p(47.025715033422216) q[7];
p(23.512857516711108) q[8];
p(11.756428758355554) q[9];
cx q[9],q[14];
p(-11.756428758355554) q[14];
cx q[9],q[14];
p(11.756428758355554) q[14];
cx q[8],q[14];
p(-23.512857516711108) q[14];
cx q[8],q[14];
p(23.512857516711108) q[14];
cx q[7],q[14];
p(-47.025715033422216) q[14];
cx q[7],q[14];
p(47.025715033422216) q[14];
cx q[6],q[14];
p(-94.05143006684443) q[14];
cx q[6],q[14];
p(94.05143006684443) q[14];
cx q[5],q[14];
p(-188.10286013368886) q[14];
cx q[5],q[14];
p(188.10286013368886) q[14];
cx q[4],q[14];
p(-376.2057202673777) q[14];
cx q[4],q[14];
p(376.2057202673777) q[14];
cx q[3],q[14];
p(-752.4114405347555) q[14];
cx q[3],q[14];
p(752.4114405347555) q[14];
cx q[2],q[14];
p(-479*pi) q[14];
cx q[2],q[14];
p(479*pi) q[14];
ry(-pi/2) q[0];
rz(-pi) q[0];
rz(-pi/4) q[1];
cx q[0],q[1];
rz(pi/4) q[1];
cx q[0],q[1];
rz(-pi/2048) q[10];
rz(-pi/4096) q[11];
rz(-pi/8192) q[12];
rz(-pi/16384) q[13];
rz(-pi/8) q[2];
cx q[0],q[2];
rz(pi/8) q[2];
cx q[0],q[2];
rz(-pi/16) q[3];
cx q[0],q[3];
rz(pi/16) q[3];
cx q[0],q[3];
rz(-pi/32) q[4];
cx q[0],q[4];
rz(pi/32) q[4];
cx q[0],q[4];
rz(-pi/64) q[5];
cx q[0],q[5];
rz(pi/64) q[5];
cx q[0],q[5];
rz(-pi/128) q[6];
cx q[0],q[6];
rz(pi/128) q[6];
cx q[0],q[6];
rz(-pi/256) q[7];
cx q[0],q[7];
rz(pi/256) q[7];
cx q[0],q[7];
rz(-pi/512) q[8];
cx q[0],q[8];
rz(pi/512) q[8];
cx q[0],q[8];
rz(-pi/1024) q[9];
cx q[0],q[9];
rz(pi/1024) q[9];
cx q[0],q[9];
cx q[0],q[10];
rz(pi/2048) q[10];
cx q[0],q[10];
cx q[0],q[11];
rz(pi/4096) q[11];
cx q[0],q[11];
cx q[0],q[12];
rz(pi/8192) q[12];
cx q[0],q[12];
cx q[0],q[13];
rz(pi/16384) q[13];
cx q[0],q[13];
ry(-pi/2) q[1];
rz(-pi) q[1];
rz(-pi/1024) q[10];
rz(-pi/2048) q[11];
rz(-pi/4096) q[12];
rz(-pi/8192) q[13];
rz(-pi/4) q[2];
cx q[1],q[2];
rz(pi/4) q[2];
cx q[1],q[2];
rz(-pi/8) q[3];
cx q[1],q[3];
rz(pi/8) q[3];
cx q[1],q[3];
rz(-pi/16) q[4];
cx q[1],q[4];
rz(pi/16) q[4];
cx q[1],q[4];
rz(-pi/32) q[5];
cx q[1],q[5];
rz(pi/32) q[5];
cx q[1],q[5];
rz(-pi/64) q[6];
cx q[1],q[6];
rz(pi/64) q[6];
cx q[1],q[6];
rz(-pi/128) q[7];
cx q[1],q[7];
rz(pi/128) q[7];
cx q[1],q[7];
rz(-pi/256) q[8];
cx q[1],q[8];
rz(pi/256) q[8];
cx q[1],q[8];
rz(-pi/512) q[9];
cx q[1],q[9];
rz(pi/512) q[9];
cx q[1],q[9];
cx q[1],q[10];
rz(pi/1024) q[10];
cx q[1],q[10];
cx q[1],q[11];
rz(pi/2048) q[11];
cx q[1],q[11];
cx q[1],q[12];
rz(pi/4096) q[12];
cx q[1],q[12];
cx q[1],q[13];
rz(pi/8192) q[13];
cx q[1],q[13];
rz(-pi/512) q[10];
rz(-pi/1024) q[11];
rz(-pi/2048) q[12];
rz(-pi/4096) q[13];
ry(-pi/2) q[2];
rz(-pi) q[2];
rz(-pi/4) q[3];
cx q[2],q[3];
rz(pi/4) q[3];
cx q[2],q[3];
rz(-pi/8) q[4];
cx q[2],q[4];
rz(pi/8) q[4];
cx q[2],q[4];
rz(-pi/16) q[5];
cx q[2],q[5];
rz(pi/16) q[5];
cx q[2],q[5];
rz(-pi/32) q[6];
cx q[2],q[6];
rz(pi/32) q[6];
cx q[2],q[6];
rz(-pi/64) q[7];
cx q[2],q[7];
rz(pi/64) q[7];
cx q[2],q[7];
rz(-pi/128) q[8];
cx q[2],q[8];
rz(pi/128) q[8];
cx q[2],q[8];
rz(-pi/256) q[9];
cx q[2],q[9];
rz(pi/256) q[9];
cx q[2],q[9];
cx q[2],q[10];
rz(pi/512) q[10];
cx q[2],q[10];
cx q[2],q[11];
rz(pi/1024) q[11];
cx q[2],q[11];
cx q[2],q[12];
rz(pi/2048) q[12];
cx q[2],q[12];
cx q[2],q[13];
rz(pi/4096) q[13];
cx q[2],q[13];
rz(-pi/256) q[10];
rz(-pi/512) q[11];
rz(-pi/1024) q[12];
rz(-pi/2048) q[13];
ry(-pi/2) q[3];
rz(-pi) q[3];
rz(-pi/4) q[4];
cx q[3],q[4];
rz(pi/4) q[4];
cx q[3],q[4];
rz(-pi/8) q[5];
cx q[3],q[5];
rz(pi/8) q[5];
cx q[3],q[5];
rz(-pi/16) q[6];
cx q[3],q[6];
rz(pi/16) q[6];
cx q[3],q[6];
rz(-pi/32) q[7];
cx q[3],q[7];
rz(pi/32) q[7];
cx q[3],q[7];
rz(-pi/64) q[8];
cx q[3],q[8];
rz(pi/64) q[8];
cx q[3],q[8];
rz(-pi/128) q[9];
cx q[3],q[9];
rz(pi/128) q[9];
cx q[3],q[9];
cx q[3],q[10];
rz(pi/256) q[10];
cx q[3],q[10];
cx q[3],q[11];
rz(pi/512) q[11];
cx q[3],q[11];
cx q[3],q[12];
rz(pi/1024) q[12];
cx q[3],q[12];
cx q[3],q[13];
rz(pi/2048) q[13];
cx q[3],q[13];
rz(-pi/128) q[10];
rz(-pi/256) q[11];
rz(-pi/512) q[12];
rz(-pi/1024) q[13];
ry(-pi/2) q[4];
rz(-pi) q[4];
rz(-pi/4) q[5];
cx q[4],q[5];
rz(pi/4) q[5];
cx q[4],q[5];
rz(-pi/8) q[6];
cx q[4],q[6];
rz(pi/8) q[6];
cx q[4],q[6];
rz(-pi/16) q[7];
cx q[4],q[7];
rz(pi/16) q[7];
cx q[4],q[7];
rz(-pi/32) q[8];
cx q[4],q[8];
rz(pi/32) q[8];
cx q[4],q[8];
rz(-pi/64) q[9];
cx q[4],q[9];
rz(pi/64) q[9];
cx q[4],q[9];
cx q[4],q[10];
rz(pi/128) q[10];
cx q[4],q[10];
cx q[4],q[11];
rz(pi/256) q[11];
cx q[4],q[11];
cx q[4],q[12];
rz(pi/512) q[12];
cx q[4],q[12];
cx q[4],q[13];
rz(pi/1024) q[13];
cx q[4],q[13];
rz(-pi/64) q[10];
rz(-pi/128) q[11];
rz(-pi/256) q[12];
rz(-pi/512) q[13];
ry(-pi/2) q[5];
rz(-pi) q[5];
rz(-pi/4) q[6];
cx q[5],q[6];
rz(pi/4) q[6];
cx q[5],q[6];
rz(-pi/8) q[7];
cx q[5],q[7];
rz(pi/8) q[7];
cx q[5],q[7];
rz(-pi/16) q[8];
cx q[5],q[8];
rz(pi/16) q[8];
cx q[5],q[8];
rz(-pi/32) q[9];
cx q[5],q[9];
rz(pi/32) q[9];
cx q[5],q[9];
cx q[5],q[10];
rz(pi/64) q[10];
cx q[5],q[10];
cx q[5],q[11];
rz(pi/128) q[11];
cx q[5],q[11];
cx q[5],q[12];
rz(pi/256) q[12];
cx q[5],q[12];
cx q[5],q[13];
rz(pi/512) q[13];
cx q[5],q[13];
rz(-pi/32) q[10];
rz(-pi/64) q[11];
rz(-pi/128) q[12];
rz(-pi/256) q[13];
ry(-pi/2) q[6];
rz(-pi) q[6];
rz(-pi/4) q[7];
cx q[6],q[7];
rz(pi/4) q[7];
cx q[6],q[7];
rz(-pi/8) q[8];
cx q[6],q[8];
rz(pi/8) q[8];
cx q[6],q[8];
rz(-pi/16) q[9];
cx q[6],q[9];
rz(pi/16) q[9];
cx q[6],q[9];
cx q[6],q[10];
rz(pi/32) q[10];
cx q[6],q[10];
cx q[6],q[11];
rz(pi/64) q[11];
cx q[6],q[11];
cx q[6],q[12];
rz(pi/128) q[12];
cx q[6],q[12];
cx q[6],q[13];
rz(pi/256) q[13];
cx q[6],q[13];
rz(-pi/16) q[10];
rz(-pi/32) q[11];
rz(-pi/64) q[12];
rz(-pi/128) q[13];
ry(-pi/2) q[7];
rz(-pi) q[7];
rz(-pi/4) q[8];
cx q[7],q[8];
rz(pi/4) q[8];
cx q[7],q[8];
rz(-pi/8) q[9];
cx q[7],q[9];
rz(pi/8) q[9];
cx q[7],q[9];
cx q[7],q[10];
rz(pi/16) q[10];
cx q[7],q[10];
cx q[7],q[11];
rz(pi/32) q[11];
cx q[7],q[11];
cx q[7],q[12];
rz(pi/64) q[12];
cx q[7],q[12];
cx q[7],q[13];
rz(pi/128) q[13];
cx q[7],q[13];
rz(-pi/8) q[10];
rz(-pi/16) q[11];
rz(-pi/32) q[12];
rz(-pi/64) q[13];
ry(-pi/2) q[8];
rz(-pi) q[8];
rz(-pi/4) q[9];
cx q[8],q[9];
rz(pi/4) q[9];
cx q[8],q[9];
cx q[8],q[10];
rz(pi/8) q[10];
cx q[8],q[10];
cx q[8],q[11];
rz(pi/16) q[11];
cx q[8],q[11];
cx q[8],q[12];
rz(pi/32) q[12];
cx q[8],q[12];
cx q[8],q[13];
rz(pi/64) q[13];
cx q[8],q[13];
rz(-pi/4) q[10];
rz(-pi/8) q[11];
rz(-pi/16) q[12];
rz(-pi/32) q[13];
ry(-pi/2) q[9];
rz(-pi) q[9];
cx q[9],q[10];
rz(pi/4) q[10];
cx q[9],q[10];
cx q[9],q[11];
rz(pi/8) q[11];
cx q[9],q[11];
cx q[9],q[12];
rz(pi/16) q[12];
cx q[9],q[12];
cx q[9],q[13];
rz(pi/32) q[13];
cx q[9],q[13];
ry(-pi/2) q[10];
rz(-pi) q[10];
rz(-pi/4) q[11];
cx q[10],q[11];
rz(pi/4) q[11];
cx q[10],q[11];
rz(-pi/8) q[12];
cx q[10],q[12];
rz(pi/8) q[12];
cx q[10],q[12];
rz(-pi/16) q[13];
cx q[10],q[13];
rz(pi/16) q[13];
cx q[10],q[13];
ry(-pi/2) q[11];
rz(-pi) q[11];
rz(-pi/4) q[12];
cx q[11],q[12];
rz(pi/4) q[12];
cx q[11],q[12];
rz(-pi/8) q[13];
cx q[11],q[13];
rz(pi/8) q[13];
cx q[11],q[13];
ry(-pi/2) q[12];
rz(-pi) q[12];
rz(-pi/4) q[13];
cx q[12],q[13];
rz(pi/4) q[13];
cx q[12],q[13];
ry(-pi/2) q[13];
rz(-pi) q[13];
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
