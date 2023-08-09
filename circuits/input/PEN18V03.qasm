OPENQASM 2.0;
include "qelib1.inc";
qreg q[18];
creg c[18];
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
ry(-pi/2) q[15];
rz(-pi) q[15];
ry(-pi/2) q[16];
rz(-pi) q[16];
rx(-pi) q[17];
p(26836*pi) q[0];
p(13418*pi) q[1];
p(82.33181684739813) q[10];
p(41.165908423699065) q[11];
p(20.582954211849533) q[12];
p(10.291477105924766) q[13];
p(5.145738552962383) q[14];
p(2.5728692764811916) q[15];
p(1.2864346382405958) q[16];
cx q[16],q[17];
p(-1.2864346382405958) q[17];
cx q[16],q[17];
p(1.2864346382405958) q[17];
cx q[15],q[17];
p(-2.5728692764811916) q[17];
cx q[15],q[17];
p(2.5728692764811916) q[17];
cx q[14],q[17];
p(-5.145738552962383) q[17];
cx q[14],q[17];
p(5.145738552962383) q[17];
cx q[13],q[17];
p(-10.291477105924766) q[17];
cx q[13],q[17];
p(10.291477105924766) q[17];
cx q[12],q[17];
p(-20.582954211849533) q[17];
cx q[12],q[17];
p(20.582954211849533) q[17];
cx q[11],q[17];
p(-41.165908423699065) q[17];
cx q[11],q[17];
p(41.165908423699065) q[17];
cx q[10],q[17];
p(-82.33181684739813) q[17];
cx q[10],q[17];
p(82.33181684739813) q[17];
p(6709*pi) q[2];
p(10538.47255646696) q[3];
p(5269.23627823348) q[4];
p(2634.61813911674) q[5];
p(1317.30906955837) q[6];
p(658.654534779185) q[7];
p(329.3272673895925) q[8];
p(164.66363369479626) q[9];
cx q[9],q[17];
p(-164.66363369479626) q[17];
cx q[9],q[17];
p(164.66363369479626) q[17];
cx q[8],q[17];
p(-329.3272673895925) q[17];
cx q[8],q[17];
p(329.3272673895925) q[17];
cx q[7],q[17];
p(-658.654534779185) q[17];
cx q[7],q[17];
p(658.654534779185) q[17];
cx q[6],q[17];
p(-1317.30906955837) q[17];
cx q[6],q[17];
p(1317.30906955837) q[17];
cx q[5],q[17];
p(-2634.61813911674) q[17];
cx q[5],q[17];
p(2634.61813911674) q[17];
cx q[4],q[17];
p(-5269.23627823348) q[17];
cx q[4],q[17];
p(5269.23627823348) q[17];
cx q[3],q[17];
p(-10538.47255646696) q[17];
cx q[3],q[17];
p(10538.47255646696) q[17];
cx q[2],q[17];
p(-6709*pi) q[17];
cx q[2],q[17];
p(6709*pi) q[17];
cx q[1],q[17];
p(-13418*pi) q[17];
cx q[1],q[17];
p(13418*pi) q[17];
cx q[0],q[17];
p(-26836*pi) q[17];
cx q[0],q[17];
p(26836*pi) q[17];
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
rz(-pi/32768) q[14];
rz(-pi/65536) q[15];
rz(-pi/131072) q[16];
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
cx q[0],q[14];
rz(pi/32768) q[14];
cx q[0],q[14];
cx q[0],q[15];
rz(pi/65536) q[15];
cx q[0],q[15];
cx q[0],q[16];
rz(pi/131072) q[16];
cx q[0],q[16];
ry(-pi/2) q[1];
rz(-pi) q[1];
rz(-pi/1024) q[10];
rz(-pi/2048) q[11];
rz(-pi/4096) q[12];
rz(-pi/8192) q[13];
rz(-pi/16384) q[14];
rz(-pi/32768) q[15];
rz(-pi/65536) q[16];
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
cx q[1],q[14];
rz(pi/16384) q[14];
cx q[1],q[14];
cx q[1],q[15];
rz(pi/32768) q[15];
cx q[1],q[15];
cx q[1],q[16];
rz(pi/65536) q[16];
cx q[1],q[16];
rz(-pi/512) q[10];
rz(-pi/1024) q[11];
rz(-pi/2048) q[12];
rz(-pi/4096) q[13];
rz(-pi/8192) q[14];
rz(-pi/16384) q[15];
rz(-pi/32768) q[16];
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
cx q[2],q[14];
rz(pi/8192) q[14];
cx q[2],q[14];
cx q[2],q[15];
rz(pi/16384) q[15];
cx q[2],q[15];
cx q[2],q[16];
rz(pi/32768) q[16];
cx q[2],q[16];
rz(-pi/256) q[10];
rz(-pi/512) q[11];
rz(-pi/1024) q[12];
rz(-pi/2048) q[13];
rz(-pi/4096) q[14];
rz(-pi/8192) q[15];
rz(-pi/16384) q[16];
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
cx q[3],q[14];
rz(pi/4096) q[14];
cx q[3],q[14];
cx q[3],q[15];
rz(pi/8192) q[15];
cx q[3],q[15];
cx q[3],q[16];
rz(pi/16384) q[16];
cx q[3],q[16];
rz(-pi/128) q[10];
rz(-pi/256) q[11];
rz(-pi/512) q[12];
rz(-pi/1024) q[13];
rz(-pi/2048) q[14];
rz(-pi/4096) q[15];
rz(-pi/8192) q[16];
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
cx q[4],q[14];
rz(pi/2048) q[14];
cx q[4],q[14];
cx q[4],q[15];
rz(pi/4096) q[15];
cx q[4],q[15];
cx q[4],q[16];
rz(pi/8192) q[16];
cx q[4],q[16];
rz(-pi/64) q[10];
rz(-pi/128) q[11];
rz(-pi/256) q[12];
rz(-pi/512) q[13];
rz(-pi/1024) q[14];
rz(-pi/2048) q[15];
rz(-pi/4096) q[16];
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
cx q[5],q[14];
rz(pi/1024) q[14];
cx q[5],q[14];
cx q[5],q[15];
rz(pi/2048) q[15];
cx q[5],q[15];
cx q[5],q[16];
rz(pi/4096) q[16];
cx q[5],q[16];
rz(-pi/32) q[10];
rz(-pi/64) q[11];
rz(-pi/128) q[12];
rz(-pi/256) q[13];
rz(-pi/512) q[14];
rz(-pi/1024) q[15];
rz(-pi/2048) q[16];
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
cx q[6],q[14];
rz(pi/512) q[14];
cx q[6],q[14];
cx q[6],q[15];
rz(pi/1024) q[15];
cx q[6],q[15];
cx q[6],q[16];
rz(pi/2048) q[16];
cx q[6],q[16];
rz(-pi/16) q[10];
rz(-pi/32) q[11];
rz(-pi/64) q[12];
rz(-pi/128) q[13];
rz(-pi/256) q[14];
rz(-pi/512) q[15];
rz(-pi/1024) q[16];
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
cx q[7],q[14];
rz(pi/256) q[14];
cx q[7],q[14];
cx q[7],q[15];
rz(pi/512) q[15];
cx q[7],q[15];
cx q[7],q[16];
rz(pi/1024) q[16];
cx q[7],q[16];
rz(-pi/8) q[10];
rz(-pi/16) q[11];
rz(-pi/32) q[12];
rz(-pi/64) q[13];
rz(-pi/128) q[14];
rz(-pi/256) q[15];
rz(-pi/512) q[16];
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
cx q[8],q[14];
rz(pi/128) q[14];
cx q[8],q[14];
cx q[8],q[15];
rz(pi/256) q[15];
cx q[8],q[15];
cx q[8],q[16];
rz(pi/512) q[16];
cx q[8],q[16];
rz(-pi/4) q[10];
rz(-pi/8) q[11];
rz(-pi/16) q[12];
rz(-pi/32) q[13];
rz(-pi/64) q[14];
rz(-pi/128) q[15];
rz(-pi/256) q[16];
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
cx q[9],q[14];
rz(pi/64) q[14];
cx q[9],q[14];
cx q[9],q[15];
rz(pi/128) q[15];
cx q[9],q[15];
cx q[9],q[16];
rz(pi/256) q[16];
cx q[9],q[16];
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
rz(-pi/32) q[14];
cx q[10],q[14];
rz(pi/32) q[14];
cx q[10],q[14];
rz(-pi/64) q[15];
cx q[10],q[15];
rz(pi/64) q[15];
cx q[10],q[15];
rz(-pi/128) q[16];
cx q[10],q[16];
rz(pi/128) q[16];
cx q[10],q[16];
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
rz(-pi/16) q[14];
cx q[11],q[14];
rz(pi/16) q[14];
cx q[11],q[14];
rz(-pi/32) q[15];
cx q[11],q[15];
rz(pi/32) q[15];
cx q[11],q[15];
rz(-pi/64) q[16];
cx q[11],q[16];
rz(pi/64) q[16];
cx q[11],q[16];
ry(-pi/2) q[12];
rz(-pi) q[12];
rz(-pi/4) q[13];
cx q[12],q[13];
rz(pi/4) q[13];
cx q[12],q[13];
rz(-pi/8) q[14];
cx q[12],q[14];
rz(pi/8) q[14];
cx q[12],q[14];
rz(-pi/16) q[15];
cx q[12],q[15];
rz(pi/16) q[15];
cx q[12],q[15];
rz(-pi/32) q[16];
cx q[12],q[16];
rz(pi/32) q[16];
cx q[12],q[16];
ry(-pi/2) q[13];
rz(-pi) q[13];
rz(-pi/4) q[14];
cx q[13],q[14];
rz(pi/4) q[14];
cx q[13],q[14];
rz(-pi/8) q[15];
cx q[13],q[15];
rz(pi/8) q[15];
cx q[13],q[15];
rz(-pi/16) q[16];
cx q[13],q[16];
rz(pi/16) q[16];
cx q[13],q[16];
ry(-pi/2) q[14];
rz(-pi) q[14];
rz(-pi/4) q[15];
cx q[14],q[15];
rz(pi/4) q[15];
cx q[14],q[15];
rz(-pi/8) q[16];
cx q[14],q[16];
rz(pi/8) q[16];
cx q[14],q[16];
ry(-pi/2) q[15];
rz(-pi) q[15];
rz(-pi/4) q[16];
cx q[15],q[16];
rz(pi/4) q[16];
cx q[15],q[16];
ry(-pi/2) q[16];
rz(-pi) q[16];
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
measure q[16] -> c[16];
measure q[17] -> c[17];
