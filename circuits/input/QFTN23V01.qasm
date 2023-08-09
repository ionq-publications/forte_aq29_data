OPENQASM 2.0;
include "qelib1.inc";
qreg q[23];
creg c[23];
rx(-pi) q[3];
rx(-pi) q[4];
rx(-pi) q[5];
rx(-pi) q[7];
rx(-pi) q[8];
rx(-pi) q[9];
rx(-pi) q[10];
rx(-pi) q[11];
rx(-pi) q[12];
rx(-pi) q[16];
rx(-pi) q[18];
rx(-pi) q[20];
rx(-pi) q[21];
rx(-pi) q[22];
h q[22];
crz(pi/2) q[21],q[22];
h q[21];
crz(pi/4) q[20],q[22];
crz(pi/2) q[20],q[21];
h q[20];
crz(pi/8) q[19],q[22];
crz(pi/4) q[19],q[21];
crz(pi/2) q[19],q[20];
h q[19];
crz(pi/16) q[18],q[22];
crz(pi/8) q[18],q[21];
crz(pi/4) q[18],q[20];
crz(pi/2) q[18],q[19];
h q[18];
crz(pi/32) q[17],q[22];
crz(pi/16) q[17],q[21];
crz(pi/8) q[17],q[20];
crz(pi/4) q[17],q[19];
crz(pi/2) q[17],q[18];
h q[17];
crz(pi/64) q[16],q[22];
crz(pi/32) q[16],q[21];
crz(pi/16) q[16],q[20];
crz(pi/8) q[16],q[19];
crz(pi/4) q[16],q[18];
crz(pi/2) q[16],q[17];
h q[16];
crz(pi/128) q[15],q[22];
crz(pi/64) q[15],q[21];
crz(pi/32) q[15],q[20];
crz(pi/16) q[15],q[19];
crz(pi/8) q[15],q[18];
crz(pi/4) q[15],q[17];
crz(pi/2) q[15],q[16];
h q[15];
crz(pi/256) q[14],q[22];
crz(pi/128) q[14],q[21];
crz(pi/64) q[14],q[20];
crz(pi/32) q[14],q[19];
crz(pi/16) q[14],q[18];
crz(pi/8) q[14],q[17];
crz(pi/4) q[14],q[16];
crz(pi/2) q[14],q[15];
h q[14];
crz(pi/512) q[13],q[22];
crz(pi/256) q[13],q[21];
crz(pi/128) q[13],q[20];
crz(pi/64) q[13],q[19];
crz(pi/32) q[13],q[18];
crz(pi/16) q[13],q[17];
crz(pi/8) q[13],q[16];
crz(pi/4) q[13],q[15];
crz(pi/2) q[13],q[14];
h q[13];
crz(pi/1024) q[12],q[22];
crz(pi/512) q[12],q[21];
crz(pi/256) q[12],q[20];
crz(pi/128) q[12],q[19];
crz(pi/64) q[12],q[18];
crz(pi/32) q[12],q[17];
crz(pi/16) q[12],q[16];
crz(pi/8) q[12],q[15];
crz(pi/4) q[12],q[14];
crz(pi/2) q[12],q[13];
h q[12];
crz(pi/2048) q[11],q[22];
crz(pi/1024) q[11],q[21];
crz(pi/512) q[11],q[20];
crz(pi/256) q[11],q[19];
crz(pi/128) q[11],q[18];
crz(pi/64) q[11],q[17];
crz(pi/32) q[11],q[16];
crz(pi/16) q[11],q[15];
crz(pi/8) q[11],q[14];
crz(pi/4) q[11],q[13];
crz(pi/2) q[11],q[12];
h q[11];
crz(pi/4096) q[10],q[22];
crz(pi/2048) q[10],q[21];
crz(pi/1024) q[10],q[20];
crz(pi/512) q[10],q[19];
crz(pi/256) q[10],q[18];
crz(pi/128) q[10],q[17];
crz(pi/64) q[10],q[16];
crz(pi/32) q[10],q[15];
crz(pi/16) q[10],q[14];
crz(pi/8) q[10],q[13];
crz(pi/4) q[10],q[12];
crz(pi/2) q[10],q[11];
h q[10];
crz(pi/8192) q[9],q[22];
crz(pi/4096) q[9],q[21];
crz(pi/2048) q[9],q[20];
crz(pi/1024) q[9],q[19];
crz(pi/512) q[9],q[18];
crz(pi/256) q[9],q[17];
crz(pi/128) q[9],q[16];
crz(pi/64) q[9],q[15];
crz(pi/32) q[9],q[14];
crz(pi/16) q[9],q[13];
crz(pi/8) q[9],q[12];
crz(pi/4) q[9],q[11];
crz(pi/2) q[9],q[10];
h q[9];
crz(pi/16384) q[8],q[22];
crz(pi/8192) q[8],q[21];
crz(pi/4096) q[8],q[20];
crz(pi/2048) q[8],q[19];
crz(pi/1024) q[8],q[18];
crz(pi/512) q[8],q[17];
crz(pi/256) q[8],q[16];
crz(pi/128) q[8],q[15];
crz(pi/64) q[8],q[14];
crz(pi/32) q[8],q[13];
crz(pi/16) q[8],q[12];
crz(pi/8) q[8],q[11];
crz(pi/4) q[8],q[10];
crz(pi/2) q[8],q[9];
h q[8];
crz(pi/32768) q[7],q[22];
crz(pi/16384) q[7],q[21];
crz(pi/8192) q[7],q[20];
crz(pi/4096) q[7],q[19];
crz(pi/2048) q[7],q[18];
crz(pi/1024) q[7],q[17];
crz(pi/512) q[7],q[16];
crz(pi/256) q[7],q[15];
crz(pi/128) q[7],q[14];
crz(pi/64) q[7],q[13];
crz(pi/32) q[7],q[12];
crz(pi/16) q[7],q[11];
crz(pi/8) q[7],q[10];
crz(pi/4) q[7],q[9];
crz(pi/2) q[7],q[8];
h q[7];
crz(pi/65536) q[6],q[22];
crz(pi/32768) q[6],q[21];
crz(pi/16384) q[6],q[20];
crz(pi/8192) q[6],q[19];
crz(pi/4096) q[6],q[18];
crz(pi/2048) q[6],q[17];
crz(pi/1024) q[6],q[16];
crz(pi/512) q[6],q[15];
crz(pi/256) q[6],q[14];
crz(pi/128) q[6],q[13];
crz(pi/64) q[6],q[12];
crz(pi/32) q[6],q[11];
crz(pi/16) q[6],q[10];
crz(pi/8) q[6],q[9];
crz(pi/4) q[6],q[8];
crz(pi/2) q[6],q[7];
h q[6];
crz(pi/131072) q[5],q[22];
crz(pi/65536) q[5],q[21];
crz(pi/32768) q[5],q[20];
crz(pi/16384) q[5],q[19];
crz(pi/8192) q[5],q[18];
crz(pi/4096) q[5],q[17];
crz(pi/2048) q[5],q[16];
crz(pi/1024) q[5],q[15];
crz(pi/512) q[5],q[14];
crz(pi/256) q[5],q[13];
crz(pi/128) q[5],q[12];
crz(pi/64) q[5],q[11];
crz(pi/32) q[5],q[10];
crz(pi/16) q[5],q[9];
crz(pi/8) q[5],q[8];
crz(pi/4) q[5],q[7];
crz(pi/2) q[5],q[6];
h q[5];
crz(pi/262144) q[4],q[22];
crz(pi/131072) q[4],q[21];
crz(pi/65536) q[4],q[20];
crz(pi/32768) q[4],q[19];
crz(pi/16384) q[4],q[18];
crz(pi/8192) q[4],q[17];
crz(pi/4096) q[4],q[16];
crz(pi/2048) q[4],q[15];
crz(pi/1024) q[4],q[14];
crz(pi/512) q[4],q[13];
crz(pi/256) q[4],q[12];
crz(pi/128) q[4],q[11];
crz(pi/64) q[4],q[10];
crz(pi/32) q[4],q[9];
crz(pi/16) q[4],q[8];
crz(pi/8) q[4],q[7];
crz(pi/4) q[4],q[6];
crz(pi/2) q[4],q[5];
h q[4];
crz(pi/524288) q[3],q[22];
crz(pi/262144) q[3],q[21];
crz(pi/131072) q[3],q[20];
crz(pi/65536) q[3],q[19];
crz(pi/32768) q[3],q[18];
crz(pi/16384) q[3],q[17];
crz(pi/8192) q[3],q[16];
crz(pi/4096) q[3],q[15];
crz(pi/2048) q[3],q[14];
crz(pi/1024) q[3],q[13];
crz(pi/512) q[3],q[12];
crz(pi/256) q[3],q[11];
crz(pi/128) q[3],q[10];
crz(pi/64) q[3],q[9];
crz(pi/32) q[3],q[8];
crz(pi/16) q[3],q[7];
crz(pi/8) q[3],q[6];
crz(pi/4) q[3],q[5];
crz(pi/2) q[3],q[4];
h q[3];
crz(pi/1048576) q[2],q[22];
crz(pi/524288) q[2],q[21];
crz(pi/262144) q[2],q[20];
crz(pi/131072) q[2],q[19];
crz(pi/65536) q[2],q[18];
crz(pi/32768) q[2],q[17];
crz(pi/16384) q[2],q[16];
crz(pi/8192) q[2],q[15];
crz(pi/4096) q[2],q[14];
crz(pi/2048) q[2],q[13];
crz(pi/1024) q[2],q[12];
crz(pi/512) q[2],q[11];
crz(pi/256) q[2],q[10];
crz(pi/128) q[2],q[9];
crz(pi/64) q[2],q[8];
crz(pi/32) q[2],q[7];
crz(pi/16) q[2],q[6];
crz(pi/8) q[2],q[5];
crz(pi/4) q[2],q[4];
crz(pi/2) q[2],q[3];
h q[2];
crz(pi/2097152) q[1],q[22];
crz(pi/1048576) q[1],q[21];
crz(pi/524288) q[1],q[20];
crz(pi/262144) q[1],q[19];
crz(pi/131072) q[1],q[18];
crz(pi/65536) q[1],q[17];
crz(pi/32768) q[1],q[16];
crz(pi/16384) q[1],q[15];
crz(pi/8192) q[1],q[14];
crz(pi/4096) q[1],q[13];
crz(pi/2048) q[1],q[12];
crz(pi/1024) q[1],q[11];
crz(pi/512) q[1],q[10];
crz(pi/256) q[1],q[9];
crz(pi/128) q[1],q[8];
crz(pi/64) q[1],q[7];
crz(pi/32) q[1],q[6];
crz(pi/16) q[1],q[5];
crz(pi/8) q[1],q[4];
crz(pi/4) q[1],q[3];
crz(pi/2) q[1],q[2];
h q[1];
crz(pi/4194304) q[0],q[22];
crz(pi/2097152) q[0],q[21];
crz(pi/1048576) q[0],q[20];
crz(pi/524288) q[0],q[19];
crz(pi/262144) q[0],q[18];
crz(pi/131072) q[0],q[17];
crz(pi/65536) q[0],q[16];
crz(pi/32768) q[0],q[15];
crz(pi/16384) q[0],q[14];
crz(pi/8192) q[0],q[13];
crz(pi/4096) q[0],q[12];
crz(pi/2048) q[0],q[11];
crz(pi/1024) q[0],q[10];
crz(pi/512) q[0],q[9];
crz(pi/256) q[0],q[8];
crz(pi/128) q[0],q[7];
crz(pi/64) q[0],q[6];
crz(pi/32) q[0],q[5];
crz(pi/16) q[0],q[4];
crz(pi/8) q[0],q[3];
crz(pi/4) q[0],q[2];
crz(pi/2) q[0],q[1];
h q[0];
p(pi) q[0];
p(pi/2) q[1];
p(pi/1024) q[10];
p(pi/2048) q[11];
p(pi/4096) q[12];
p(pi/8192) q[13];
p(pi/16384) q[14];
p(pi/32768) q[15];
p(pi/65536) q[16];
p(pi/131072) q[17];
p(pi/262144) q[18];
p(pi/524288) q[19];
p(pi/4) q[2];
p(pi/1048576) q[20];
p(pi/2097152) q[21];
p(pi/4194304) q[22];
p(pi/8) q[3];
p(pi/16) q[4];
p(pi/32) q[5];
p(pi/64) q[6];
p(pi/128) q[7];
p(pi/256) q[8];
p(pi/512) q[9];
h q[0];
crz(-pi/2) q[0],q[1];
crz(-pi/4) q[0],q[2];
crz(-pi/8) q[0],q[3];
crz(-pi/16) q[0],q[4];
crz(-pi/32) q[0],q[5];
crz(-pi/64) q[0],q[6];
crz(-pi/128) q[0],q[7];
crz(-pi/256) q[0],q[8];
crz(-pi/512) q[0],q[9];
crz(-pi/1024) q[0],q[10];
crz(-pi/2048) q[0],q[11];
crz(-pi/4096) q[0],q[12];
crz(-pi/8192) q[0],q[13];
crz(-pi/16384) q[0],q[14];
crz(-pi/32768) q[0],q[15];
crz(-pi/65536) q[0],q[16];
crz(-pi/131072) q[0],q[17];
crz(-pi/262144) q[0],q[18];
crz(-pi/524288) q[0],q[19];
crz(-pi/1048576) q[0],q[20];
crz(-pi/2097152) q[0],q[21];
crz(-pi/4194304) q[0],q[22];
h q[1];
crz(-pi/2) q[1],q[2];
crz(-pi/4) q[1],q[3];
crz(-pi/8) q[1],q[4];
crz(-pi/16) q[1],q[5];
crz(-pi/32) q[1],q[6];
crz(-pi/64) q[1],q[7];
crz(-pi/128) q[1],q[8];
crz(-pi/256) q[1],q[9];
crz(-pi/512) q[1],q[10];
crz(-pi/1024) q[1],q[11];
crz(-pi/2048) q[1],q[12];
crz(-pi/4096) q[1],q[13];
crz(-pi/8192) q[1],q[14];
crz(-pi/16384) q[1],q[15];
crz(-pi/32768) q[1],q[16];
crz(-pi/65536) q[1],q[17];
crz(-pi/131072) q[1],q[18];
crz(-pi/262144) q[1],q[19];
crz(-pi/524288) q[1],q[20];
crz(-pi/1048576) q[1],q[21];
crz(-pi/2097152) q[1],q[22];
h q[2];
crz(-pi/2) q[2],q[3];
crz(-pi/4) q[2],q[4];
crz(-pi/8) q[2],q[5];
crz(-pi/16) q[2],q[6];
crz(-pi/32) q[2],q[7];
crz(-pi/64) q[2],q[8];
crz(-pi/128) q[2],q[9];
crz(-pi/256) q[2],q[10];
crz(-pi/512) q[2],q[11];
crz(-pi/1024) q[2],q[12];
crz(-pi/2048) q[2],q[13];
crz(-pi/4096) q[2],q[14];
crz(-pi/8192) q[2],q[15];
crz(-pi/16384) q[2],q[16];
crz(-pi/32768) q[2],q[17];
crz(-pi/65536) q[2],q[18];
crz(-pi/131072) q[2],q[19];
crz(-pi/262144) q[2],q[20];
crz(-pi/524288) q[2],q[21];
crz(-pi/1048576) q[2],q[22];
h q[3];
crz(-pi/2) q[3],q[4];
crz(-pi/4) q[3],q[5];
crz(-pi/8) q[3],q[6];
crz(-pi/16) q[3],q[7];
crz(-pi/32) q[3],q[8];
crz(-pi/64) q[3],q[9];
crz(-pi/128) q[3],q[10];
crz(-pi/256) q[3],q[11];
crz(-pi/512) q[3],q[12];
crz(-pi/1024) q[3],q[13];
crz(-pi/2048) q[3],q[14];
crz(-pi/4096) q[3],q[15];
crz(-pi/8192) q[3],q[16];
crz(-pi/16384) q[3],q[17];
crz(-pi/32768) q[3],q[18];
crz(-pi/65536) q[3],q[19];
crz(-pi/131072) q[3],q[20];
crz(-pi/262144) q[3],q[21];
crz(-pi/524288) q[3],q[22];
h q[4];
crz(-pi/2) q[4],q[5];
crz(-pi/4) q[4],q[6];
crz(-pi/8) q[4],q[7];
crz(-pi/16) q[4],q[8];
crz(-pi/32) q[4],q[9];
crz(-pi/64) q[4],q[10];
crz(-pi/128) q[4],q[11];
crz(-pi/256) q[4],q[12];
crz(-pi/512) q[4],q[13];
crz(-pi/1024) q[4],q[14];
crz(-pi/2048) q[4],q[15];
crz(-pi/4096) q[4],q[16];
crz(-pi/8192) q[4],q[17];
crz(-pi/16384) q[4],q[18];
crz(-pi/32768) q[4],q[19];
crz(-pi/65536) q[4],q[20];
crz(-pi/131072) q[4],q[21];
crz(-pi/262144) q[4],q[22];
h q[5];
crz(-pi/2) q[5],q[6];
crz(-pi/4) q[5],q[7];
crz(-pi/8) q[5],q[8];
crz(-pi/16) q[5],q[9];
crz(-pi/32) q[5],q[10];
crz(-pi/64) q[5],q[11];
crz(-pi/128) q[5],q[12];
crz(-pi/256) q[5],q[13];
crz(-pi/512) q[5],q[14];
crz(-pi/1024) q[5],q[15];
crz(-pi/2048) q[5],q[16];
crz(-pi/4096) q[5],q[17];
crz(-pi/8192) q[5],q[18];
crz(-pi/16384) q[5],q[19];
crz(-pi/32768) q[5],q[20];
crz(-pi/65536) q[5],q[21];
crz(-pi/131072) q[5],q[22];
h q[6];
crz(-pi/2) q[6],q[7];
crz(-pi/4) q[6],q[8];
crz(-pi/8) q[6],q[9];
crz(-pi/16) q[6],q[10];
crz(-pi/32) q[6],q[11];
crz(-pi/64) q[6],q[12];
crz(-pi/128) q[6],q[13];
crz(-pi/256) q[6],q[14];
crz(-pi/512) q[6],q[15];
crz(-pi/1024) q[6],q[16];
crz(-pi/2048) q[6],q[17];
crz(-pi/4096) q[6],q[18];
crz(-pi/8192) q[6],q[19];
crz(-pi/16384) q[6],q[20];
crz(-pi/32768) q[6],q[21];
crz(-pi/65536) q[6],q[22];
h q[7];
crz(-pi/2) q[7],q[8];
crz(-pi/4) q[7],q[9];
crz(-pi/8) q[7],q[10];
crz(-pi/16) q[7],q[11];
crz(-pi/32) q[7],q[12];
crz(-pi/64) q[7],q[13];
crz(-pi/128) q[7],q[14];
crz(-pi/256) q[7],q[15];
crz(-pi/512) q[7],q[16];
crz(-pi/1024) q[7],q[17];
crz(-pi/2048) q[7],q[18];
crz(-pi/4096) q[7],q[19];
crz(-pi/8192) q[7],q[20];
crz(-pi/16384) q[7],q[21];
crz(-pi/32768) q[7],q[22];
h q[8];
crz(-pi/2) q[8],q[9];
crz(-pi/4) q[8],q[10];
crz(-pi/8) q[8],q[11];
crz(-pi/16) q[8],q[12];
crz(-pi/32) q[8],q[13];
crz(-pi/64) q[8],q[14];
crz(-pi/128) q[8],q[15];
crz(-pi/256) q[8],q[16];
crz(-pi/512) q[8],q[17];
crz(-pi/1024) q[8],q[18];
crz(-pi/2048) q[8],q[19];
crz(-pi/4096) q[8],q[20];
crz(-pi/8192) q[8],q[21];
crz(-pi/16384) q[8],q[22];
h q[9];
crz(-pi/2) q[9],q[10];
crz(-pi/4) q[9],q[11];
crz(-pi/8) q[9],q[12];
crz(-pi/16) q[9],q[13];
crz(-pi/32) q[9],q[14];
crz(-pi/64) q[9],q[15];
crz(-pi/128) q[9],q[16];
crz(-pi/256) q[9],q[17];
crz(-pi/512) q[9],q[18];
crz(-pi/1024) q[9],q[19];
crz(-pi/2048) q[9],q[20];
crz(-pi/4096) q[9],q[21];
crz(-pi/8192) q[9],q[22];
h q[10];
crz(-pi/2) q[10],q[11];
crz(-pi/4) q[10],q[12];
crz(-pi/8) q[10],q[13];
crz(-pi/16) q[10],q[14];
crz(-pi/32) q[10],q[15];
crz(-pi/64) q[10],q[16];
crz(-pi/128) q[10],q[17];
crz(-pi/256) q[10],q[18];
crz(-pi/512) q[10],q[19];
crz(-pi/1024) q[10],q[20];
crz(-pi/2048) q[10],q[21];
crz(-pi/4096) q[10],q[22];
h q[11];
crz(-pi/2) q[11],q[12];
crz(-pi/4) q[11],q[13];
crz(-pi/8) q[11],q[14];
crz(-pi/16) q[11],q[15];
crz(-pi/32) q[11],q[16];
crz(-pi/64) q[11],q[17];
crz(-pi/128) q[11],q[18];
crz(-pi/256) q[11],q[19];
crz(-pi/512) q[11],q[20];
crz(-pi/1024) q[11],q[21];
crz(-pi/2048) q[11],q[22];
h q[12];
crz(-pi/2) q[12],q[13];
crz(-pi/4) q[12],q[14];
crz(-pi/8) q[12],q[15];
crz(-pi/16) q[12],q[16];
crz(-pi/32) q[12],q[17];
crz(-pi/64) q[12],q[18];
crz(-pi/128) q[12],q[19];
crz(-pi/256) q[12],q[20];
crz(-pi/512) q[12],q[21];
crz(-pi/1024) q[12],q[22];
h q[13];
crz(-pi/2) q[13],q[14];
crz(-pi/4) q[13],q[15];
crz(-pi/8) q[13],q[16];
crz(-pi/16) q[13],q[17];
crz(-pi/32) q[13],q[18];
crz(-pi/64) q[13],q[19];
crz(-pi/128) q[13],q[20];
crz(-pi/256) q[13],q[21];
crz(-pi/512) q[13],q[22];
h q[14];
crz(-pi/2) q[14],q[15];
crz(-pi/4) q[14],q[16];
crz(-pi/8) q[14],q[17];
crz(-pi/16) q[14],q[18];
crz(-pi/32) q[14],q[19];
crz(-pi/64) q[14],q[20];
crz(-pi/128) q[14],q[21];
crz(-pi/256) q[14],q[22];
h q[15];
crz(-pi/2) q[15],q[16];
crz(-pi/4) q[15],q[17];
crz(-pi/8) q[15],q[18];
crz(-pi/16) q[15],q[19];
crz(-pi/32) q[15],q[20];
crz(-pi/64) q[15],q[21];
crz(-pi/128) q[15],q[22];
h q[16];
crz(-pi/2) q[16],q[17];
crz(-pi/4) q[16],q[18];
crz(-pi/8) q[16],q[19];
crz(-pi/16) q[16],q[20];
crz(-pi/32) q[16],q[21];
crz(-pi/64) q[16],q[22];
h q[17];
crz(-pi/2) q[17],q[18];
crz(-pi/4) q[17],q[19];
crz(-pi/8) q[17],q[20];
crz(-pi/16) q[17],q[21];
crz(-pi/32) q[17],q[22];
h q[18];
crz(-pi/2) q[18],q[19];
crz(-pi/4) q[18],q[20];
crz(-pi/8) q[18],q[21];
crz(-pi/16) q[18],q[22];
h q[19];
crz(-pi/2) q[19],q[20];
crz(-pi/4) q[19],q[21];
crz(-pi/8) q[19],q[22];
h q[20];
crz(-pi/2) q[20],q[21];
crz(-pi/4) q[20],q[22];
h q[21];
crz(-pi/2) q[21],q[22];
h q[22];
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
measure q[18] -> c[18];
measure q[19] -> c[19];
measure q[20] -> c[20];
measure q[21] -> c[21];
measure q[22] -> c[22];
