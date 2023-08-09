OPENQASM 2.0;
include "qelib1.inc";
qreg q[21];
creg c[21];
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
ry(-pi/2) q[17];
rz(-pi) q[17];
ry(-pi/2) q[18];
rz(-pi) q[18];
ry(-pi/2) q[19];
rz(-pi) q[19];
rx(-pi) q[20];
p(218779*pi) q[0];
p(343657.24957986065) q[1];
p(671.2055655856653) q[10];
p(335.60278279283267) q[11];
p(167.80139139641634) q[12];
p(83.90069569820817) q[13];
p(41.950347849104084) q[14];
p(20.975173924552042) q[15];
p(10.487586962276021) q[16];
p(5.2437934811380105) q[17];
p(2.6218967405690052) q[18];
p(1.3109483702845026) q[19];
cx q[19],q[20];
p(171828.62478993033) q[2];
p(-1.3109483702845026) q[20];
cx q[19],q[20];
p(1.3109483702845026) q[20];
cx q[18],q[20];
p(-2.6218967405690052) q[20];
cx q[18],q[20];
p(2.6218967405690052) q[20];
cx q[17],q[20];
p(-5.2437934811380105) q[20];
cx q[17],q[20];
p(5.2437934811380105) q[20];
cx q[16],q[20];
p(-10.487586962276021) q[20];
cx q[16],q[20];
p(10.487586962276021) q[20];
cx q[15],q[20];
p(-20.975173924552042) q[20];
cx q[15],q[20];
p(20.975173924552042) q[20];
cx q[14],q[20];
p(-41.950347849104084) q[20];
cx q[14],q[20];
p(41.950347849104084) q[20];
cx q[13],q[20];
p(-83.90069569820817) q[20];
cx q[13],q[20];
p(83.90069569820817) q[20];
cx q[12],q[20];
p(-167.80139139641634) q[20];
cx q[12],q[20];
p(167.80139139641634) q[20];
cx q[11],q[20];
p(-335.60278279283267) q[20];
cx q[11],q[20];
p(335.60278279283267) q[20];
cx q[10],q[20];
p(-671.2055655856653) q[20];
cx q[10],q[20];
p(671.2055655856653) q[20];
p(85914.31239496516) q[3];
p(42957.15619748258) q[4];
p(21478.57809874129) q[5];
p(10739.289049370645) q[6];
p(5369.644524685323) q[7];
p(2684.8222623426614) q[8];
p(1342.4111311713307) q[9];
cx q[9],q[20];
p(-1342.4111311713307) q[20];
cx q[9],q[20];
p(1342.4111311713307) q[20];
cx q[8],q[20];
p(-2684.8222623426614) q[20];
cx q[8],q[20];
p(2684.8222623426614) q[20];
cx q[7],q[20];
p(-5369.644524685323) q[20];
cx q[7],q[20];
p(5369.644524685323) q[20];
cx q[6],q[20];
p(-10739.289049370645) q[20];
cx q[6],q[20];
p(10739.289049370645) q[20];
cx q[5],q[20];
p(-21478.57809874129) q[20];
cx q[5],q[20];
p(21478.57809874129) q[20];
cx q[4],q[20];
p(-42957.15619748258) q[20];
cx q[4],q[20];
p(42957.15619748258) q[20];
cx q[3],q[20];
p(-85914.31239496516) q[20];
cx q[3],q[20];
p(85914.31239496516) q[20];
cx q[2],q[20];
p(-171828.62478993033) q[20];
cx q[2],q[20];
p(171828.62478993033) q[20];
cx q[1],q[20];
p(-343657.24957986065) q[20];
cx q[1],q[20];
p(343657.24957986065) q[20];
cx q[0],q[20];
p(-218779*pi) q[20];
cx q[0],q[20];
p(218779*pi) q[20];
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
rz(-pi/262144) q[17];
rz(-pi/524288) q[18];
rz(-pi/1048576) q[19];
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
cx q[0],q[17];
rz(pi/262144) q[17];
cx q[0],q[17];
cx q[0],q[18];
rz(pi/524288) q[18];
cx q[0],q[18];
cx q[0],q[19];
rz(pi/1048576) q[19];
cx q[0],q[19];
ry(-pi/2) q[1];
rz(-pi) q[1];
rz(-pi/1024) q[10];
rz(-pi/2048) q[11];
rz(-pi/4096) q[12];
rz(-pi/8192) q[13];
rz(-pi/16384) q[14];
rz(-pi/32768) q[15];
rz(-pi/65536) q[16];
rz(-pi/131072) q[17];
rz(-pi/262144) q[18];
rz(-pi/524288) q[19];
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
cx q[1],q[17];
rz(pi/131072) q[17];
cx q[1],q[17];
cx q[1],q[18];
rz(pi/262144) q[18];
cx q[1],q[18];
cx q[1],q[19];
rz(pi/524288) q[19];
cx q[1],q[19];
rz(-pi/512) q[10];
rz(-pi/1024) q[11];
rz(-pi/2048) q[12];
rz(-pi/4096) q[13];
rz(-pi/8192) q[14];
rz(-pi/16384) q[15];
rz(-pi/32768) q[16];
rz(-pi/65536) q[17];
rz(-pi/131072) q[18];
rz(-pi/262144) q[19];
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
cx q[2],q[17];
rz(pi/65536) q[17];
cx q[2],q[17];
cx q[2],q[18];
rz(pi/131072) q[18];
cx q[2],q[18];
cx q[2],q[19];
rz(pi/262144) q[19];
cx q[2],q[19];
rz(-pi/256) q[10];
rz(-pi/512) q[11];
rz(-pi/1024) q[12];
rz(-pi/2048) q[13];
rz(-pi/4096) q[14];
rz(-pi/8192) q[15];
rz(-pi/16384) q[16];
rz(-pi/32768) q[17];
rz(-pi/65536) q[18];
rz(-pi/131072) q[19];
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
cx q[3],q[17];
rz(pi/32768) q[17];
cx q[3],q[17];
cx q[3],q[18];
rz(pi/65536) q[18];
cx q[3],q[18];
cx q[3],q[19];
rz(pi/131072) q[19];
cx q[3],q[19];
rz(-pi/128) q[10];
rz(-pi/256) q[11];
rz(-pi/512) q[12];
rz(-pi/1024) q[13];
rz(-pi/2048) q[14];
rz(-pi/4096) q[15];
rz(-pi/8192) q[16];
rz(-pi/16384) q[17];
rz(-pi/32768) q[18];
rz(-pi/65536) q[19];
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
cx q[4],q[17];
rz(pi/16384) q[17];
cx q[4],q[17];
cx q[4],q[18];
rz(pi/32768) q[18];
cx q[4],q[18];
cx q[4],q[19];
rz(pi/65536) q[19];
cx q[4],q[19];
rz(-pi/64) q[10];
rz(-pi/128) q[11];
rz(-pi/256) q[12];
rz(-pi/512) q[13];
rz(-pi/1024) q[14];
rz(-pi/2048) q[15];
rz(-pi/4096) q[16];
rz(-pi/8192) q[17];
rz(-pi/16384) q[18];
rz(-pi/32768) q[19];
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
cx q[5],q[17];
rz(pi/8192) q[17];
cx q[5],q[17];
cx q[5],q[18];
rz(pi/16384) q[18];
cx q[5],q[18];
cx q[5],q[19];
rz(pi/32768) q[19];
cx q[5],q[19];
rz(-pi/32) q[10];
rz(-pi/64) q[11];
rz(-pi/128) q[12];
rz(-pi/256) q[13];
rz(-pi/512) q[14];
rz(-pi/1024) q[15];
rz(-pi/2048) q[16];
rz(-pi/4096) q[17];
rz(-pi/8192) q[18];
rz(-pi/16384) q[19];
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
cx q[6],q[17];
rz(pi/4096) q[17];
cx q[6],q[17];
cx q[6],q[18];
rz(pi/8192) q[18];
cx q[6],q[18];
cx q[6],q[19];
rz(pi/16384) q[19];
cx q[6],q[19];
rz(-pi/16) q[10];
rz(-pi/32) q[11];
rz(-pi/64) q[12];
rz(-pi/128) q[13];
rz(-pi/256) q[14];
rz(-pi/512) q[15];
rz(-pi/1024) q[16];
rz(-pi/2048) q[17];
rz(-pi/4096) q[18];
rz(-pi/8192) q[19];
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
cx q[7],q[17];
rz(pi/2048) q[17];
cx q[7],q[17];
cx q[7],q[18];
rz(pi/4096) q[18];
cx q[7],q[18];
cx q[7],q[19];
rz(pi/8192) q[19];
cx q[7],q[19];
rz(-pi/8) q[10];
rz(-pi/16) q[11];
rz(-pi/32) q[12];
rz(-pi/64) q[13];
rz(-pi/128) q[14];
rz(-pi/256) q[15];
rz(-pi/512) q[16];
rz(-pi/1024) q[17];
rz(-pi/2048) q[18];
rz(-pi/4096) q[19];
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
cx q[8],q[17];
rz(pi/1024) q[17];
cx q[8],q[17];
cx q[8],q[18];
rz(pi/2048) q[18];
cx q[8],q[18];
cx q[8],q[19];
rz(pi/4096) q[19];
cx q[8],q[19];
rz(-pi/4) q[10];
rz(-pi/8) q[11];
rz(-pi/16) q[12];
rz(-pi/32) q[13];
rz(-pi/64) q[14];
rz(-pi/128) q[15];
rz(-pi/256) q[16];
rz(-pi/512) q[17];
rz(-pi/1024) q[18];
rz(-pi/2048) q[19];
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
cx q[9],q[17];
rz(pi/512) q[17];
cx q[9],q[17];
cx q[9],q[18];
rz(pi/1024) q[18];
cx q[9],q[18];
cx q[9],q[19];
rz(pi/2048) q[19];
cx q[9],q[19];
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
rz(-pi/256) q[17];
cx q[10],q[17];
rz(pi/256) q[17];
cx q[10],q[17];
rz(-pi/512) q[18];
cx q[10],q[18];
rz(pi/512) q[18];
cx q[10],q[18];
rz(-pi/1024) q[19];
cx q[10],q[19];
rz(pi/1024) q[19];
cx q[10],q[19];
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
rz(-pi/128) q[17];
cx q[11],q[17];
rz(pi/128) q[17];
cx q[11],q[17];
rz(-pi/256) q[18];
cx q[11],q[18];
rz(pi/256) q[18];
cx q[11],q[18];
rz(-pi/512) q[19];
cx q[11],q[19];
rz(pi/512) q[19];
cx q[11],q[19];
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
rz(-pi/64) q[17];
cx q[12],q[17];
rz(pi/64) q[17];
cx q[12],q[17];
rz(-pi/128) q[18];
cx q[12],q[18];
rz(pi/128) q[18];
cx q[12],q[18];
rz(-pi/256) q[19];
cx q[12],q[19];
rz(pi/256) q[19];
cx q[12],q[19];
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
rz(-pi/32) q[17];
cx q[13],q[17];
rz(pi/32) q[17];
cx q[13],q[17];
rz(-pi/64) q[18];
cx q[13],q[18];
rz(pi/64) q[18];
cx q[13],q[18];
rz(-pi/128) q[19];
cx q[13],q[19];
rz(pi/128) q[19];
cx q[13],q[19];
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
rz(-pi/16) q[17];
cx q[14],q[17];
rz(pi/16) q[17];
cx q[14],q[17];
rz(-pi/32) q[18];
cx q[14],q[18];
rz(pi/32) q[18];
cx q[14],q[18];
rz(-pi/64) q[19];
cx q[14],q[19];
rz(pi/64) q[19];
cx q[14],q[19];
ry(-pi/2) q[15];
rz(-pi) q[15];
rz(-pi/4) q[16];
cx q[15],q[16];
rz(pi/4) q[16];
cx q[15],q[16];
rz(-pi/8) q[17];
cx q[15],q[17];
rz(pi/8) q[17];
cx q[15],q[17];
rz(-pi/16) q[18];
cx q[15],q[18];
rz(pi/16) q[18];
cx q[15],q[18];
rz(-pi/32) q[19];
cx q[15],q[19];
rz(pi/32) q[19];
cx q[15],q[19];
ry(-pi/2) q[16];
rz(-pi) q[16];
rz(-pi/4) q[17];
cx q[16],q[17];
rz(pi/4) q[17];
cx q[16],q[17];
rz(-pi/8) q[18];
cx q[16],q[18];
rz(pi/8) q[18];
cx q[16],q[18];
rz(-pi/16) q[19];
cx q[16],q[19];
rz(pi/16) q[19];
cx q[16],q[19];
ry(-pi/2) q[17];
rz(-pi) q[17];
rz(-pi/4) q[18];
cx q[17],q[18];
rz(pi/4) q[18];
cx q[17],q[18];
rz(-pi/8) q[19];
cx q[17],q[19];
rz(pi/8) q[19];
cx q[17],q[19];
ry(-pi/2) q[18];
rz(-pi) q[18];
rz(-pi/4) q[19];
cx q[18],q[19];
rz(pi/4) q[19];
cx q[18],q[19];
ry(-pi/2) q[19];
rz(-pi) q[19];
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
