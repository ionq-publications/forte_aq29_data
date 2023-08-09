OPENQASM 2.0;
include "qelib1.inc";
qreg q[9];
creg c[9];
u3(-pi,-pi/2,pi/2) q[0];
u3(-pi,-pi/2,pi/2) q[2];
u3(-pi,-pi/2,pi/2) q[5];
u3(pi/2,pi/4,-pi) q[8];
cx q[7],q[8];
u3(0,0,-pi/4) q[8];
cx q[7],q[8];
u3(pi/2,pi/4,-pi) q[7];
u3(0,0,pi/8) q[8];
cx q[6],q[8];
u3(0,0,-pi/8) q[8];
cx q[6],q[8];
cx q[6],q[7];
u3(0,0,-pi/4) q[7];
cx q[6],q[7];
u3(pi/2,pi/4,-pi) q[6];
u3(0,0,pi/8) q[7];
u3(0,0,pi/16) q[8];
cx q[5],q[8];
u3(0,0,-pi/16) q[8];
cx q[5],q[8];
cx q[5],q[7];
u3(0,0,-pi/8) q[7];
cx q[5],q[7];
cx q[5],q[6];
u3(0,0,-pi/4) q[6];
cx q[5],q[6];
u3(pi/2,pi/4,-pi) q[5];
u3(0,0,pi/8) q[6];
u3(0,0,pi/16) q[7];
u3(0,0,pi/32) q[8];
cx q[4],q[8];
u3(0,0,-pi/32) q[8];
cx q[4],q[8];
cx q[4],q[7];
u3(0,0,-pi/16) q[7];
cx q[4],q[7];
cx q[4],q[6];
u3(0,0,-pi/8) q[6];
cx q[4],q[6];
cx q[4],q[5];
u3(0,0,-pi/4) q[5];
cx q[4],q[5];
u3(pi/2,pi/4,-pi) q[4];
u3(0,0,pi/8) q[5];
u3(0,0,pi/16) q[6];
u3(0,0,pi/32) q[7];
u3(0,0,pi/64) q[8];
cx q[3],q[8];
u3(0,0,-pi/64) q[8];
cx q[3],q[8];
cx q[3],q[7];
u3(0,0,-pi/32) q[7];
cx q[3],q[7];
cx q[3],q[6];
u3(0,0,-pi/16) q[6];
cx q[3],q[6];
cx q[3],q[5];
u3(0,0,-pi/8) q[5];
cx q[3],q[5];
cx q[3],q[4];
u3(0,0,-pi/4) q[4];
cx q[3],q[4];
u3(pi/2,pi/4,-pi) q[3];
u3(0,0,pi/8) q[4];
u3(0,0,pi/16) q[5];
u3(0,0,pi/32) q[6];
u3(0,0,pi/64) q[7];
u3(0,0,pi/128) q[8];
cx q[2],q[8];
u3(0,0,-pi/128) q[8];
cx q[2],q[8];
cx q[2],q[7];
u3(0,0,-pi/64) q[7];
cx q[2],q[7];
cx q[2],q[6];
u3(0,0,-pi/32) q[6];
cx q[2],q[6];
cx q[2],q[5];
u3(0,0,-pi/16) q[5];
cx q[2],q[5];
cx q[2],q[4];
u3(0,0,-pi/8) q[4];
cx q[2],q[4];
cx q[2],q[3];
u3(0,0,-pi/4) q[3];
cx q[2],q[3];
u3(pi/2,pi/4,-pi) q[2];
u3(0,0,pi/8) q[3];
u3(0,0,pi/16) q[4];
u3(0,0,pi/32) q[5];
u3(0,0,pi/64) q[6];
u3(0,0,pi/128) q[7];
u3(0,0,pi/256) q[8];
cx q[1],q[8];
u3(0,0,-pi/256) q[8];
cx q[1],q[8];
cx q[1],q[7];
u3(0,0,-pi/128) q[7];
cx q[1],q[7];
cx q[1],q[6];
u3(0,0,-pi/64) q[6];
cx q[1],q[6];
cx q[1],q[5];
u3(0,0,-pi/32) q[5];
cx q[1],q[5];
cx q[1],q[4];
u3(0,0,-pi/16) q[4];
cx q[1],q[4];
cx q[1],q[3];
u3(0,0,-pi/8) q[3];
cx q[1],q[3];
cx q[1],q[2];
u3(0,0,-pi/4) q[2];
cx q[1],q[2];
u3(pi/2,pi/4,-pi) q[1];
u3(0,0,pi/8) q[2];
u3(0,0,pi/16) q[3];
u3(0,0,pi/32) q[4];
u3(0,0,pi/64) q[5];
u3(0,0,pi/128) q[6];
u3(0,0,pi/256) q[7];
u3(0,0,pi/512) q[8];
cx q[0],q[8];
u3(0,0,-0.00613592315154255) q[8];
cx q[0],q[8];
cx q[0],q[7];
u3(0,0,-pi/256) q[7];
cx q[0],q[7];
cx q[0],q[6];
u3(0,0,-pi/128) q[6];
cx q[0],q[6];
cx q[0],q[5];
u3(0,0,-pi/64) q[5];
cx q[0],q[5];
cx q[0],q[4];
u3(0,0,-pi/32) q[4];
cx q[0],q[4];
cx q[0],q[3];
u3(0,0,-pi/16) q[3];
cx q[0],q[3];
cx q[0],q[2];
u3(0,0,-pi/8) q[2];
cx q[0],q[2];
cx q[0],q[1];
u3(0,0,-pi/4) q[1];
cx q[0],q[1];
u3(pi,-0.5990270408223095,2.5425656127674836) q[0];
u3(0,pi/8,pi/8) q[1];
cx q[0],q[1];
u3(0,0,pi/4) q[1];
cx q[0],q[1];
u3(pi/2,0,pi) q[1];
u3(0,pi/16,pi/16) q[2];
cx q[0],q[2];
u3(0,0,pi/8) q[2];
cx q[0],q[2];
u3(0,0,-pi/4) q[2];
cx q[1],q[2];
u3(0,0,pi/4) q[2];
cx q[1],q[2];
u3(pi/2,0,pi) q[2];
u3(0,pi/32,pi/32) q[3];
cx q[0],q[3];
u3(0,0,pi/16) q[3];
cx q[0],q[3];
u3(0,0,-pi/8) q[3];
cx q[1],q[3];
u3(0,0,pi/8) q[3];
cx q[1],q[3];
u3(0,0,-pi/4) q[3];
cx q[2],q[3];
u3(0,0,pi/4) q[3];
cx q[2],q[3];
u3(pi/2,0,pi) q[3];
u3(0,pi/64,pi/64) q[4];
cx q[0],q[4];
u3(0,0,pi/32) q[4];
cx q[0],q[4];
u3(0,0,-pi/16) q[4];
cx q[1],q[4];
u3(0,0,pi/16) q[4];
cx q[1],q[4];
u3(0,0,-pi/8) q[4];
cx q[2],q[4];
u3(0,0,pi/8) q[4];
cx q[2],q[4];
u3(0,0,-pi/4) q[4];
cx q[3],q[4];
u3(0,0,pi/4) q[4];
cx q[3],q[4];
u3(pi/2,0,pi) q[4];
u3(0,pi/128,pi/128) q[5];
cx q[0],q[5];
u3(0,0,pi/64) q[5];
cx q[0],q[5];
u3(0,0,-pi/32) q[5];
cx q[1],q[5];
u3(0,0,pi/32) q[5];
cx q[1],q[5];
u3(0,0,-pi/16) q[5];
cx q[2],q[5];
u3(0,0,pi/16) q[5];
cx q[2],q[5];
u3(0,0,-pi/8) q[5];
cx q[3],q[5];
u3(0,0,pi/8) q[5];
cx q[3],q[5];
u3(0,0,-pi/4) q[5];
cx q[4],q[5];
u3(0,0,pi/4) q[5];
cx q[4],q[5];
u3(pi/2,0,pi) q[5];
u3(0,pi/256,pi/256) q[6];
cx q[0],q[6];
u3(0,0,pi/128) q[6];
cx q[0],q[6];
u3(0,0,-pi/64) q[6];
cx q[1],q[6];
u3(0,0,pi/64) q[6];
cx q[1],q[6];
u3(0,0,-pi/32) q[6];
cx q[2],q[6];
u3(0,0,pi/32) q[6];
cx q[2],q[6];
u3(0,0,-pi/16) q[6];
cx q[3],q[6];
u3(0,0,pi/16) q[6];
cx q[3],q[6];
u3(0,0,-pi/8) q[6];
cx q[4],q[6];
u3(0,0,pi/8) q[6];
cx q[4],q[6];
u3(0,0,-pi/4) q[6];
cx q[5],q[6];
u3(0,0,pi/4) q[6];
cx q[5],q[6];
u3(pi/2,0,pi) q[6];
u3(0,0.006135923151542766,0.006135923151542766) q[7];
cx q[0],q[7];
u3(0,0,pi/256) q[7];
cx q[0],q[7];
u3(0,0,-pi/128) q[7];
cx q[1],q[7];
u3(0,0,pi/128) q[7];
cx q[1],q[7];
u3(0,0,-pi/64) q[7];
cx q[2],q[7];
u3(0,0,pi/64) q[7];
cx q[2],q[7];
u3(0,0,-pi/32) q[7];
cx q[3],q[7];
u3(0,0,pi/32) q[7];
cx q[3],q[7];
u3(0,0,-pi/16) q[7];
cx q[4],q[7];
u3(0,0,pi/16) q[7];
cx q[4],q[7];
u3(0,0,-pi/8) q[7];
cx q[5],q[7];
u3(0,0,pi/8) q[7];
cx q[5],q[7];
u3(0,0,-pi/4) q[7];
cx q[6],q[7];
u3(0,0,pi/4) q[7];
cx q[6],q[7];
u3(pi/2,0,pi) q[7];
u3(0,0.003067961575771161,0.003067961575771161) q[8];
cx q[0],q[8];
u3(0,0,0.00613592315154255) q[8];
cx q[0],q[8];
u3(0,0,-pi/256) q[8];
cx q[1],q[8];
u3(0,0,pi/256) q[8];
cx q[1],q[8];
u3(0,0,-pi/128) q[8];
cx q[2],q[8];
u3(0,0,pi/128) q[8];
cx q[2],q[8];
u3(0,0,-pi/64) q[8];
cx q[3],q[8];
u3(0,0,pi/64) q[8];
cx q[3],q[8];
u3(0,0,-pi/32) q[8];
cx q[4],q[8];
u3(0,0,pi/32) q[8];
cx q[4],q[8];
u3(0,0,-pi/16) q[8];
cx q[5],q[8];
u3(0,0,pi/16) q[8];
cx q[5],q[8];
u3(0,0,-pi/8) q[8];
cx q[6],q[8];
u3(0,0,pi/8) q[8];
cx q[6],q[8];
u3(0,0,-pi/4) q[8];
cx q[7],q[8];
u3(0,0,pi/4) q[8];
cx q[7],q[8];
u3(pi/2,0,pi) q[8];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
measure q[8] -> c[8];
