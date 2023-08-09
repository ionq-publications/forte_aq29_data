OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(0,-1.40658297059163,1.4065829705916304) q[0];
u3(pi/2,-pi,-pi) q[1];
u3(pi/2,pi/2,-pi) q[2];
u3(pi/2,-3*pi/4,-pi) q[3];
u3(pi/2,5*pi/8,-pi) q[4];
u3(pi/2,5*pi/16,-pi) q[5];
u3(-pi,-pi/2,pi/2) q[6];
cx q[5],q[6];
u3(0,0,-5*pi/16) q[6];
cx q[5],q[6];
u3(0,0,-pi/64) q[5];
u3(0,0,5*pi/16) q[6];
cx q[4],q[6];
u3(0,0,-5*pi/8) q[6];
cx q[4],q[6];
u3(0,0,-pi/32) q[4];
u3(0,0,5*pi/8) q[6];
cx q[3],q[6];
u3(0,0,-5*pi/4) q[6];
cx q[3],q[6];
u3(0,0,-pi/16) q[3];
u3(0,0,5*pi/4) q[6];
cx q[2],q[6];
u3(0,0,-5*pi/2) q[6];
cx q[2],q[6];
u3(0,0,-pi/8) q[2];
u3(0,0,5*pi/2) q[6];
cx q[1],q[6];
u3(0,0,-5*pi) q[6];
cx q[1],q[6];
u3(0,0,-pi/4) q[1];
cx q[0],q[1];
u3(0,0,pi/4) q[1];
cx q[0],q[1];
cx q[0],q[2];
u3(pi/2,0,-pi) q[1];
u3(0,0,pi/8) q[2];
cx q[0],q[2];
cx q[0],q[3];
u3(0,0,-pi/4) q[2];
cx q[1],q[2];
u3(0,0,pi/4) q[2];
cx q[1],q[2];
u3(pi/2,0,-pi) q[2];
u3(0,0,pi/16) q[3];
cx q[0],q[3];
cx q[0],q[4];
u3(0,0,-pi/8) q[3];
cx q[1],q[3];
u3(0,0,pi/8) q[3];
cx q[1],q[3];
u3(0,0,-pi/4) q[3];
cx q[2],q[3];
u3(0,0,pi/4) q[3];
cx q[2],q[3];
u3(pi/2,0,-pi) q[3];
u3(0,0,pi/32) q[4];
cx q[0],q[4];
cx q[0],q[5];
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
u3(pi/2,0,-pi) q[4];
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
u3(pi/2,0,-pi) q[5];
u3(0,0,5*pi) q[6];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
