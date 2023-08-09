OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(0,-1.40658297059163,1.4065829705916304) q[0];
u3(pi/2,-pi,-pi) q[1];
u3(pi/2,-pi/2,-pi) q[2];
u3(pi/2,3*pi/4,-pi) q[3];
u3(-pi,-pi/2,pi/2) q[4];
cx q[3],q[4];
u3(0,0,-3*pi/4) q[4];
cx q[3],q[4];
u3(0,0,-pi/16) q[3];
u3(0,0,3*pi/4) q[4];
cx q[2],q[4];
u3(0,0,-3*pi/2) q[4];
cx q[2],q[4];
u3(0,0,-pi/8) q[2];
u3(0,0,3*pi/2) q[4];
cx q[1],q[4];
u3(0,0,-3*pi) q[4];
cx q[1],q[4];
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
u3(0,0,-pi/8) q[3];
cx q[1],q[3];
u3(0,0,pi/8) q[3];
cx q[1],q[3];
u3(0,0,-pi/4) q[3];
cx q[2],q[3];
u3(0,0,pi/4) q[3];
cx q[2],q[3];
u3(pi/2,0,-pi) q[3];
u3(0,0,3*pi) q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
