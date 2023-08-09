OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[2],q[3];
rzz(-pi/2) q[1],q[3];
u3(pi/2,pi,-pi) q[3];
u3(pi/2,15*pi/8,-15*pi/8) q[3];
rzz(-pi/2) q[1],q[3];
rzz(pi/2) q[1],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(-pi/2) q[1],q[2];
rzz(pi/2) q[0],q[3];
rzz(-pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi/2,0,0) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,15*pi/8,-15*pi/8) q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
u3(pi/2,13*pi/8,-13*pi/8) q[3];
rzz(-pi/2) q[1],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,5*pi/8,-5*pi/8) q[3];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
rzz(pi/2) q[2],q[3];
u3(pi,0,0) q[0];
u3(pi,pi/2,-pi/2) q[1];
u3(pi,0,0) q[2];
u3(pi/2,7*pi/8,-7*pi/8) q[3];
u3(pi/2,pi,-pi) q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
