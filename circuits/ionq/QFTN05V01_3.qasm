OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
u3(pi/2,0,0) q[4];
rzz(-pi/2) q[3],q[4];
rzz(-pi/2) q[2],q[4];
u3(pi/2,0,0) q[4];
u3(pi/2,7*pi/8,-7*pi/8) q[4];
rzz(-pi/2) q[2],q[4];
rzz(-pi/2) q[2],q[3];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(-pi/2) q[2],q[3];
rzz(-pi/2) q[1],q[4];
u3(pi/2,15*pi/8,-15*pi/8) q[4];
u3(pi/2,2.552229871776348,-2.552229871776348) q[4];
rzz(-pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
u3(pi/2,13*pi/8,-13*pi/8) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[2];
u3(pi/2,0,0) q[2];
u3(pi/2,pi/4,-pi/4) q[2];
rzz(-pi/2) q[1],q[2];
rzz(-pi/2) q[0],q[4];
rzz(pi/2) q[0],q[3];
rzz(-pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi/2,0,0) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,5*pi/8,-5*pi/8) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[4];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
rzz(pi/2) q[1],q[2];
rzz(-pi/2) q[1],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
u3(pi/2,5.301751762198135,-5.301751762198135) q[4];
rzz(pi/2) q[1],q[4];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,pi,-pi) q[2];
rzz(-pi/2) q[2],q[3];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
u3(pi/2,pi/8,-pi/8) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,5.301751762198135,-5.301751762198135) q[4];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[4];
rzz(pi/2) q[2],q[4];
u3(pi/2,13*pi/8,-13*pi/8) q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(-pi/2) q[3],q[4];
u3(pi/2,4.9084243619686925,-4.9084243619686925) q[4];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[4];
rzz(pi/2) q[3],q[4];
u3(pi,0,0) q[0];
u3(pi,0,0) q[1];
u3(pi,pi,-pi) q[2];
u3(pi,pi/4,-pi/4) q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[4];
u3(pi/2,15*pi/8,-15*pi/8) q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
