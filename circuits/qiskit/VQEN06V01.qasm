OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(-pi/2,0.0,0.0) q[0];
u3(pi/2,0,pi/2) q[1];
cx q[1],q[0];
u3(0,0,1.76406) q[0];
cx q[1],q[0];
u3(-pi/2,-pi/2,pi/2) q[0];
u3(pi/2,-pi/2,pi/2) q[1];
cx q[1],q[0];
u3(0,0,-1.76406) q[0];
cx q[1],q[0];
u3(pi/2,-pi/2,pi/2) q[0];
u3(pi/2,0,pi) q[1];
u3(pi/2,0,pi/2) q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,0.40016) q[0];
cx q[1],q[0];
u3(-pi/2,-pi/2,pi/2) q[0];
cx q[2],q[1];
u3(pi/2,-pi/2,pi/2) q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,-0.40016) q[0];
cx q[1],q[0];
u3(pi/2,-pi/2,pi/2) q[0];
cx q[2],q[1];
u3(pi/2,0,pi) q[1];
u3(-pi/2,0.0,0.0) q[3];
u3(pi/2,0,pi/2) q[4];
cx q[4],q[3];
u3(0,0,0.97874) q[3];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
u3(pi/2,-pi/2,pi/2) q[4];
cx q[4],q[3];
u3(0,0,-0.97874) q[3];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
u3(pi/2,0,pi) q[4];
u3(pi/2,0,pi/2) q[5];
cx q[5],q[4];
cx q[4],q[3];
u3(0,0,2.2409) q[3];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
cx q[5],q[4];
u3(pi/2,-pi/2,pi/2) q[5];
cx q[5],q[4];
cx q[4],q[3];
u3(0,0,-2.2409) q[3];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
cx q[5],q[4];
u3(pi/2,0,pi/2) q[4];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,-0.46688) q[0];
cx q[1],q[0];
cx q[3],q[1];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
u3(pi/2,-pi/2,pi/2) q[4];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,0.46688) q[0];
cx q[1],q[0];
cx q[3],q[1];
u3(-pi/2,-pi/2,pi/2) q[1];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,-0.46688) q[0];
cx q[1],q[0];
cx q[3],q[1];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
u3(-pi/2,-pi/2,pi/2) q[4];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,-0.46688) q[0];
cx q[1],q[0];
u3(-pi/2,-pi/2,pi/2) q[0];
cx q[3],q[1];
u3(pi/2,-pi/2,pi/2) q[1];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
u3(pi/2,-pi/2,pi/2) q[4];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,0.46688) q[0];
cx q[1],q[0];
cx q[3],q[1];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
u3(-pi/2,-pi/2,pi/2) q[4];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,0.46688) q[0];
cx q[1],q[0];
cx q[3],q[1];
u3(-pi/2,-pi/2,pi/2) q[1];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,-0.46688) q[0];
cx q[1],q[0];
cx q[3],q[1];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
u3(pi/2,-pi/2,pi/2) q[4];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,0.46688) q[0];
cx q[1],q[0];
u3(pi/2,-pi/2,pi/2) q[0];
cx q[3],q[1];
u3(pi/2,-pi/2,pi/2) q[1];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
u3(pi/2,0,pi) q[4];
u3(-pi/2,-pi/2,pi/2) q[5];
cx q[5],q[4];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,0.24432) q[0];
cx q[1],q[0];
cx q[3],q[1];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
cx q[5],q[4];
u3(pi/2,-pi/2,pi/2) q[5];
cx q[5],q[4];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,-0.24432) q[0];
cx q[1],q[0];
cx q[3],q[1];
u3(-pi/2,-pi/2,pi/2) q[1];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,0.24432) q[0];
cx q[1],q[0];
cx q[3],q[1];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
cx q[5],q[4];
u3(-pi/2,-pi/2,pi/2) q[5];
cx q[5],q[4];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,0.24432) q[0];
cx q[1],q[0];
u3(-pi/2,-pi/2,pi/2) q[0];
cx q[3],q[1];
u3(pi/2,-pi/2,pi/2) q[1];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
cx q[5],q[4];
u3(pi/2,-pi/2,pi/2) q[5];
cx q[5],q[4];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,-0.24432) q[0];
cx q[1],q[0];
cx q[3],q[1];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
cx q[5],q[4];
u3(-pi/2,-pi/2,pi/2) q[5];
cx q[5],q[4];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,-0.24432) q[0];
cx q[1],q[0];
cx q[3],q[1];
u3(-pi/2,-pi/2,pi/2) q[1];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,0.24432) q[0];
cx q[1],q[0];
cx q[3],q[1];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
cx q[5],q[4];
u3(pi/2,-pi/2,pi/2) q[5];
cx q[5],q[4];
cx q[4],q[3];
cx q[3],q[1];
cx q[1],q[0];
u3(0,0,-0.24432) q[0];
cx q[1],q[0];
u3(pi/2,-pi/2,pi/2) q[0];
cx q[3],q[1];
u3(pi/2,pi/2,-pi) q[1];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
cx q[5],q[4];
u3(pi/2,0,pi/2) q[4];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,-0.23752) q[0];
cx q[1],q[0];
cx q[2],q[1];
cx q[3],q[2];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
u3(pi/2,-pi/2,pi/2) q[4];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,0.23752) q[0];
cx q[1],q[0];
cx q[2],q[1];
cx q[3],q[2];
u3(-pi/2,-pi/2,pi/2) q[2];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,-0.23752) q[0];
cx q[1],q[0];
cx q[2],q[1];
cx q[3],q[2];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
u3(-pi/2,-pi/2,pi/2) q[4];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,-0.23752) q[0];
cx q[1],q[0];
u3(-pi/2,-pi/2,pi/2) q[0];
cx q[2],q[1];
cx q[3],q[2];
u3(pi/2,-pi/2,pi/2) q[2];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
u3(pi/2,-pi/2,pi/2) q[4];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,0.23752) q[0];
cx q[1],q[0];
cx q[2],q[1];
cx q[3],q[2];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
u3(-pi/2,-pi/2,pi/2) q[4];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,0.23752) q[0];
cx q[1],q[0];
cx q[2],q[1];
cx q[3],q[2];
u3(-pi/2,-pi/2,pi/2) q[2];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,-0.23752) q[0];
cx q[1],q[0];
cx q[2],q[1];
cx q[3],q[2];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
u3(pi/2,-pi/2,pi/2) q[4];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,0.23752) q[0];
cx q[1],q[0];
u3(pi/2,-pi/2,pi/2) q[0];
cx q[2],q[1];
cx q[3],q[2];
u3(pi/2,-pi/2,pi/2) q[2];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
u3(pi/2,0,pi) q[4];
u3(-pi/2,-pi/2,pi/2) q[5];
cx q[5],q[4];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,0.03784) q[0];
cx q[1],q[0];
cx q[2],q[1];
cx q[3],q[2];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
cx q[5],q[4];
u3(pi/2,-pi/2,pi/2) q[5];
cx q[5],q[4];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,-0.03784) q[0];
cx q[1],q[0];
cx q[2],q[1];
cx q[3],q[2];
u3(-pi/2,-pi/2,pi/2) q[2];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,0.03784) q[0];
cx q[1],q[0];
cx q[2],q[1];
cx q[3],q[2];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
cx q[5],q[4];
u3(-pi/2,-pi/2,pi/2) q[5];
cx q[5],q[4];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,0.03784) q[0];
cx q[1],q[0];
u3(-pi/2,-pi/2,pi/2) q[0];
cx q[2],q[1];
cx q[3],q[2];
u3(pi/2,-pi/2,pi/2) q[2];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
cx q[5],q[4];
u3(pi/2,-pi/2,pi/2) q[5];
cx q[5],q[4];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,-0.03784) q[0];
cx q[1],q[0];
cx q[2],q[1];
cx q[3],q[2];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
cx q[5],q[4];
u3(-pi/2,-pi/2,pi/2) q[5];
cx q[5],q[4];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,-0.03784) q[0];
cx q[1],q[0];
cx q[2],q[1];
cx q[3],q[2];
u3(-pi/2,-pi/2,pi/2) q[2];
cx q[4],q[3];
u3(pi/2,-pi/2,pi/2) q[3];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,0.03784) q[0];
cx q[1],q[0];
cx q[2],q[1];
cx q[3],q[2];
cx q[4],q[3];
u3(-pi/2,-pi/2,pi/2) q[3];
cx q[5],q[4];
u3(pi/2,-pi/2,pi/2) q[5];
cx q[5],q[4];
cx q[4],q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
u3(0,0,-0.03784) q[0];
cx q[1],q[0];
u3(pi/2,pi/2,-pi) q[0];
cx q[2],q[1];
cx q[3],q[2];
u3(pi/2,pi/2,-pi) q[2];
cx q[4],q[3];
u3(pi/2,pi/2,-pi) q[3];
cx q[5],q[4];
u3(pi/2,0,pi) q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
