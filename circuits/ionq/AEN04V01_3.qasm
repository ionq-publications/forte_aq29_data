OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[3],q[1];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[3],q[1];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(-pi/2) q[1],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
u3(pi/2,0,0) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(-pi/2) q[1],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(-pi/2) q[1],q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(pi/2) q[1],q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[1],q[3];
u3(pi/2,pi,-pi) q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(-pi/2) q[1],q[3];
u3(pi/2,0,0) q[2];
u3(pi/2,4.882034983678539,-4.882034983678539) q[3];
u3(pi/2,5.929441974385376,-5.929441974385376) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,0.26200882730938874,-0.26200882730938874) q[3];
u3(pi/2,2.6182033175017336,-2.6182033175017336) q[3];
rzz(pi/2) q[1],q[3];
u3(pi/2,1.0474069907068368,-1.0474069907068368) q[3];
u3(pi/2,3.4036014808991815,-3.4036014808991815) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,1.8328051541042854,-1.8328051541042854) q[3];
u3(pi/2,4.18899964429663,-4.18899964429663) q[3];
rzz(pi/2) q[1],q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(-pi/2) q[1],q[2];
u3(pi/2,0.26200882730938874,-0.26200882730938874) q[3];
rzz(pi/2) q[1],q[3];
u3(pi/2,1.8328051541042854,-1.8328051541042854) q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
rzz(-pi/2) q[3],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
u3(pi/2,0,0) q[2];
rzz(-pi/2) q[3],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,5.759795971091527,-5.759795971091527) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[0];
rzz(pi/2) q[3],q[0];
u3(pi/2,4.18899964429663,-4.18899964429663) q[3];
u3(pi/2,4.9743978076940785,-4.9743978076940785) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3.4036014808991815,-3.4036014808991815) q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(-pi/2) q[3],q[2];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,1.0474069907068368,-1.0474069907068368) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,1.8328051541042854,-1.8328051541042854) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3.4036014808991815,-3.4036014808991815) q[3];
u3(pi/2,4.18899964429663,-4.18899964429663) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,2.6182033175017336,-2.6182033175017336) q[3];
u3(pi/2,4.9743978076940785,-4.9743978076940785) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
u3(pi/2,1.217052994000686,-1.217052994000686) q[3];
u3(pi/2,2.264459984707523,-2.264459984707523) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,2.8795838262804043,-2.8795838262804043) q[3];
u3(pi/2,5.235778316472749,-5.235778316472749) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3.664981989677853,-3.664981989677853) q[3];
u3(pi/2,6.021176479870197,-6.021176479870197) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,pi/4,-pi/4) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,1.308787499485508,-1.308787499485508) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,0.5233893360880595,-0.5233893360880595) q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[0],q[2];
rzz(pi/2) q[0],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
u3(pi/2,pi/4,-pi/4) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,0.5233893360880595,-0.5233893360880595) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,1.308787499485508,-1.308787499485508) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,2.8795838262804043,-2.8795838262804043) q[3];
u3(pi/2,3.664981989677853,-3.664981989677853) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,2.094185662882956,-2.094185662882956) q[3];
u3(pi/2,4.450380153075301,-4.450380153075301) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,0.6936636579126263,-0.6936636579126263) q[3];
u3(pi/2,1.7404423300887455,-1.7404423300887455) q[3];
rzz(-pi/2) q[2],q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,0,0) q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(-pi/2) q[3],q[2];
u3(pi/2,0,0) q[2];
u3(pi/2,pi/4,-pi/4) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(-pi/2) q[3],q[2];
u3(pi/2,0,0) q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,0,0) q[0];
u3(pi/2,0,0) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
u3(pi/2,pi/4,-pi/4) q[1];
rzz(pi/2) q[0],q[1];
u3(pi,0,0) q[0];
u3(pi/2,pi/4,-pi/4) q[1];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
u3(pi,pi/2,-pi/2) q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
