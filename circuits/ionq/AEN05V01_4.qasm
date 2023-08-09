OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[3],q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[4];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[4],q[2];
u3(pi/2,0,0) q[4];
rzz(pi/2) q[4],q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[2],q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(-pi/2) q[2],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,0,0) q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[4];
rzz(pi/2) q[2],q[4];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[2],q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[4];
rzz(pi/2) q[2],q[4];
u3(pi/2,pi,-pi) q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
rzz(pi/2) q[2],q[4];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,1.7404423300887455,-1.7404423300887455) q[4];
u3(pi/2,2.787849320795582,-2.787849320795582) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.4036014808991815,-3.4036014808991815) q[4];
u3(pi/2,5.759795971091527,-5.759795971091527) q[4];
rzz(-pi/2) q[2],q[4];
u3(pi/2,1.0474069907068368,-1.0474069907068368) q[4];
u3(pi/2,3.4036014808991815,-3.4036014808991815) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,1.8328051541042854,-1.8328051541042854) q[4];
u3(pi/2,4.18899964429663,-4.18899964429663) q[4];
rzz(pi/2) q[2],q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,0.26200882730938874,-0.26200882730938874) q[4];
rzz(-pi/2) q[2],q[4];
u3(pi/2,4.9743978076940785,-4.9743978076940785) q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,0,0) q[3];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,pi,-pi) q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(-pi/2) q[2],q[3];
u3(pi/2,5.759795971091527,-5.759795971091527) q[4];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[4],q[1];
u3(pi/2,4.18899964429663,-4.18899964429663) q[4];
u3(pi/2,4.9743978076940785,-4.9743978076940785) q[4];
rzz(pi/2) q[1],q[4];
u3(pi/2,3.4036014808991815,-3.4036014808991815) q[4];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(-pi/2) q[1],q[3];
u3(pi/2,0,0) q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,4.18899964429663,-4.18899964429663) q[4];
rzz(pi/2) q[1],q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[1],q[3];
u3(pi/2,4.9743978076940785,-4.9743978076940785) q[4];
rzz(pi/2) q[1],q[4];
u3(pi/2,0.26200882730938874,-0.26200882730938874) q[4];
u3(pi/2,1.0474069907068368,-1.0474069907068368) q[4];
rzz(pi/2) q[1],q[4];
u3(pi/2,5.759795971091527,-5.759795971091527) q[4];
u3(pi/2,1.8328051541042854,-1.8328051541042854) q[4];
rzz(pi/2) q[1],q[4];
u3(pi/2,0,0) q[3];
u3(pi/2,4.358645647590479,-4.358645647590479) q[4];
u3(pi/2,5.406052638297316,-5.406052638297316) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,6.021176479870197,-6.021176479870197) q[4];
u3(pi/2,2.094185662882956,-2.094185662882956) q[4];
rzz(pi/2) q[1],q[4];
u3(pi/2,0.5233893360880595,-0.5233893360880595) q[4];
u3(pi/2,2.8795838262804043,-2.8795838262804043) q[4];
rzz(-pi/2) q[3],q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(pi/2) q[1],q[3];
u3(pi/2,1.308787499485508,-1.308787499485508) q[4];
rzz(-pi/2) q[1],q[4];
u3(pi/2,0.5233893360880595,-0.5233893360880595) q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,0,0) q[3];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(-pi/2) q[1],q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[1],q[3];
rzz(-pi/2) q[1],q[4];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,pi,-pi) q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(-pi/2) q[1],q[3];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,0.5233893360880595,-0.5233893360880595) q[4];
rzz(pi/2) q[1],q[4];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[1],q[3];
u3(pi/2,1.308787499485508,-1.308787499485508) q[4];
rzz(-pi/2) q[1],q[4];
u3(pi/2,6.021176479870197,-6.021176479870197) q[4];
u3(pi/2,0.5233893360880595,-0.5233893360880595) q[4];
rzz(pi/2) q[1],q[4];
u3(pi/2,5.235778316472749,-5.235778316472749) q[4];
u3(pi/2,1.308787499485508,-1.308787499485508) q[4];
rzz(pi/2) q[1],q[4];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,3.8352563115024196,-3.8352563115024196) q[4];
u3(pi/2,4.882034983678539,-4.882034983678539) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
u3(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[1],q[4];
u3(pi/2,pi,-pi) q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
u3(pi/2,0,0) q[4];
rzz(pi/2) q[1],q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(pi/2) q[1],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[1],q[4];
u3(pi/2,pi/4,-pi/4) q[4];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(pi/2) q[1],q[3];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(pi/2) q[1],q[3];
u3(pi/2,pi/2,-pi/2) q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[0];
rzz(-pi/2) q[4],q[0];
u3(pi/2,pi,-pi) q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,0,0) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
u3(pi/2,pi,-pi) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,pi/2,-pi/2) q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,3.311238656883642,-3.311238656883642) q[4];
u3(pi/2,4.358645647590479,-4.358645647590479) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,4.9743978076940785,-4.9743978076940785) q[4];
u3(pi/2,1.0474069907068368,-1.0474069907068368) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,5.759795971091527,-5.759795971091527) q[4];
u3(pi/2,1.8328051541042854,-1.8328051541042854) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,3.4036014808991815,-3.4036014808991815) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,2.6182033175017336,-2.6182033175017336) q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[0],q[3];
rzz(-pi/2) q[0],q[4];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,2.6182033175017336,-2.6182033175017336) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3.4036014808991815,-3.4036014808991815) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,4.9743978076940785,-4.9743978076940785) q[4];
u3(pi/2,5.759795971091527,-5.759795971091527) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,4.18899964429663,-4.18899964429663) q[4];
u3(pi/2,0.26200882730938874,-0.26200882730938874) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,0,0) q[3];
u3(pi/2,2.787849320795582,-2.787849320795582) q[4];
u3(pi/2,3.8352563115024196,-3.8352563115024196) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,4.450380153075301,-4.450380153075301) q[4];
u3(pi/2,0.5233893360880595,-0.5233893360880595) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,5.235778316472749,-5.235778316472749) q[4];
u3(pi/2,1.308787499485508,-1.308787499485508) q[4];
rzz(-pi/2) q[3],q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,6.021176479870197,-6.021176479870197) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,5.235778316472749,-5.235778316472749) q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,0,0) q[3];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[0],q[3];
rzz(-pi/2) q[0],q[4];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,pi,-pi) q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,5.235778316472749,-5.235778316472749) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,6.021176479870197,-6.021176479870197) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,4.450380153075301,-4.450380153075301) q[4];
u3(pi/2,5.235778316472749,-5.235778316472749) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,3.664981989677853,-3.664981989677853) q[4];
u3(pi/2,6.021176479870197,-6.021176479870197) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,2.264459984707523,-2.264459984707523) q[4];
u3(pi/2,3.311238656883642,-3.311238656883642) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
u3(pi/2,0,0) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[4];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[3],q[4];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[0],q[3];
rzz(-pi/2) q[0],q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,pi/4,-pi/4) q[4];
u3(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,0,0) q[4];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,0,0) q[3];
u3(pi/2,4.882034983678539,-4.882034983678539) q[4];
u3(pi/2,5.929441974385376,-5.929441974385376) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,0.26200882730938874,-0.26200882730938874) q[4];
u3(pi/2,2.6182033175017336,-2.6182033175017336) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,4.18899964429663,-4.18899964429663) q[4];
u3(pi/2,0.26200882730938874,-0.26200882730938874) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,4.9743978076940785,-4.9743978076940785) q[4];
u3(pi/2,1.0474069907068368,-1.0474069907068368) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3.4036014808991815,-3.4036014808991815) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,1.8328051541042854,-1.8328051541042854) q[4];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,pi,-pi) q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,0,0) q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,pi/2,-pi/2) q[0];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi,-pi) q[1];
u3(pi/2,7*pi/4,-7*pi/4) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
u3(pi/2,11*pi/8,-11*pi/8) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
rzz(-pi/2) q[1],q[2];
u3(pi/2,3*pi/8,-3*pi/8) q[2];
u3(pi/2,9*pi/8,-9*pi/8) q[2];
rzz(pi/2) q[1],q[2];
u3(pi,0,0) q[0];
u3(pi,pi,-pi) q[1];
u3(pi/2,5*pi/8,-5*pi/8) q[2];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi,1.0474069907068368,-1.0474069907068368) q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
