OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
u3(pi/2,0.19352210746113127,-0.19352210746113127) q[0];
u3(pi/2,5.600831382819883,-5.600831382819883) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi,-pi) q[0];
u3(pi/2,4.030035056024986,-4.030035056024986) q[1];
u3(pi/2,0.6949202949740623,-0.6949202949740623) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,5.112627884452029,-5.112627884452029) q[2];
u3(pi/2,3.5418315576571326,-3.5418315576571326) q[3];
rzz(-pi/2) q[2],q[3];
u3(pi/2,5.112627884452029,-5.112627884452029) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,0,0) q[0];
rzz(-pi/2) q[2],q[0];
u3(pi/2,pi/2,-pi/2) q[0];
u3(pi/2,2.265716621768959,-2.265716621768959) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,3.8365129485638554,-3.8365129485638554) q[1];
u3(pi/2,0.4505043865247763,-0.4505043865247763) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[0];
rzz(-pi/2) q[2],q[0];
u3(pi/2,pi,-pi) q[2];
rzz(-pi/2) q[3],q[2];
rzz(pi/2) q[3],q[0];
u3(pi/2,3*pi/2,-3*pi/2) q[0];
rzz(pi/2) q[0],q[1];
u3(pi/2,0.4505043865247763,-0.4505043865247763) q[1];
u3(pi/2,3.347052813134565,-3.347052813134565) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,pi/2,-pi/2) q[0];
rzz(pi/2) q[3],q[0];
u3(pi/2,pi,-pi) q[0];
rzz(pi/2) q[0],q[1];
u3(pi/2,3.347052813134565,-3.347052813134565) q[1];
u3(pi/2,6.244229558275073,-6.244229558275073) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,pi,-pi) q[0];
rzz(-pi/2) q[3],q[0];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[3],q[0];
u3(pi/2,0,0) q[0];
rzz(-pi/2) q[0],q[1];
u3(pi/2,3.10263690468528,-3.10263690468528) q[1];
u3(pi/2,5.9991853312950685,-5.9991853312950685) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi,-pi) q[0];
rzz(pi/2) q[3],q[0];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[0];
u3(pi/2,pi/2,-pi/2) q[0];
u3(pi/2,2.857592677705276,-2.857592677705276) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,1.2867963509103792,-1.2867963509103792) q[1];
u3(pi/2,4.183973096050887,-4.183973096050887) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi/2,-pi/2) q[0];
rzz(-pi/2) q[2],q[0];
u3(pi/2,0,0) q[2];
rzz(-pi/2) q[2],q[3];
rzz(pi/2) q[3],q[0];
u3(pi/2,pi/2,-pi/2) q[0];
rzz(pi/2) q[0],q[1];
u3(pi/2,4.183973096050887,-4.183973096050887) q[1];
u3(pi/2,0.7973362154810896,-0.7973362154810896) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,0,0) q[3];
u3(pi/2,pi/2,-pi/2) q[0];
rzz(pi/2) q[3],q[0];
u3(pi/2,pi,-pi) q[0];
rzz(-pi/2) q[0],q[1];
u3(pi/2,0.7973362154810896,-0.7973362154810896) q[1];
u3(pi/2,3.6945129606215965,-3.6945129606215965) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,0,0) q[0];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,0,0) q[3];
rzz(-pi/2) q[3],q[0];
u3(pi/2,pi,-pi) q[0];
rzz(pi/2) q[0],q[1];
u3(pi/2,0.5529203070318035,-0.5529203070318035) q[1];
u3(pi/2,3.450097052172311,-3.450097052172311) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi,-pi) q[0];
rzz(pi/2) q[3],q[0];
u3(pi/2,0,0) q[2];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[2],q[3];
u3(pi,pi,-pi) q[0];
u3(pi,0.8224689567098078,-0.8224689567098078) q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
u3(pi,1.6600175581568466,-1.6600175581568466) q[3];
measure q[0] -> c[1];
measure q[1] -> c[0];
measure q[2] -> c[3];
measure q[3] -> c[2];
