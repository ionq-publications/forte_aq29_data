OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(pi/2,pi/4,-pi/4) q[4];
u3(pi/2,11*pi/8,-11*pi/8) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,15*pi/8,-15*pi/8) q[5];
u3(pi/2,2.5528581903070657,-2.5528581903070657) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(pi/2) q[3],q[5];
u3(pi/2,2.552229871776348,-2.552229871776348) q[5];
u3(pi/2,5.301751762198135,-5.301751762198135) q[5];
rzz(-pi/2) q[3],q[5];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(pi/2) q[2],q[5];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[5];
u3(pi/2,4.516353598800687,-4.516353598800687) q[5];
rzz(pi/2) q[2],q[5];
u3(pi/2,pi/2,-pi/2) q[1];
u3(pi/2,6.087149925595583,-6.087149925595583) q[5];
rzz(pi/2) q[1],q[5];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,pi,-pi) q[1];
u3(pi/2,7*pi/4,-7*pi/4) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
u3(pi/2,pi/8,-pi/8) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,15*pi/8,-15*pi/8) q[3];
u3(pi/2,2.552229871776348,-2.552229871776348) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
u3(pi/2,5.399769452990137,-5.399769452990137) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[1];
u3(pi/2,0,0) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,pi/8,-pi/8) q[2];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(-pi/2) q[1],q[3];
u3(pi/2,2.552229871776348,-2.552229871776348) q[3];
u3(pi/2,5.301751762198135,-5.301751762198135) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,5.399769452990137,-5.399769452990137) q[4];
u3(pi/2,2.061513099285622,-2.061513099285622) q[4];
rzz(pi/2) q[1],q[4];
u3(pi/2,3*pi/8,-3*pi/8) q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(-pi/2) q[2],q[3];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[3];
u3(pi/2,4.516353598800687,-4.516353598800687) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,2.061513099285622,-2.061513099285622) q[4];
u3(pi/2,4.810406671176691,-4.810406671176691) q[4];
rzz(pi/2) q[2],q[4];
u3(pi/2,2.94555727200579,-2.94555727200579) q[3];
u3(pi/2,13*pi/8,-13*pi/8) q[3];
rzz(-pi/2) q[3],q[4];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[4];
u3(pi/2,4.025008507779242,-4.025008507779242) q[4];
rzz(pi/2) q[3],q[4];
u3(pi,3*pi/2,-3*pi/2) q[1];
u3(pi,3*pi/4,-3*pi/4) q[2];
u3(pi,pi/8,-pi/8) q[3];
u3(pi/2,5.595804834574139,-5.595804834574139) q[4];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[4];
u3(pi,3*pi/2,-3*pi/2) q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
