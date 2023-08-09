OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
u3(pi/2,7*pi/4,-7*pi/4) q[6];
rzz(-pi/2) q[5],q[6];
u3(pi/2,5*pi/4,-5*pi/4) q[6];
u3(pi/2,5.9885039162728635,-5.9885039162728635) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[4],q[6];
u3(pi/2,5.9885039162728635,-5.9885039162728635) q[6];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[6];
rzz(pi/2) q[4],q[6];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[3],q[6];
u3(pi/2,3.82897312619524,-3.82897312619524) q[6];
u3(pi/2,5.007070371291412,-5.007070371291412) q[6];
rzz(-pi/2) q[3],q[6];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[2],q[6];
u3(pi/2,5.007070371291412,-5.007070371291412) q[6];
u3(pi/2,5.792468534688861,-5.792468534688861) q[6];
rzz(pi/2) q[2],q[6];
u3(pi/2,0,0) q[1];
u3(pi/2,4.221672207893964,-4.221672207893964) q[6];
rzz(pi/2) q[1],q[6];
u3(pi/2,pi,-pi) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
u3(pi/2,pi/4,-pi/4) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,15*pi/8,-15*pi/8) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,0,0) q[3];
u3(pi/2,2.944928953475072,-2.944928953475072) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
u3(pi/2,5.399769452990137,-5.399769452990137) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,5.399769452990137,-5.399769452990137) q[5];
u3(pi/2,2.2091679540043425,-2.2091679540043425) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,7*pi/4,-7*pi/4) q[1];
u3(pi/2,0,0) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,15*pi/8,-15*pi/8) q[2];
u3(pi/2,5*pi/8,-5*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,6.087149925595583,-6.087149925595583) q[3];
u3(pi/2,2.552229871776348,-2.552229871776348) q[3];
rzz(pi/2) q[1],q[3];
rzz(-pi/2) q[1],q[4];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[4];
u3(pi/2,5.203105752875415,-5.203105752875415) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
u3(pi/2,2.110521944681623,-2.110521944681623) q[5];
rzz(pi/2) q[1],q[5];
u3(pi/2,9*pi/8,-9*pi/8) q[2];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
rzz(-pi/2) q[2],q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,5.203105752875415,-5.203105752875415) q[4];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[4];
rzz(pi/2) q[2],q[4];
rzz(-pi/2) q[2],q[5];
u3(pi/2,5.252114598271416,-5.252114598271416) q[5];
u3(pi/2,1.91448656309762,-1.91448656309762) q[5];
rzz(pi/2) q[2],q[5];
u3(pi/2,3.3376280351737964,-3.3376280351737964) q[3];
u3(pi/2,9*pi/8,-9*pi/8) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[4];
u3(pi/2,4.025008507779242,-4.025008507779242) q[4];
rzz(pi/2) q[3],q[4];
rzz(-pi/2) q[3],q[5];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[5];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[5];
rzz(pi/2) q[3],q[5];
u3(pi/2,5.595804834574139,-5.595804834574139) q[4];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[5];
u3(pi/2,3.8779819715912405,-3.8779819715912405) q[5];
rzz(pi/2) q[4],q[5];
u3(pi,3*pi/2,-3*pi/2) q[1];
u3(pi,7*pi/4,-7*pi/4) q[2];
u3(pi,7*pi/8,-7*pi/8) q[3];
u3(pi,0.5893627818134451,-0.5893627818134451) q[4];
u3(pi/2,5.448778298386137,-5.448778298386137) q[5];
u3(pi/2,15*pi/8,-15*pi/8) q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
