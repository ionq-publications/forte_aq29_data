OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
u3(pi/2,7*pi/4,-7*pi/4) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,pi/4,-pi/4) q[6];
u3(pi/2,2.8469112626830704,-2.8469112626830704) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[4],q[6];
u3(pi/2,2.8469112626830704,-2.8469112626830704) q[6];
u3(pi/2,3.82897312619524,-3.82897312619524) q[6];
rzz(-pi/2) q[4],q[6];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[3],q[6];
u3(pi/2,3.82897312619524,-3.82897312619524) q[6];
u3(pi/2,5.007070371291412,-5.007070371291412) q[6];
rzz(pi/2) q[3],q[6];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[2],q[6];
u3(pi/2,1.865477717701619,-1.865477717701619) q[6];
u3(pi/2,2.6508758810990676,-2.6508758810990676) q[6];
rzz(pi/2) q[2],q[6];
u3(pi/2,0,0) q[1];
u3(pi/2,1.0800795543041708,-1.0800795543041708) q[6];
rzz(-pi/2) q[1],q[6];
u3(pi/2,0,0) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
u3(pi/2,pi/4,-pi/4) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,15*pi/8,-15*pi/8) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,pi,-pi) q[3];
u3(pi/2,6.086521607064865,-6.086521607064865) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,3.82897312619524,-3.82897312619524) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,5.399769452990137,-5.399769452990137) q[5];
u3(pi/2,2.2091679540043425,-2.2091679540043425) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,7*pi/4,-7*pi/4) q[1];
u3(pi/2,0,0) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
u3(pi/2,13*pi/8,-13*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(-pi/2) q[1],q[3];
u3(pi/2,2.944928953475072,-2.944928953475072) q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[4];
u3(pi/2,5.203105752875415,-5.203105752875415) q[4];
rzz(pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
u3(pi/2,2.110521944681623,-2.110521944681623) q[5];
rzz(pi/2) q[1],q[5];
u3(pi/2,pi/8,-pi/8) q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[3];
rzz(pi/2) q[2],q[3];
rzz(-pi/2) q[2],q[4];
u3(pi/2,2.061513099285622,-2.061513099285622) q[4];
u3(pi/2,4.810406671176691,-4.810406671176691) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,2.110521944681623,-2.110521944681623) q[5];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[5];
rzz(pi/2) q[2],q[5];
u3(pi/2,3.3376280351737964,-3.3376280351737964) q[3];
u3(pi/2,9*pi/8,-9*pi/8) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,4.810406671176691,-4.810406671176691) q[4];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[5];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[5];
rzz(pi/2) q[3],q[5];
u3(pi/2,2.4542121809843462,-2.4542121809843462) q[4];
u3(pi/2,4.9084243619686925,-4.9084243619686925) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[5];
u3(pi/2,3.8779819715912405,-3.8779819715912405) q[5];
rzz(-pi/2) q[4],q[5];
u3(pi,pi/2,-pi/2) q[1];
u3(pi,7*pi/4,-7*pi/4) q[2];
u3(pi,15*pi/8,-15*pi/8) q[3];
u3(pi,5.301751762198135,-5.301751762198135) q[4];
u3(pi/2,2.3071856447963444,-2.3071856447963444) q[5];
u3(pi/2,7*pi/8,-7*pi/8) q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
