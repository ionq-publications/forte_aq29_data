OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
u3(pi/2,4.417707589477967,-4.417707589477967) q[7];
rzz(pi/2) q[6],q[7];
u3(pi/2,5.9885039162728635,-5.9885039162728635) q[7];
u3(pi/2,0.7363893180014475,-0.7363893180014475) q[7];
rzz(pi/2) q[6],q[7];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(pi/2) q[5],q[7];
u3(pi/2,3.8779819715912405,-3.8779819715912405) q[7];
u3(pi/2,4.9580615258954115,-4.9580615258954115) q[7];
rzz(-pi/2) q[5],q[7];
u3(pi/2,5.693822525366141,-5.693822525366141) q[4];
rzz(pi/2) q[4],q[7];
u3(pi/2,4.9580615258954115,-4.9580615258954115) q[7];
u3(pi/2,5.9394950708768635,-5.9394950708768635) q[7];
rzz(pi/2) q[4],q[7];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(pi/2) q[3],q[7];
u3(pi/2,2.7979024172870695,-2.7979024172870695) q[7];
u3(pi/2,3.9759996623832423,-3.9759996623832423) q[7];
rzz(pi/2) q[3],q[7];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
rzz(-pi/2) q[2],q[7];
u3(pi/2,3.9759996623832423,-3.9759996623832423) q[7];
u3(pi/2,4.761397825780691,-4.761397825780691) q[7];
rzz(pi/2) q[2],q[7];
u3(pi/2,pi/2,-pi/2) q[1];
u3(pi/2,3.190601498985794,-3.190601498985794) q[7];
rzz(pi/2) q[1],q[7];
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
u3(pi/2,11*pi/8,-11*pi/8) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,7*pi/8,-7*pi/8) q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,2.552229871776348,-2.552229871776348) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,4.123026198571244,-4.123026198571244) q[4];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,5.399769452990137,-5.399769452990137) q[5];
u3(pi/2,2.2091679540043425,-2.2091679540043425) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,5.301123443667417,-5.301123443667417) q[6];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,pi/4,-pi/4) q[1];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,pi/8,-pi/8) q[2];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(-pi/2) q[1],q[3];
u3(pi/2,2.5528581903070657,-2.5528581903070657) q[3];
u3(pi/2,5.301751762198135,-5.301751762198135) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[4];
u3(pi/2,3.82897312619524,-3.82897312619524) q[4];
rzz(pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
u3(pi/2,2.110521944681623,-2.110521944681623) q[5];
rzz(pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[6];
u3(pi/2,5.227610175573416,-5.227610175573416) q[6];
rzz(pi/2) q[1],q[6];
u3(pi/2,11*pi/8,-11*pi/8) q[2];
u3(pi/2,0,0) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,5.301751762198135,-5.301751762198135) q[3];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,3.82897312619524,-3.82897312619524) q[4];
u3(pi/2,0.2946813909067226,-0.2946813909067226) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,2.110521944681623,-2.110521944681623) q[5];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[5];
rzz(-pi/2) q[2],q[5];
rzz(pi/2) q[2],q[6];
u3(pi/2,5.227610175573416,-5.227610175573416) q[6];
u3(pi/2,1.9879998311916212,-1.9879998311916212) q[6];
rzz(pi/2) q[2],q[6];
u3(pi/2,2.944928953475072,-2.944928953475072) q[3];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,0.2946813909067226,-0.2946813909067226) q[4];
u3(pi/2,2.6508758810990676,-2.6508758810990676) q[4];
rzz(-pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,1.91448656309762,-1.91448656309762) q[5];
u3(pi/2,4.663380134988689,-4.663380134988689) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,1.9879998311916212,-1.9879998311916212) q[6];
u3(pi/2,4.933557103197411,-4.933557103197411) q[6];
rzz(pi/2) q[3],q[6];
u3(pi/2,1.0800795543041708,-1.0800795543041708) q[4];
u3(pi/2,9*pi/8,-9*pi/8) q[4];
rzz(-pi/2) q[4],q[5];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[5];
u3(pi/2,3.8779819715912405,-3.8779819715912405) q[5];
rzz(pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,4.933557103197411,-4.933557103197411) q[6];
u3(pi/2,1.399265367908894,-1.399265367908894) q[6];
rzz(pi/2) q[4],q[6];
u3(pi/2,5.448778298386137,-5.448778298386137) q[5];
u3(pi/2,15*pi/8,-15*pi/8) q[5];
rzz(-pi/2) q[5],q[6];
u3(pi/2,4.540858021498687,-4.540858021498687) q[6];
u3(pi/2,0.6138672045114455,-0.6138672045114455) q[6];
rzz(pi/2) q[5],q[6];
u3(pi,pi,-pi) q[1];
u3(pi,pi,-pi) q[2];
u3(pi,3*pi/4,-3*pi/4) q[3];
u3(pi,7*pi/4,-7*pi/4) q[4];
u3(pi,4.516353598800687,-4.516353598800687) q[5];
u3(pi/2,5.326256184896136,-5.326256184896136) q[6];
u3(pi/2,15*pi/8,-15*pi/8) q[6];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
