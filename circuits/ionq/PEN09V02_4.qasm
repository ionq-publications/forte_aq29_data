OPENQASM 2.0;
include "qelib1.inc";
qreg q[9];
creg c[9];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
u3(pi/2,9*pi/8,-9*pi/8) q[8];
rzz(-pi/2) q[7],q[8];
u3(pi/2,5*pi/8,-5*pi/8) q[8];
u3(pi/2,4.516353598800687,-4.516353598800687) q[8];
rzz(pi/2) q[7],q[8];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
rzz(pi/2) q[6],q[8];
u3(pi/2,4.516353598800687,-4.516353598800687) q[8];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[8];
rzz(pi/2) q[6],q[8];
u3(pi/2,5.301751762198135,-5.301751762198135) q[5];
rzz(pi/2) q[5],q[8];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[8];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[8];
rzz(-pi/2) q[5],q[8];
u3(pi/2,pi/4,-pi/4) q[4];
u3(pi/2,2.552229871776348,-2.552229871776348) q[8];
rzz(pi/2) q[4],q[8];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
u3(pi/2,7*pi/4,-7*pi/4) q[1];
rzz(-pi/2) q[0],q[1];
rzz(pi/2) q[0],q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
rzz(pi/2) q[0],q[2];
rzz(-pi/2) q[0],q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
u3(pi/2,2.159530790077624,-2.159530790077624) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
u3(pi/2,5.399769452990137,-5.399769452990137) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,5.301751762198135,-5.301751762198135) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[5];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,2.159530790077624,-2.159530790077624) q[6];
u3(pi/2,5.276619020969417,-5.276619020969417) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,1.3251237812841747,-1.3251237812841747) q[7];
u3(pi/2,4.454778382790327,-4.454778382790327) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,pi/4,-pi/4) q[1];
u3(pi/2,pi,-pi) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
u3(pi/2,13*pi/8,-13*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(-pi/2) q[1],q[3];
u3(pi/2,5.301751762198135,-5.301751762198135) q[3];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,5.399769452990137,-5.399769452990137) q[4];
u3(pi/2,2.061513099285622,-2.061513099285622) q[4];
rzz(pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[5];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[5];
u3(pi/2,3.583300580684518,-3.583300580684518) q[5];
rzz(pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,5.276619020969417,-5.276619020969417) q[6];
u3(pi/2,2.0860175219836226,-2.0860175219836226) q[6];
rzz(pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u3(pi/2,4.454778382790327,-4.454778382790327) q[7];
u3(pi/2,1.288681306502533,-1.288681306502533) q[7];
rzz(pi/2) q[1],q[7];
u3(pi/2,9*pi/8,-9*pi/8) q[2];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[3];
u3(pi/2,4.123026198571244,-4.123026198571244) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,2.061513099285622,-2.061513099285622) q[4];
u3(pi/2,4.810406671176691,-4.810406671176691) q[4];
rzz(-pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,3.583300580684518,-3.583300580684518) q[5];
u3(pi/2,0.24567254551072185,-0.24567254551072185) q[5];
rzz(pi/2) q[2],q[5];
rzz(pi/2) q[2],q[6];
u3(pi/2,2.0860175219836226,-2.0860175219836226) q[6];
u3(pi/2,5.129592484781415,-5.129592484781415) q[6];
rzz(-pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,1.288681306502533,-1.288681306502533) q[7];
u3(pi/2,4.381265114696325,-4.381265114696325) q[7];
rzz(pi/2) q[2],q[7];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[4];
u3(pi/2,4.025008507779242,-4.025008507779242) q[4];
rzz(pi/2) q[3],q[4];
rzz(-pi/2) q[3],q[5];
u3(pi/2,3.387265199100515,-3.387265199100515) q[5];
u3(pi/2,6.136158770991584,-6.136158770991584) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,1.9879998311916212,-1.9879998311916212) q[6];
u3(pi/2,4.933557103197411,-4.933557103197411) q[6];
rzz(pi/2) q[3],q[6];
rzz(-pi/2) q[3],q[7];
u3(pi/2,1.2396724611065324,-1.2396724611065324) q[7];
u3(pi/2,4.282619105373606,-4.282619105373606) q[7];
rzz(pi/2) q[3],q[7];
u3(pi/2,2.4542121809843462,-2.4542121809843462) q[4];
u3(pi/2,2.5528581903070657,-2.5528581903070657) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,6.136158770991584,-6.136158770991584) q[5];
u3(pi/2,2.2091679540043425,-2.2091679540043425) q[5];
rzz(pi/2) q[4],q[5];
rzz(-pi/2) q[4],q[6];
u3(pi/2,1.7919644496076181,-1.7919644496076181) q[6];
u3(pi/2,4.540858021498687,-4.540858021498687) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,4.282619105373606,-4.282619105373606) q[7];
u3(pi/2,0.9449910701998098,-0.9449910701998098) q[7];
rzz(pi/2) q[4],q[7];
u3(pi/2,3.779964280799239,-3.779964280799239) q[5];
u3(pi/2,6.087149925595583,-6.087149925595583) q[5];
rzz(-pi/2) q[5],q[6];
u3(pi/2,1.399265367908894,-1.399265367908894) q[6];
u3(pi/2,3.7554598581012386,-3.7554598581012386) q[6];
rzz(pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,0.9449910701998098,-0.9449910701998098) q[7];
u3(pi/2,3.6938846420908784,-3.6938846420908784) q[7];
rzz(pi/2) q[5],q[7];
u3(pi/2,2.184663531306342,-2.184663531306342) q[6];
u3(pi/2,4.9580615258954115,-4.9580615258954115) q[6];
rzz(pi/2) q[6],q[7];
u3(pi/2,3.6938846420908784,-3.6938846420908784) q[7];
u3(pi/2,6.050079132283224,-6.050079132283224) q[7];
rzz(-pi/2) q[6],q[7];
u3(pi,pi/2,-pi/2) q[1];
u3(pi,pi/4,-pi/4) q[2];
u3(pi,pi/8,-pi/8) q[3];
u3(pi,0,0) q[4];
u3(pi,4.025008507779242,-4.025008507779242) q[5];
u3(pi,3*pi/2,-3*pi/2) q[6];
u3(pi/2,4.479282805488327,-4.479282805488327) q[7];
u3(pi/2,5.448778298386137,-5.448778298386137) q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
measure q[8] -> c[8];
