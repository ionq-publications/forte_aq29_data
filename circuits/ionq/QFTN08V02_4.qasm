OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(pi/2,2.110521944681623,-2.110521944681623) q[7];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[7];
rzz(pi/2) q[6],q[7];
rzz(-pi/2) q[5],q[7];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
u3(pi/2,2.503221026380347,-2.503221026380347) q[7];
rzz(-pi/2) q[5],q[7];
rzz(pi/2) q[5],q[6];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
u3(pi/2,1.3741326266801754,-1.3741326266801754) q[6];
rzz(-pi/2) q[5],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,2.503221026380347,-2.503221026380347) q[7];
u3(pi/2,5.448778298386137,-5.448778298386137) q[7];
rzz(-pi/2) q[4],q[7];
rzz(-pi/2) q[4],q[6];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[6];
u3(pi/2,4.123026198571244,-4.123026198571244) q[6];
rzz(-pi/2) q[4],q[6];
rzz(-pi/2) q[4],q[5];
u3(pi/2,9*pi/8,-9*pi/8) q[5];
u3(pi/2,11*pi/8,-11*pi/8) q[5];
rzz(-pi/2) q[4],q[5];
rzz(-pi/2) q[3],q[7];
u3(pi/2,2.3071856447963444,-2.3071856447963444) q[7];
u3(pi/2,5.350760607594136,-5.350760607594136) q[7];
rzz(-pi/2) q[3],q[7];
rzz(-pi/2) q[3],q[6];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[6];
u3(pi/2,5*pi/4,-5*pi/4) q[6];
rzz(pi/2) q[3],q[6];
rzz(-pi/2) q[3],q[5];
u3(pi/2,11*pi/8,-11*pi/8) q[5];
u3(pi/2,pi/4,-pi/4) q[5];
rzz(-pi/2) q[3],q[5];
rzz(pi/2) q[3],q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[4];
rzz(pi/2) q[3],q[4];
rzz(-pi/2) q[2],q[7];
u3(pi/2,2.2091679540043425,-2.2091679540043425) q[7];
u3(pi/2,5.301751762198135,-5.301751762198135) q[7];
rzz(-pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[6];
u3(pi/2,5*pi/4,-5*pi/4) q[6];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[6];
rzz(pi/2) q[2],q[6];
rzz(-pi/2) q[2],q[5];
u3(pi/2,5*pi/4,-5*pi/4) q[5];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[5];
rzz(-pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[4];
u3(pi/2,pi/2,-pi/2) q[4];
u3(pi/2,11*pi/8,-11*pi/8) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(-pi/2) q[2],q[3];
rzz(-pi/2) q[1],q[7];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[7];
u3(pi/2,5.276619020969417,-5.276619020969417) q[7];
rzz(-pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[6];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[6];
u3(pi/2,3.779964280799239,-3.779964280799239) q[6];
rzz(-pi/2) q[1],q[6];
rzz(-pi/2) q[1],q[5];
u3(pi/2,3.730955435403238,-3.730955435403238) q[5];
u3(pi/2,0.4907167724907257,-0.4907167724907257) q[5];
rzz(-pi/2) q[1],q[5];
rzz(pi/2) q[1],q[4];
u3(pi/2,3*pi/8,-3*pi/8) q[4];
u3(pi/2,4.123026198571244,-4.123026198571244) q[4];
rzz(-pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(-pi/2) q[1],q[3];
rzz(-pi/2) q[1],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(-pi/2) q[1],q[2];
rzz(-pi/2) q[0],q[7];
rzz(-pi/2) q[0],q[6];
rzz(-pi/2) q[0],q[5];
rzz(pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[3];
rzz(-pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi/2,0,0) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,4.123026198571244,-4.123026198571244) q[4];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,3.6323094260805187,-3.6323094260805187) q[5];
u3(pi/2,pi/8,-pi/8) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,0.638371627209446,-0.638371627209446) q[6];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[7];
u3(pi/2,5.252114598271416,-5.252114598271416) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
u3(pi/2,pi/4,-pi/4) q[2];
rzz(-pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,9*pi/8,-9*pi/8) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,pi/4,-pi/4) q[4];
u3(pi/2,3.730955435403238,-3.730955435403238) q[4];
rzz(pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[5];
u3(pi/2,pi/8,-pi/8) q[5];
u3(pi/2,3.436274044496516,-3.436274044496516) q[5];
rzz(pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[6];
rzz(pi/2) q[1],q[6];
rzz(-pi/2) q[1],q[7];
u3(pi/2,2.110521944681623,-2.110521944681623) q[7];
u3(pi/2,5.227610175573416,-5.227610175573416) q[7];
rzz(pi/2) q[1],q[7];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,9*pi/8,-9*pi/8) q[3];
u3(pi/2,15*pi/8,-15*pi/8) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,3.730955435403238,-3.730955435403238) q[4];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[4];
rzz(-pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,3.436274044496516,-3.436274044496516) q[5];
u3(pi/2,0.09801769079200154,-0.09801769079200154) q[5];
rzz(pi/2) q[2],q[5];
rzz(pi/2) q[2],q[6];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[6];
u3(pi/2,3.583300580684518,-3.583300580684518) q[6];
rzz(-pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,5.227610175573416,-5.227610175573416) q[7];
u3(pi/2,2.0370086765876216,-2.0370086765876216) q[7];
rzz(pi/2) q[2],q[7];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.3376280351737964,-3.3376280351737964) q[4];
u3(pi/2,5.693822525366141,-5.693822525366141) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,0.09801769079200154,-0.09801769079200154) q[5];
u3(pi/2,2.8469112626830704,-2.8469112626830704) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,0.4417079270947249,-0.4417079270947249) q[6];
u3(pi/2,3.387265199100515,-3.387265199100515) q[6];
rzz(-pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,2.0370086765876216,-2.0370086765876216) q[7];
u3(pi/2,5.080583639385413,-5.080583639385413) q[7];
rzz(pi/2) q[3],q[7];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,2.8469112626830704,-2.8469112626830704) q[5];
u3(pi/2,5.203105752875415,-5.203105752875415) q[5];
rzz(pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,0.24567254551072185,-0.24567254551072185) q[6];
u3(pi/2,2.994566117401791,-2.994566117401791) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,5.080583639385413,-5.080583639385413) q[7];
u3(pi/2,1.7423272856808991,-1.7423272856808991) q[7];
rzz(pi/2) q[4],q[7];
u3(pi/2,3.6323094260805187,-3.6323094260805187) q[5];
u3(pi/2,3.7303271168725205,-3.7303271168725205) q[5];
rzz(-pi/2) q[5],q[6];
u3(pi/2,6.136158770991584,-6.136158770991584) q[6];
u3(pi/2,2.2091679540043425,-2.2091679540043425) q[6];
rzz(pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,1.7423272856808991,-1.7423272856808991) q[7];
u3(pi/2,4.491220857571968,-4.491220857571968) q[7];
rzz(pi/2) q[5],q[7];
u3(pi/2,3.779964280799239,-3.779964280799239) q[6];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
rzz(-pi/2) q[6],q[7];
u3(pi/2,1.349628203982175,-1.349628203982175) q[7];
u3(pi/2,3.70582269417452,-3.70582269417452) q[7];
rzz(pi/2) q[6],q[7];
u3(pi,pi,-pi) q[0];
u3(pi,0,0) q[1];
u3(pi,3*pi/2,-3*pi/2) q[2];
u3(pi,0,0) q[3];
u3(pi,pi/4,-pi/4) q[4];
u3(pi,0.5893627818134451,-0.5893627818134451) q[5];
u3(pi,0.5893627818134451,-0.5893627818134451) q[6];
u3(pi/2,5.276619020969417,-5.276619020969417) q[7];
u3(pi/2,2.110521944681623,-2.110521944681623) q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
