OPENQASM 2.0;
include "qelib1.inc";
qreg q[14];
creg c[14];
u3(pi/2,0.9116901880717581,-0.9116901880717581) q[13];
u3(pi/2,1.6970883514692063,-1.6970883514692063) q[13];
rzz(-pi/2) q[12],q[13];
rzz(-pi/2) q[11],q[13];
u3(pi/2,1.6970883514692063,-1.6970883514692063) q[13];
u3(pi/2,4.445981923360275,-4.445981923360275) q[13];
rzz(pi/2) q[11],q[13];
rzz(pi/2) q[11],q[12];
u3(pi/2,3.960920017646011,-3.960920017646011) q[12];
u3(pi/2,4.746318181043459,-4.746318181043459) q[12];
rzz(-pi/2) q[11],q[12];
rzz(-pi/2) q[10],q[13];
u3(pi/2,1.3043892697704822,-1.3043892697704822) q[13];
u3(pi/2,4.249946541776272,-4.249946541776272) q[13];
rzz(pi/2) q[10],q[13];
rzz(-pi/2) q[10],q[12];
u3(pi/2,1.6047255274536665,-1.6047255274536665) q[12];
u3(pi/2,4.353619099344735,-4.353619099344735) q[12];
rzz(-pi/2) q[10],q[12];
rzz(pi/2) q[10],q[11];
u3(pi/2,2.0740794698999814,-2.0740794698999814) q[11];
u3(pi/2,2.8594776332974297,-2.8594776332974297) q[11];
rzz(-pi/2) q[10],q[11];
rzz(pi/2) q[9],q[13];
u3(pi/2,1.108353888186479,-1.108353888186479) q[13];
u3(pi/2,4.151928850984271,-4.151928850984271) q[13];
rzz(-pi/2) q[9],q[13];
rzz(-pi/2) q[9],q[12];
u3(pi/2,1.2120264457549421,-1.2120264457549421) q[12];
u3(pi/2,4.156955399230014,-4.156955399230014) q[12];
rzz(-pi/2) q[9],q[12];
rzz(pi/2) q[9],q[11];
u3(pi/2,6.001070286887223,-6.001070286887223) q[11];
u3(pi/2,2.4667785515987055,-2.4667785515987055) q[11];
rzz(-pi/2) q[9],q[11];
rzz(pi/2) q[9],q[10];
u3(pi/2,4.84747746448905,-4.84747746448905) q[10];
u3(pi/2,5.632875627886499,-5.632875627886499) q[10];
rzz(-pi/2) q[9],q[10];
rzz(pi/2) q[8],q[13];
u3(pi/2,1.0103361973944776,-1.0103361973944776) q[13];
u3(pi/2,4.1029200055882695,-4.1029200055882695) q[13];
rzz(-pi/2) q[8],q[13];
rzz(-pi/2) q[8],q[12];
u3(pi/2,4.156955399230014,-4.156955399230014) q[12];
u3(pi/2,0.9173450548482195,-0.9173450548482195) q[12];
rzz(-pi/2) q[8],q[12];
rzz(pi/2) q[8],q[11];
u3(pi/2,5.608371205188498,-5.608371205188498) q[11];
u3(pi/2,2.2701148514839846,-2.2701148514839846) q[11];
rzz(-pi/2) q[8],q[11];
rzz(-pi/2) q[8],q[10];
u3(pi/2,2.491282974296706,-2.491282974296706) q[10];
u3(pi/2,5.240176546187775,-5.240176546187775) q[10];
rzz(-pi/2) q[8],q[10];
rzz(pi/2) q[8],q[9];
u3(pi/2,6.025574709585223,-6.025574709585223) q[9];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[9];
rzz(-pi/2) q[8],q[9];
rzz(-pi/2) q[7],q[13];
u3(pi/2,4.1029200055882695,-4.1029200055882695) q[13];
u3(pi/2,0.9361946107697583,-0.9361946107697583) q[13];
rzz(-pi/2) q[7],q[13];
rzz(-pi/2) q[7],q[12];
u3(pi/2,0.9173450548482195,-0.9173450548482195) q[12];
u3(pi/2,4.009928863042012,-4.009928863042012) q[12];
rzz(pi/2) q[7],q[12];
rzz(-pi/2) q[7],q[11];
u3(pi/2,2.2701148514839846,-2.2701148514839846) q[11];
u3(pi/2,5.313689814281776,-5.313689814281776) q[11];
rzz(-pi/2) q[7],q[11];
rzz(-pi/2) q[7],q[10];
u3(pi/2,5.240176546187775,-5.240176546187775) q[10];
u3(pi/2,1.901920192483261,-1.901920192483261) q[10];
rzz(pi/2) q[7],q[10];
rzz(-pi/2) q[7],q[9];
u3(pi/2,3.669380219392878,-3.669380219392878) q[9];
u3(pi/2,0.13508848410436108,-0.13508848410436108) q[9];
rzz(-pi/2) q[7],q[9];
rzz(pi/2) q[7],q[8];
u3(pi/2,5.252114598271416,-5.252114598271416) q[8];
u3(pi/2,6.037512761668864,-6.037512761668864) q[8];
rzz(-pi/2) q[7],q[8];
rzz(pi/2) q[6],q[13];
u3(pi/2,0.9361946107697583,-0.9361946107697583) q[13];
u3(pi/2,4.065849212275911,-4.065849212275911) q[13];
rzz(-pi/2) q[6],q[13];
rzz(-pi/2) q[6],q[12];
u3(pi/2,4.009928863042012,-4.009928863042012) q[12];
u3(pi/2,0.8438317867542184,-0.8438317867542184) q[12];
rzz(-pi/2) q[6],q[12];
rzz(pi/2) q[6],q[11];
u3(pi/2,5.313689814281776,-5.313689814281776) q[11];
u3(pi/2,2.123088315295982,-2.123088315295982) q[11];
rzz(-pi/2) q[6],q[11];
rzz(-pi/2) q[6],q[10];
u3(pi/2,1.901920192483261,-1.901920192483261) q[10];
u3(pi/2,4.945495155281052,-4.945495155281052) q[10];
rzz(-pi/2) q[6],q[10];
rzz(pi/2) q[6],q[9];
u3(pi/2,0.13508848410436108,-0.13508848410436108) q[9];
u3(pi/2,3.0800174375794334,-3.0800174375794334) q[9];
rzz(-pi/2) q[6],q[9];
rzz(-pi/2) q[6],q[8];
u3(pi/2,6.037512761668864,-6.037512761668864) q[8];
u3(pi/2,2.503221026380347,-2.503221026380347) q[8];
rzz(-pi/2) q[6],q[8];
rzz(-pi/2) q[6],q[7];
u3(pi/2,5.841477380084861,-5.841477380084861) q[7];
u3(pi/2,0.3436902363027234,-0.3436902363027234) q[7];
rzz(-pi/2) q[6],q[7];
rzz(-pi/2) q[5],q[13];
rzz(-pi/2) q[5],q[13];
rzz(-pi/2) q[5],q[12];
u3(pi/2,0.8438317867542184,-0.8438317867542184) q[12];
u3(pi/2,3.9728580697296523,-3.9728580697296523) q[12];
rzz(pi/2) q[5],q[12];
rzz(-pi/2) q[5],q[11];
u3(pi/2,2.123088315295982,-2.123088315295982) q[11];
u3(pi/2,5.240176546187775,-5.240176546187775) q[11];
rzz(-pi/2) q[5],q[11];
rzz(-pi/2) q[5],q[10];
u3(pi/2,4.945495155281052,-4.945495155281052) q[10];
u3(pi/2,1.7548936562952584,-1.7548936562952584) q[10];
rzz(pi/2) q[5],q[10];
rzz(-pi/2) q[5],q[9];
u3(pi/2,3.0800174375794334,-3.0800174375794334) q[9];
u3(pi/2,6.123592400377225,-6.123592400377225) q[9];
rzz(-pi/2) q[5],q[9];
rzz(-pi/2) q[5],q[8];
u3(pi/2,2.503221026380347,-2.503221026380347) q[8];
u3(pi/2,5.448778298386137,-5.448778298386137) q[8];
rzz(-pi/2) q[5],q[8];
rzz(pi/2) q[5],q[7];
u3(pi/2,0.3436902363027234,-0.3436902363027234) q[7];
u3(pi/2,3.0925838081937926,-3.0925838081937926) q[7];
rzz(-pi/2) q[5],q[7];
rzz(pi/2) q[5],q[6];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
u3(pi/2,4.516353598800687,-4.516353598800687) q[6];
rzz(-pi/2) q[5],q[6];
rzz(pi/2) q[4],q[13];
rzz(-pi/2) q[4],q[13];
rzz(-pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[11];
u3(pi/2,2.098583892597982,-2.098583892597982) q[11];
u3(pi/2,5.227610175573416,-5.227610175573416) q[11];
rzz(-pi/2) q[4],q[11];
rzz(-pi/2) q[4],q[10];
u3(pi/2,1.7548936562952584,-1.7548936562952584) q[10];
u3(pi/2,4.871981887187051,-4.871981887187051) q[10];
rzz(-pi/2) q[4],q[10];
rzz(pi/2) q[4],q[9];
u3(pi/2,6.123592400377225,-6.123592400377225) q[9];
u3(pi/2,2.9329909013914306,-2.9329909013914306) q[9];
rzz(-pi/2) q[4],q[9];
rzz(-pi/2) q[4],q[8];
u3(pi/2,2.3071856447963444,-2.3071856447963444) q[8];
u3(pi/2,5.350760607594136,-5.350760607594136) q[8];
rzz(-pi/2) q[4],q[8];
rzz(-pi/2) q[4],q[7];
u3(pi/2,6.234176461783585,-6.234176461783585) q[7];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[7];
rzz(-pi/2) q[4],q[7];
rzz(-pi/2) q[4],q[6];
u3(pi/2,4.516353598800687,-4.516353598800687) q[6];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[6];
rzz(-pi/2) q[4],q[6];
rzz(pi/2) q[4],q[5];
u3(pi/2,3*pi/4,-3*pi/4) q[5];
u3(pi/2,pi,-pi) q[5];
rzz(pi/2) q[4],q[5];
rzz(-pi/2) q[3],q[13];
rzz(-pi/2) q[3],q[13];
rzz(-pi/2) q[3],q[12];
rzz(pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[11];
rzz(-pi/2) q[3],q[11];
rzz(-pi/2) q[3],q[10];
u3(pi/2,1.730389233597258,-1.730389233597258) q[10];
u3(pi/2,4.8594155165726916,-4.8594155165726916) q[10];
rzz(-pi/2) q[3],q[10];
rzz(pi/2) q[3],q[9];
u3(pi/2,2.9329909013914306,-2.9329909013914306) q[9];
u3(pi/2,6.050079132283224,-6.050079132283224) q[9];
rzz(-pi/2) q[3],q[9];
rzz(-pi/2) q[3],q[8];
u3(pi/2,2.2091679540043425,-2.2091679540043425) q[8];
u3(pi/2,5.301751762198135,-5.301751762198135) q[8];
rzz(-pi/2) q[3],q[8];
rzz(pi/2) q[3],q[7];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[7];
u3(pi/2,5.9394950708768635,-5.9394950708768635) q[7];
rzz(-pi/2) q[3],q[7];
rzz(-pi/2) q[3],q[6];
u3(pi/2,4.123026198571244,-4.123026198571244) q[6];
u3(pi/2,pi/4,-pi/4) q[6];
rzz(-pi/2) q[3],q[6];
rzz(-pi/2) q[3],q[5];
u3(pi/2,0,0) q[5];
u3(pi/2,7*pi/8,-7*pi/8) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[4];
u3(pi/2,pi/4,-pi/4) q[4];
u3(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[3],q[4];
rzz(pi/2) q[2],q[13];
rzz(-pi/2) q[2],q[13];
rzz(-pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[11];
rzz(pi/2) q[2],q[11];
rzz(-pi/2) q[2],q[10];
rzz(-pi/2) q[2],q[10];
rzz(-pi/2) q[2],q[9];
u3(pi/2,6.050079132283224,-6.050079132283224) q[9];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[9];
rzz(pi/2) q[2],q[9];
rzz(-pi/2) q[2],q[8];
u3(pi/2,5.301751762198135,-5.301751762198135) q[8];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[8];
rzz(-pi/2) q[2],q[8];
rzz(-pi/2) q[2],q[7];
u3(pi/2,5.9394950708768635,-5.9394950708768635) q[7];
u3(pi/2,7*pi/8,-7*pi/8) q[7];
rzz(-pi/2) q[2],q[7];
rzz(pi/2) q[2],q[6];
u3(pi/2,5*pi/4,-5*pi/4) q[6];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[6];
rzz(-pi/2) q[2],q[6];
rzz(-pi/2) q[2],q[5];
u3(pi/2,15*pi/8,-15*pi/8) q[5];
u3(pi/2,2.552229871776348,-2.552229871776348) q[5];
rzz(-pi/2) q[2],q[5];
rzz(pi/2) q[2],q[4];
u3(pi/2,pi/2,-pi/2) q[4];
u3(pi/2,11*pi/8,-11*pi/8) q[4];
rzz(-pi/2) q[2],q[4];
rzz(pi/2) q[2],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(-pi/2) q[2],q[3];
rzz(-pi/2) q[1],q[13];
rzz(-pi/2) q[1],q[13];
rzz(-pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[11];
rzz(pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[10];
rzz(-pi/2) q[1],q[10];
rzz(-pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[8];
u3(pi/2,5.276619020969417,-5.276619020969417) q[8];
u3(pi/2,2.123088315295982,-2.123088315295982) q[8];
rzz(-pi/2) q[1],q[8];
rzz(-pi/2) q[1],q[7];
u3(pi/2,15*pi/8,-15*pi/8) q[7];
u3(pi/2,2.7243891491930685,-2.7243891491930685) q[7];
rzz(pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[6];
u3(pi/2,3.82897312619524,-3.82897312619524) q[6];
u3(pi/2,0.638371627209446,-0.638371627209446) q[6];
rzz(-pi/2) q[1],q[6];
rzz(-pi/2) q[1],q[5];
u3(pi/2,5.693822525366141,-5.693822525366141) q[5];
u3(pi/2,2.4542121809843462,-2.4542121809843462) q[5];
rzz(pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[4];
u3(pi/2,3*pi/8,-3*pi/8) q[4];
u3(pi/2,4.123026198571244,-4.123026198571244) q[4];
rzz(-pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
u3(pi/2,5*pi/8,-5*pi/8) q[3];
rzz(-pi/2) q[1],q[3];
rzz(-pi/2) q[1],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(-pi/2) q[1],q[2];
rzz(-pi/2) q[0],q[13];
rzz(-pi/2) q[0],q[12];
rzz(pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[7];
rzz(pi/2) q[0],q[6];
rzz(-pi/2) q[0],q[5];
rzz(-pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[3];
rzz(pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi/2,0,0) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,5*pi/8,-5*pi/8) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,4.123026198571244,-4.123026198571244) q[4];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,5.595804834574139,-5.595804834574139) q[5];
u3(pi/2,3*pi/4,-3*pi/4) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,3.779964280799239,-3.779964280799239) q[6];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,5.865981802782861,-5.865981802782861) q[7];
u3(pi/2,2.6998847264950685,-2.6998847264950685) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,2.123088315295982,-2.123088315295982) q[8];
u3(pi/2,5.252114598271416,-5.252114598271416) q[8];
rzz(pi/2) q[0],q[8];
rzz(pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[12];
rzz(pi/2) q[0],q[13];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
u3(pi/2,pi/4,-pi/4) q[2];
rzz(pi/2) q[1],q[2];
rzz(-pi/2) q[1],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,pi/4,-pi/4) q[4];
u3(pi/2,3.730955435403238,-3.730955435403238) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,3*pi/4,-3*pi/4) q[5];
u3(pi/2,5.399769452990137,-5.399769452990137) q[5];
rzz(pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[6];
rzz(pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u3(pi/2,5.841477380084861,-5.841477380084861) q[7];
u3(pi/2,2.675380303797068,-2.675380303797068) q[7];
rzz(-pi/2) q[1],q[7];
rzz(pi/2) q[1],q[8];
u3(pi/2,5.252114598271416,-5.252114598271416) q[8];
u3(pi/2,2.098583892597982,-2.098583892597982) q[8];
rzz(pi/2) q[1],q[8];
rzz(pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[9];
rzz(pi/2) q[1],q[10];
rzz(pi/2) q[1],q[10];
rzz(pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[11];
rzz(pi/2) q[1],q[12];
rzz(pi/2) q[1],q[12];
rzz(pi/2) q[1],q[13];
rzz(-pi/2) q[1],q[13];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
u3(pi/2,0,0) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
u3(pi/2,9*pi/8,-9*pi/8) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,3.730955435403238,-3.730955435403238) q[4];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[4];
rzz(pi/2) q[2],q[4];
rzz(-pi/2) q[2],q[5];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[5];
u3(pi/2,5.203105752875415,-5.203105752875415) q[5];
rzz(pi/2) q[2],q[5];
rzz(pi/2) q[2],q[6];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[6];
u3(pi/2,0.4417079270947249,-0.4417079270947249) q[6];
rzz(pi/2) q[2],q[6];
rzz(-pi/2) q[2],q[7];
u3(pi/2,2.675380303797068,-2.675380303797068) q[7];
u3(pi/2,5.76796411199086,-5.76796411199086) q[7];
rzz(pi/2) q[2],q[7];
rzz(pi/2) q[2],q[8];
u3(pi/2,2.098583892597982,-2.098583892597982) q[8];
u3(pi/2,5.2156721234897745,-5.2156721234897745) q[8];
rzz(pi/2) q[2],q[8];
rzz(pi/2) q[2],q[9];
u3(pi/2,2.8776988706882505,-2.8776988706882505) q[9];
u3(pi/2,6.007353472194402,-6.007353472194402) q[9];
rzz(-pi/2) q[2],q[9];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[11];
rzz(-pi/2) q[2],q[11];
rzz(pi/2) q[2],q[12];
rzz(pi/2) q[2],q[12];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[13];
u3(pi/2,5*pi/8,-5*pi/8) q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(-pi/2) q[3],q[4];
u3(pi/2,3.3376280351737964,-3.3376280351737964) q[4];
u3(pi/2,5.693822525366141,-5.693822525366141) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,5.203105752875415,-5.203105752875415) q[5];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[5];
rzz(-pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,0.4417079270947249,-0.4417079270947249) q[6];
u3(pi/2,3.387265199100515,-3.387265199100515) q[6];
rzz(pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,5.76796411199086,-5.76796411199086) q[7];
u3(pi/2,2.5277254490783476,-2.5277254490783476) q[7];
rzz(pi/2) q[3],q[7];
rzz(-pi/2) q[3],q[8];
u3(pi/2,2.0740794698999814,-2.0740794698999814) q[8];
u3(pi/2,5.166663278093774,-5.166663278093774) q[8];
rzz(pi/2) q[3],q[8];
rzz(pi/2) q[3],q[9];
u3(pi/2,2.8657608186046093,-2.8657608186046093) q[9];
u3(pi/2,5.982220730965683,-5.982220730965683) q[9];
rzz(pi/2) q[3],q[9];
rzz(-pi/2) q[3],q[10];
u3(pi/2,1.7027432182456679,-1.7027432182456679) q[10];
u3(pi/2,4.8317695012211015,-4.8317695012211015) q[10];
rzz(pi/2) q[3],q[10];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[12];
rzz(pi/2) q[3],q[12];
rzz(pi/2) q[3],q[13];
rzz(pi/2) q[3],q[13];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,4.810406671176691,-4.810406671176691) q[5];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[5];
rzz(-pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,3.387265199100515,-3.387265199100515) q[6];
u3(pi/2,6.136158770991584,-6.136158770991584) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,2.5277254490783476,-2.5277254490783476) q[7];
u3(pi/2,5.4732827210841375,-5.4732827210841375) q[7];
rzz(-pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,5.166663278093774,-5.166663278093774) q[8];
u3(pi/2,1.926424615181261,-1.926424615181261) q[8];
rzz(pi/2) q[4],q[8];
rzz(pi/2) q[4],q[9];
u3(pi/2,5.982220730965683,-5.982220730965683) q[9];
u3(pi/2,2.79161923197989,-2.79161923197989) q[9];
rzz(pi/2) q[4],q[9];
rzz(-pi/2) q[4],q[10];
u3(pi/2,1.6901768476313088,-1.6901768476313088) q[10];
u3(pi/2,4.807265078523101,-4.807265078523101) q[10];
rzz(pi/2) q[4],q[10];
rzz(pi/2) q[4],q[11];
u3(pi/2,5.216928760551211,-5.216928760551211) q[11];
u3(pi/2,2.0633980548777764,-2.0633980548777764) q[11];
rzz(-pi/2) q[4],q[11];
rzz(pi/2) q[4],q[12];
rzz(pi/2) q[4],q[12];
rzz(pi/2) q[4],q[13];
rzz(pi/2) q[4],q[13];
u3(pi/2,2.4542121809843462,-2.4542121809843462) q[5];
u3(pi/2,2.552229871776348,-2.552229871776348) q[5];
rzz(-pi/2) q[5],q[6];
u3(pi/2,2.994566117401791,-2.994566117401791) q[6];
u3(pi/2,5.350760607594136,-5.350760607594136) q[6];
rzz(pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,2.3316900674943444,-2.3316900674943444) q[7];
u3(pi/2,5.080583639385413,-5.080583639385413) q[7];
rzz(pi/2) q[5],q[7];
rzz(-pi/2) q[5],q[8];
u3(pi/2,5.0680172687710545,-5.0680172687710545) q[8];
u3(pi/2,1.730389233597258,-1.730389233597258) q[8];
rzz(pi/2) q[5],q[8];
rzz(pi/2) q[5],q[9];
u3(pi/2,2.79161923197989,-2.79161923197989) q[9];
u3(pi/2,5.8351941947776815,-5.8351941947776815) q[9];
rzz(pi/2) q[5],q[9];
rzz(pi/2) q[5],q[10];
u3(pi/2,4.807265078523101,-4.807265078523101) q[10];
u3(pi/2,1.6166635795373074,-1.6166635795373074) q[10];
rzz(-pi/2) q[5],q[10];
rzz(pi/2) q[5],q[11];
u3(pi/2,5.2049907084675695,-5.2049907084675695) q[11];
u3(pi/2,2.038893632179776,-2.038893632179776) q[11];
rzz(pi/2) q[5],q[11];
rzz(pi/2) q[5],q[12];
u3(pi/2,0.8180707269947822,-0.8180707269947822) q[12];
u3(pi/2,3.9477253285009337,-3.9477253285009337) q[12];
rzz(-pi/2) q[5],q[12];
rzz(pi/2) q[5],q[13];
rzz(pi/2) q[5],q[13];
u3(pi/2,3.779964280799239,-3.779964280799239) q[6];
u3(pi/2,3.82897312619524,-3.82897312619524) q[6];
rzz(pi/2) q[6],q[7];
u3(pi/2,5.080583639385413,-5.080583639385413) q[7];
u3(pi/2,1.1535928223981722,-1.1535928223981722) q[7];
rzz(pi/2) q[6],q[7];
rzz(-pi/2) q[6],q[8];
u3(pi/2,4.871981887187051,-4.871981887187051) q[8];
u3(pi/2,1.337690151898534,-1.337690151898534) q[8];
rzz(pi/2) q[6],q[8];
rzz(pi/2) q[6],q[9];
u3(pi/2,5.8351941947776815,-5.8351941947776815) q[9];
u3(pi/2,2.4975661596038856,-2.4975661596038856) q[9];
rzz(pi/2) q[6],q[9];
rzz(-pi/2) q[6],q[10];
u3(pi/2,1.6166635795373074,-1.6166635795373074) q[10];
u3(pi/2,4.660238542335099,-4.660238542335099) q[10];
rzz(pi/2) q[6],q[10];
rzz(pi/2) q[6],q[11];
u3(pi/2,2.038893632179776,-2.038893632179776) q[11];
u3(pi/2,5.131477440373568,-5.131477440373568) q[11];
rzz(pi/2) q[6],q[11];
rzz(-pi/2) q[6],q[12];
u3(pi/2,3.9477253285009337,-3.9477253285009337) q[12];
u3(pi/2,0.7816282522131405,-0.7816282522131405) q[12];
rzz(pi/2) q[6],q[12];
rzz(pi/2) q[6],q[13];
u3(pi/2,0.9047786842338603,-0.9047786842338603) q[13];
u3(pi/2,4.033804967209295,-4.033804967209295) q[13];
rzz(pi/2) q[6],q[13];
u3(pi/2,2.7243891491930685,-2.7243891491930685) q[7];
u3(pi/2,5.841477380084861,-5.841477380084861) q[7];
rzz(-pi/2) q[7],q[8];
u3(pi/2,4.479282805488327,-4.479282805488327) q[8];
u3(pi/2,0.5522919885010856,-0.5522919885010856) q[8];
rzz(pi/2) q[7],q[8];
rzz(pi/2) q[7],q[9];
u3(pi/2,2.4975661596038856,-2.4975661596038856) q[9];
u3(pi/2,5.246459731494954,-5.246459731494954) q[9];
rzz(pi/2) q[7],q[9];
rzz(pi/2) q[7],q[10];
u3(pi/2,4.660238542335099,-4.660238542335099) q[10];
u3(pi/2,1.321982188630585,-1.321982188630585) q[10];
rzz(pi/2) q[7],q[10];
rzz(pi/2) q[7],q[11];
u3(pi/2,5.131477440373568,-5.131477440373568) q[11];
u3(pi/2,1.8912387774610553,-1.8912387774610553) q[11];
rzz(pi/2) q[7],q[11];
rzz(pi/2) q[7],q[12];
u3(pi/2,0.7816282522131405,-0.7816282522131405) q[12];
u3(pi/2,3.8742120604069332,-3.8742120604069332) q[12];
rzz(-pi/2) q[7],q[12];
rzz(pi/2) q[7],q[13];
u3(pi/2,4.033804967209295,-4.033804967209295) q[13];
u3(pi/2,0.8677078909215009,-0.8677078909215009) q[13];
rzz(pi/2) q[7],q[13];
u3(pi/2,2.123088315295982,-2.123088315295982) q[8];
u3(pi/2,5.252114598271416,-5.252114598271416) q[8];
rzz(pi/2) q[8],q[9];
u3(pi/2,5.246459731494954,-5.246459731494954) q[9];
u3(pi/2,1.319468914507713,-1.319468914507713) q[9];
rzz(-pi/2) q[8],q[9];
rzz(pi/2) q[8],q[10];
u3(pi/2,1.321982188630585,-1.321982188630585) q[10];
u3(pi/2,4.070875760521654,-4.070875760521654) q[10];
rzz(pi/2) q[8],q[10];
rzz(pi/2) q[8],q[11];
u3(pi/2,1.8912387774610553,-1.8912387774610553) q[11];
u3(pi/2,4.836796049466845,-4.836796049466845) q[11];
rzz(pi/2) q[8],q[11];
rzz(-pi/2) q[8],q[12];
u3(pi/2,3.8742120604069332,-3.8742120604069332) q[12];
u3(pi/2,0.6346017160251383,-0.6346017160251383) q[12];
rzz(pi/2) q[8],q[12];
rzz(pi/2) q[8],q[13];
u3(pi/2,0.8677078909215009,-0.8677078909215009) q[13];
u3(pi/2,3.960291699115293,-3.960291699115293) q[13];
rzz(pi/2) q[8],q[13];
rzz(-pi/2) q[9],q[10];
u3(pi/2,0.9292831069318608,-0.9292831069318608) q[10];
u3(pi/2,3.285477597124206,-3.285477597124206) q[10];
rzz(pi/2) q[9],q[10];
rzz(pi/2) q[9],q[11];
u3(pi/2,4.836796049466845,-4.836796049466845) q[11];
u3(pi/2,1.3025043141783283,-1.3025043141783283) q[11];
rzz(pi/2) q[9],q[11];
rzz(-pi/2) q[9],q[12];
u3(pi/2,3.7761943696149314,-3.7761943696149314) q[12];
u3(pi/2,0.43793801591041714,-0.43793801591041714) q[12];
rzz(pi/2) q[9],q[12];
rzz(pi/2) q[9],q[13];
u3(pi/2,3.960291699115293,-3.960291699115293) q[13];
u3(pi/2,0.7206813547334985,-0.7206813547334985) q[13];
rzz(pi/2) q[9],q[13];
u3(pi,1.7146812703293088,-1.7146812703293088) q[10];
rzz(-pi/2) q[10],q[11];
u3(pi/2,4.444096967768122,-4.444096967768122) q[11];
u3(pi/2,0.51710615078088,-0.51710615078088) q[11];
rzz(pi/2) q[10],q[11];
rzz(pi/2) q[10],q[12];
u3(pi/2,0.43793801591041714,-0.43793801591041714) q[12];
u3(pi/2,3.186831587801486,-3.186831587801486) q[12];
rzz(pi/2) q[10],q[12];
rzz(pi/2) q[10],q[13];
u3(pi/2,0.7206813547334985,-0.7206813547334985) q[13];
u3(pi/2,3.6656103082085707,-3.6656103082085707) q[13];
rzz(-pi/2) q[10],q[13];
u3(pi,5.22949513116557,-5.22949513116557) q[11];
rzz(pi/2) q[11],q[12];
u3(pi/2,3.186831587801486,-3.186831587801486) q[12];
u3(pi/2,5.543026077993831,-5.543026077993831) q[12];
rzz(pi/2) q[11],q[12];
rzz(pi/2) q[11],q[13];
u3(pi/2,0.5240176546187775,-0.5240176546187775) q[13];
u3(pi/2,3.2729112265098466,-3.2729112265098466) q[13];
rzz(-pi/2) q[11],q[13];
rzz(pi/2) q[12],q[13];
u3(pi/2,0.13131857292005333,-0.13131857292005333) q[13];
u3(pi/2,2.487513063112398,-2.487513063112398) q[13];
rzz(pi/2) q[12],q[13];
u3(pi,0,0) q[0];
u3(pi,pi/2,-pi/2) q[1];
u3(pi,3*pi/2,-3*pi/2) q[2];
u3(pi,pi/2,-pi/2) q[3];
u3(pi,9*pi/8,-9*pi/8) q[4];
u3(pi,pi/4,-pi/4) q[5];
u3(pi,2.94555727200579,-2.94555727200579) q[6];
u3(pi,4.098521775873244,-4.098521775873244) q[7];
u3(pi,2.675380303797068,-2.675380303797068) q[8];
u3(pi,4.73689340308269,-4.73689340308269) q[9];
u3(pi,0.28211502029236346,-0.28211502029236346) q[10];
u3(pi,0.5893627818134451,-0.5893627818134451) q[11];
u3(pi,4.006787270388423,-4.006787270388423) q[12];
u3(pi,0.9167167363175016,-0.9167167363175016) q[13];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
measure q[8] -> c[8];
measure q[9] -> c[9];
measure q[10] -> c[10];
measure q[11] -> c[11];
measure q[12] -> c[12];
measure q[13] -> c[13];
