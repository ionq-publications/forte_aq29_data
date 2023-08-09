OPENQASM 2.0;
include "qelib1.inc";
qreg q[9];
creg c[9];
u3(pi/2,2.503221026380347,-2.503221026380347) q[8];
u3(pi/2,3.2886191897777954,-3.2886191897777954) q[8];
rzz(pi/2) q[7],q[8];
rzz(-pi/2) q[6],q[8];
u3(pi/2,3.2886191897777954,-3.2886191897777954) q[8];
u3(pi/2,6.037512761668864,-6.037512761668864) q[8];
rzz(-pi/2) q[6],q[8];
rzz(pi/2) q[6],q[7];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[7];
u3(pi/2,4.466716434873968,-4.466716434873968) q[7];
rzz(-pi/2) q[6],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,6.037512761668864,-6.037512761668864) q[8];
u3(pi/2,2.6998847264950685,-2.6998847264950685) q[8];
rzz(-pi/2) q[5],q[8];
rzz(-pi/2) q[5],q[7];
u3(pi/2,4.466716434873968,-4.466716434873968) q[7];
u3(pi/2,0.9324246995854506,-0.9324246995854506) q[7];
rzz(-pi/2) q[5],q[7];
rzz(-pi/2) q[5],q[6];
u3(pi/2,2.061513099285622,-2.061513099285622) q[6];
u3(pi/2,2.8469112626830704,-2.8469112626830704) q[6];
rzz(-pi/2) q[5],q[6];
rzz(-pi/2) q[4],q[8];
u3(pi/2,5.841477380084861,-5.841477380084861) q[8];
u3(pi/2,2.6018670357030667,-2.6018670357030667) q[8];
rzz(-pi/2) q[4],q[8];
rzz(-pi/2) q[4],q[7];
u3(pi/2,4.074017353175243,-4.074017353175243) q[7];
u3(pi/2,0.7363893180014475,-0.7363893180014475) q[7];
rzz(pi/2) q[4],q[7];
rzz(-pi/2) q[4],q[6];
u3(pi/2,2.8469112626830704,-2.8469112626830704) q[6];
u3(pi/2,5.595804834574139,-5.595804834574139) q[6];
rzz(-pi/2) q[4],q[6];
rzz(pi/2) q[4],q[5];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[5];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[5];
rzz(pi/2) q[4],q[5];
rzz(-pi/2) q[3],q[8];
u3(pi/2,2.6018670357030667,-2.6018670357030667) q[8];
u3(pi/2,5.694450843896859,-5.694450843896859) q[8];
rzz(-pi/2) q[3],q[8];
rzz(-pi/2) q[3],q[7];
u3(pi/2,3.8779819715912405,-3.8779819715912405) q[7];
u3(pi/2,0.638371627209446,-0.638371627209446) q[7];
rzz(-pi/2) q[3],q[7];
rzz(pi/2) q[3],q[6];
u3(pi/2,2.4542121809843462,-2.4542121809843462) q[6];
u3(pi/2,5.399769452990137,-5.399769452990137) q[6];
rzz(-pi/2) q[3],q[6];
rzz(-pi/2) q[3],q[5];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[5];
u3(pi/2,4.123026198571244,-4.123026198571244) q[5];
rzz(-pi/2) q[3],q[5];
rzz(-pi/2) q[3],q[4];
u3(pi/2,pi/4,-pi/4) q[4];
u3(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[3],q[4];
rzz(-pi/2) q[2],q[8];
u3(pi/2,2.5528581903070657,-2.5528581903070657) q[8];
u3(pi/2,5.66931810266814,-5.66931810266814) q[8];
rzz(-pi/2) q[2],q[8];
rzz(pi/2) q[2],q[7];
u3(pi/2,0.638371627209446,-0.638371627209446) q[7];
u3(pi/2,3.730955435403238,-3.730955435403238) q[7];
rzz(-pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[6];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[6];
u3(pi/2,5.301751762198135,-5.301751762198135) q[6];
rzz(-pi/2) q[2],q[6];
rzz(pi/2) q[2],q[5];
u3(pi/2,4.123026198571244,-4.123026198571244) q[5];
u3(pi/2,pi/4,-pi/4) q[5];
rzz(-pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[4];
u3(pi/2,pi/2,-pi/2) q[4];
u3(pi/2,11*pi/8,-11*pi/8) q[4];
rzz(-pi/2) q[2],q[4];
rzz(pi/2) q[2],q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(-pi/2) q[2],q[3];
rzz(-pi/2) q[1],q[8];
u3(pi/2,2.5277254490783476,-2.5277254490783476) q[8];
u3(pi/2,5.657380050584499,-5.657380050584499) q[8];
rzz(-pi/2) q[1],q[8];
rzz(-pi/2) q[1],q[7];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[7];
u3(pi/2,3.70582269417452,-3.70582269417452) q[7];
rzz(pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[6];
u3(pi/2,2.159530790077624,-2.159530790077624) q[6];
u3(pi/2,5.252114598271416,-5.252114598271416) q[6];
rzz(-pi/2) q[1],q[6];
rzz(-pi/2) q[1],q[5];
u3(pi/2,5*pi/4,-5*pi/4) q[5];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[4];
u3(pi/2,3*pi/8,-3*pi/8) q[4];
u3(pi/2,4.123026198571244,-4.123026198571244) q[4];
rzz(-pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(-pi/2) q[1],q[3];
rzz(-pi/2) q[1],q[2];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(-pi/2) q[1],q[2];
rzz(-pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[7];
rzz(pi/2) q[0],q[6];
rzz(-pi/2) q[0],q[5];
rzz(-pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[3];
rzz(pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi,-pi) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
u3(pi/2,3.730955435403238,-3.730955435403238) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,5.252114598271416,-5.252114598271416) q[6];
u3(pi/2,2.061513099285622,-2.061513099285622) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,3.70582269417452,-3.70582269417452) q[7];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,2.515787396994706,-2.515787396994706) q[8];
u3(pi/2,5.64481367997014,-5.64481367997014) q[8];
rzz(-pi/2) q[0],q[8];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,0,0) q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,9*pi/8,-9*pi/8) q[3];
rzz(-pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,pi/4,-pi/4) q[4];
u3(pi/2,3.730955435403238,-3.730955435403238) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[5];
u3(pi/2,3.6323094260805187,-3.6323094260805187) q[5];
rzz(pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[6];
u3(pi/2,2.061513099285622,-2.061513099285622) q[6];
u3(pi/2,5.154096907479415,-5.154096907479415) q[6];
rzz(pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[7];
u3(pi/2,0.5152211951887261,-0.5152211951887261) q[7];
rzz(pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[8];
u3(pi/2,2.503221026380347,-2.503221026380347) q[8];
u3(pi/2,5.632875627886499,-5.632875627886499) q[8];
rzz(pi/2) q[1],q[8];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,pi/8,-pi/8) q[3];
u3(pi/2,7*pi/8,-7*pi/8) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,3.730955435403238,-3.730955435403238) q[4];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,3.6323094260805187,-3.6323094260805187) q[5];
u3(pi/2,0.2946813909067226,-0.2946813909067226) q[5];
rzz(pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[6];
u3(pi/2,2.012504253889621,-2.012504253889621) q[6];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[6];
rzz(pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,0.5152211951887261,-0.5152211951887261) q[7];
u3(pi/2,3.6078050033825186,-3.6078050033825186) q[7];
rzz(pi/2) q[2],q[7];
rzz(pi/2) q[2],q[8];
u3(pi/2,5.632875627886499,-5.632875627886499) q[8];
u3(pi/2,2.4667785515987055,-2.4667785515987055) q[8];
rzz(pi/2) q[2],q[8];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[4];
u3(pi/2,2.552229871776348,-2.552229871776348) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,0.2946813909067226,-0.2946813909067226) q[5];
u3(pi/2,3.0435749627977917,-3.0435749627977917) q[5];
rzz(-pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[6];
u3(pi/2,1.7178228629828987,-1.7178228629828987) q[6];
rzz(pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,3.6078050033825186,-3.6078050033825186) q[7];
u3(pi/2,0.36819465900072373,-0.36819465900072373) q[7];
rzz(-pi/2) q[3],q[7];
rzz(pi/2) q[3],q[8];
u3(pi/2,2.4667785515987055,-2.4667785515987055) q[8];
u3(pi/2,5.559362359792498,-5.559362359792498) q[8];
rzz(pi/2) q[3],q[8];
u3(pi/2,4.123026198571244,-4.123026198571244) q[4];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,6.185167616387585,-6.185167616387585) q[5];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[5];
rzz(pi/2) q[4],q[5];
rzz(-pi/2) q[4],q[6];
u3(pi/2,4.8594155165726916,-4.8594155165726916) q[6];
u3(pi/2,1.3251237812841747,-1.3251237812841747) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,3.509787312590517,-3.509787312590517) q[7];
u3(pi/2,0.1715309588860027,-0.1715309588860027) q[7];
rzz(pi/2) q[4],q[7];
rzz(-pi/2) q[4],q[8];
u3(pi/2,2.4177697062027046,-2.4177697062027046) q[8];
u3(pi/2,5.460716350469778,-5.460716350469778) q[8];
rzz(pi/2) q[4],q[8];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
u3(pi/2,pi/4,-pi/4) q[5];
rzz(pi/2) q[5],q[6];
u3(pi/2,1.3251237812841747,-1.3251237812841747) q[6];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[6];
rzz(pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,0.1715309588860027,-0.1715309588860027) q[7];
u3(pi/2,2.9204245307770718,-2.9204245307770718) q[7];
rzz(-pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,5.460716350469778,-5.460716350469778) q[8];
u3(pi/2,2.123088315295982,-2.123088315295982) q[8];
rzz(pi/2) q[5],q[8];
u3(pi/2,5.252114598271416,-5.252114598271416) q[6];
u3(pi/2,2.061513099285622,-2.061513099285622) q[6];
rzz(-pi/2) q[6],q[7];
u3(pi/2,2.9204245307770718,-2.9204245307770718) q[7];
u3(pi/2,5.276619020969417,-5.276619020969417) q[7];
rzz(pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,2.123088315295982,-2.123088315295982) q[8];
u3(pi/2,4.871981887187051,-4.871981887187051) q[8];
rzz(pi/2) q[6],q[8];
u3(pi/2,0.5642300405847269,-0.5642300405847269) q[7];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,4.871981887187051,-4.871981887187051) q[8];
u3(pi/2,0.9449910701998098,-0.9449910701998098) q[8];
rzz(-pi/2) q[7],q[8];
u3(pi,0,0) q[1];
u3(pi,3*pi/2,-3*pi/2) q[2];
u3(pi,5*pi/4,-5*pi/4) q[3];
u3(pi,7*pi/4,-7*pi/4) q[4];
u3(pi,3*pi/4,-3*pi/4) q[5];
u3(pi,4.417707589477967,-4.417707589477967) q[6];
u3(pi,0.14702653618800232,-0.14702653618800232) q[7];
u3(pi/2,2.515787396994706,-2.515787396994706) q[8];
u3(pi/2,2.5277254490783476,-2.5277254490783476) q[8];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
measure q[8] -> c[8];
