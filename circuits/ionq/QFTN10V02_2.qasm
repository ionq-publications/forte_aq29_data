OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(pi/2,2.3932652835047046,-2.3932652835047046) q[9];
u3(pi/2,3.178663446902153,-3.178663446902153) q[9];
rzz(-pi/2) q[8],q[9];
rzz(-pi/2) q[7],q[9];
u3(pi/2,0.03707079331235956,-0.03707079331235956) q[9];
u3(pi/2,2.7859643652034287,-2.7859643652034287) q[9];
rzz(pi/2) q[7],q[9];
rzz(pi/2) q[7],q[8];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[8];
u3(pi/2,1.3251237812841747,-1.3251237812841747) q[8];
rzz(-pi/2) q[7],q[8];
rzz(-pi/2) q[6],q[9];
u3(pi/2,2.7859643652034287,-2.7859643652034287) q[9];
u3(pi/2,5.7308933186785005,-5.7308933186785005) q[9];
rzz(pi/2) q[6],q[9];
rzz(-pi/2) q[6],q[8];
u3(pi/2,1.3251237812841747,-1.3251237812841747) q[8];
u3(pi/2,4.074017353175243,-4.074017353175243) q[8];
rzz(-pi/2) q[6],q[8];
rzz(pi/2) q[6],q[7];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[7];
u3(pi/2,2.4542121809843462,-2.4542121809843462) q[7];
rzz(-pi/2) q[6],q[7];
rzz(pi/2) q[5],q[9];
u3(pi/2,2.589300665088708,-2.589300665088708) q[9];
u3(pi/2,5.632875627886499,-5.632875627886499) q[9];
rzz(-pi/2) q[5],q[9];
rzz(-pi/2) q[5],q[8];
u3(pi/2,0.9324246995854506,-0.9324246995854506) q[8];
u3(pi/2,3.8779819715912405,-3.8779819715912405) q[8];
rzz(-pi/2) q[5],q[8];
rzz(pi/2) q[5],q[7];
u3(pi/2,5.595804834574139,-5.595804834574139) q[7];
u3(pi/2,2.061513099285622,-2.061513099285622) q[7];
rzz(-pi/2) q[5],q[7];
rzz(pi/2) q[5],q[6];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
u3(pi/2,4.516353598800687,-4.516353598800687) q[6];
rzz(-pi/2) q[5],q[6];
rzz(-pi/2) q[4],q[9];
u3(pi/2,5.632875627886499,-5.632875627886499) q[9];
u3(pi/2,2.442274128900705,-2.442274128900705) q[9];
rzz(pi/2) q[4],q[9];
rzz(-pi/2) q[4],q[8];
u3(pi/2,3.8779819715912405,-3.8779819715912405) q[8];
u3(pi/2,0.638371627209446,-0.638371627209446) q[8];
rzz(-pi/2) q[4],q[8];
rzz(-pi/2) q[4],q[7];
u3(pi/2,2.061513099285622,-2.061513099285622) q[7];
u3(pi/2,5.007070371291412,-5.007070371291412) q[7];
rzz(pi/2) q[4],q[7];
rzz(-pi/2) q[4],q[6];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[6];
u3(pi/2,4.123026198571244,-4.123026198571244) q[6];
rzz(-pi/2) q[4],q[6];
rzz(pi/2) q[4],q[5];
u3(pi/2,13*pi/8,-13*pi/8) q[5];
u3(pi/2,15*pi/8,-15*pi/8) q[5];
rzz(pi/2) q[4],q[5];
rzz(-pi/2) q[3],q[9];
u3(pi/2,2.442274128900705,-2.442274128900705) q[9];
u3(pi/2,5.559362359792498,-5.559362359792498) q[9];
rzz(-pi/2) q[3],q[9];
rzz(-pi/2) q[3],q[8];
u3(pi/2,3.779964280799239,-3.779964280799239) q[8];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[8];
rzz(pi/2) q[3],q[8];
rzz(-pi/2) q[3],q[7];
u3(pi/2,5.007070371291412,-5.007070371291412) q[7];
u3(pi/2,1.7674600269096177,-1.7674600269096177) q[7];
rzz(-pi/2) q[3],q[7];
rzz(-pi/2) q[3],q[6];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[6];
u3(pi/2,5*pi/4,-5*pi/4) q[6];
rzz(-pi/2) q[3],q[6];
rzz(-pi/2) q[3],q[5];
u3(pi/2,7*pi/8,-7*pi/8) q[5];
u3(pi/2,7*pi/4,-7*pi/4) q[5];
rzz(-pi/2) q[3],q[5];
rzz(pi/2) q[3],q[4];
u3(pi/2,pi/4,-pi/4) q[4];
u3(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[3],q[4];
rzz(pi/2) q[2],q[9];
u3(pi/2,2.4177697062027046,-2.4177697062027046) q[9];
u3(pi/2,5.546795989178139,-5.546795989178139) q[9];
rzz(-pi/2) q[2],q[9];
rzz(-pi/2) q[2],q[8];
u3(pi/2,3.730955435403238,-3.730955435403238) q[8];
u3(pi/2,0.5642300405847269,-0.5642300405847269) q[8];
rzz(-pi/2) q[2],q[8];
rzz(pi/2) q[2],q[7];
u3(pi/2,4.909052680499411,-4.909052680499411) q[7];
u3(pi/2,1.7178228629828987,-1.7178228629828987) q[7];
rzz(-pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[6];
u3(pi/2,5*pi/4,-5*pi/4) q[6];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[6];
rzz(-pi/2) q[2],q[6];
rzz(-pi/2) q[2],q[5];
u3(pi/2,7*pi/4,-7*pi/4) q[5];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[5];
rzz(pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[4];
u3(pi/2,3*pi/8,-3*pi/8) q[4];
rzz(-pi/2) q[2],q[4];
rzz(pi/2) q[2],q[3];
u3(pi/2,0,0) q[3];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(pi/2) q[2],q[3];
rzz(-pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[8];
u3(pi/2,3.70582269417452,-3.70582269417452) q[8];
u3(pi/2,0.5522919885010856,-0.5522919885010856) q[8];
rzz(pi/2) q[1],q[8];
rzz(-pi/2) q[1],q[7];
u3(pi/2,4.8594155165726916,-4.8594155165726916) q[7];
u3(pi/2,1.6933184402848986,-1.6933184402848986) q[7];
rzz(-pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[6];
u3(pi/2,3.82897312619524,-3.82897312619524) q[6];
u3(pi/2,0.638371627209446,-0.638371627209446) q[6];
rzz(pi/2) q[1],q[6];
rzz(-pi/2) q[1],q[5];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[5];
u3(pi/2,5.203105752875415,-5.203105752875415) q[5];
rzz(-pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[4];
u3(pi/2,11*pi/8,-11*pi/8) q[4];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[4];
rzz(-pi/2) q[1],q[4];
rzz(pi/2) q[1],q[3];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,9*pi/8,-9*pi/8) q[3];
rzz(-pi/2) q[1],q[3];
rzz(pi/2) q[1],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(-pi/2) q[1],q[2];
rzz(pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[8];
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
u3(pi/2,pi/8,-pi/8) q[3];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,2.061513099285622,-2.061513099285622) q[5];
u3(pi/2,13*pi/8,-13*pi/8) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,0.638371627209446,-0.638371627209446) q[6];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,4.834911093874691,-4.834911093874691) q[7];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,0.5522919885010856,-0.5522919885010856) q[8];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[8];
rzz(pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[9];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,0,0) q[3];
u3(pi/2,7*pi/8,-7*pi/8) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,pi/4,-pi/4) q[4];
u3(pi/2,3.730955435403238,-3.730955435403238) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,5*pi/8,-5*pi/8) q[5];
u3(pi/2,5.007070371291412,-5.007070371291412) q[5];
rzz(pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[6];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[6];
rzz(pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u3(pi/2,4.810406671176691,-4.810406671176691) q[7];
u3(pi/2,1.6443095948888977,-1.6443095948888977) q[7];
rzz(pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[8];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[8];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[8];
rzz(pi/2) q[1],q[8];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[9];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,7*pi/8,-7*pi/8) q[3];
u3(pi/2,13*pi/8,-13*pi/8) q[3];
rzz(-pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,3.730955435403238,-3.730955435403238) q[4];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,5.007070371291412,-5.007070371291412) q[5];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[5];
rzz(-pi/2) q[2],q[5];
rzz(pi/2) q[2],q[6];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[6];
u3(pi/2,3.583300580684518,-3.583300580684518) q[6];
rzz(pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,1.6443095948888977,-1.6443095948888977) q[7];
u3(pi/2,4.73689340308269,-4.73689340308269) q[7];
rzz(pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[8];
u3(pi/2,3.669380219392878,-3.669380219392878) q[8];
u3(pi/2,0.5032831431050849,-0.5032831431050849) q[8];
rzz(pi/2) q[2],q[8];
rzz(pi/2) q[2],q[9];
u3(pi/2,5.528574751787318,-5.528574751787318) q[9];
u3(pi/2,2.3744157275831657,-2.3744157275831657) q[9];
rzz(-pi/2) q[2],q[9];
u3(pi/2,pi/8,-pi/8) q[3];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[4];
u3(pi/2,2.552229871776348,-2.552229871776348) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,4.810406671176691,-4.810406671176691) q[5];
u3(pi/2,1.276114935888174,-1.276114935888174) q[5];
rzz(pi/2) q[3],q[5];
rzz(-pi/2) q[3],q[6];
u3(pi/2,0.4417079270947249,-0.4417079270947249) q[6];
u3(pi/2,3.387265199100515,-3.387265199100515) q[6];
rzz(pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,4.73689340308269,-4.73689340308269) q[7];
u3(pi/2,1.4972830587008954,-1.4972830587008954) q[7];
rzz(pi/2) q[3],q[7];
rzz(-pi/2) q[3],q[8];
u3(pi/2,3.6448757966948775,-3.6448757966948775) q[8];
u3(pi/2,0.4542742977090841,-0.4542742977090841) q[8];
rzz(pi/2) q[3],q[8];
rzz(pi/2) q[3],q[9];
u3(pi/2,5.516008381172959,-5.516008381172959) q[9];
u3(pi/2,2.349911304885165,-2.349911304885165) q[9];
rzz(pi/2) q[3],q[9];
u3(pi/2,4.123026198571244,-4.123026198571244) q[4];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,1.276114935888174,-1.276114935888174) q[5];
u3(pi/2,3.6323094260805187,-3.6323094260805187) q[5];
rzz(-pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,3.387265199100515,-3.387265199100515) q[6];
u3(pi/2,6.136158770991584,-6.136158770991584) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,1.4972830587008954,-1.4972830587008954) q[7];
u3(pi/2,4.442212012175967,-4.442212012175967) q[7];
rzz(-pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,0.4542742977090841,-0.4542742977090841) q[8];
u3(pi/2,3.4972209419761575,-3.4972209419761575) q[8];
rzz(pi/2) q[4],q[8];
rzz(pi/2) q[4],q[9];
u3(pi/2,2.349911304885165,-2.349911304885165) q[9];
u3(pi/2,5.442495113078958,-5.442495113078958) q[9];
rzz(pi/2) q[4],q[9];
u3(pi/2,2.061513099285622,-2.061513099285622) q[5];
u3(pi/2,13*pi/8,-13*pi/8) q[5];
rzz(pi/2) q[5],q[6];
u3(pi/2,6.136158770991584,-6.136158770991584) q[6];
u3(pi/2,2.2091679540043425,-2.2091679540043425) q[6];
rzz(pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,1.3006193585861743,-1.3006193585861743) q[7];
u3(pi/2,4.049512930477243,-4.049512930477243) q[7];
rzz(pi/2) q[5],q[7];
rzz(-pi/2) q[5],q[8];
u3(pi/2,0.35562828838636457,-0.35562828838636457) q[8];
u3(pi/2,3.3011855603921543,-3.3011855603921543) q[8];
rzz(pi/2) q[5],q[8];
rzz(pi/2) q[5],q[9];
u3(pi/2,5.442495113078958,-5.442495113078958) q[9];
u3(pi/2,2.202884768697163,-2.202884768697163) q[9];
rzz(pi/2) q[5],q[9];
u3(pi/2,0.638371627209446,-0.638371627209446) q[6];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[6];
rzz(pi/2) q[6],q[7];
u3(pi/2,4.049512930477243,-4.049512930477243) q[7];
u3(pi/2,0.12252211349000193,-0.12252211349000193) q[7];
rzz(pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,3.3011855603921543,-3.3011855603921543) q[8];
u3(pi/2,6.050079132283224,-6.050079132283224) q[8];
rzz(pi/2) q[6],q[8];
rzz(-pi/2) q[6],q[9];
u3(pi/2,5.344477422286956,-5.344477422286956) q[9];
u3(pi/2,2.0062210685824415,-2.0062210685824415) q[9];
rzz(pi/2) q[6],q[9];
u3(pi/2,1.6933184402848986,-1.6933184402848986) q[7];
u3(pi/2,4.810406671176691,-4.810406671176691) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,6.050079132283224,-6.050079132283224) q[8];
u3(pi/2,2.123088315295982,-2.123088315295982) q[8];
rzz(pi/2) q[7],q[8];
rzz(-pi/2) q[7],q[9];
u3(pi/2,5.147813722172235,-5.147813722172235) q[9];
u3(pi/2,1.6135219868837176,-1.6135219868837176) q[9];
rzz(pi/2) q[7],q[9];
u3(pi/2,3.6938846420908784,-3.6938846420908784) q[8];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[8];
rzz(pi/2) q[8],q[9];
u3(pi/2,1.6135219868837176,-1.6135219868837176) q[9];
u3(pi/2,3.9697164770760627,-3.9697164770760627) q[9];
rzz(pi/2) q[8],q[9];
u3(pi,0,0) q[1];
u3(pi,pi/2,-pi/2) q[2];
u3(pi,5*pi/4,-5*pi/4) q[3];
u3(pi,7*pi/4,-7*pi/4) q[4];
u3(pi,15*pi/8,-15*pi/8) q[5];
u3(pi,1.0800795543041708,-1.0800795543041708) q[6];
u3(pi,0.2946813909067226,-0.2946813909067226) q[7];
u3(pi,1.4237697906068942,-1.4237697906068942) q[8];
u3(pi,2.3989201502811657,-2.3989201502811657) q[9];
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
