OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(pi/2,5.488362365821369,-5.488362365821369) q[12];
u3(pi/2,6.273760529218817,-6.273760529218817) q[12];
rzz(-pi/2) q[11],q[12];
rzz(-pi/2) q[10],q[12];
u3(pi/2,3.132167875629024,-3.132167875629024) q[12];
u3(pi/2,5.881061447520093,-5.881061447520093) q[12];
rzz(-pi/2) q[10],q[12];
rzz(-pi/2) q[10],q[11];
u3(pi/2,1.0398671683382217,-1.0398671683382217) q[11];
u3(pi/2,1.8252653317356697,-1.8252653317356697) q[11];
rzz(-pi/2) q[10],q[11];
rzz(-pi/2) q[9],q[12];
u3(pi/2,5.881061447520093,-5.881061447520093) q[12];
u3(pi/2,2.5434334123462965,-2.5434334123462965) q[12];
rzz(-pi/2) q[9],q[12];
rzz(-pi/2) q[9],q[11];
u3(pi/2,4.966857985325463,-4.966857985325463) q[11];
u3(pi/2,1.4325662500369458,-1.4325662500369458) q[11];
rzz(pi/2) q[9],q[11];
rzz(pi/2) q[9],q[10];
u3(pi/2,4.69982260977033,-4.69982260977033) q[10];
u3(pi/2,5.485220773167779,-5.485220773167779) q[10];
rzz(-pi/2) q[9],q[10];
rzz(-pi/2) q[8],q[12];
u3(pi/2,5.68502606593609,-5.68502606593609) q[12];
u3(pi/2,2.445415721554295,-2.445415721554295) q[12];
rzz(pi/2) q[8],q[12];
rzz(-pi/2) q[8],q[11];
u3(pi/2,1.4325662500369458,-1.4325662500369458) q[11];
u3(pi/2,4.378123522042736,-4.378123522042736) q[11];
rzz(-pi/2) q[8],q[11];
rzz(-pi/2) q[8],q[10];
u3(pi/2,5.485220773167779,-5.485220773167779) q[10];
u3(pi/2,1.9509290378792614,-1.9509290378792614) q[10];
rzz(-pi/2) q[8],q[10];
rzz(-pi/2) q[8],q[9];
u3(pi/2,2.589300665088708,-2.589300665088708) q[9];
u3(pi/2,3.374698828486156,-3.374698828486156) q[9];
rzz(-pi/2) q[8],q[9];
rzz(-pi/2) q[7],q[12];
u3(pi/2,5.587008375144088,-5.587008375144088) q[12];
u3(pi/2,2.395778557627576,-2.395778557627576) q[12];
rzz(-pi/2) q[7],q[12];
rzz(pi/2) q[7],q[11];
u3(pi/2,1.2365308684529426,-1.2365308684529426) q[11];
u3(pi/2,4.280105831250735,-4.280105831250735) q[11];
rzz(-pi/2) q[7],q[11];
rzz(-pi/2) q[7],q[10];
u3(pi/2,1.9509290378792614,-1.9509290378792614) q[10];
u3(pi/2,4.896486309885051,-4.896486309885051) q[10];
rzz(-pi/2) q[7],q[10];
rzz(pi/2) q[7],q[9];
u3(pi/2,3.374698828486156,-3.374698828486156) q[9];
u3(pi/2,6.123592400377225,-6.123592400377225) q[9];
rzz(-pi/2) q[7],q[9];
rzz(pi/2) q[7],q[8];
u3(pi/2,2.577362613005066,-2.577362613005066) q[8];
u3(pi/2,3.3627607764025145,-3.3627607764025145) q[8];
rzz(-pi/2) q[7],q[8];
rzz(pi/2) q[6],q[12];
u3(pi/2,5.537371211217369,-5.537371211217369) q[12];
u3(pi/2,2.371274134929576,-2.371274134929576) q[12];
rzz(-pi/2) q[6],q[12];
rzz(-pi/2) q[6],q[11];
u3(pi/2,4.280105831250735,-4.280105831250735) q[11];
u3(pi/2,1.0888760137342224,-1.0888760137342224) q[11];
rzz(-pi/2) q[6],q[11];
rzz(-pi/2) q[6],q[10];
u3(pi/2,4.896486309885051,-4.896486309885051) q[10];
u3(pi/2,1.6568759655032568,-1.6568759655032568) q[10];
rzz(-pi/2) q[6],q[10];
rzz(-pi/2) q[6],q[9];
u3(pi/2,6.123592400377225,-6.123592400377225) q[9];
u3(pi/2,2.7859643652034287,-2.7859643652034287) q[9];
rzz(-pi/2) q[6],q[9];
rzz(-pi/2) q[6],q[8];
u3(pi/2,0.22116812281272144,-0.22116812281272144) q[8];
u3(pi/2,2.9700616947037903,-2.9700616947037903) q[8];
rzz(pi/2) q[6],q[8];
rzz(pi/2) q[6],q[7];
u3(pi/2,2.503221026380347,-2.503221026380347) q[7];
u3(pi/2,3.2886191897777954,-3.2886191897777954) q[7];
rzz(-pi/2) q[6],q[7];
rzz(-pi/2) q[5],q[12];
u3(pi/2,2.371274134929576,-2.371274134929576) q[12];
u3(pi/2,5.5009287364357276,-5.5009287364357276) q[12];
rzz(pi/2) q[5],q[12];
rzz(-pi/2) q[5],q[11];
u3(pi/2,1.0888760137342224,-1.0888760137342224) q[11];
u3(pi/2,4.205964244626015,-4.205964244626015) q[11];
rzz(-pi/2) q[5],q[11];
rzz(-pi/2) q[5],q[10];
u3(pi/2,1.6568759655032568,-1.6568759655032568) q[10];
u3(pi/2,4.749459773697049,-4.749459773697049) q[10];
rzz(-pi/2) q[5],q[10];
rzz(pi/2) q[5],q[9];
u3(pi/2,5.927557018793221,-5.927557018793221) q[9];
u3(pi/2,2.687318355880709,-2.687318355880709) q[9];
rzz(-pi/2) q[5],q[9];
rzz(-pi/2) q[5],q[8];
u3(pi/2,6.111654348293584,-6.111654348293584) q[8];
u3(pi/2,2.7733979945890694,-2.7733979945890694) q[8];
rzz(-pi/2) q[5],q[8];
rzz(pi/2) q[5],q[7];
u3(pi/2,3.2886191897777954,-3.2886191897777954) q[7];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
rzz(-pi/2) q[5],q[7];
rzz(pi/2) q[5],q[6];
u3(pi/2,2.4542121809843462,-2.4542121809843462) q[6];
u3(pi/2,3.2396103443817945,-3.2396103443817945) q[6];
rzz(-pi/2) q[5],q[6];
rzz(pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[11];
u3(pi/2,4.205964244626015,-4.205964244626015) q[11];
u3(pi/2,1.0524335389525807,-1.0524335389525807) q[11];
rzz(-pi/2) q[4],q[11];
rzz(pi/2) q[4],q[10];
u3(pi/2,1.6078671201072563,-1.6078671201072563) q[10];
u3(pi/2,4.724955350999049,-4.724955350999049) q[10];
rzz(-pi/2) q[4],q[10];
rzz(-pi/2) q[4],q[9];
u3(pi/2,2.687318355880709,-2.687318355880709) q[9];
u3(pi/2,5.779902164074501,-5.779902164074501) q[9];
rzz(-pi/2) q[4],q[9];
rzz(-pi/2) q[4],q[8];
u3(pi/2,2.7733979945890694,-2.7733979945890694) q[8];
u3(pi/2,5.81697295738686,-5.81697295738686) q[8];
rzz(pi/2) q[4],q[8];
rzz(-pi/2) q[4],q[7];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
u3(pi/2,2.6998847264950685,-2.6998847264950685) q[7];
rzz(-pi/2) q[4],q[7];
rzz(-pi/2) q[4],q[6];
u3(pi/2,0.09801769079200154,-0.09801769079200154) q[6];
u3(pi/2,2.8469112626830704,-2.8469112626830704) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[5];
u3(pi/2,5.301751762198135,-5.301751762198135) q[5];
u3(pi/2,6.087149925595583,-6.087149925595583) q[5];
rzz(-pi/2) q[4],q[5];
rzz(-pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[12];
rzz(pi/2) q[3],q[11];
rzz(-pi/2) q[3],q[11];
rzz(-pi/2) q[3],q[10];
u3(pi/2,4.724955350999049,-4.724955350999049) q[10];
u3(pi/2,pi/2,-pi/2) q[10];
rzz(-pi/2) q[3],q[10];
rzz(pi/2) q[3],q[9];
u3(pi/2,2.6383095104847083,-2.6383095104847083) q[9];
u3(pi/2,5.755397741376501,-5.755397741376501) q[9];
rzz(-pi/2) q[3],q[9];
rzz(-pi/2) q[3],q[8];
u3(pi/2,2.675380303797068,-2.675380303797068) q[8];
u3(pi/2,5.76796411199086,-5.76796411199086) q[8];
rzz(-pi/2) q[3],q[8];
rzz(-pi/2) q[3],q[7];
u3(pi/2,2.6998847264950685,-2.6998847264950685) q[7];
u3(pi/2,5.74345968929286,-5.74345968929286) q[7];
rzz(-pi/2) q[3],q[7];
rzz(-pi/2) q[3],q[6];
u3(pi/2,5.9885039162728635,-5.9885039162728635) q[6];
u3(pi/2,2.6508758810990676,-2.6508758810990676) q[6];
rzz(-pi/2) q[3],q[6];
rzz(pi/2) q[3],q[5];
u3(pi/2,6.087149925595583,-6.087149925595583) q[5];
u3(pi/2,2.552229871776348,-2.552229871776348) q[5];
rzz(-pi/2) q[3],q[5];
rzz(pi/2) q[3],q[4];
u3(pi/2,5*pi/8,-5*pi/8) q[4];
u3(pi/2,7*pi/8,-7*pi/8) q[4];
rzz(-pi/2) q[3],q[4];
rzz(-pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[11];
rzz(-pi/2) q[2],q[11];
rzz(-pi/2) q[2],q[10];
rzz(pi/2) q[2],q[10];
rzz(-pi/2) q[2],q[9];
u3(pi/2,5.755397741376501,-5.755397741376501) q[9];
u3(pi/2,2.6018670357030667,-2.6018670357030667) q[9];
rzz(-pi/2) q[2],q[9];
rzz(-pi/2) q[2],q[8];
u3(pi/2,5.76796411199086,-5.76796411199086) q[8];
u3(pi/2,2.6018670357030667,-2.6018670357030667) q[8];
rzz(pi/2) q[2],q[8];
rzz(-pi/2) q[2],q[7];
u3(pi/2,5.74345968929286,-5.74345968929286) q[7];
u3(pi/2,2.5528581903070657,-2.5528581903070657) q[7];
rzz(-pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[6];
u3(pi/2,2.6508758810990676,-2.6508758810990676) q[6];
u3(pi/2,5.694450843896859,-5.694450843896859) q[6];
rzz(-pi/2) q[2],q[6];
rzz(pi/2) q[2],q[5];
u3(pi/2,5.693822525366141,-5.693822525366141) q[5];
u3(pi/2,3*pi/4,-3*pi/4) q[5];
rzz(-pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[4];
u3(pi/2,15*pi/8,-15*pi/8) q[4];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[2],q[4];
rzz(-pi/2) q[2],q[3];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(-pi/2) q[2],q[3];
rzz(-pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[11];
rzz(pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[10];
rzz(-pi/2) q[1],q[10];
rzz(pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[8];
u3(pi/2,5.74345968929286,-5.74345968929286) q[8];
u3(pi/2,2.589300665088708,-2.589300665088708) q[8];
rzz(-pi/2) q[1],q[8];
rzz(-pi/2) q[1],q[7];
u3(pi/2,2.5528581903070657,-2.5528581903070657) q[7];
u3(pi/2,5.66931810266814,-5.66931810266814) q[7];
rzz(pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[6];
u3(pi/2,5.693822525366141,-5.693822525366141) q[6];
u3(pi/2,2.503221026380347,-2.503221026380347) q[6];
rzz(-pi/2) q[1],q[6];
rzz(-pi/2) q[1],q[5];
u3(pi/2,3*pi/4,-3*pi/4) q[5];
u3(pi/2,5.399769452990137,-5.399769452990137) q[5];
rzz(pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[4];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
u3(pi/2,5.301751762198135,-5.301751762198135) q[4];
rzz(-pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(-pi/2) q[1],q[3];
rzz(-pi/2) q[1],q[2];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(-pi/2) q[1],q[2];
rzz(-pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[11];
rzz(pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[7];
rzz(-pi/2) q[0],q[6];
rzz(-pi/2) q[0],q[5];
rzz(-pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[3];
rzz(-pi/2) q[0],q[2];
rzz(-pi/2) q[0],q[1];
u3(pi/2,0,0) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
u3(pi/2,0,0) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[4];
u3(pi/2,13*pi/8,-13*pi/8) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,5.399769452990137,-5.399769452990137) q[5];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,5.64481367997014,-5.64481367997014) q[6];
u3(pi/2,2.4542121809843462,-2.4542121809843462) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,5.66931810266814,-5.66931810266814) q[7];
u3(pi/2,2.503221026380347,-2.503221026380347) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,5.7308933186785005,-5.7308933186785005) q[8];
u3(pi/2,2.577362613005066,-2.577362613005066) q[8];
rzz(-pi/2) q[0],q[8];
rzz(pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[12];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
u3(pi/2,pi/8,-pi/8) q[3];
rzz(pi/2) q[1],q[3];
rzz(-pi/2) q[1],q[4];
u3(pi/2,13*pi/8,-13*pi/8) q[4];
u3(pi/2,1.7674600269096177,-1.7674600269096177) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,5.301751762198135,-5.301751762198135) q[5];
u3(pi/2,2.061513099285622,-2.061513099285622) q[5];
rzz(pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[6];
u3(pi/2,2.4542121809843462,-2.4542121809843462) q[6];
u3(pi/2,5.546795989178139,-5.546795989178139) q[6];
rzz(pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u3(pi/2,5.64481367997014,-5.64481367997014) q[7];
u3(pi/2,2.4787166036823467,-2.4787166036823467) q[7];
rzz(pi/2) q[1],q[7];
rzz(pi/2) q[1],q[8];
u3(pi/2,2.577362613005066,-2.577362613005066) q[8];
u3(pi/2,5.7063888959805,-5.7063888959805) q[8];
rzz(-pi/2) q[1],q[8];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[10];
rzz(-pi/2) q[1],q[10];
rzz(pi/2) q[1],q[11];
rzz(pi/2) q[1],q[11];
rzz(pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[12];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,pi/8,-pi/8) q[3];
u3(pi/2,7*pi/8,-7*pi/8) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[4];
u3(pi/2,4.516353598800687,-4.516353598800687) q[4];
rzz(-pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,2.061513099285622,-2.061513099285622) q[5];
u3(pi/2,5.007070371291412,-5.007070371291412) q[5];
rzz(pi/2) q[2],q[5];
rzz(pi/2) q[2],q[6];
u3(pi/2,5.546795989178139,-5.546795989178139) q[6];
u3(pi/2,2.3071856447963444,-2.3071856447963444) q[6];
rzz(pi/2) q[2],q[6];
rzz(-pi/2) q[2],q[7];
u3(pi/2,5.620309257272139,-5.620309257272139) q[7];
u3(pi/2,2.4297077582863458,-2.4297077582863458) q[7];
rzz(pi/2) q[2],q[7];
rzz(pi/2) q[2],q[8];
u3(pi/2,2.5647962423907074,-2.5647962423907074) q[8];
u3(pi/2,5.6818844732825,-5.6818844732825) q[8];
rzz(pi/2) q[2],q[8];
rzz(-pi/2) q[2],q[9];
u3(pi/2,5.724610133371321,-5.724610133371321) q[9];
u3(pi/2,2.5710794276978866,-2.5710794276978866) q[9];
rzz(pi/2) q[2],q[9];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[12];
rzz(pi/2) q[2],q[12];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[4];
u3(pi/2,3.730955435403238,-3.730955435403238) q[4];
rzz(-pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,5.007070371291412,-5.007070371291412) q[5];
u3(pi/2,1.472778636002895,-1.472778636002895) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,2.3071856447963444,-2.3071856447963444) q[6];
u3(pi/2,5.252114598271416,-5.252114598271416) q[6];
rzz(-pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,2.4297077582863458,-2.4297077582863458) q[7];
u3(pi/2,5.4732827210841375,-5.4732827210841375) q[7];
rzz(pi/2) q[3],q[7];
rzz(pi/2) q[3],q[8];
u3(pi/2,5.6818844732825,-5.6818844732825) q[8];
u3(pi/2,2.491282974296706,-2.491282974296706) q[8];
rzz(-pi/2) q[3],q[8];
rzz(pi/2) q[3],q[9];
u3(pi/2,2.5710794276978866,-2.5710794276978866) q[9];
u3(pi/2,5.68816765858968,-5.68816765858968) q[9];
rzz(pi/2) q[3],q[9];
rzz(pi/2) q[3],q[10];
u3(pi/2,1.561371548834127,-1.561371548834127) q[10];
u3(pi/2,4.69102615034028,-4.69102615034028) q[10];
rzz(-pi/2) q[3],q[10];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[12];
rzz(pi/2) q[3],q[12];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[4];
u3(pi/2,13*pi/8,-13*pi/8) q[4];
rzz(-pi/2) q[4],q[5];
u3(pi/2,4.614371289592689,-4.614371289592689) q[5];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,2.110521944681623,-2.110521944681623) q[6];
u3(pi/2,4.8594155165726916,-4.8594155165726916) q[6];
rzz(pi/2) q[4],q[6];
rzz(-pi/2) q[4],q[7];
u3(pi/2,2.3316900674943444,-2.3316900674943444) q[7];
u3(pi/2,5.276619020969417,-5.276619020969417) q[7];
rzz(pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,5.632875627886499,-5.632875627886499) q[8];
u3(pi/2,2.3932652835047046,-2.3932652835047046) q[8];
rzz(pi/2) q[4],q[8];
rzz(pi/2) q[4],q[9];
u3(pi/2,5.68816765858968,-5.68816765858968) q[9];
u3(pi/2,2.4975661596038856,-2.4975661596038856) q[9];
rzz(-pi/2) q[4],q[9];
rzz(pi/2) q[4],q[10];
u3(pi/2,1.5494334967504861,-1.5494334967504861) q[10];
u3(pi/2,4.666521727642279,-4.666521727642279) q[10];
rzz(pi/2) q[4],q[10];
rzz(pi/2) q[4],q[11];
u3(pi/2,4.180203184866579,-4.180203184866579) q[11];
u3(pi/2,1.0260441606624264,-1.0260441606624264) q[11];
rzz(-pi/2) q[4],q[11];
rzz(pi/2) q[4],q[12];
rzz(pi/2) q[4],q[12];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[5];
u3(pi/2,5.301751762198135,-5.301751762198135) q[5];
rzz(pi/2) q[5],q[6];
u3(pi/2,4.8594155165726916,-4.8594155165726916) q[6];
u3(pi/2,0.9324246995854506,-0.9324246995854506) q[6];
rzz(pi/2) q[5],q[6];
rzz(-pi/2) q[5],q[7];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[7];
u3(pi/2,4.883919939270692,-4.883919939270692) q[7];
rzz(pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,2.3932652835047046,-2.3932652835047046) q[8];
u3(pi/2,5.338194236979777,-5.338194236979777) q[8];
rzz(-pi/2) q[5],q[8];
rzz(pi/2) q[5],q[9];
u3(pi/2,5.639158813193679,-5.639158813193679) q[9];
u3(pi/2,2.3989201502811657,-2.3989201502811657) q[9];
rzz(pi/2) q[5],q[9];
rzz(pi/2) q[5],q[10];
u3(pi/2,4.666521727642279,-4.666521727642279) q[10];
u3(pi/2,1.4759202286564848,-1.4759202286564848) q[10];
rzz(pi/2) q[5],q[10];
rzz(-pi/2) q[5],q[11];
u3(pi/2,1.0260441606624264,-1.0260441606624264) q[11];
u3(pi/2,4.143132391554219,-4.143132391554219) q[11];
rzz(pi/2) q[5],q[11];
rzz(pi/2) q[5],q[12];
u3(pi/2,2.3417431639858317,-2.3417431639858317) q[12];
u3(pi/2,5.470769446961266,-5.470769446961266) q[12];
rzz(pi/2) q[5],q[12];
u3(pi/2,2.503221026380347,-2.503221026380347) q[6];
u3(pi/2,5.595804834574139,-5.595804834574139) q[6];
rzz(-pi/2) q[6],q[7];
u3(pi/2,1.7423272856808991,-1.7423272856808991) q[7];
u3(pi/2,4.098521775873244,-4.098521775873244) q[7];
rzz(pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,2.1966015833899837,-2.1966015833899837) q[8];
u3(pi/2,4.945495155281052,-4.945495155281052) q[8];
rzz(pi/2) q[6],q[8];
rzz(pi/2) q[6],q[9];
u3(pi/2,2.3989201502811657,-2.3989201502811657) q[9];
u3(pi/2,5.344477422286956,-5.344477422286956) q[9];
rzz(pi/2) q[6],q[9];
rzz(pi/2) q[6],q[10];
u3(pi/2,1.4759202286564848,-1.4759202286564848) q[10];
u3(pi/2,4.518866872923558,-4.518866872923558) q[10];
rzz(pi/2) q[6],q[10];
rzz(pi/2) q[6],q[11];
u3(pi/2,4.143132391554219,-4.143132391554219) q[11];
u3(pi/2,0.9525308925684254,-0.9525308925684254) q[11];
rzz(-pi/2) q[6],q[11];
rzz(pi/2) q[6],q[12];
u3(pi/2,5.470769446961266,-5.470769446961266) q[12];
u3(pi/2,2.3046723706734724,-2.3046723706734724) q[12];
rzz(pi/2) q[6],q[12];
u3(pi/2,5.66931810266814,-5.66931810266814) q[7];
u3(pi/2,2.503221026380347,-2.503221026380347) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,4.945495155281052,-4.945495155281052) q[8];
u3(pi/2,1.0185043382938108,-1.0185043382938108) q[8];
rzz(-pi/2) q[7],q[8];
rzz(pi/2) q[7],q[9];
u3(pi/2,5.344477422286956,-5.344477422286956) q[9];
u3(pi/2,1.810185686998439,-1.810185686998439) q[9];
rzz(pi/2) q[7],q[9];
rzz(pi/2) q[7],q[10];
u3(pi/2,4.518866872923558,-4.518866872923558) q[10];
u3(pi/2,1.1812388377497622,-1.1812388377497622) q[10];
rzz(pi/2) q[7],q[10];
rzz(-pi/2) q[7],q[11];
u3(pi/2,0.9525308925684254,-0.9525308925684254) q[11];
u3(pi/2,3.996105855366217,-3.996105855366217) q[11];
rzz(pi/2) q[7],q[11];
rzz(pi/2) q[7],q[12];
u3(pi/2,2.3046723706734724,-2.3046723706734724) q[12];
u3(pi/2,5.397256178867265,-5.397256178867265) q[12];
rzz(-pi/2) q[7],q[12];
u3(pi/2,5.7308933186785005,-5.7308933186785005) q[8];
u3(pi/2,2.577362613005066,-2.577362613005066) q[8];
rzz(pi/2) q[8],q[9];
u3(pi/2,1.810185686998439,-1.810185686998439) q[9];
u3(pi/2,4.166380177190784,-4.166380177190784) q[9];
rzz(pi/2) q[8],q[9];
rzz(pi/2) q[8],q[10];
u3(pi/2,1.1812388377497622,-1.1812388377497622) q[10];
u3(pi/2,3.930132409640831,-3.930132409640831) q[10];
rzz(pi/2) q[8],q[10];
rzz(-pi/2) q[8],q[11];
u3(pi/2,0.8545132017764238,-0.8545132017764238) q[11];
u3(pi/2,3.799442155251496,-3.799442155251496) q[11];
rzz(pi/2) q[8],q[11];
rzz(pi/2) q[8],q[12];
u3(pi/2,2.2556635252774715,-2.2556635252774715) q[12];
u3(pi/2,5.299238488075264,-5.299238488075264) q[12];
rzz(pi/2) q[8],q[12];
u3(pi,2.595583850395887,-2.595583850395887) q[9];
rzz(-pi/2) q[9],q[10];
u3(pi/2,0.7885397560510381,-0.7885397560510381) q[10];
u3(pi/2,3.1447342462433827,-3.1447342462433827) q[10];
rzz(pi/2) q[9],q[10];
rzz(pi/2) q[9],q[11];
u3(pi/2,3.799442155251496,-3.799442155251496) q[11];
u3(pi/2,0.26515041996297856,-0.26515041996297856) q[11];
rzz(pi/2) q[9],q[11];
rzz(pi/2) q[9],q[12];
u3(pi/2,5.299238488075264,-5.299238488075264) q[12];
u3(pi/2,1.9609821343707488,-1.9609821343707488) q[12];
rzz(-pi/2) q[9],q[12];
u3(pi,1.5739379194484864,-1.5739379194484864) q[10];
rzz(pi/2) q[10],q[11];
u3(pi/2,0.26515041996297856,-0.26515041996297856) q[11];
u3(pi/2,2.6213449101553237,-2.6213449101553237) q[11];
rzz(pi/2) q[10],q[11];
rzz(pi/2) q[10],q[12];
u3(pi/2,5.102574787960542,-5.102574787960542) q[12];
u3(pi/2,1.5682830526720246,-1.5682830526720246) q[12];
rzz(-pi/2) q[10],q[12];
u3(pi,1.0505485833604267,-1.0505485833604267) q[11];
rzz(pi/2) q[11],q[12];
u3(pi/2,4.709875706261818,-4.709875706261818) q[12];
u3(pi/2,0.7828848892745764,-0.7828848892745764) q[12];
rzz(pi/2) q[11],q[12];
u3(pi,pi/2,-pi/2) q[1];
u3(pi,0,0) q[2];
u3(pi,7*pi/4,-7*pi/4) q[3];
u3(pi,15*pi/8,-15*pi/8) q[4];
u3(pi,4.123026198571244,-4.123026198571244) q[5];
u3(pi,1.865477717701619,-1.865477717701619) q[6];
u3(pi,0.638371627209446,-0.638371627209446) q[7];
u3(pi,2.4297077582863458,-2.4297077582863458) q[8];
u3(pi,4.098521775873244,-4.098521775873244) q[9];
u3(pi,4.68160137237951,-4.68160137237951) q[10];
u3(pi,1.8164688723056186,-1.8164688723056186) q[11];
u3(pi,5.4952738696592665,-5.4952738696592665) q[12];
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
