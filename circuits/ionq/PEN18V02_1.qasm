OPENQASM 2.0;
include "qelib1.inc";
qreg q[18];
creg c[18];
u3(pi/2,3.666238626739289,-3.666238626739289) q[16];
u3(pi/2,5.518521655295831,-5.518521655295831) q[17];
rzz(-pi/2) q[16],q[17];
u3(pi/2,3.9477253285009337,-3.9477253285009337) q[17];
u3(pi/2,4.430902278623044,-4.430902278623044) q[17];
rzz(pi/2) q[16],q[17];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[15];
rzz(pi/2) q[15],q[17];
u3(pi/2,1.289309625033251,-1.289309625033251) q[17];
u3(pi/2,3.464548378378824,-3.464548378378824) q[17];
rzz(pi/2) q[15],q[17];
u3(pi/2,4.389433255595659,-4.389433255595659) q[14];
rzz(pi/2) q[14],q[17];
u3(pi/2,3.464548378378824,-3.464548378378824) q[17];
u3(pi/2,4.674061550010895,-4.674061550010895) q[17];
rzz(-pi/2) q[14],q[17];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[13];
rzz(pi/2) q[13],q[17];
u3(pi/2,4.674061550010895,-4.674061550010895) q[17];
u3(pi/2,5.397256178867265,-5.397256178867265) q[17];
rzz(pi/2) q[13],q[17];
u3(pi/2,3.416796170044259,-3.416796170044259) q[12];
rzz(pi/2) q[12],q[17];
u3(pi/2,2.2556635252774715,-2.2556635252774715) q[17];
u3(pi/2,3.950238602623806,-3.950238602623806) q[17];
rzz(-pi/2) q[12],q[17];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[11];
rzz(pi/2) q[11],q[17];
u3(pi/2,0.8086459490340128,-0.8086459490340128) q[17];
u3(pi/2,1.0574600871983244,-1.0574600871983244) q[17];
rzz(pi/2) q[11],q[17];
u3(pi/2,2.672238711143478,-2.672238711143478) q[10];
rzz(pi/2) q[10],q[17];
u3(pi/2,4.199052740788117,-4.199052740788117) q[17];
u3(pi/2,0.5604601294004191,-0.5604601294004191) q[17];
rzz(pi/2) q[10],q[17];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
rzz(-pi/2) q[9],q[17];
u3(pi/2,3.702052782990212,-3.702052782990212) q[17];
u3(pi/2,5.849645520984195,-5.849645520984195) q[17];
rzz(pi/2) q[9],q[17];
u3(pi/2,5.976565864189222,-5.976565864189222) q[8];
rzz(pi/2) q[8],q[17];
u3(pi/2,5.849645520984195,-5.849645520984195) q[17];
u3(pi/2,0.7200530362027805,-0.7200530362027805) q[17];
rzz(pi/2) q[8],q[17];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
rzz(-pi/2) q[7],q[17];
u3(pi/2,0.7200530362027805,-0.7200530362027805) q[17];
u3(pi/2,1.5544600449962296,-1.5544600449962296) q[17];
rzz(pi/2) q[7],q[17];
u3(pi/2,0.3436902363027234,-0.3436902363027234) q[6];
rzz(pi/2) q[6],q[17];
u3(pi/2,4.696052698586023,-4.696052698586023) q[17];
u3(pi/2,6.168831334588917,-6.168831334588917) q[17];
rzz(pi/2) q[6],q[17];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(pi/2) q[5],q[17];
u3(pi/2,3.0272386809991247,-3.0272386809991247) q[17];
u3(pi/2,3.223274062583128,-3.223274062583128) q[17];
rzz(-pi/2) q[5],q[17];
u3(pi/2,3.3376280351737964,-3.3376280351737964) q[4];
rzz(pi/2) q[4],q[17];
u3(pi/2,3.223274062583128,-3.223274062583128) q[17];
u3(pi/2,5.9721676344741965,-5.9721676344741965) q[17];
rzz(pi/2) q[4],q[17];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(pi/2) q[3],q[17];
u3(pi/2,5.9721676344741965,-5.9721676344741965) q[17];
u3(pi/2,2.0451768174869556,-2.0451768174869556) q[17];
rzz(-pi/2) q[3],q[17];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,0.47438049069205873,-0.47438049069205873) q[17];
rzz(pi/2) q[2],q[17];
u3(pi/2,pi,-pi) q[1];
u3(pi/2,7*pi/4,-7*pi/4) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
u3(pi/2,13*pi/8,-13*pi/8) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,15*pi/8,-15*pi/8) q[3];
u3(pi/2,2.552229871776348,-2.552229871776348) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,4.9084243619686925,-4.9084243619686925) q[4];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,5.399769452990137,-5.399769452990137) q[5];
u3(pi/2,2.2091679540043425,-2.2091679540043425) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,3.4852828898925163,-3.4852828898925163) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,1.91448656309762,-1.91448656309762) q[6];
u3(pi/2,5.031574793989412,-5.031574793989412) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,4.466716434873968,-4.466716434873968) q[7];
u3(pi/2,1.3131857292005336,-1.3131857292005336) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,2.834973210599429,-2.834973210599429) q[8];
rzz(-pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[8];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
rzz(pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[9];
u3(pi/2,5.813831364733272,-5.813831364733272) q[10];
rzz(-pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[10];
u3(pi/2,3.669380219392878,-3.669380219392878) q[11];
rzz(pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[11];
u3(pi/2,3.416796170044259,-3.416796170044259) q[12];
rzz(-pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[12];
u3(pi/2,3.669380219392878,-3.669380219392878) q[13];
rzz(-pi/2) q[0],q[13];
rzz(pi/2) q[0],q[13];
u3(pi/2,4.389433255595659,-4.389433255595659) q[14];
rzz(-pi/2) q[0],q[14];
rzz(-pi/2) q[0],q[14];
u3(pi/2,3.669380219392878,-3.669380219392878) q[15];
rzz(-pi/2) q[0],q[15];
rzz(pi/2) q[0],q[15];
u3(pi/2,3.666238626739289,-3.666238626739289) q[16];
rzz(-pi/2) q[0],q[16];
rzz(-pi/2) q[0],q[16];
u3(pi/2,pi/4,-pi/4) q[1];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,5*pi/8,-5*pi/8) q[2];
u3(pi/2,11*pi/8,-11*pi/8) q[2];
rzz(-pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,2.552229871776348,-2.552229871776348) q[3];
u3(pi/2,5.301751762198135,-5.301751762198135) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,4.810406671176691,-4.810406671176691) q[4];
u3(pi/2,1.472778636002895,-1.472778636002895) q[4];
rzz(-pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,2.2091679540043425,-2.2091679540043425) q[5];
u3(pi/2,5.252114598271416,-5.252114598271416) q[5];
rzz(pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,5.031574793989412,-5.031574793989412) q[6];
u3(pi/2,1.8409732950036186,-1.8409732950036186) q[6];
rzz(pi/2) q[1],q[6];
rzz(-pi/2) q[1],q[7];
u3(pi/2,4.454778382790327,-4.454778382790327) q[7];
u3(pi/2,1.288681306502533,-1.288681306502533) q[7];
rzz(pi/2) q[1],q[7];
rzz(pi/2) q[1],q[8];
u3(pi/2,4.399486352087147,-4.399486352087147) q[8];
u3(pi/2,1.245327327882994,-1.245327327882994) q[8];
rzz(pi/2) q[1],q[8];
rzz(-pi/2) q[1],q[9];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[10];
rzz(pi/2) q[1],q[10];
rzz(pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[11];
rzz(pi/2) q[1],q[12];
rzz(pi/2) q[1],q[12];
rzz(pi/2) q[1],q[13];
rzz(-pi/2) q[1],q[13];
rzz(pi/2) q[1],q[14];
rzz(pi/2) q[1],q[14];
rzz(pi/2) q[1],q[15];
rzz(pi/2) q[1],q[15];
rzz(pi/2) q[1],q[16];
rzz(pi/2) q[1],q[16];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,5.301751762198135,-5.301751762198135) q[3];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[3];
rzz(-pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,4.614371289592689,-4.614371289592689) q[4];
u3(pi/2,1.0800795543041708,-1.0800795543041708) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,5.252114598271416,-5.252114598271416) q[5];
u3(pi/2,1.91448656309762,-1.91448656309762) q[5];
rzz(pi/2) q[2],q[5];
rzz(pi/2) q[2],q[6];
u3(pi/2,1.8409732950036186,-1.8409732950036186) q[6];
u3(pi/2,4.883919939270692,-4.883919939270692) q[6];
rzz(pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,1.288681306502533,-1.288681306502533) q[7];
u3(pi/2,4.381265114696325,-4.381265114696325) q[7];
rzz(pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[8];
u3(pi/2,4.386919981472787,-4.386919981472787) q[8];
u3(pi/2,1.2208229051849937,-1.2208229051849937) q[8];
rzz(pi/2) q[2],q[8];
rzz(pi/2) q[2],q[9];
u3(pi/2,5.237034953534185,-5.237034953534185) q[9];
u3(pi/2,2.082875929330033,-2.082875929330033) q[9];
rzz(pi/2) q[2],q[9];
rzz(-pi/2) q[2],q[10];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[12];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[14];
rzz(-pi/2) q[2],q[14];
rzz(pi/2) q[2],q[15];
rzz(pi/2) q[2],q[15];
rzz(pi/2) q[2],q[16];
rzz(pi/2) q[2],q[16];
u3(pi/2,2.94555727200579,-2.94555727200579) q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(-pi/2) q[3],q[4];
u3(pi/2,4.221672207893964,-4.221672207893964) q[4];
u3(pi/2,0.2946813909067226,-0.2946813909067226) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,1.91448656309762,-1.91448656309762) q[5];
u3(pi/2,4.663380134988689,-4.663380134988689) q[5];
rzz(-pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,4.883919939270692,-4.883919939270692) q[6];
u3(pi/2,1.5462919040968963,-1.5462919040968963) q[6];
rzz(pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,4.381265114696325,-4.381265114696325) q[7];
u3(pi/2,1.141026451783813,-1.141026451783813) q[7];
rzz(pi/2) q[3],q[7];
rzz(-pi/2) q[3],q[8];
u3(pi/2,4.362415558774787,-4.362415558774787) q[8];
u3(pi/2,1.1718140597889928,-1.1718140597889928) q[8];
rzz(pi/2) q[3],q[8];
rzz(pi/2) q[3],q[9];
u3(pi/2,2.082875929330033,-2.082875929330033) q[9];
u3(pi/2,5.199964160221826,-5.199964160221826) q[9];
rzz(pi/2) q[3],q[9];
rzz(-pi/2) q[3],q[10];
u3(pi/2,1.093902561979966,-1.093902561979966) q[10];
u3(pi/2,4.2229288449554,-4.2229288449554) q[10];
rzz(pi/2) q[3],q[10];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[12];
rzz(pi/2) q[3],q[13];
rzz(pi/2) q[3],q[13];
rzz(pi/2) q[3],q[14];
rzz(-pi/2) q[3],q[14];
rzz(pi/2) q[3],q[15];
rzz(pi/2) q[3],q[15];
rzz(pi/2) q[3],q[16];
rzz(pi/2) q[3],q[16];
u3(pi/2,5.007070371291412,-5.007070371291412) q[4];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[5];
u3(pi/2,3.8779819715912405,-3.8779819715912405) q[5];
rzz(pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,1.5462919040968963,-1.5462919040968963) q[6];
u3(pi/2,4.295185475987965,-4.295185475987965) q[6];
rzz(pi/2) q[4],q[6];
rzz(-pi/2) q[4],q[7];
u3(pi/2,4.282619105373606,-4.282619105373606) q[7];
u3(pi/2,0.9449910701998098,-0.9449910701998098) q[7];
rzz(pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,1.1718140597889928,-1.1718140597889928) q[8];
u3(pi/2,4.215389022586785,-4.215389022586785) q[8];
rzz(pi/2) q[4],q[8];
rzz(pi/2) q[4],q[9];
u3(pi/2,5.199964160221826,-5.199964160221826) q[9];
u3(pi/2,2.0093626612360316,-2.0093626612360316) q[9];
rzz(pi/2) q[4],q[9];
rzz(pi/2) q[4],q[10];
u3(pi/2,4.2229288449554,-4.2229288449554) q[10];
u3(pi/2,1.0568317686676063,-1.0568317686676063) q[10];
rzz(pi/2) q[4],q[10];
rzz(-pi/2) q[4],q[11];
u3(pi/2,5.231380086757723,-5.231380086757723) q[11];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[11];
rzz(pi/2) q[4],q[11];
rzz(pi/2) q[4],q[12];
rzz(pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[13];
rzz(pi/2) q[4],q[13];
rzz(pi/2) q[4],q[14];
rzz(pi/2) q[4],q[14];
rzz(pi/2) q[4],q[15];
rzz(-pi/2) q[4],q[15];
rzz(pi/2) q[4],q[16];
rzz(pi/2) q[4],q[16];
u3(pi/2,5.448778298386137,-5.448778298386137) q[5];
u3(pi/2,0.4907167724907257,-0.4907167724907257) q[5];
rzz(pi/2) q[5],q[6];
u3(pi/2,4.295185475987965,-4.295185475987965) q[6];
u3(pi/2,0.36819465900072373,-0.36819465900072373) q[6];
rzz(-pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,0.9449910701998098,-0.9449910701998098) q[7];
u3(pi/2,3.6938846420908784,-3.6938846420908784) q[7];
rzz(pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,4.215389022586785,-4.215389022586785) q[8];
u3(pi/2,0.8771326688822703,-0.8771326688822703) q[8];
rzz(pi/2) q[5],q[8];
rzz(-pi/2) q[5],q[9];
u3(pi/2,5.150955314825825,-5.150955314825825) q[9];
u3(pi/2,1.9113449704440304,-1.9113449704440304) q[9];
rzz(pi/2) q[5],q[9];
rzz(pi/2) q[5],q[10];
u3(pi/2,1.0568317686676063,-1.0568317686676063) q[10];
u3(pi/2,4.149415576861399,-4.149415576861399) q[10];
rzz(pi/2) q[5],q[10];
rzz(-pi/2) q[5],q[11];
u3(pi/2,5.218813716143364,-5.218813716143364) q[11];
u3(pi/2,2.052716639855571,-2.052716639855571) q[11];
rzz(pi/2) q[5],q[11];
rzz(pi/2) q[5],q[12];
u3(pi/2,4.979424355939822,-4.979424355939822) q[12];
u3(pi/2,1.8258936502663878,-1.8258936502663878) q[12];
rzz(pi/2) q[5],q[12];
rzz(-pi/2) q[5],q[13];
rzz(pi/2) q[5],q[13];
rzz(pi/2) q[5],q[14];
rzz(pi/2) q[5],q[14];
rzz(-pi/2) q[5],q[15];
rzz(pi/2) q[5],q[15];
rzz(pi/2) q[5],q[16];
rzz(pi/2) q[5],q[16];
u3(pi/2,5.080583639385413,-5.080583639385413) q[6];
u3(pi/2,1.8164688723056186,-1.8164688723056186) q[6];
rzz(pi/2) q[6],q[7];
u3(pi/2,3.6938846420908784,-3.6938846420908784) q[7];
u3(pi/2,6.050079132283224,-6.050079132283224) q[7];
rzz(-pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,0.8771326688822703,-0.8771326688822703) q[8];
u3(pi/2,3.626026240773339,-3.626026240773339) q[8];
rzz(pi/2) q[6],q[8];
rzz(pi/2) q[6],q[9];
u3(pi/2,1.9113449704440304,-1.9113449704440304) q[9];
u3(pi/2,4.856273923919102,-4.856273923919102) q[9];
rzz(-pi/2) q[6],q[9];
rzz(pi/2) q[6],q[10];
u3(pi/2,4.149415576861399,-4.149415576861399) q[10];
u3(pi/2,0.9098052324796042,-0.9098052324796042) q[10];
rzz(pi/2) q[6],q[10];
rzz(pi/2) q[6],q[11];
u3(pi/2,2.052716639855571,-2.052716639855571) q[11];
u3(pi/2,5.145300448049363,-5.145300448049363) q[11];
rzz(pi/2) q[6],q[11];
rzz(pi/2) q[6],q[12];
u3(pi/2,1.8258936502663878,-1.8258936502663878) q[12];
u3(pi/2,4.94298188115818,-4.94298188115818) q[12];
rzz(pi/2) q[6],q[12];
rzz(pi/2) q[6],q[13];
u3(pi/2,2.0891591146372126,-2.0891591146372126) q[13];
u3(pi/2,5.218813716143364,-5.218813716143364) q[13];
rzz(pi/2) q[6],q[13];
rzz(-pi/2) q[6],q[14];
rzz(pi/2) q[6],q[14];
rzz(pi/2) q[6],q[15];
rzz(pi/2) q[6],q[15];
rzz(-pi/2) q[6],q[16];
rzz(pi/2) q[6],q[16];
u3(pi/2,4.479282805488327,-4.479282805488327) q[7];
u3(pi/2,5.203105752875415,-5.203105752875415) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,3.626026240773339,-3.626026240773339) q[8];
u3(pi/2,5.982220730965683,-5.982220730965683) q[8];
rzz(pi/2) q[7],q[8];
rzz(-pi/2) q[7],q[9];
u3(pi/2,4.856273923919102,-4.856273923919102) q[9];
u3(pi/2,1.321982188630585,-1.321982188630585) q[9];
rzz(pi/2) q[7],q[9];
rzz(pi/2) q[7],q[10];
u3(pi/2,0.9098052324796042,-0.9098052324796042) q[10];
u3(pi/2,3.8547341859546767,-3.8547341859546767) q[10];
rzz(pi/2) q[7],q[10];
rzz(-pi/2) q[7],q[11];
u3(pi/2,2.00370779445957,-2.00370779445957) q[11];
u3(pi/2,5.047282757257362,-5.047282757257362) q[11];
rzz(pi/2) q[7],q[11];
rzz(pi/2) q[7],q[12];
u3(pi/2,4.94298188115818,-4.94298188115818) q[12];
u3(pi/2,1.7517520636416686,-1.7517520636416686) q[12];
rzz(pi/2) q[7],q[12];
rzz(pi/2) q[7],q[13];
u3(pi/2,5.218813716143364,-5.218813716143364) q[13];
u3(pi/2,2.052716639855571,-2.052716639855571) q[13];
rzz(-pi/2) q[7],q[13];
rzz(pi/2) q[7],q[14];
u3(pi/2,2.808583832309275,-2.808583832309275) q[14];
u3(pi/2,5.938238433815427,-5.938238433815427) q[14];
rzz(pi/2) q[7],q[14];
rzz(pi/2) q[7],q[15];
rzz(-pi/2) q[7],q[15];
rzz(pi/2) q[7],q[16];
rzz(pi/2) q[7],q[16];
u3(pi/2,1.2698317505809944,-1.2698317505809944) q[8];
u3(pi/2,3.9885660329976016,-3.9885660329976016) q[8];
rzz(pi/2) q[8],q[9];
u3(pi/2,1.321982188630585,-1.321982188630585) q[9];
u3(pi/2,3.67817667882293,-3.67817667882293) q[9];
rzz(pi/2) q[8],q[9];
rzz(-pi/2) q[8],q[10];
u3(pi/2,0.7131415323648831,-0.7131415323648831) q[10];
u3(pi/2,3.4620351042559525,-3.4620351042559525) q[10];
rzz(pi/2) q[8],q[10];
rzz(pi/2) q[8],q[11];
u3(pi/2,5.047282757257362,-5.047282757257362) q[11];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[11];
rzz(pi/2) q[8],q[11];
rzz(-pi/2) q[8],q[12];
u3(pi/2,4.893344717231462,-4.893344717231462) q[12];
u3(pi/2,1.653734372849667,-1.653734372849667) q[12];
rzz(pi/2) q[8],q[12];
rzz(pi/2) q[8],q[13];
u3(pi/2,5.194309293445364,-5.194309293445364) q[13];
u3(pi/2,2.00370779445957,-2.00370779445957) q[13];
rzz(pi/2) q[8],q[13];
rzz(pi/2) q[8],q[14];
u3(pi/2,5.938238433815427,-5.938238433815427) q[14];
u3(pi/2,2.7715130389969156,-2.7715130389969156) q[14];
rzz(-pi/2) q[8],q[14];
rzz(pi/2) q[8],q[15];
u3(pi/2,2.0860175219836226,-2.0860175219836226) q[15];
u3(pi/2,5.2156721234897745,-5.2156721234897745) q[15];
rzz(pi/2) q[8],q[15];
rzz(-pi/2) q[8],q[16];
rzz(pi/2) q[8],q[16];
u3(pi/2,2.107380352028033,-2.107380352028033) q[9];
u3(pi/2,4.67531818707233,-4.67531818707233) q[9];
rzz(pi/2) q[9],q[10];
u3(pi/2,3.4620351042559525,-3.4620351042559525) q[10];
u3(pi/2,5.818229594448297,-5.818229594448297) q[10];
rzz(pi/2) q[9],q[10];
rzz(pi/2) q[9],q[11];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[11];
u3(pi/2,4.457919975443916,-4.457919975443916) q[11];
rzz(-pi/2) q[9],q[11];
rzz(pi/2) q[9],q[12];
u3(pi/2,1.653734372849667,-1.653734372849667) q[12];
u3(pi/2,4.599291644855457,-4.599291644855457) q[12];
rzz(pi/2) q[9],q[12];
rzz(pi/2) q[9],q[13];
u3(pi/2,2.00370779445957,-2.00370779445957) q[13];
u3(pi/2,5.047282757257362,-5.047282757257362) q[13];
rzz(-pi/2) q[9],q[13];
rzz(pi/2) q[9],q[14];
u3(pi/2,5.913105692586709,-5.913105692586709) q[14];
u3(pi/2,2.7225041936009147,-2.7225041936009147) q[14];
rzz(pi/2) q[9],q[14];
rzz(pi/2) q[9],q[15];
u3(pi/2,5.2156721234897745,-5.2156721234897745) q[15];
u3(pi/2,2.049575047201981,-2.049575047201981) q[15];
rzz(pi/2) q[9],q[15];
rzz(pi/2) q[9],q[16];
u3(pi/2,5.225725219981262,-5.225725219981262) q[16];
u3(pi/2,2.0715661957771094,-2.0715661957771094) q[16];
rzz(pi/2) q[9],q[16];
u3(pi/2,4.247433267653401,-4.247433267653401) q[10];
u3(pi/2,0.03392920065876977,-0.03392920065876977) q[10];
rzz(pi/2) q[10],q[11];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[11];
u3(pi/2,3.6725218120464684,-3.6725218120464684) q[11];
rzz(pi/2) q[10],q[11];
rzz(-pi/2) q[10],q[12];
u3(pi/2,1.4576989912656642,-1.4576989912656642) q[12];
u3(pi/2,4.2065925631567325,-4.2065925631567325) q[12];
rzz(pi/2) q[10],q[12];
rzz(pi/2) q[10],q[13];
u3(pi/2,1.9056901036675686,-1.9056901036675686) q[13];
u3(pi/2,4.850619057142641,-4.850619057142641) q[13];
rzz(pi/2) q[10],q[13];
rzz(-pi/2) q[10],q[14];
u3(pi/2,5.864096847190708,-5.864096847190708) q[14];
u3(pi/2,2.6244865028089133,-2.6244865028089133) q[14];
rzz(pi/2) q[10],q[14];
rzz(pi/2) q[10],q[15];
u3(pi/2,2.049575047201981,-2.049575047201981) q[15];
u3(pi/2,5.142158855395773,-5.142158855395773) q[15];
rzz(pi/2) q[10],q[15];
rzz(pi/2) q[10],q[16];
u3(pi/2,2.0715661957771094,-2.0715661957771094) q[16];
u3(pi/2,5.1886544266689025,-5.1886544266689025) q[16];
rzz(-pi/2) q[10],q[16];
u3(pi/2,5.2433181388413646,-5.2433181388413646) q[11];
u3(pi/2,0.28211502029236346,-0.28211502029236346) q[11];
rzz(pi/2) q[11],q[12];
u3(pi/2,4.2065925631567325,-4.2065925631567325) q[12];
u3(pi/2,0.27960174616949157,-0.27960174616949157) q[12];
rzz(pi/2) q[11],q[12];
rzz(-pi/2) q[11],q[13];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[13];
u3(pi/2,4.457919975443916,-4.457919975443916) q[13];
rzz(pi/2) q[11],q[13];
rzz(pi/2) q[11],q[14];
u3(pi/2,2.6244865028089133,-2.6244865028089133) q[14];
u3(pi/2,5.570043774814703,-5.570043774814703) q[14];
rzz(pi/2) q[11],q[14];
rzz(pi/2) q[11],q[15];
u3(pi/2,5.142158855395773,-5.142158855395773) q[15];
u3(pi/2,1.9025485110139788,-1.9025485110139788) q[15];
rzz(-pi/2) q[11],q[15];
rzz(pi/2) q[11],q[16];
u3(pi/2,2.047061773079109,-2.047061773079109) q[16];
u3(pi/2,5.139645581272902,-5.139645581272902) q[16];
rzz(pi/2) q[11],q[16];
u3(pi/2,4.991990726554181,-4.991990726554181) q[12];
u3(pi/2,5.116397795636337,-5.116397795636337) q[12];
rzz(pi/2) q[12],q[13];
u3(pi/2,4.457919975443916,-4.457919975443916) q[13];
u3(pi/2,0.5309291584566751,-0.5309291584566751) q[13];
rzz(-pi/2) q[12],q[13];
rzz(pi/2) q[12],q[14];
u3(pi/2,5.570043774814703,-5.570043774814703) q[14];
u3(pi/2,2.035752039526186,-2.035752039526186) q[14];
rzz(pi/2) q[12],q[14];
rzz(pi/2) q[12],q[15];
u3(pi/2,5.044141164603771,-5.044141164603771) q[15];
u3(pi/2,1.7058848108992577,-1.7058848108992577) q[15];
rzz(pi/2) q[12],q[15];
rzz(-pi/2) q[12],q[16];
u3(pi/2,1.9980529276831085,-1.9980529276831085) q[16];
u3(pi/2,5.0416278904809,-5.0416278904809) q[16];
rzz(pi/2) q[12],q[16];
u3(pi/2,5.2433181388413646,-5.2433181388413646) q[13];
u3(pi/2,1.2541237873130453,-1.2541237873130453) q[13];
rzz(pi/2) q[13],q[14];
u3(pi/2,2.035752039526186,-2.035752039526186) q[14];
u3(pi/2,4.39194652971853,-4.39194652971853) q[14];
rzz(pi/2) q[13],q[14];
rzz(-pi/2) q[13],q[15];
u3(pi/2,4.84747746448905,-4.84747746448905) q[15];
u3(pi/2,1.3131857292005336,-1.3131857292005336) q[15];
rzz(pi/2) q[13],q[15];
rzz(pi/2) q[13],q[16];
u3(pi/2,5.0416278904809,-5.0416278904809) q[16];
u3(pi/2,1.703371536776386,-1.703371536776386) q[16];
rzz(pi/2) q[13],q[16];
u3(pi/2,5.962742856513427,-5.962742856513427) q[14];
u3(pi/2,0.04084070449666731,-0.04084070449666731) q[14];
rzz(pi/2) q[14],q[15];
u3(pi/2,1.3131857292005336,-1.3131857292005336) q[15];
u3(pi/2,3.669380219392878,-3.669380219392878) q[15];
rzz(-pi/2) q[14],q[15];
rzz(pi/2) q[14],q[16];
u3(pi/2,1.703371536776386,-1.703371536776386) q[16];
u3(pi/2,4.452265108667455,-4.452265108667455) q[16];
rzz(pi/2) q[14],q[16];
u3(pi/2,5.240176546187775,-5.240176546187775) q[15];
u3(pi/2,5.844618972738451,-5.844618972738451) q[15];
rzz(pi/2) q[15],q[16];
u3(pi/2,4.452265108667455,-4.452265108667455) q[16];
u3(pi/2,0.5252742916802133,-0.5252742916802133) q[16];
rzz(-pi/2) q[15],q[16];
u3(pi,3*pi/2,-3*pi/2) q[1];
u3(pi,3*pi/2,-3*pi/2) q[2];
u3(pi,7*pi/8,-7*pi/8) q[3];
u3(pi,2.94555727200579,-2.94555727200579) q[4];
u3(pi,1.7668317083788996,-1.7668317083788996) q[5];
u3(pi,1.6688140175868982,-1.6688140175868982) q[6];
u3(pi,2.7733979945890694,-2.7733979945890694) q[7];
u3(pi,2.7733979945890694,-2.7733979945890694) q[8];
u3(pi,0.1413716694115407,-0.1413716694115407) q[9];
u3(pi,1.6933184402848986,-1.6933184402848986) q[10];
u3(pi,5.8358225133084,-5.8358225133084) q[11];
u3(pi,3.178663446902153,-3.178663446902153) q[12];
u3(pi,5.008955326883567,-5.008955326883567) q[13];
u3(pi,4.2719376903514,-4.2719376903514) q[14];
u3(pi,0.9035220471724246,-0.9035220471724246) q[15];
u3(pi/2,2.09607061847511,-2.09607061847511) q[16];
u3(pi/2,3.1836899951478967,-3.1836899951478967) q[16];
u3(pi,2.4096015653033716,-2.4096015653033716) q[17];
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
measure q[14] -> c[14];
measure q[15] -> c[15];
measure q[16] -> c[16];
measure q[17] -> c[17];
