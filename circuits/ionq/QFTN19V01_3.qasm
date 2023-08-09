OPENQASM 2.0;
include "qelib1.inc";
qreg q[19];
creg c[19];
u3(pi/2,2.070937877246392,-2.070937877246392) q[18];
u3(pi/2,4.427132367438737,-4.427132367438737) q[18];
rzz(pi/2) q[17],q[18];
u3(pi/2,2.85633604064384,-2.85633604064384) q[18];
rzz(pi/8) q[16],q[18];
u3(pi/2,5.235778316472749,-5.235778316472749) q[17];
rzz(0.7853830994606742) q[16],q[17];
u3(pi,0.14199998794225863,-0.14199998794225863) q[18];
rzz(0.19640131429629326) q[15],q[18];
u3(pi,3.9037430313506767,-3.9037430313506767) q[17];
rzz(pi/8) q[15],q[17];
u3(pi/2,2.0005662018059804,-2.0005662018059804) q[16];
rzz(0.7853830994606742) q[15],q[16];
u3(pi,3.232070522013179,-3.232070522013179) q[18];
rzz(pi/32) q[14],q[18];
u3(pi,5.086238506161875,-5.086238506161875) q[17];
rzz(0.19640131429629326) q[14],q[17];
u3(pi,5.042256209011618,-5.042256209011618) q[16];
rzz(pi/8) q[14],q[16];
u3(pi/2,2.094185662882956,-2.094185662882956) q[15];
rzz(0.7853830994606742) q[14],q[15];
u3(pi,0.6716725093374978,-0.6716725093374978) q[18];
rzz(0.049100328574073315) q[13],q[18];
u3(pi,0.3612831551628262,-0.3612831551628262) q[17];
rzz(pi/32) q[13],q[17];
u3(pi,0.7087433026498573,-0.7087433026498573) q[16];
rzz(0.19640131429629326) q[13],q[16];
u3(pi,1.1353715850073511,-1.1353715850073511) q[15];
rzz(pi/8) q[13],q[15];
u3(pi/2,1.7197078185750527,-1.7197078185750527) q[14];
rzz(0.7853830994606742) q[13],q[14];
u3(pi,4.444096967768122,-4.444096967768122) q[18];
rzz(pi/128) q[12],q[18];
u3(pi,1.0813361913656068,-1.0813361913656068) q[17];
rzz(0.049100328574073315) q[12],q[17];
u3(pi,3.1849466322093325,-3.1849466322093325) q[16];
rzz(pi/32) q[12],q[16];
u3(pi,2.358707764315217,-2.358707764315217) q[15];
rzz(0.19640131429629326) q[12],q[15];
u3(pi,4.4321589156844805,-4.4321589156844805) q[14];
rzz(pi/8) q[12],q[14];
u3(pi/2,2.0929290258215203,-2.0929290258215203) q[13];
rzz(0.7853830994606742) q[12],q[13];
u3(pi,0.8997521359881168,-0.8997521359881168) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi,1.3961237752553042,-1.3961237752553042) q[17];
rzz(pi/128) q[11],q[17];
u3(pi,1.45581403567351,-1.45581403567351) q[16];
rzz(0.049100328574073315) q[11],q[16];
u3(pi,3.2283006108288714,-3.2283006108288714) q[15];
rzz(pi/32) q[11],q[15];
u3(pi,0.24567254551072185,-0.24567254551072185) q[14];
rzz(0.19640131429629326) q[11],q[14];
u3(pi,0.6848671984825749,-0.6848671984825749) q[13];
rzz(pi/8) q[11],q[13];
u3(pi/2,0.5969026041820606,-0.5969026041820606) q[12];
rzz(0.7853830994606742) q[11],q[12];
u3(pi/2,2.386982098197525,-2.386982098197525) q[10];
u3(pi/2,4.068362486398782,-4.068362486398782) q[18];
rzz(0) q[10],q[18];
u3(pi/2,5.528574751787318,-5.528574751787318) q[10];
u3(pi,1.8434865691264906,-1.8434865691264906) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi,5.9175039223017345,-5.9175039223017345) q[16];
rzz(pi/128) q[10],q[16];
u3(pi,3.8032120664358033,-3.8032120664358033) q[15];
rzz(0.049100328574073315) q[10],q[15];
u3(pi,3.353335998441745,-3.353335998441745) q[14];
rzz(pi/32) q[10],q[14];
u3(pi,1.1649025559510953,-1.1649025559510953) q[13];
rzz(0.19640131429629326) q[10],q[13];
u3(pi,3.1987696398851275,-3.1987696398851275) q[12];
rzz(pi/8) q[10],q[12];
u3(pi/2,2.0891591146372126,-2.0891591146372126) q[11];
rzz(0.7853830994606742) q[10],q[11];
u3(pi/2,2.0740794698999814,-2.0740794698999814) q[9];
rzz(0) q[9],q[18];
u3(pi,5.642300405847268,-5.642300405847268) q[9];
u3(pi/2,0.5114512840044183,-0.5114512840044183) q[17];
rzz(0) q[9],q[17];
u3(pi/2,2.927964353145687,-2.927964353145687) q[9];
rzz(0.012275082143518329) q[9],q[16];
rzz(pi/128) q[9],q[15];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/32) q[9],q[13];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/8) q[9],q[11];
u3(pi/2,5.528574751787318,-5.528574751787318) q[10];
rzz(0.7853830994606742) q[9],q[10];
u3(pi/2,2.0370086765876216,-2.0370086765876216) q[8];
rzz(0) q[8],q[18];
u3(pi,6.164433104873892,-6.164433104873892) q[8];
rzz(0) q[8],q[17];
u3(pi,4.993247363615617,-4.993247363615617) q[8];
u3(pi/2,3.1101767270538954,-3.1101767270538954) q[16];
rzz(0) q[8],q[16];
u3(pi/2,2.837486484722301,-2.837486484722301) q[8];
u3(pi,2.426566165632756,-2.426566165632756) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi,2.5270971305476295,-2.5270971305476295) q[14];
rzz(pi/128) q[8],q[14];
u3(pi,0.5095663284122645,-0.5095663284122645) q[13];
rzz(0.049100328574073315) q[8],q[13];
u3(pi,1.3464866113285852,-1.3464866113285852) q[12];
rzz(pi/32) q[8],q[12];
u3(pi,4.428389004500172,-4.428389004500172) q[11];
rzz(0.19640131429629326) q[8],q[11];
u3(pi,4.094751864688936,-4.094751864688936) q[10];
rzz(pi/8) q[8],q[10];
u3(pi/2,6.06955700673548,-6.06955700673548) q[9];
rzz(0.7853830994606742) q[8],q[9];
u3(pi/2,2.012504253889621,-2.012504253889621) q[7];
rzz(0) q[7],q[18];
u3(pi,4.351734143752581,-4.351734143752581) q[7];
rzz(0) q[7],q[17];
u3(pi,5.88797295135799,-5.88797295135799) q[7];
rzz(0) q[7],q[16];
u3(pi,1.141026451783813,-1.141026451783813) q[7];
u3(pi/2,5.444380068671112,-5.444380068671112) q[15];
rzz(0) q[7],q[15];
u3(pi/2,3.479628023116055,-3.479628023116055) q[7];
u3(pi,3.303070515984308,-3.303070515984308) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi,3.7095926053588277,-3.7095926053588277) q[13];
rzz(pi/128) q[7],q[13];
u3(pi,1.8126989611213105,-1.8126989611213105) q[12];
rzz(0.049100328574073315) q[7],q[12];
u3(pi,5.580725189836908,-5.580725189836908) q[11];
rzz(pi/32) q[7],q[11];
u3(pi,4.606203148693354,-4.606203148693354) q[10];
rzz(0.19640131429629326) q[7],q[10];
u3(pi,2.3115838745113697,-2.3115838745113697) q[9];
rzz(pi/8) q[7],q[9];
u3(pi/2,5.979079138312095,-5.979079138312095) q[8];
rzz(0.7853830994606742) q[7],q[8];
u3(pi/2,0.4907167724907257,-0.4907167724907257) q[6];
rzz(0) q[6],q[18];
u3(pi,3.659327122901391,-3.659327122901391) q[6];
rzz(0) q[6],q[17];
u3(pi,0.5711415444226243,-0.5711415444226243) q[6];
rzz(0) q[6],q[16];
u3(pi,3.766141273123444,-3.766141273123444) q[6];
rzz(0) q[6],q[15];
u3(pi,0.6779556946446773,-0.6779556946446773) q[6];
u3(pi/2,2.3461413937008575,-2.3461413937008575) q[14];
rzz(0) q[6],q[14];
u3(pi/2,3.8459377265246246,-3.8459377265246246) q[6];
u3(pi,0.45553093477052,-0.45553093477052) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi,2.3329467045557806,-2.3329467045557806) q[12];
rzz(pi/128) q[6],q[12];
u3(pi,1.1297167182308896,-1.1297167182308896) q[11];
rzz(0.049100328574073315) q[6],q[11];
u3(pi,5.595176516043422,-5.595176516043422) q[10];
rzz(pi/32) q[6],q[10];
u3(pi,4.408282811517198,-4.408282811517198) q[9];
rzz(0.19640131429629326) q[6],q[9];
u3(pi,3.419309444167131,-3.419309444167131) q[8];
rzz(pi/8) q[6],q[8];
u3(pi/2,0.33803536952626173,-0.33803536952626173) q[7];
rzz(0.7853830994606742) q[6],q[7];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[5];
rzz(0) q[5],q[18];
u3(pi,4.292672201865093,-4.292672201865093) q[5];
rzz(0) q[5],q[17];
u3(pi,6.2021322167169695,-6.2021322167169695) q[5];
rzz(0) q[5],q[16];
u3(pi,1.8284069243892596,-1.8284069243892596) q[5];
rzz(0) q[5],q[15];
u3(pi,3.737866939241136,-3.737866939241136) q[5];
rzz(0) q[5],q[14];
u3(pi,5.647326954093012,-5.647326954093012) q[5];
u3(pi/2,4.023751870717807,-4.023751870717807) q[13];
rzz(0) q[5],q[13];
u3(pi/2,1.8899821403996195,-1.8899821403996195) q[5];
rzz(0.012275082143518329) q[5],q[12];
rzz(pi/128) q[5],q[11];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/32) q[5],q[9];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/8) q[5],q[7];
u3(pi/2,0.7043450729348316,-0.7043450729348316) q[6];
rzz(0.7853830994606742) q[5],q[6];
u3(pi/2,5*pi/8,-5*pi/8) q[4];
rzz(0) q[4],q[18];
u3(pi,4.302096979825863,-4.302096979825863) q[4];
rzz(0) q[4],q[17];
u3(pi,5.838335787431272,-5.838335787431272) q[4];
rzz(0) q[4],q[16];
u3(pi,1.091389287857094,-1.091389287857094) q[4];
rzz(0) q[4],q[15];
u3(pi,2.627628095462503,-2.627628095462503) q[4];
rzz(0) q[4],q[14];
u3(pi,4.163866903067912,-4.163866903067912) q[4];
rzz(0) q[4],q[13];
u3(pi,5.700105710673321,-5.700105710673321) q[4];
u3(pi/2,1.0744246875277093,-1.0744246875277093) q[12];
rzz(0) q[4],q[12];
u3(pi/2,1.7561502933566946,-1.7561502933566946) q[4];
u3(pi,5.857813661883529,-5.857813661883529) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi,4.209734155810323,-4.209734155810323) q[10];
rzz(pi/128) q[4],q[10];
u3(pi,2.5453183679384503,-2.5453183679384503) q[9];
rzz(0.049100328574073315) q[4],q[9];
u3(pi,3.042318325736356,-3.042318325736356) q[8];
rzz(pi/32) q[4],q[8];
u3(pi,5.246459731494954,-5.246459731494954) q[7];
rzz(0.19640131429629326) q[4],q[7];
u3(pi,3.3068404271686163,-3.3068404271686163) q[6];
rzz(pi/8) q[4],q[6];
u3(pi/2,1.8899821403996195,-1.8899821403996195) q[5];
rzz(0.7853830994606742) q[4],q[5];
u3(pi/2,pi,-pi) q[3];
rzz(0) q[3],q[18];
u3(pi,5.74345968929286,-5.74345968929286) q[3];
rzz(0) q[3],q[17];
u3(pi,1.5230441184603318,-1.5230441184603318) q[3];
rzz(0) q[3],q[16];
u3(pi,3.5858138548073897,-3.5858138548073897) q[3];
rzz(0) q[3],q[15];
u3(pi,5.648583591154448,-5.648583591154448) q[3];
rzz(0) q[3],q[14];
u3(pi,1.4275397017912022,-1.4275397017912022) q[3];
rzz(0) q[3],q[13];
u3(pi,3.4903094381382602,-3.4903094381382602) q[3];
rzz(0) q[3],q[12];
u3(pi,5.553079174485318,-5.553079174485318) q[3];
u3(pi/2,4.4258757303773,-4.4258757303773) q[11];
rzz(0) q[3],q[11];
u3(pi/2,1.8717609030087987,-1.8717609030087987) q[3];
rzz(0.012275082143518329) q[3],q[10];
rzz(pi/128) q[3],q[9];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/32) q[3],q[7];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/8) q[3],q[5];
u3(pi/2,1.7561502933566946,-1.7561502933566946) q[4];
rzz(0.7853830994606742) q[3],q[4];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(0) q[2],q[18];
u3(pi,5.1390172627421835,-5.1390172627421835) q[2];
rzz(0) q[2],q[17];
u3(pi,2.851309492398096,-2.851309492398096) q[2];
rzz(0) q[2],q[16];
u3(pi,0.5636017220540089,-0.5636017220540089) q[2];
rzz(0) q[2],q[15];
u3(pi,4.559079258889508,-4.559079258889508) q[2];
rzz(0) q[2],q[14];
u3(pi,2.2713714885454204,-2.2713714885454204) q[2];
rzz(0) q[2],q[13];
u3(pi,6.266849025380919,-6.266849025380919) q[2];
rzz(0) q[2],q[12];
u3(pi,3.979141255036832,-3.979141255036832) q[2];
rzz(0) q[2],q[11];
u3(pi,1.6914334846927446,-1.6914334846927446) q[2];
u3(pi/2,0.6389999457401639,-0.6389999457401639) q[10];
rzz(0) q[2],q[10];
u3(pi/2,5.26028273917075,-5.26028273917075) q[2];
u3(pi,0.4951150022057514,-0.4951150022057514) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi,1.4671237692264334,-1.4671237692264334) q[8];
rzz(pi/128) q[2],q[8];
u3(pi,3.7485483542633413,-3.7485483542633413) q[7];
rzz(0.049100328574073315) q[2],q[7];
u3(pi,1.5136193404995624,-1.5136193404995624) q[6];
rzz(pi/32) q[2],q[6];
u3(pi,4.228583711731861,-4.228583711731861) q[5];
rzz(0.19640131429629326) q[2],q[5];
u3(pi,0.4090353634973911,-0.4090353634973911) q[4];
rzz(pi/8) q[2],q[4];
u3(pi/2,1.8717609030087987,-1.8717609030087987) q[3];
rzz(0.7853830994606742) q[2],q[3];
u3(pi/2,pi,-pi) q[1];
rzz(0) q[1],q[18];
u3(pi,5.74345968929286,-5.74345968929286) q[1];
rzz(0) q[1],q[17];
u3(pi,1.5230441184603318,-1.5230441184603318) q[1];
rzz(0) q[1],q[16];
u3(pi,3.5858138548073897,-3.5858138548073897) q[1];
rzz(0) q[1],q[15];
u3(pi,5.648583591154448,-5.648583591154448) q[1];
rzz(0) q[1],q[14];
u3(pi,1.4275397017912022,-1.4275397017912022) q[1];
rzz(0) q[1],q[13];
u3(pi,3.4903094381382602,-3.4903094381382602) q[1];
rzz(0) q[1],q[12];
u3(pi,5.553079174485318,-5.553079174485318) q[1];
rzz(0) q[1],q[11];
u3(pi,1.3326636036527904,-1.3326636036527904) q[1];
rzz(0) q[1],q[10];
u3(pi,3.3948050214691303,-3.3948050214691303) q[1];
u3(pi/2,3.020955495691945,-3.020955495691945) q[9];
rzz(0) q[1],q[9];
u3(pi/2,5.997300375702915,-5.997300375702915) q[1];
u3(pi,4.990734089492745,-4.990734089492745) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi,0.054663712172462395,-0.054663712172462395) q[7];
rzz(pi/128) q[1],q[7];
u3(pi,1.9371060302034666,-1.9371060302034666) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi,1.5123627034381264,-1.5123627034381264) q[5];
rzz(pi/32) q[1],q[5];
u3(pi,0.9506459369762713,-0.9506459369762713) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi,0.5585751738082653,-0.5585751738082653) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,5.26028273917075,-5.26028273917075) q[2];
rzz(0.7853830994606742) q[1],q[2];
rzz(pi/2) q[0],q[18];
rzz(pi/2) q[0],q[17];
rzz(-pi/2) q[0],q[16];
rzz(pi/2) q[0],q[15];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[11];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[9];
u3(pi/2,2.8469112626830704,-2.8469112626830704) q[8];
rzz(pi/2) q[0],q[8];
u3(pi/2,2.76711480928189,-2.76711480928189) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,0.4724955350999049,-0.4724955350999049) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,5.881689766050811,-5.881689766050811) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,5.98033577537353,-5.98033577537353) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,5.529203070318036,-5.529203070318036) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,2.1186900855809565,-2.1186900855809565) q[2];
rzz(-pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi,0,0) q[0];
u3(pi/2,2.855707722113122,-2.855707722113122) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,3.6894864123758535,-3.6894864123758535) q[2];
u3(pi/2,6.045680902568198,-6.045680902568198) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,3.958406743523139,-3.958406743523139) q[3];
u3(pi/2,0.4241150082346221,-0.4241150082346221) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,4.4095394485786334,-4.4095394485786334) q[4];
u3(pi/2,1.0712830948741194,-1.0712830948741194) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,4.310893439255914,-4.310893439255914) q[5];
u3(pi/2,1.0712830948741194,-1.0712830948741194) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,5.184884515484595,-5.184884515484595) q[6];
u3(pi/2,1.9942830164988008,-1.9942830164988008) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,4.337911136076786,-4.337911136076786) q[7];
u3(pi/2,1.1718140597889928,-1.1718140597889928) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,1.276114935888174,-1.276114935888174) q[8];
u3(pi/2,4.405769537394326,-4.405769537394326) q[8];
rzz(pi/2) q[0],q[8];
rzz(pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[11];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[15];
rzz(pi/2) q[0],q[16];
rzz(pi/2) q[0],q[17];
rzz(pi/2) q[0],q[18];
u3(pi/2,4.4265040489080185,-4.4265040489080185) q[1];
u3(pi/2,1.3332919221835082,-1.3332919221835082) q[2];
rzz(0.7853830994606742) q[1],q[2];
u3(pi/2,5.136503988619312,-5.136503988619312) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,2.642079421669016,-2.642079421669016) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi/2,2.642079421669016,-2.642079421669016) q[5];
rzz(pi/32) q[1],q[5];
u3(pi/2,3.5650793432936974,-3.5650793432936974) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi/2,5.8842030401736825,-5.8842030401736825) q[7];
rzz(pi/128) q[1],q[7];
u3(pi/2,5.976565864189222,-5.976565864189222) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi/2,1.2849113953182254,-1.2849113953182254) q[1];
rzz(0) q[1],q[9];
u3(pi,4.853132331265512,-4.853132331265512) q[1];
rzz(0) q[1],q[10];
u3(pi,2.565424560921425,-2.565424560921425) q[1];
rzz(0) q[1],q[11];
u3(pi,0.2777167905773377,-0.2777167905773377) q[1];
rzz(0) q[1],q[12];
u3(pi,4.273194327412837,-4.273194327412837) q[1];
rzz(0) q[1],q[13];
u3(pi,1.9854865570687492,-1.9854865570687492) q[1];
rzz(0) q[1],q[14];
u3(pi,5.980964093904248,-5.980964093904248) q[1];
rzz(0) q[1],q[15];
u3(pi,3.6932563235601608,-3.6932563235601608) q[1];
rzz(0) q[1],q[16];
u3(pi,1.4055485532160734,-1.4055485532160734) q[1];
rzz(0) q[1],q[17];
u3(pi,5.401026090051572,-5.401026090051572) q[1];
rzz(0) q[1],q[18];
u3(pi/2,3.6894864123758535,-3.6894864123758535) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,3.008389125077586,-3.008389125077586) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi/2,5.26028273917075,-5.26028273917075) q[2];
rzz(0) q[2],q[10];
u3(pi,2.5453183679384503,-2.5453183679384503) q[2];
rzz(0) q[2],q[11];
u3(pi,0.25761059759436306,-0.25761059759436306) q[2];
rzz(0) q[2],q[12];
u3(pi,4.253088134429862,-4.253088134429862) q[2];
rzz(0) q[2],q[13];
u3(pi,1.9653803640857748,-1.9653803640857748) q[2];
rzz(0) q[2],q[14];
u3(pi,5.960857900921273,-5.960857900921273) q[2];
rzz(0) q[2],q[15];
u3(pi,3.673150130577186,-3.673150130577186) q[2];
rzz(0) q[2],q[16];
u3(pi,1.3854423602330987,-1.3854423602330987) q[2];
rzz(0) q[2],q[17];
u3(pi,5.380919897068598,-5.380919897068598) q[2];
rzz(0) q[2],q[18];
u3(pi/2,4.743804906920587,-4.743804906920587) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,3.764884636062008,-3.764884636062008) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,3.173008580125691,-3.173008580125691) q[3];
rzz(0) q[3],q[11];
u3(pi,0.45804420889339187,-0.45804420889339187) q[3];
rzz(0) q[3],q[12];
u3(pi,4.4535217457288905,-4.4535217457288905) q[3];
rzz(0) q[3],q[13];
u3(pi,2.1658139753848036,-2.1658139753848036) q[3];
rzz(0) q[3],q[14];
u3(pi,6.161291512220302,-6.161291512220302) q[3];
rzz(0) q[3],q[15];
u3(pi,3.873583741876215,-3.873583741876215) q[3];
rzz(0) q[3],q[16];
u3(pi,1.5858759715321276,-1.5858759715321276) q[3];
rzz(0) q[3],q[17];
u3(pi,5.581353508367626,-5.581353508367626) q[3];
rzz(0) q[3],q[18];
u3(pi/2,2.4460440400850127,-2.4460440400850127) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,1.2704600691117123,-1.2704600691117123) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,4.01684036687991,-4.01684036687991) q[4];
rzz(0) q[4],q[12];
u3(pi,0.07225663103256524,-0.07225663103256524) q[4];
rzz(0) q[4],q[13];
u3(pi,1.6084954386379742,-1.6084954386379742) q[4];
rzz(0) q[4],q[14];
u3(pi,3.1447342462433827,-3.1447342462433827) q[4];
rzz(0) q[4],q[15];
u3(pi,4.680973053848792,-4.680973053848792) q[4];
rzz(0) q[4],q[16];
u3(pi,6.217211861454201,-6.217211861454201) q[4];
rzz(0) q[4],q[17];
u3(pi,1.4702653618800232,-1.4702653618800232) q[4];
rzz(0) q[4],q[18];
u3(pi/2,5.68502606593609,-5.68502606593609) q[5];
rzz(0.7853830994606742) q[5],q[6];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,1.0568317686676063,-1.0568317686676063) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,4.1142297391411935,-4.1142297391411935) q[5];
rzz(0) q[5],q[13];
u3(pi,1.3998936864396119,-1.3998936864396119) q[5];
rzz(0) q[5],q[14];
u3(pi,5.395371223275111,-5.395371223275111) q[5];
rzz(0) q[5],q[15];
u3(pi,3.1076634529310234,-3.1076634529310234) q[5];
rzz(0) q[5],q[16];
u3(pi,0.8199556825869361,-0.8199556825869361) q[5];
rzz(0) q[5],q[17];
u3(pi,4.8154332194224345,-4.8154332194224345) q[5];
rzz(0) q[5],q[18];
u3(pi/2,0.3738495257771854,-0.3738495257771854) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,4.006787270388423,-4.006787270388423) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,1.944645852572082,-1.944645852572082) q[6];
rzz(0) q[6],q[14];
u3(pi,5.684397747405371,-5.684397747405371) q[6];
rzz(0) q[6],q[15];
u3(pi,3.7384952577718535,-3.7384952577718535) q[6];
rzz(0) q[6],q[16];
u3(pi,1.7932210866690539,-1.7932210866690539) q[6];
rzz(0) q[6],q[17];
u3(pi,6.130503904215122,-6.130503904215122) q[6];
rzz(0) q[6],q[18];
u3(pi/2,5.859698617475682,-5.859698617475682) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,5.46951280989983,-5.46951280989983) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,4.288902290680785,-4.288902290680785) q[7];
rzz(0) q[7],q[15];
u3(pi,0.26138050877867075,-0.26138050877867075) q[7];
rzz(0) q[7],q[16];
u3(pi,1.6317432242745384,-1.6317432242745384) q[7];
rzz(0) q[7],q[17];
u3(pi,3.0014776212396885,-3.0014776212396885) q[7];
rzz(0) q[7],q[18];
u3(pi/2,2.82240683998507,-2.82240683998507) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,5.426787149811009,-5.426787149811009) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,4.393203166779967,-4.393203166779967) q[8];
rzz(0) q[8],q[16];
u3(pi,0.8224689567098078,-0.8224689567098078) q[8];
rzz(0) q[8],q[17];
u3(pi,3.1057784973388696,-3.1057784973388696) q[8];
rzz(0) q[8],q[18];
u3(pi/2,2.9964510729939446,-2.9964510729939446) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,3.091955489663074,-3.091955489663074) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,1.425654746199048,-1.425654746199048) q[9];
rzz(0) q[9],q[17];
u3(pi,4.993875682146335,-4.993875682146335) q[9];
rzz(0) q[9],q[18];
u3(pi/2,3.7560881766319567,-3.7560881766319567) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,3.6348227002033906,-3.6348227002033906) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi/2,5.329397777549724,-5.329397777549724) q[10];
rzz(0) q[10],q[18];
u3(pi/2,4.40451290033289,-4.40451290033289) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,4.049512930477243,-4.049512930477243) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,1.0499202648297088,-1.0499202648297088) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/128) q[12],q[18];
u3(pi/2,0.8582831129607315,-0.8582831129607315) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
u3(pi/2,5.463229624592651,-5.463229624592651) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
u3(pi/2,2.278911310914036,-2.278911310914036) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
u3(pi/2,3.0856723043558945,-3.0856723043558945) q[16];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/8) q[16],q[18];
u3(pi/2,3.628539514896211,-3.628539514896211) q[17];
rzz(0.7853830994606742) q[17],q[18];
u3(pi,pi,-pi) q[0];
u3(pi/2,2.686061718819273,-2.686061718819273) q[1];
u3(pi/2,2.6659555258362984,-2.6659555258362984) q[2];
u3(pi/2,2.8663891371353274,-2.8663891371353274) q[3];
u3(pi/2,3.809495251742983,-3.809495251742983) q[4];
u3(pi/2,2.1004688481901357,-2.1004688481901357) q[5];
u3(pi/2,3.5870704918688254,-3.5870704918688254) q[6];
u3(pi/2,5.25714114651716,-5.25714114651716) q[7];
u3(pi/2,5.818229594448297,-5.818229594448297) q[8];
u3(pi/2,2.279539629444754,-2.279539629444754) q[9];
u3(pi/2,2.187805123959932,-2.187805123959932) q[10];
u3(pi/2,0.9022654101109886,-0.9022654101109886) q[18];
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
measure q[18] -> c[18];
