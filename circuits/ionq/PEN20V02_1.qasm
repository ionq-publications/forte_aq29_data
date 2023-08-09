OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
creg c[20];
u3(pi/2,3.666238626739289,-3.666238626739289) q[18];
u3(pi/2,5.979707456842812,-5.979707456842812) q[19];
rzz(-pi/2) q[18],q[19];
u3(pi/2,4.408911130047915,-4.408911130047915) q[19];
u3(pi/2,4.971256215040489,-4.971256215040489) q[19];
rzz(pi/2) q[18],q[19];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[17];
rzz(pi/2) q[17],q[19];
u3(pi/2,1.8296635614506955,-1.8296635614506955) q[19];
u3(pi/2,3.8465660450553423,-3.8465660450553423) q[19];
rzz(pi/2) q[17],q[19];
u3(pi/2,5.022778334559361,-5.022778334559361) q[16];
rzz(pi/2) q[16],q[19];
u3(pi/2,3.8465660450553423,-3.8465660450553423) q[19];
u3(pi/2,4.739406677205562,-4.739406677205562) q[19];
rzz(-pi/2) q[16],q[19];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[15];
rzz(pi/2) q[15],q[19];
u3(pi/2,4.739406677205562,-4.739406677205562) q[19];
u3(pi/2,6.095946385025634,-6.095946385025634) q[19];
rzz(pi/2) q[15],q[19];
u3(pi/2,5.9508048044297865,-5.9508048044297865) q[14];
rzz(pi/2) q[14],q[19];
u3(pi/2,2.9543537314358415,-2.9543537314358415) q[19];
u3(pi/2,3.382866969385489,-3.382866969385489) q[19];
rzz(-pi/2) q[14],q[19];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[13];
rzz(pi/2) q[13],q[19];
u3(pi/2,3.382866969385489,-3.382866969385489) q[19];
u3(pi/2,5.667433147075987,-5.667433147075987) q[19];
rzz(pi/2) q[13],q[19];
u3(pi/2,3.3784687396704633,-3.3784687396704633) q[12];
rzz(pi/2) q[12],q[19];
u3(pi/2,5.667433147075987,-5.667433147075987) q[19];
u3(pi/2,0.8111592231568845,-0.8111592231568845) q[19];
rzz(pi/2) q[12],q[19];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[11];
rzz(-pi/2) q[11],q[19];
u3(pi/2,0.8111592231568845,-0.8111592231568845) q[19];
u3(pi/2,1.0983007916949918,-1.0983007916949918) q[19];
rzz(pi/2) q[11],q[19];
u3(pi/2,2.518928989648296,-2.518928989648296) q[10];
rzz(pi/2) q[10],q[19];
u3(pi/2,4.2398934452847845,-4.2398934452847845) q[19];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[19];
rzz(pi/2) q[10],q[19];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
rzz(-pi/2) q[9],q[19];
u3(pi/2,3.666238626739289,-3.666238626739289) q[19];
u3(pi/2,5.66052164323809,-5.66052164323809) q[19];
rzz(pi/2) q[9],q[19];
u3(pi/2,5.362698659677777,-5.362698659677777) q[8];
rzz(pi/2) q[8],q[19];
u3(pi/2,5.66052164323809,-5.66052164323809) q[19];
u3(pi/2,0.22368139693559327,-0.22368139693559327) q[19];
rzz(pi/2) q[8],q[19];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
rzz(pi/2) q[7],q[19];
u3(pi/2,3.365274050525386,-3.365274050525386) q[19];
u3(pi/2,4.813548263830281,-4.813548263830281) q[19];
rzz(-pi/2) q[7],q[19];
u3(pi/2,4.172663362497963,-4.172663362497963) q[6];
rzz(pi/2) q[6],q[19];
u3(pi/2,4.813548263830281,-4.813548263830281) q[19];
u3(pi/2,5.059220809341003,-5.059220809341003) q[19];
rzz(pi/2) q[6],q[19];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(pi/2) q[5],q[19];
u3(pi/2,1.9176281557512098,-1.9176281557512098) q[19];
u3(pi/2,4.567875718319559,-4.567875718319559) q[19];
rzz(-pi/2) q[5],q[19];
u3(pi/2,5.693822525366141,-5.693822525366141) q[4];
rzz(pi/2) q[4],q[19];
u3(pi/2,1.4262830647297662,-1.4262830647297662) q[19];
u3(pi/2,3.5864421733381078,-3.5864421733381078) q[19];
rzz(pi/2) q[4],q[19];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(pi/2) q[3],q[19];
u3(pi/2,3.5864421733381078,-3.5864421733381078) q[19];
u3(pi/2,4.76453941843428,-4.76453941843428) q[19];
rzz(-pi/2) q[3],q[19];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
rzz(pi/2) q[2],q[19];
u3(pi/2,4.76453941843428,-4.76453941843428) q[19];
u3(pi/2,5.5499375818317285,-5.5499375818317285) q[19];
rzz(pi/2) q[2],q[19];
u3(pi/2,pi/2,-pi/2) q[1];
u3(pi/2,3.979141255036832,-3.979141255036832) q[19];
rzz(pi/2) q[1],q[19];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,0,0) q[1];
u3(pi/2,3*pi/4,-3*pi/4) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,9*pi/8,-9*pi/8) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,15*pi/8,-15*pi/8) q[3];
u3(pi/2,2.552229871776348,-2.552229871776348) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,2.552229871776348,-2.552229871776348) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[4];
u3(pi/2,4.025008507779242,-4.025008507779242) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,5.399769452990137,-5.399769452990137) q[5];
u3(pi/2,2.2091679540043425,-2.2091679540043425) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,1.03107070890817,-1.03107070890817) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,5.74345968929286,-5.74345968929286) q[6];
u3(pi/2,2.577362613005066,-2.577362613005066) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,4.466716434873968,-4.466716434873968) q[7];
u3(pi/2,1.3131857292005336,-1.3131857292005336) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,2.2211060060879837,-2.2211060060879837) q[8];
rzz(-pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[8];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
rzz(-pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[9];
u3(pi/2,5.66052164323809,-5.66052164323809) q[10];
rzz(pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[10];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[11];
rzz(-pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[11];
u3(pi/2,0.23687608608067037,-0.23687608608067037) q[12];
rzz(pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[12];
u3(pi/2,3.669380219392878,-3.669380219392878) q[13];
rzz(-pi/2) q[0],q[13];
rzz(-pi/2) q[0],q[13];
u3(pi/2,5.9508048044297865,-5.9508048044297865) q[14];
rzz(pi/2) q[0],q[14];
rzz(-pi/2) q[0],q[14];
u3(pi/2,3.669380219392878,-3.669380219392878) q[15];
rzz(-pi/2) q[0],q[15];
rzz(-pi/2) q[0],q[15];
u3(pi/2,5.022778334559361,-5.022778334559361) q[16];
rzz(pi/2) q[0],q[16];
rzz(-pi/2) q[0],q[16];
u3(pi/2,3.669380219392878,-3.669380219392878) q[17];
rzz(-pi/2) q[0],q[17];
rzz(-pi/2) q[0],q[17];
u3(pi/2,3.666238626739289,-3.666238626739289) q[18];
rzz(-pi/2) q[0],q[18];
rzz(pi/2) q[0],q[18];
u3(pi/2,pi/4,-pi/4) q[1];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,9*pi/8,-9*pi/8) q[2];
u3(pi/2,15*pi/8,-15*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,2.552229871776348,-2.552229871776348) q[3];
u3(pi/2,5.301751762198135,-5.301751762198135) q[3];
rzz(-pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,4.025008507779242,-4.025008507779242) q[4];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
u3(pi/2,2.110521944681623,-2.110521944681623) q[5];
rzz(pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[6];
u3(pi/2,5.718955266594859,-5.718955266594859) q[6];
u3(pi/2,2.5277254490783476,-2.5277254490783476) q[6];
rzz(pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u3(pi/2,4.454778382790327,-4.454778382790327) q[7];
u3(pi/2,1.288681306502533,-1.288681306502533) q[7];
rzz(pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[8];
u3(pi/2,0.6440264939859075,-0.6440264939859075) q[8];
u3(pi/2,3.7736810954920594,-3.7736810954920594) q[8];
rzz(pi/2) q[1],q[8];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[10];
rzz(pi/2) q[1],q[10];
rzz(pi/2) q[1],q[11];
rzz(pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[12];
rzz(pi/2) q[1],q[12];
rzz(pi/2) q[1],q[13];
rzz(pi/2) q[1],q[13];
rzz(pi/2) q[1],q[14];
rzz(pi/2) q[1],q[14];
rzz(pi/2) q[1],q[15];
rzz(pi/2) q[1],q[15];
rzz(pi/2) q[1],q[16];
rzz(-pi/2) q[1],q[16];
rzz(pi/2) q[1],q[17];
rzz(pi/2) q[1],q[17];
rzz(pi/2) q[1],q[18];
rzz(-pi/2) q[1],q[18];
u3(pi/2,3*pi/8,-3*pi/8) q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[3];
u3(pi/2,4.516353598800687,-4.516353598800687) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[4];
u3(pi/2,3.436274044496516,-3.436274044496516) q[4];
rzz(pi/2) q[2],q[4];
rzz(-pi/2) q[2],q[5];
u3(pi/2,5.252114598271416,-5.252114598271416) q[5];
u3(pi/2,1.91448656309762,-1.91448656309762) q[5];
rzz(pi/2) q[2],q[5];
rzz(pi/2) q[2],q[6];
u3(pi/2,2.5277254490783476,-2.5277254490783476) q[6];
u3(pi/2,5.571300411876139,-5.571300411876139) q[6];
rzz(pi/2) q[2],q[6];
rzz(-pi/2) q[2],q[7];
u3(pi/2,4.430273960092326,-4.430273960092326) q[7];
u3(pi/2,1.2396724611065324,-1.2396724611065324) q[7];
rzz(pi/2) q[2],q[7];
rzz(pi/2) q[2],q[8];
u3(pi/2,3.7736810954920594,-3.7736810954920594) q[8];
u3(pi/2,0.607584019204266,-0.607584019204266) q[8];
rzz(pi/2) q[2],q[8];
rzz(pi/2) q[2],q[9];
u3(pi/2,5.237034953534185,-5.237034953534185) q[9];
u3(pi/2,2.082875929330033,-2.082875929330033) q[9];
rzz(-pi/2) q[2],q[9];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[10];
rzz(-pi/2) q[2],q[11];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[12];
rzz(pi/2) q[2],q[12];
rzz(pi/2) q[2],q[13];
rzz(-pi/2) q[2],q[13];
rzz(pi/2) q[2],q[14];
rzz(pi/2) q[2],q[14];
rzz(pi/2) q[2],q[15];
rzz(-pi/2) q[2],q[15];
rzz(pi/2) q[2],q[16];
rzz(pi/2) q[2],q[16];
rzz(pi/2) q[2],q[17];
rzz(pi/2) q[2],q[17];
rzz(-pi/2) q[2],q[18];
rzz(pi/2) q[2],q[18];
u3(pi/2,6.086521607064865,-6.086521607064865) q[3];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.436274044496516,-3.436274044496516) q[4];
u3(pi/2,5.792468534688861,-5.792468534688861) q[4];
rzz(pi/2) q[3],q[4];
rzz(-pi/2) q[3],q[5];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[5];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,5.571300411876139,-5.571300411876139) q[6];
u3(pi/2,2.233672376702343,-2.233672376702343) q[6];
rzz(pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,1.2396724611065324,-1.2396724611065324) q[7];
u3(pi/2,4.282619105373606,-4.282619105373606) q[7];
rzz(pi/2) q[3],q[7];
rzz(pi/2) q[3],q[8];
u3(pi/2,0.607584019204266,-0.607584019204266) q[8];
u3(pi/2,3.700167827398058,-3.700167827398058) q[8];
rzz(pi/2) q[3],q[8];
rzz(pi/2) q[3],q[9];
u3(pi/2,5.224468582919826,-5.224468582919826) q[9];
u3(pi/2,2.0583715066320325,-2.0583715066320325) q[9];
rzz(-pi/2) q[3],q[9];
rzz(pi/2) q[3],q[10];
u3(pi/2,4.079043901420987,-4.079043901420987) q[10];
u3(pi/2,0.9248848772168351,-0.9248848772168351) q[10];
rzz(pi/2) q[3],q[10];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[12];
rzz(pi/2) q[3],q[12];
rzz(pi/2) q[3],q[13];
rzz(-pi/2) q[3],q[13];
rzz(pi/2) q[3],q[14];
rzz(pi/2) q[3],q[14];
rzz(pi/2) q[3],q[15];
rzz(-pi/2) q[3],q[15];
rzz(pi/2) q[3],q[16];
rzz(pi/2) q[3],q[16];
rzz(pi/2) q[3],q[17];
rzz(pi/2) q[3],q[17];
rzz(-pi/2) q[3],q[18];
rzz(pi/2) q[3],q[18];
u3(pi/2,4.221672207893964,-4.221672207893964) q[4];
u3(pi/2,4.9084243619686925,-4.9084243619686925) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[5];
u3(pi/2,3.8779819715912405,-3.8779819715912405) q[5];
rzz(pi/2) q[4],q[5];
rzz(-pi/2) q[4],q[6];
u3(pi/2,5.3752650302921365,-5.3752650302921365) q[6];
u3(pi/2,1.8409732950036186,-1.8409732950036186) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,4.282619105373606,-4.282619105373606) q[7];
u3(pi/2,0.9449910701998098,-0.9449910701998098) q[7];
rzz(pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,3.700167827398058,-3.700167827398058) q[8];
u3(pi/2,0.4599291644855457,-0.4599291644855457) q[8];
rzz(-pi/2) q[4],q[8];
rzz(pi/2) q[4],q[9];
u3(pi/2,5.199964160221826,-5.199964160221826) q[9];
u3(pi/2,2.0093626612360316,-2.0093626612360316) q[9];
rzz(pi/2) q[4],q[9];
rzz(pi/2) q[4],q[10];
u3(pi/2,0.9248848772168351,-0.9248848772168351) q[10];
u3(pi/2,4.041973108108627,-4.041973108108627) q[10];
rzz(-pi/2) q[4],q[10];
rzz(pi/2) q[4],q[11];
u3(pi/2,2.0897874331679303,-2.0897874331679303) q[11];
u3(pi/2,5.218813716143364,-5.218813716143364) q[11];
rzz(pi/2) q[4],q[11];
rzz(pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[12];
rzz(pi/2) q[4],q[13];
rzz(pi/2) q[4],q[13];
rzz(pi/2) q[4],q[14];
rzz(-pi/2) q[4],q[14];
rzz(pi/2) q[4],q[15];
rzz(pi/2) q[4],q[15];
rzz(pi/2) q[4],q[16];
rzz(pi/2) q[4],q[16];
rzz(-pi/2) q[4],q[17];
rzz(pi/2) q[4],q[17];
rzz(pi/2) q[4],q[18];
rzz(pi/2) q[4],q[18];
u3(pi/2,2.3071856447963444,-2.3071856447963444) q[5];
u3(pi/2,3.436274044496516,-3.436274044496516) q[5];
rzz(-pi/2) q[5],q[6];
u3(pi/2,4.982565948593412,-4.982565948593412) q[6];
u3(pi/2,1.0555751316061706,-1.0555751316061706) q[6];
rzz(pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,0.9449910701998098,-0.9449910701998098) q[7];
u3(pi/2,3.6938846420908784,-3.6938846420908784) q[7];
rzz(pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,3.601521818075339,-3.601521818075339) q[8];
u3(pi/2,0.26389378290154264,-0.26389378290154264) q[8];
rzz(pi/2) q[5],q[8];
rzz(pi/2) q[5],q[9];
u3(pi/2,2.0093626612360316,-2.0093626612360316) q[9];
u3(pi/2,5.052937624033824,-5.052937624033824) q[9];
rzz(pi/2) q[5],q[9];
rzz(pi/2) q[5],q[10];
u3(pi/2,0.9003804545188347,-0.9003804545188347) q[10];
u3(pi/2,3.992964262712627,-3.992964262712627) q[10];
rzz(-pi/2) q[5],q[10];
rzz(pi/2) q[5],q[11];
u3(pi/2,5.218813716143364,-5.218813716143364) q[11];
u3(pi/2,2.052716639855571,-2.052716639855571) q[11];
rzz(pi/2) q[5],q[11];
rzz(pi/2) q[5],q[12];
u3(pi/2,4.941096925566026,-4.941096925566026) q[12];
u3(pi/2,1.7875662198925921,-1.7875662198925921) q[12];
rzz(-pi/2) q[5],q[12];
rzz(pi/2) q[5],q[13];
rzz(pi/2) q[5],q[13];
rzz(pi/2) q[5],q[14];
rzz(-pi/2) q[5],q[14];
rzz(pi/2) q[5],q[15];
rzz(pi/2) q[5],q[15];
rzz(pi/2) q[5],q[16];
rzz(-pi/2) q[5],q[16];
rzz(pi/2) q[5],q[17];
rzz(pi/2) q[5],q[17];
rzz(pi/2) q[5],q[18];
rzz(pi/2) q[5],q[18];
u3(pi/2,2.626371458401067,-2.626371458401067) q[6];
u3(pi/2,4.417707589477967,-4.417707589477967) q[6];
rzz(-pi/2) q[6],q[7];
u3(pi/2,0.5522919885010856,-0.5522919885010856) q[7];
u3(pi/2,2.90848647869343,-2.90848647869343) q[7];
rzz(pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,0.26389378290154264,-0.26389378290154264) q[8];
u3(pi/2,3.0127873547926116,-3.0127873547926116) q[8];
rzz(pi/2) q[6],q[8];
rzz(-pi/2) q[6],q[9];
u3(pi/2,1.9113449704440304,-1.9113449704440304) q[9];
u3(pi/2,4.856273923919102,-4.856273923919102) q[9];
rzz(pi/2) q[6],q[9];
rzz(pi/2) q[6],q[10];
u3(pi/2,0.8513716091228339,-0.8513716091228339) q[10];
u3(pi/2,3.8949465719206255,-3.8949465719206255) q[10];
rzz(pi/2) q[6],q[10];
rzz(pi/2) q[6],q[11];
u3(pi/2,2.052716639855571,-2.052716639855571) q[11];
u3(pi/2,5.145300448049363,-5.145300448049363) q[11];
rzz(-pi/2) q[6],q[11];
rzz(pi/2) q[6],q[12];
u3(pi/2,4.929158873482385,-4.929158873482385) q[12];
u3(pi/2,1.763061797194592,-1.763061797194592) q[12];
rzz(pi/2) q[6],q[12];
rzz(pi/2) q[6],q[13];
u3(pi/2,2.0891591146372126,-2.0891591146372126) q[13];
u3(pi/2,5.218813716143364,-5.218813716143364) q[13];
rzz(-pi/2) q[6],q[13];
rzz(pi/2) q[6],q[14];
rzz(pi/2) q[6],q[14];
rzz(pi/2) q[6],q[15];
rzz(pi/2) q[6],q[15];
rzz(-pi/2) q[6],q[16];
rzz(pi/2) q[6],q[16];
rzz(pi/2) q[6],q[17];
rzz(-pi/2) q[6],q[17];
rzz(pi/2) q[6],q[18];
rzz(pi/2) q[6],q[18];
u3(pi/2,4.479282805488327,-4.479282805488327) q[7];
u3(pi/2,4.589866866894688,-4.589866866894688) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,3.0127873547926116,-3.0127873547926116) q[8];
u3(pi/2,5.3689818449849565,-5.3689818449849565) q[8];
rzz(pi/2) q[7],q[8];
rzz(-pi/2) q[7],q[9];
u3(pi/2,1.7146812703293088,-1.7146812703293088) q[9];
u3(pi/2,4.463574842220378,-4.463574842220378) q[9];
rzz(pi/2) q[7],q[9];
rzz(pi/2) q[7],q[10];
u3(pi/2,3.8949465719206255,-3.8949465719206255) q[10];
u3(pi/2,0.5566902182161113,-0.5566902182161113) q[10];
rzz(pi/2) q[7],q[10];
rzz(-pi/2) q[7],q[11];
u3(pi/2,5.145300448049363,-5.145300448049363) q[11];
u3(pi/2,1.9056901036675686,-1.9056901036675686) q[11];
rzz(pi/2) q[7],q[11];
rzz(pi/2) q[7],q[12];
u3(pi/2,1.763061797194592,-1.763061797194592) q[12];
u3(pi/2,4.855017286857667,-4.855017286857667) q[12];
rzz(pi/2) q[7],q[12];
rzz(pi/2) q[7],q[13];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[13];
u3(pi/2,5.194309293445364,-5.194309293445364) q[13];
rzz(pi/2) q[7],q[13];
rzz(pi/2) q[7],q[14];
u3(pi/2,4.371212018204838,-4.371212018204838) q[14];
u3(pi/2,1.2176813125314039,-1.2176813125314039) q[14];
rzz(pi/2) q[7],q[14];
rzz(pi/2) q[7],q[15];
rzz(-pi/2) q[7],q[15];
rzz(pi/2) q[7],q[16];
rzz(pi/2) q[7],q[16];
rzz(pi/2) q[7],q[17];
rzz(-pi/2) q[7],q[17];
rzz(pi/2) q[7],q[18];
rzz(pi/2) q[7],q[18];
u3(pi/2,0.6565928646002668,-0.6565928646002668) q[8];
u3(pi/2,3.068079385495792,-3.068079385495792) q[8];
rzz(pi/2) q[8],q[9];
u3(pi/2,4.463574842220378,-4.463574842220378) q[9];
u3(pi/2,0.5365840252331366,-0.5365840252331366) q[9];
rzz(pi/2) q[8],q[9];
rzz(-pi/2) q[8],q[10];
u3(pi/2,3.6982828718059046,-3.6982828718059046) q[10];
u3(pi/2,0.1639911365173872,-0.1639911365173872) q[10];
rzz(pi/2) q[8],q[10];
rzz(pi/2) q[8],q[11];
u3(pi/2,1.9056901036675686,-1.9056901036675686) q[11];
u3(pi/2,4.850619057142641,-4.850619057142641) q[11];
rzz(-pi/2) q[8],q[11];
rzz(pi/2) q[8],q[12];
u3(pi/2,4.855017286857667,-4.855017286857667) q[12];
u3(pi/2,1.6154069424758717,-1.6154069424758717) q[12];
rzz(pi/2) q[8],q[12];
rzz(pi/2) q[8],q[13];
u3(pi/2,5.194309293445364,-5.194309293445364) q[13];
u3(pi/2,2.00370779445957,-2.00370779445957) q[13];
rzz(pi/2) q[8],q[13];
rzz(-pi/2) q[8],q[14];
u3(pi/2,4.359273966121196,-4.359273966121196) q[14];
u3(pi/2,1.1931768898334034,-1.1931768898334034) q[14];
rzz(pi/2) q[8],q[14];
rzz(pi/2) q[8],q[15];
u3(pi/2,5.230751768227005,-5.230751768227005) q[15];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[15];
rzz(pi/2) q[8],q[15];
rzz(-pi/2) q[8],q[16];
rzz(pi/2) q[8],q[16];
rzz(pi/2) q[8],q[17];
rzz(pi/2) q[8],q[17];
rzz(pi/2) q[8],q[18];
rzz(-pi/2) q[8],q[18];
u3(pi/2,5.2489730056178265,-5.2489730056178265) q[9];
u3(pi/2,1.6876635735084369,-1.6876635735084369) q[9];
rzz(pi/2) q[9],q[10];
u3(pi/2,0.1639911365173872,-0.1639911365173872) q[10];
u3(pi/2,2.5201856267097322,-2.5201856267097322) q[10];
rzz(pi/2) q[9],q[10];
rzz(pi/2) q[9],q[11];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[11];
u3(pi/2,4.457919975443916,-4.457919975443916) q[11];
rzz(-pi/2) q[9],q[11];
rzz(pi/2) q[9],q[12];
u3(pi/2,1.6154069424758717,-1.6154069424758717) q[12];
u3(pi/2,4.560964214481662,-4.560964214481662) q[12];
rzz(pi/2) q[9],q[12];
rzz(pi/2) q[9],q[13];
u3(pi/2,2.00370779445957,-2.00370779445957) q[13];
u3(pi/2,5.047282757257362,-5.047282757257362) q[13];
rzz(pi/2) q[9],q[13];
rzz(-pi/2) q[9],q[14];
u3(pi/2,4.334769543423196,-4.334769543423196) q[14];
u3(pi/2,1.1441680444374027,-1.1441680444374027) q[14];
rzz(pi/2) q[9],q[14];
rzz(pi/2) q[9],q[15];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[15];
u3(pi/2,5.194309293445364,-5.194309293445364) q[15];
rzz(pi/2) q[9],q[15];
rzz(-pi/2) q[9],q[16];
u3(pi/2,0.29719466502959446,-0.29719466502959446) q[16];
u3(pi/2,3.426849266535746,-3.426849266535746) q[16];
rzz(pi/2) q[9],q[16];
rzz(pi/2) q[9],q[17];
rzz(pi/2) q[9],q[17];
rzz(-pi/2) q[9],q[18];
rzz(pi/2) q[9],q[18];
u3(pi/2,0.9493892999148356,-0.9493892999148356) q[10];
u3(pi/2,3.0957254008473822,-3.0957254008473822) q[10];
rzz(pi/2) q[10],q[11];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[11];
u3(pi/2,3.6725218120464684,-3.6725218120464684) q[11];
rzz(pi/2) q[10],q[11];
rzz(-pi/2) q[10],q[12];
u3(pi/2,1.4193715608918684,-1.4193715608918684) q[12];
u3(pi/2,4.168265132782937,-4.168265132782937) q[12];
rzz(pi/2) q[10],q[12];
rzz(pi/2) q[10],q[13];
u3(pi/2,5.047282757257362,-5.047282757257362) q[13];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[13];
rzz(pi/2) q[10],q[13];
rzz(pi/2) q[10],q[14];
u3(pi/2,1.1441680444374027,-1.1441680444374027) q[14];
u3(pi/2,4.187114688704476,-4.187114688704476) q[14];
rzz(pi/2) q[10],q[14];
rzz(pi/2) q[10],q[15];
u3(pi/2,5.194309293445364,-5.194309293445364) q[15];
u3(pi/2,2.00370779445957,-2.00370779445957) q[15];
rzz(pi/2) q[10],q[15];
rzz(pi/2) q[10],q[16];
u3(pi/2,3.426849266535746,-3.426849266535746) q[16];
u3(pi/2,0.2607521902479528,-0.2607521902479528) q[16];
rzz(-pi/2) q[10],q[16];
rzz(pi/2) q[10],q[17];
u3(pi/2,5.230751768227005,-5.230751768227005) q[17];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[17];
rzz(pi/2) q[10],q[17];
rzz(pi/2) q[10],q[18];
rzz(-pi/2) q[10],q[18];
u3(pi/2,5.2433181388413646,-5.2433181388413646) q[11];
u3(pi/2,0.24378758991856794,-0.24378758991856794) q[11];
rzz(pi/2) q[11],q[12];
u3(pi/2,4.168265132782937,-4.168265132782937) q[12];
u3(pi/2,0.24127431579569608,-0.24127431579569608) q[12];
rzz(pi/2) q[11],q[12];
rzz(pi/2) q[11],q[13];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[13];
u3(pi/2,4.457919975443916,-4.457919975443916) q[13];
rzz(pi/2) q[11],q[13];
rzz(-pi/2) q[11],q[14];
u3(pi/2,1.045522035114683,-1.045522035114683) q[14];
u3(pi/2,3.991079307120473,-3.991079307120473) q[14];
rzz(pi/2) q[11],q[14];
rzz(pi/2) q[11],q[15];
u3(pi/2,2.00370779445957,-2.00370779445957) q[15];
u3(pi/2,5.046654438726644,-5.046654438726644) q[15];
rzz(pi/2) q[11],q[15];
rzz(-pi/2) q[11],q[16];
u3(pi/2,0.2607521902479528,-0.2607521902479528) q[16];
u3(pi/2,3.352707679911027,-3.352707679911027) q[16];
rzz(pi/2) q[11],q[16];
rzz(pi/2) q[11],q[17];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[17];
u3(pi/2,5.194309293445364,-5.194309293445364) q[17];
rzz(pi/2) q[11],q[17];
rzz(-pi/2) q[11],q[18];
u3(pi/2,5.228238494104134,-5.228238494104134) q[18];
u3(pi/2,2.0747077884306995,-2.0747077884306995) q[18];
rzz(pi/2) q[11],q[18];
u3(pi/2,1.8120706425905926,-1.8120706425905926) q[12];
u3(pi/2,4.809778352645973,-4.809778352645973) q[12];
rzz(pi/2) q[12],q[13];
u3(pi/2,4.457919975443916,-4.457919975443916) q[13];
u3(pi/2,0.5309291584566751,-0.5309291584566751) q[13];
rzz(pi/2) q[12],q[13];
rzz(-pi/2) q[12],q[14];
u3(pi/2,0.84948665353068,-0.84948665353068) q[14];
u3(pi/2,3.598380225421749,-3.598380225421749) q[14];
rzz(pi/2) q[12],q[14];
rzz(pi/2) q[12],q[15];
u3(pi/2,5.046654438726644,-5.046654438726644) q[15];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[15];
rzz(pi/2) q[12],q[15];
rzz(pi/2) q[12],q[16];
u3(pi/2,3.352707679911027,-3.352707679911027) q[16];
u3(pi/2,0.11309733552923254,-0.11309733552923254) q[16];
rzz(-pi/2) q[12],q[16];
rzz(pi/2) q[12],q[17];
u3(pi/2,5.194309293445364,-5.194309293445364) q[17];
u3(pi/2,2.00370779445957,-2.00370779445957) q[17];
rzz(pi/2) q[12],q[17];
rzz(pi/2) q[12],q[18];
u3(pi/2,2.0747077884306995,-2.0747077884306995) q[18];
u3(pi/2,5.191796019322492,-5.191796019322492) q[18];
rzz(-pi/2) q[12],q[18];
u3(pi/2,5.2433181388413646,-5.2433181388413646) q[13];
u3(pi/2,1.3879556343559707,-1.3879556343559707) q[13];
rzz(pi/2) q[13],q[14];
u3(pi/2,3.598380225421749,-3.598380225421749) q[14];
u3(pi/2,5.954574715614093,-5.954574715614093) q[14];
rzz(pi/2) q[13],q[14];
rzz(pi/2) q[13],q[15];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[15];
u3(pi/2,4.457919975443916,-4.457919975443916) q[15];
rzz(pi/2) q[13],q[15];
rzz(-pi/2) q[13],q[16];
u3(pi/2,0.11309733552923254,-0.11309733552923254) q[16];
u3(pi/2,3.0586546075350225,-3.0586546075350225) q[16];
rzz(pi/2) q[13],q[16];
rzz(pi/2) q[13],q[17];
u3(pi/2,2.00370779445957,-2.00370779445957) q[17];
u3(pi/2,5.046654438726644,-5.046654438726644) q[17];
rzz(pi/2) q[13],q[17];
rzz(-pi/2) q[13],q[18];
u3(pi/2,5.191796019322492,-5.191796019322492) q[18];
u3(pi/2,2.001194520336698,-2.001194520336698) q[18];
rzz(pi/2) q[13],q[18];
u3(pi/2,1.2421857352294041,-1.2421857352294041) q[14];
u3(pi/2,2.384468824074653,-2.384468824074653) q[14];
rzz(pi/2) q[14],q[15];
u3(pi/2,4.457919975443916,-4.457919975443916) q[15];
u3(pi/2,0.5309291584566751,-0.5309291584566751) q[15];
rzz(pi/2) q[14],q[15];
rzz(pi/2) q[14],q[16];
u3(pi/2,3.0586546075350225,-3.0586546075350225) q[16];
u3(pi/2,5.8075481794260915,-5.8075481794260915) q[16];
rzz(pi/2) q[14],q[16];
rzz(pi/2) q[14],q[17];
u3(pi/2,5.046654438726644,-5.046654438726644) q[17];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[17];
rzz(pi/2) q[14],q[17];
rzz(-pi/2) q[14],q[18];
u3(pi/2,5.142787173926491,-5.142787173926491) q[18];
u3(pi/2,1.9025485110139788,-1.9025485110139788) q[18];
rzz(pi/2) q[14],q[18];
u3(pi/2,2.101725485251572,-2.101725485251572) q[15];
u3(pi/2,5.0290615198665405,-5.0290615198665405) q[15];
rzz(pi/2) q[15],q[16];
u3(pi/2,5.8075481794260915,-5.8075481794260915) q[16];
u3(pi/2,1.8805573624388503,-1.8805573624388503) q[16];
rzz(pi/2) q[15],q[16];
rzz(pi/2) q[15],q[17];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[17];
u3(pi/2,4.457919975443916,-4.457919975443916) q[17];
rzz(pi/2) q[15],q[17];
rzz(pi/2) q[15],q[18];
u3(pi/2,1.9025485110139788,-1.9025485110139788) q[18];
u3(pi/2,4.8481057830197685,-4.8481057830197685) q[18];
rzz(pi/2) q[15],q[18];
u3(pi/2,3.4513536892337466,-3.4513536892337466) q[16];
u3(pi/2,4.129309383878424,-4.129309383878424) q[16];
rzz(pi/2) q[16],q[17];
u3(pi/2,4.457919975443916,-4.457919975443916) q[17];
u3(pi/2,0.5309291584566751,-0.5309291584566751) q[17];
rzz(-pi/2) q[16],q[17];
rzz(pi/2) q[16],q[18];
u3(pi/2,4.8481057830197685,-4.8481057830197685) q[18];
u3(pi/2,1.3138140477312514,-1.3138140477312514) q[18];
rzz(pi/2) q[16],q[18];
u3(pi/2,2.101725485251572,-2.101725485251572) q[17];
u3(pi/2,2.5478316420613223,-2.5478316420613223) q[17];
rzz(pi/2) q[17],q[18];
u3(pi/2,1.3138140477312514,-1.3138140477312514) q[18];
u3(pi/2,3.670008537923596,-3.670008537923596) q[18];
rzz(-pi/2) q[17],q[18];
u3(pi,pi/2,-pi/2) q[1];
u3(pi,0,0) q[2];
u3(pi,7*pi/4,-7*pi/4) q[3];
u3(pi,3.730955435403238,-3.730955435403238) q[4];
u3(pi,0.4907167724907257,-0.4907167724907257) q[5];
u3(pi,0.5893627818134451,-0.5893627818134451) q[6];
u3(pi,0.319185813604723,-0.319185813604723) q[7];
u3(pi,2.503221026380347,-2.503221026380347) q[8];
u3(pi,2.9757165614802523,-2.9757165614802523) q[9];
u3(pi,4.525150058230738,-4.525150058230738) q[10];
u3(pi,2.5409201382234246,-2.5409201382234246) q[11];
u3(pi,4.385035025880633,-4.385035025880633) q[12];
u3(pi,3.5437165132492865,-3.5437165132492865) q[13];
u3(pi,2.674123666735632,-2.674123666735632) q[14];
u3(pi,5.970910997412761,-5.970910997412761) q[15];
u3(pi,1.4564423542042282,-1.4564423542042282) q[16];
u3(pi,2.793504187572044,-2.793504187572044) q[17];
u3(pi/2,5.240804864718492,-5.240804864718492) q[18];
u3(pi/2,6.249256106520817,-6.249256106520817) q[18];
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
measure q[19] -> c[19];
