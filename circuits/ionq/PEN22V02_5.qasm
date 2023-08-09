OPENQASM 2.0;
include "qelib1.inc";
qreg q[22];
creg c[22];
u3(pi/2,3.666238626739289,-3.666238626739289) q[20];
u3(pi/2,5.742831370762142,-5.742831370762142) q[21];
rzz(pi/2) q[20],q[21];
u3(pi,1.0304423903774522,-1.0304423903774522) q[21];
rzz(pi/2) q[20],q[21];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[19];
rzz(-pi/2) q[19],q[21];
u3(pi/2,1.035468938623196,-1.035468938623196) q[21];
u3(pi/2,4.167008495721501,-4.167008495721501) q[21];
rzz(pi/2) q[19],q[21];
u3(pi/2,0.5648583591154448,-0.5648583591154448) q[18];
rzz(pi/2) q[18],q[21];
u3(pi/2,4.167008495721501,-4.167008495721501) q[21];
u3(pi/2,1.0053096491487339,-1.0053096491487339) q[21];
rzz(pi/2) q[18],q[21];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[17];
rzz(-pi/2) q[17],q[21];
u3(pi/2,4.1469023027385274,-4.1469023027385274) q[21];
u3(pi/2,0.9657255817135024,-0.9657255817135024) q[21];
rzz(pi/2) q[17],q[21];
u3(pi/2,0.6848671984825749,-0.6848671984825749) q[16];
rzz(pi/2) q[16],q[21];
u3(pi/2,0.9657255817135024,-0.9657255817135024) q[21];
u3(pi/2,4.026893463371397,-4.026893463371397) q[21];
rzz(pi/2) q[16],q[21];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[15];
rzz(pi/2) q[15],q[21];
u3(pi/2,4.026893463371397,-4.026893463371397) q[21];
u3(pi/2,0.7257079029792423,-0.7257079029792423) q[21];
rzz(-pi/2) q[15],q[21];
u3(pi/2,1.1642742374203774,-1.1642742374203774) q[14];
rzz(pi/2) q[14],q[21];
u3(pi/2,3.8673005565690355,-3.8673005565690355) q[21];
u3(pi/2,0.4052654523130833,-0.4052654523130833) q[21];
rzz(pi/2) q[14],q[21];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[13];
rzz(pi/2) q[13],q[21];
u3(pi/2,0.4052654523130833,-0.4052654523130833) q[21];
u3(pi/2,2.9072298416319944,-2.9072298416319944) q[21];
rzz(-pi/2) q[13],q[21];
u3(pi/2,3.083787348763741,-3.083787348763741) q[12];
rzz(pi/2) q[12],q[21];
u3(pi/2,6.0488224952217875,-6.0488224952217875) q[21];
u3(pi/2,1.626716676028795,-1.626716676028795) q[21];
rzz(pi/2) q[12],q[21];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[11];
rzz(pi/2) q[11],q[21];
u3(pi/2,1.626716676028795,-1.626716676028795) q[21];
u3(pi/2,2.2085396354736244,-2.2085396354736244) q[21];
rzz(-pi/2) q[11],q[21];
u3(pi/2,1.3408317445521238,-1.3408317445521238) q[10];
rzz(pi/2) q[10],q[21];
u3(pi/2,2.2085396354736244,-2.2085396354736244) q[21];
u3(pi/2,4.187114688704476,-4.187114688704476) q[21];
rzz(pi/2) q[10],q[21];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
rzz(pi/2) q[9],q[21];
u3(pi/2,4.187114688704476,-4.187114688704476) q[21];
u3(pi/2,5.003300460107105,-5.003300460107105) q[21];
rzz(-pi/2) q[9],q[21];
u3(pi/2,0.6503096792930871,-0.6503096792930871) q[8];
rzz(pi/2) q[8],q[21];
u3(pi/2,5.003300460107105,-5.003300460107105) q[21];
u3(pi/2,0.22933626371205487,-0.22933626371205487) q[21];
rzz(pi/2) q[8],q[21];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
rzz(pi/2) q[7],q[21];
u3(pi/2,3.370928917301848,-3.370928917301848) q[21];
u3(pi/2,3.4934510307918503,-3.4934510307918503) q[21];
rzz(pi/2) q[7],q[21];
u3(pi/2,4.172663362497963,-4.172663362497963) q[6];
rzz(pi/2) q[6],q[21];
u3(pi/2,0.3518583772020568,-0.3518583772020568) q[21];
u3(pi/2,3.248406803811846,-3.248406803811846) q[21];
rzz(pi/2) q[6],q[21];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(pi/2) q[5],q[21];
u3(pi/2,3.248406803811846,-3.248406803811846) q[21];
u3(pi/2,5.899282684910913,-5.899282684910913) q[21];
rzz(pi/2) q[5],q[21];
u3(pi/2,5.693822525366141,-5.693822525366141) q[4];
rzz(-pi/2) q[4],q[21];
u3(pi/2,2.7576900313211206,-2.7576900313211206) q[21];
u3(pi/2,4.917220821398744,-4.917220821398744) q[21];
rzz(pi/2) q[4],q[21];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(pi/2) q[3],q[21];
u3(pi/2,4.917220821398744,-4.917220821398744) q[21];
u3(pi/2,6.095318066494916,-6.095318066494916) q[21];
rzz(pi/2) q[3],q[21];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
rzz(-pi/2) q[2],q[21];
u3(pi/2,6.095318066494916,-6.095318066494916) q[21];
u3(pi/2,0.5975309227127786,-0.5975309227127786) q[21];
rzz(pi/2) q[2],q[21];
u3(pi/2,pi/2,-pi/2) q[1];
u3(pi/2,5.309919903097468,-5.309919903097468) q[21];
rzz(pi/2) q[1],q[21];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,pi,-pi) q[1];
u3(pi/2,7*pi/4,-7*pi/4) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,9*pi/8,-9*pi/8) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,7*pi/8,-7*pi/8) q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[4];
u3(pi/2,4.025008507779242,-4.025008507779242) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,3.82897312619524,-3.82897312619524) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,1.03107070890817,-1.03107070890817) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,2.6018670357030667,-2.6018670357030667) q[6];
u3(pi/2,5.718955266594859,-5.718955266594859) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,1.3251237812841747,-1.3251237812841747) q[7];
u3(pi/2,4.454778382790327,-4.454778382790327) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,3.7919023328828807,-3.7919023328828807) q[8];
rzz(pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[8];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
rzz(-pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[9];
u3(pi/2,4.482424398141917,-4.482424398141917) q[10];
rzz(-pi/2) q[0],q[10];
rzz(pi/2) q[0],q[10];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[11];
rzz(-pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[11];
u3(pi/2,6.225380002353534,-6.225380002353534) q[12];
rzz(-pi/2) q[0],q[12];
rzz(pi/2) q[0],q[12];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[13];
rzz(-pi/2) q[0],q[13];
rzz(-pi/2) q[0],q[13];
u3(pi/2,4.3058668910101705,-4.3058668910101705) q[14];
rzz(-pi/2) q[0],q[14];
rzz(-pi/2) q[0],q[14];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[15];
rzz(pi/2) q[0],q[15];
rzz(-pi/2) q[0],q[15];
u3(pi/2,3.826459852072368,-3.826459852072368) q[16];
rzz(-pi/2) q[0],q[16];
rzz(-pi/2) q[0],q[16];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[17];
rzz(pi/2) q[0],q[17];
rzz(-pi/2) q[0],q[17];
u3(pi/2,3.7064510127052377,-3.7064510127052377) q[18];
rzz(-pi/2) q[0],q[18];
rzz(-pi/2) q[0],q[18];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[19];
rzz(-pi/2) q[0],q[19];
rzz(-pi/2) q[0],q[19];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[20];
rzz(-pi/2) q[0],q[20];
rzz(-pi/2) q[0],q[20];
u3(pi/2,5*pi/4,-5*pi/4) q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
rzz(-pi/2) q[1],q[2];
u3(pi/2,pi/8,-pi/8) q[2];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,4.025008507779242,-4.025008507779242) q[4];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
u3(pi/2,2.110521944681623,-2.110521944681623) q[5];
rzz(pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,5.718955266594859,-5.718955266594859) q[6];
u3(pi/2,2.5277254490783476,-2.5277254490783476) q[6];
rzz(-pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u3(pi/2,4.454778382790327,-4.454778382790327) q[7];
u3(pi/2,1.288681306502533,-1.288681306502533) q[7];
rzz(pi/2) q[1],q[7];
rzz(pi/2) q[1],q[8];
u3(pi/2,2.214822820780804,-2.214822820780804) q[8];
u3(pi/2,5.344477422286956,-5.344477422286956) q[8];
rzz(-pi/2) q[1],q[8];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[10];
rzz(pi/2) q[1],q[10];
rzz(-pi/2) q[1],q[11];
rzz(pi/2) q[1],q[11];
rzz(pi/2) q[1],q[12];
rzz(pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[13];
rzz(pi/2) q[1],q[13];
rzz(pi/2) q[1],q[14];
rzz(pi/2) q[1],q[14];
rzz(pi/2) q[1],q[15];
rzz(-pi/2) q[1],q[15];
rzz(pi/2) q[1],q[16];
rzz(pi/2) q[1],q[16];
rzz(-pi/2) q[1],q[17];
rzz(pi/2) q[1],q[17];
rzz(pi/2) q[1],q[18];
rzz(pi/2) q[1],q[18];
rzz(pi/2) q[1],q[19];
rzz(-pi/2) q[1],q[19];
rzz(pi/2) q[1],q[20];
rzz(pi/2) q[1],q[20];
u3(pi/2,11*pi/8,-11*pi/8) q[2];
u3(pi/2,0,0) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[3];
u3(pi/2,4.516353598800687,-4.516353598800687) q[3];
rzz(-pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[4];
u3(pi/2,3.436274044496516,-3.436274044496516) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,2.110521944681623,-2.110521944681623) q[5];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[5];
rzz(pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[6];
u3(pi/2,2.5277254490783476,-2.5277254490783476) q[6];
u3(pi/2,5.571300411876139,-5.571300411876139) q[6];
rzz(pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,1.288681306502533,-1.288681306502533) q[7];
u3(pi/2,4.381265114696325,-4.381265114696325) q[7];
rzz(pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[8];
u3(pi/2,5.344477422286956,-5.344477422286956) q[8];
u3(pi/2,2.1783803459991624,-2.1783803459991624) q[8];
rzz(pi/2) q[2],q[8];
rzz(pi/2) q[2],q[9];
u3(pi/2,5.237034953534185,-5.237034953534185) q[9];
u3(pi/2,2.082875929330033,-2.082875929330033) q[9];
rzz(pi/2) q[2],q[9];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[12];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[14];
rzz(pi/2) q[2],q[14];
rzz(pi/2) q[2],q[15];
rzz(pi/2) q[2],q[15];
rzz(pi/2) q[2],q[16];
rzz(-pi/2) q[2],q[16];
rzz(pi/2) q[2],q[17];
rzz(pi/2) q[2],q[17];
rzz(pi/2) q[2],q[18];
rzz(-pi/2) q[2],q[18];
rzz(pi/2) q[2],q[19];
rzz(pi/2) q[2],q[19];
rzz(pi/2) q[2],q[20];
rzz(pi/2) q[2],q[20];
u3(pi/2,2.944928953475072,-2.944928953475072) q[3];
u3(pi/2,pi,-pi) q[3];
rzz(-pi/2) q[3],q[4];
u3(pi/2,0.2946813909067226,-0.2946813909067226) q[4];
u3(pi/2,2.6508758810990676,-2.6508758810990676) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[5];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[5];
rzz(pi/2) q[3],q[5];
rzz(-pi/2) q[3],q[6];
u3(pi/2,2.4297077582863458,-2.4297077582863458) q[6];
u3(pi/2,5.3752650302921365,-5.3752650302921365) q[6];
rzz(pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,4.381265114696325,-4.381265114696325) q[7];
u3(pi/2,1.141026451783813,-1.141026451783813) q[7];
rzz(pi/2) q[3],q[7];
rzz(pi/2) q[3],q[8];
u3(pi/2,2.1783803459991624,-2.1783803459991624) q[8];
u3(pi/2,5.270964154192955,-5.270964154192955) q[8];
rzz(-pi/2) q[3],q[8];
rzz(pi/2) q[3],q[9];
u3(pi/2,2.082875929330033,-2.082875929330033) q[9];
u3(pi/2,5.199964160221826,-5.199964160221826) q[9];
rzz(pi/2) q[3],q[9];
rzz(pi/2) q[3],q[10];
u3(pi/2,2.904088248978405,-2.904088248978405) q[10];
u3(pi/2,6.033114531953839,-6.033114531953839) q[10];
rzz(-pi/2) q[3],q[10];
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
rzz(-pi/2) q[3],q[17];
rzz(pi/2) q[3],q[17];
rzz(pi/2) q[3],q[18];
rzz(pi/2) q[3],q[18];
rzz(-pi/2) q[3],q[19];
rzz(pi/2) q[3],q[19];
rzz(pi/2) q[3],q[20];
rzz(pi/2) q[3],q[20];
u3(pi/2,1.0800795543041708,-1.0800795543041708) q[4];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[5];
u3(pi/2,3.8779819715912405,-3.8779819715912405) q[5];
rzz(pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,5.3752650302921365,-5.3752650302921365) q[6];
u3(pi/2,1.8409732950036186,-1.8409732950036186) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,1.141026451783813,-1.141026451783813) q[7];
u3(pi/2,4.0865837237896026,-4.0865837237896026) q[7];
rzz(-pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,2.1293715006031615,-2.1293715006031615) q[8];
u3(pi/2,5.172318144870236,-5.172318144870236) q[8];
rzz(pi/2) q[4],q[8];
rzz(pi/2) q[4],q[9];
u3(pi/2,5.199964160221826,-5.199964160221826) q[9];
u3(pi/2,2.0093626612360316,-2.0093626612360316) q[9];
rzz(-pi/2) q[4],q[9];
rzz(pi/2) q[4],q[10];
u3(pi/2,2.8915218783640455,-2.8915218783640455) q[10];
u3(pi/2,6.008610109255838,-6.008610109255838) q[10];
rzz(pi/2) q[4],q[10];
rzz(pi/2) q[4],q[11];
u3(pi/2,2.0866458405143407,-2.0866458405143407) q[11];
u3(pi/2,5.2156721234897745,-5.2156721234897745) q[11];
rzz(-pi/2) q[4],q[11];
rzz(pi/2) q[4],q[12];
rzz(pi/2) q[4],q[12];
rzz(pi/2) q[4],q[13];
rzz(-pi/2) q[4],q[13];
rzz(pi/2) q[4],q[14];
rzz(pi/2) q[4],q[14];
rzz(pi/2) q[4],q[15];
rzz(pi/2) q[4],q[15];
rzz(-pi/2) q[4],q[16];
rzz(pi/2) q[4],q[16];
rzz(pi/2) q[4],q[17];
rzz(pi/2) q[4],q[17];
rzz(-pi/2) q[4],q[18];
rzz(pi/2) q[4],q[18];
rzz(pi/2) q[4],q[19];
rzz(pi/2) q[4],q[19];
rzz(pi/2) q[4],q[20];
rzz(-pi/2) q[4],q[20];
u3(pi/2,2.3071856447963444,-2.3071856447963444) q[5];
u3(pi/2,3.436274044496516,-3.436274044496516) q[5];
rzz(pi/2) q[5],q[6];
u3(pi/2,1.8409732950036186,-1.8409732950036186) q[6];
u3(pi/2,4.197167785195964,-4.197167785195964) q[6];
rzz(pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,0.9449910701998098,-0.9449910701998098) q[7];
u3(pi/2,3.6938846420908784,-3.6938846420908784) q[7];
rzz(-pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,5.172318144870236,-5.172318144870236) q[8];
u3(pi/2,1.834690109696439,-1.834690109696439) q[8];
rzz(pi/2) q[5],q[8];
rzz(pi/2) q[5],q[9];
u3(pi/2,5.150955314825825,-5.150955314825825) q[9];
u3(pi/2,1.9113449704440304,-1.9113449704440304) q[9];
rzz(pi/2) q[5],q[9];
rzz(-pi/2) q[5],q[10];
u3(pi/2,2.867017455666045,-2.867017455666045) q[10];
u3(pi/2,5.959601263859837,-5.959601263859837) q[10];
rzz(pi/2) q[5],q[10];
rzz(pi/2) q[5],q[11];
u3(pi/2,2.0740794698999814,-2.0740794698999814) q[11];
u3(pi/2,5.191167700791774,-5.191167700791774) q[11];
rzz(-pi/2) q[5],q[11];
rzz(pi/2) q[5],q[12];
u3(pi/2,4.646415534659305,-4.646415534659305) q[12];
u3(pi/2,1.4928848289858698,-1.4928848289858698) q[12];
rzz(pi/2) q[5],q[12];
rzz(pi/2) q[5],q[13];
rzz(pi/2) q[5],q[13];
rzz(-pi/2) q[5],q[14];
rzz(pi/2) q[5],q[14];
rzz(pi/2) q[5],q[15];
rzz(pi/2) q[5],q[15];
rzz(-pi/2) q[5],q[16];
rzz(pi/2) q[5],q[16];
rzz(pi/2) q[5],q[17];
rzz(pi/2) q[5],q[17];
rzz(pi/2) q[5],q[18];
rzz(pi/2) q[5],q[18];
rzz(pi/2) q[5],q[19];
rzz(pi/2) q[5],q[19];
rzz(pi/2) q[5],q[20];
rzz(-pi/2) q[5],q[20];
u3(pi/2,2.626371458401067,-2.626371458401067) q[6];
u3(pi/2,3.9759996623832423,-3.9759996623832423) q[6];
rzz(pi/2) q[6],q[7];
u3(pi/2,0.5522919885010856,-0.5522919885010856) q[7];
u3(pi/2,2.90848647869343,-2.90848647869343) q[7];
rzz(pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,1.834690109696439,-1.834690109696439) q[8];
u3(pi/2,4.583583681587508,-4.583583681587508) q[8];
rzz(-pi/2) q[6],q[8];
rzz(pi/2) q[6],q[9];
u3(pi/2,1.9113449704440304,-1.9113449704440304) q[9];
u3(pi/2,4.856273923919102,-4.856273923919102) q[9];
rzz(pi/2) q[6],q[9];
rzz(pi/2) q[6],q[10];
u3(pi/2,5.959601263859837,-5.959601263859837) q[10];
u3(pi/2,2.7199909194780427,-2.7199909194780427) q[10];
rzz(pi/2) q[6],q[10];
rzz(-pi/2) q[6],q[11];
u3(pi/2,5.191167700791774,-5.191167700791774) q[11];
u3(pi/2,2.0005662018059804,-2.0005662018059804) q[11];
rzz(pi/2) q[6],q[11];
rzz(pi/2) q[6],q[12];
u3(pi/2,1.4928848289858698,-1.4928848289858698) q[12];
u3(pi/2,4.609973059877663,-4.609973059877663) q[12];
rzz(-pi/2) q[6],q[12];
rzz(pi/2) q[6],q[13];
u3(pi/2,5.230123449696288,-5.230123449696288) q[13];
u3(pi/2,2.0765927440228533,-2.0765927440228533) q[13];
rzz(pi/2) q[6],q[13];
rzz(pi/2) q[6],q[14];
rzz(pi/2) q[6],q[14];
rzz(-pi/2) q[6],q[15];
rzz(pi/2) q[6],q[15];
rzz(pi/2) q[6],q[16];
rzz(pi/2) q[6],q[16];
rzz(-pi/2) q[6],q[17];
rzz(pi/2) q[6],q[17];
rzz(pi/2) q[6],q[18];
rzz(pi/2) q[6],q[18];
rzz(pi/2) q[6],q[19];
rzz(-pi/2) q[6],q[19];
rzz(pi/2) q[6],q[20];
rzz(pi/2) q[6],q[20];
u3(pi/2,4.479282805488327,-4.479282805488327) q[7];
u3(pi/2,6.160663193689585,-6.160663193689585) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,1.441991027997715,-1.441991027997715) q[8];
u3(pi/2,3.79818551819006,-3.79818551819006) q[8];
rzz(-pi/2) q[7],q[8];
rzz(pi/2) q[7],q[9];
u3(pi/2,4.856273923919102,-4.856273923919102) q[9];
u3(pi/2,1.321982188630585,-1.321982188630585) q[9];
rzz(pi/2) q[7],q[9];
rzz(pi/2) q[7],q[10];
u3(pi/2,2.7199909194780427,-2.7199909194780427) q[10];
u3(pi/2,5.6649198729531145,-5.6649198729531145) q[10];
rzz(pi/2) q[7],q[10];
rzz(-pi/2) q[7],q[11];
u3(pi/2,5.142158855395773,-5.142158855395773) q[11];
u3(pi/2,1.9025485110139788,-1.9025485110139788) q[11];
rzz(pi/2) q[7],q[11];
rzz(pi/2) q[7],q[12];
u3(pi/2,1.4683804062878691,-1.4683804062878691) q[12];
u3(pi/2,4.560964214481662,-4.560964214481662) q[12];
rzz(pi/2) q[7],q[12];
rzz(-pi/2) q[7],q[13];
u3(pi/2,5.218185397612647,-5.218185397612647) q[13];
u3(pi/2,2.052088321324853,-2.052088321324853) q[13];
rzz(pi/2) q[7],q[13];
rzz(pi/2) q[7],q[14];
u3(pi/2,5.867238439844297,-5.867238439844297) q[14];
u3(pi/2,2.7130794156401454,-2.7130794156401454) q[14];
rzz(pi/2) q[7],q[14];
rzz(-pi/2) q[7],q[15];
rzz(pi/2) q[7],q[15];
rzz(pi/2) q[7],q[16];
rzz(pi/2) q[7],q[16];
rzz(-pi/2) q[7],q[17];
rzz(pi/2) q[7],q[17];
rzz(pi/2) q[7],q[18];
rzz(pi/2) q[7],q[18];
rzz(pi/2) q[7],q[19];
rzz(pi/2) q[7],q[19];
rzz(pi/2) q[7],q[20];
rzz(pi/2) q[7],q[20];
u3(pi/2,2.2273891913951633,-2.2273891913951633) q[8];
u3(pi/2,2.282681222098344,-2.282681222098344) q[8];
rzz(pi/2) q[8],q[9];
u3(pi/2,1.321982188630585,-1.321982188630585) q[9];
u3(pi/2,3.67817667882293,-3.67817667882293) q[9];
rzz(-pi/2) q[8],q[9];
rzz(pi/2) q[8],q[10];
u3(pi/2,5.6649198729531145,-5.6649198729531145) q[10];
u3(pi/2,2.1306281376645977,-2.1306281376645977) q[10];
rzz(pi/2) q[8],q[10];
rzz(pi/2) q[8],q[11];
u3(pi/2,1.9025485110139788,-1.9025485110139788) q[11];
u3(pi/2,4.8481057830197685,-4.8481057830197685) q[11];
rzz(-pi/2) q[8],q[11];
rzz(pi/2) q[8],q[12];
u3(pi/2,4.560964214481662,-4.560964214481662) q[12];
u3(pi/2,1.320725551569149,-1.320725551569149) q[12];
rzz(pi/2) q[8],q[12];
rzz(pi/2) q[8],q[13];
u3(pi/2,2.052088321324853,-2.052088321324853) q[13];
u3(pi/2,5.144672129518645,-5.144672129518645) q[13];
rzz(pi/2) q[8],q[13];
rzz(-pi/2) q[8],q[14];
u3(pi/2,5.854672069229938,-5.854672069229938) q[14];
u3(pi/2,2.688574992942145,-2.688574992942145) q[14];
rzz(pi/2) q[8],q[14];
rzz(pi/2) q[8],q[15];
u3(pi/2,5.230751768227005,-5.230751768227005) q[15];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[15];
rzz(pi/2) q[8],q[15];
rzz(-pi/2) q[8],q[16];
rzz(pi/2) q[8],q[16];
rzz(pi/2) q[8],q[17];
rzz(pi/2) q[8],q[17];
rzz(-pi/2) q[8],q[18];
rzz(pi/2) q[8],q[18];
rzz(pi/2) q[8],q[19];
rzz(pi/2) q[8],q[19];
rzz(-pi/2) q[8],q[20];
rzz(pi/2) q[8],q[20];
u3(pi/2,2.107380352028033,-2.107380352028033) q[9];
u3(pi/2,4.491220857571968,-4.491220857571968) q[9];
rzz(pi/2) q[9],q[10];
u3(pi/2,2.1306281376645977,-2.1306281376645977) q[10];
u3(pi/2,4.4868226278569425,-4.4868226278569425) q[10];
rzz(pi/2) q[9],q[10];
rzz(pi/2) q[9],q[11];
u3(pi/2,1.7065131294299756,-1.7065131294299756) q[11];
u3(pi/2,4.455406701321044,-4.455406701321044) q[11];
rzz(-pi/2) q[9],q[11];
rzz(pi/2) q[9],q[12];
u3(pi/2,1.320725551569149,-1.320725551569149) q[12];
u3(pi/2,4.266282823574939,-4.266282823574939) q[12];
rzz(pi/2) q[9],q[12];
rzz(pi/2) q[9],q[13];
u3(pi/2,5.144672129518645,-5.144672129518645) q[13];
u3(pi/2,1.9044334666061324,-1.9044334666061324) q[13];
rzz(-pi/2) q[9],q[13];
rzz(pi/2) q[9],q[14];
u3(pi/2,2.688574992942145,-2.688574992942145) q[14];
u3(pi/2,5.781158801135938,-5.781158801135938) q[14];
rzz(pi/2) q[9],q[14];
rzz(pi/2) q[9],q[15];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[15];
u3(pi/2,5.194309293445364,-5.194309293445364) q[15];
rzz(pi/2) q[9],q[15];
rzz(-pi/2) q[9],q[16];
u3(pi/2,2.24686706584742,-2.24686706584742) q[16];
u3(pi/2,5.376521667353572,-5.376521667353572) q[16];
rzz(pi/2) q[9],q[16];
rzz(pi/2) q[9],q[17];
rzz(pi/2) q[9],q[17];
rzz(-pi/2) q[9],q[18];
rzz(pi/2) q[9],q[18];
rzz(pi/2) q[9],q[19];
rzz(pi/2) q[9],q[19];
rzz(pi/2) q[9],q[20];
rzz(pi/2) q[9],q[20];
u3(pi/2,2.916026301062046,-2.916026301062046) q[10];
u3(pi/2,5.65109686527732,-5.65109686527732) q[10];
rzz(pi/2) q[10],q[11];
u3(pi/2,1.3138140477312514,-1.3138140477312514) q[11];
u3(pi/2,3.670008537923596,-3.670008537923596) q[11];
rzz(pi/2) q[10],q[11];
rzz(-pi/2) q[10],q[12];
u3(pi/2,1.124690169985146,-1.124690169985146) q[12];
u3(pi/2,3.873583741876215,-3.873583741876215) q[12];
rzz(pi/2) q[10],q[12];
rzz(pi/2) q[10],q[13];
u3(pi/2,5.046026120195926,-5.046026120195926) q[13];
u3(pi/2,1.7083980850221294,-1.7083980850221294) q[13];
rzz(pi/2) q[10],q[13];
rzz(pi/2) q[10],q[14];
u3(pi/2,5.781158801135938,-5.781158801135938) q[14];
u3(pi/2,2.5415484567541426,-2.5415484567541426) q[14];
rzz(pi/2) q[10],q[14];
rzz(pi/2) q[10],q[15];
u3(pi/2,5.194309293445364,-5.194309293445364) q[15];
u3(pi/2,2.00370779445957,-2.00370779445957) q[15];
rzz(pi/2) q[10],q[15];
rzz(pi/2) q[10],q[16];
u3(pi/2,5.376521667353572,-5.376521667353572) q[16];
u3(pi/2,2.2104245910657783,-2.2104245910657783) q[16];
rzz(-pi/2) q[10],q[16];
rzz(pi/2) q[10],q[17];
u3(pi/2,2.0860175219836226,-2.0860175219836226) q[17];
u3(pi/2,5.2156721234897745,-5.2156721234897745) q[17];
rzz(pi/2) q[10],q[17];
rzz(pi/2) q[10],q[18];
rzz(-pi/2) q[10],q[18];
rzz(pi/2) q[10],q[19];
rzz(pi/2) q[10],q[19];
rzz(pi/2) q[10],q[20];
rzz(pi/2) q[10],q[20];
u3(pi/2,5.240804864718492,-5.240804864718492) q[11];
u3(pi/2,6.229149913537841,-6.229149913537841) q[11];
rzz(-pi/2) q[11],q[12];
u3(pi/2,0.7319910882864218,-0.7319910882864218) q[12];
u3(pi/2,3.0881855784787664,-3.0881855784787664) q[12];
rzz(pi/2) q[11],q[12];
rzz(pi/2) q[11],q[13];
u3(pi/2,1.7083980850221294,-1.7083980850221294) q[13];
u3(pi/2,4.457291656913199,-4.457291656913199) q[13];
rzz(pi/2) q[11],q[13];
rzz(-pi/2) q[11],q[14];
u3(pi/2,5.683141110343936,-5.683141110343936) q[14];
u3(pi/2,2.3448847566394213,-2.3448847566394213) q[14];
rzz(pi/2) q[11],q[14];
rzz(pi/2) q[11],q[15];
u3(pi/2,2.00370779445957,-2.00370779445957) q[15];
u3(pi/2,5.046654438726644,-5.046654438726644) q[15];
rzz(pi/2) q[11],q[15];
rzz(pi/2) q[11],q[16];
u3(pi/2,5.352017244655571,-5.352017244655571) q[16];
u3(pi/2,2.1607874271390597,-2.1607874271390597) q[16];
rzz(pi/2) q[11],q[16];
rzz(pi/2) q[11],q[17];
u3(pi/2,5.2156721234897745,-5.2156721234897745) q[17];
u3(pi/2,2.0489467286712633,-2.0489467286712633) q[17];
rzz(pi/2) q[11],q[17];
rzz(-pi/2) q[11],q[18];
u3(pi/2,5.265309287416493,-5.265309287416493) q[18];
u3(pi/2,2.111778581743059,-2.111778581743059) q[18];
rzz(pi/2) q[11],q[18];
rzz(pi/2) q[11],q[19];
rzz(pi/2) q[11],q[19];
rzz(pi/2) q[11],q[20];
rzz(-pi/2) q[11],q[20];
u3(pi/2,1.51738925168387,-1.51738925168387) q[12];
u3(pi/2,1.808300731406285,-1.808300731406285) q[12];
rzz(pi/2) q[12],q[13];
u3(pi/2,4.457291656913199,-4.457291656913199) q[13];
u3(pi/2,0.5303008399259571,-0.5303008399259571) q[13];
rzz(pi/2) q[12],q[13];
rzz(pi/2) q[12],q[14];
u3(pi/2,2.3448847566394213,-2.3448847566394213) q[14];
u3(pi/2,5.09377832853049,-5.09377832853049) q[14];
rzz(-pi/2) q[12],q[14];
rzz(pi/2) q[12],q[15];
u3(pi/2,5.046654438726644,-5.046654438726644) q[15];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[15];
rzz(pi/2) q[12],q[15];
rzz(pi/2) q[12],q[16];
u3(pi/2,2.1607874271390597,-2.1607874271390597) q[16];
u3(pi/2,5.204362389936851,-5.204362389936851) q[16];
rzz(pi/2) q[12],q[16];
rzz(-pi/2) q[12],q[17];
u3(pi/2,5.190539382261056,-5.190539382261056) q[17];
u3(pi/2,1.9999378832752626,-1.9999378832752626) q[17];
rzz(pi/2) q[12],q[17];
rzz(pi/2) q[12],q[18];
u3(pi/2,2.111778581743059,-2.111778581743059) q[18];
u3(pi/2,5.228866812634852,-5.228866812634852) q[18];
rzz(pi/2) q[12],q[18];
rzz(-pi/2) q[12],q[19];
u3(pi/2,5.230751768227005,-5.230751768227005) q[19];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[19];
rzz(pi/2) q[12],q[19];
rzz(pi/2) q[12],q[20];
rzz(pi/2) q[12],q[20];
u3(pi/2,5.242689820310647,-5.242689820310647) q[13];
u3(pi/2,6.173229564303944,-6.173229564303944) q[13];
rzz(pi/2) q[13],q[14];
u3(pi/2,1.9521856749406974,-1.9521856749406974) q[14];
u3(pi/2,4.308380165133042,-4.308380165133042) q[14];
rzz(pi/2) q[13],q[14];
rzz(pi/2) q[13],q[15];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[15];
u3(pi/2,4.457919975443916,-4.457919975443916) q[15];
rzz(pi/2) q[13],q[15];
rzz(pi/2) q[13],q[16];
u3(pi/2,5.204362389936851,-5.204362389936851) q[16];
u3(pi/2,1.8667343547630548,-1.8667343547630548) q[16];
rzz(pi/2) q[13],q[16];
rzz(pi/2) q[13],q[17];
u3(pi/2,1.9999378832752626,-1.9999378832752626) q[17];
u3(pi/2,5.043512846073054,-5.043512846073054) q[17];
rzz(pi/2) q[13],q[17];
rzz(pi/2) q[13],q[18];
u3(pi/2,5.228866812634852,-5.228866812634852) q[18];
u3(pi/2,2.038265313649058,-2.038265313649058) q[18];
rzz(pi/2) q[13],q[18];
rzz(pi/2) q[13],q[19];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[19];
u3(pi/2,5.194309293445364,-5.194309293445364) q[19];
rzz(pi/2) q[13],q[19];
rzz(pi/2) q[13],q[20];
u3(pi/2,2.0872741590450583,-2.0872741590450583) q[20];
u3(pi/2,5.216300442020493,-5.216300442020493) q[20];
rzz(-pi/2) q[13],q[20];
u3(pi/2,2.7375838383381454,-2.7375838383381454) q[14];
u3(pi/2,3.9885660329976016,-3.9885660329976016) q[14];
rzz(pi/2) q[14],q[15];
u3(pi/2,4.457919975443916,-4.457919975443916) q[15];
u3(pi/2,0.5309291584566751,-0.5309291584566751) q[15];
rzz(pi/2) q[14],q[15];
rzz(pi/2) q[14],q[16];
u3(pi/2,1.8667343547630548,-1.8667343547630548) q[16];
u3(pi/2,4.615627926654124,-4.615627926654124) q[16];
rzz(-pi/2) q[14],q[16];
rzz(pi/2) q[14],q[17];
u3(pi/2,5.043512846073054,-5.043512846073054) q[17];
u3(pi/2,1.7058848108992577,-1.7058848108992577) q[17];
rzz(pi/2) q[14],q[17];
rzz(pi/2) q[14],q[18];
u3(pi/2,2.038265313649058,-2.038265313649058) q[18];
u3(pi/2,5.081211957916131,-5.081211957916131) q[18];
rzz(pi/2) q[14],q[18];
rzz(-pi/2) q[14],q[19];
u3(pi/2,2.052716639855571,-2.052716639855571) q[19];
u3(pi/2,5.145300448049363,-5.145300448049363) q[19];
rzz(pi/2) q[14],q[19];
rzz(pi/2) q[14],q[20];
u3(pi/2,2.0747077884306995,-2.0747077884306995) q[20];
u3(pi/2,5.191796019322492,-5.191796019322492) q[20];
rzz(pi/2) q[14],q[20];
u3(pi/2,5.2433181388413646,-5.2433181388413646) q[15];
u3(pi/2,0.37070793312359557,-0.37070793312359557) q[15];
rzz(-pi/2) q[15],q[16];
u3(pi/2,4.615627926654124,-4.615627926654124) q[16];
u3(pi/2,0.6886371096668826,-0.6886371096668826) q[16];
rzz(pi/2) q[15],q[16];
rzz(pi/2) q[15],q[17];
u3(pi/2,1.7058848108992577,-1.7058848108992577) q[17];
u3(pi/2,4.454778382790327,-4.454778382790327) q[17];
rzz(pi/2) q[15],q[17];
rzz(pi/2) q[15],q[18];
u3(pi/2,5.081211957916131,-5.081211957916131) q[18];
u3(pi/2,1.7435839227423353,-1.7435839227423353) q[18];
rzz(-pi/2) q[15],q[18];
rzz(pi/2) q[15],q[19];
u3(pi/2,5.145300448049363,-5.145300448049363) q[19];
u3(pi/2,1.9050617851368508,-1.9050617851368508) q[19];
rzz(pi/2) q[15],q[19];
rzz(pi/2) q[15],q[20];
u3(pi/2,5.191796019322492,-5.191796019322492) q[20];
u3(pi/2,2.001194520336698,-2.001194520336698) q[20];
rzz(pi/2) q[15],q[20];
u3(pi/2,5.401026090051572,-5.401026090051572) q[16];
u3(pi/2,0.608212337734984,-0.608212337734984) q[16];
rzz(pi/2) q[16],q[17];
u3(pi/2,4.454778382790327,-4.454778382790327) q[17];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[17];
rzz(pi/2) q[16],q[17];
rzz(pi/2) q[16],q[18];
u3(pi/2,4.885176576332128,-4.885176576332128) q[18];
u3(pi/2,1.350884841043611,-1.350884841043611) q[18];
rzz(-pi/2) q[16],q[18];
rzz(pi/2) q[16],q[19];
u3(pi/2,1.9050617851368508,-1.9050617851368508) q[19];
u3(pi/2,4.850619057142641,-4.850619057142641) q[19];
rzz(pi/2) q[16],q[19];
rzz(pi/2) q[16],q[20];
u3(pi/2,2.001194520336698,-2.001194520336698) q[20];
u3(pi/2,5.044769483134489,-5.044769483134489) q[20];
rzz(-pi/2) q[16],q[20];
u3(pi/2,5.240176546187775,-5.240176546187775) q[17];
u3(pi/2,0.4875751798371359,-0.4875751798371359) q[17];
rzz(pi/2) q[17],q[18];
u3(pi/2,4.492477494633404,-4.492477494633404) q[18];
u3(pi/2,0.5654866776461628,-0.5654866776461628) q[18];
rzz(pi/2) q[17],q[18];
rzz(pi/2) q[17],q[19];
u3(pi/2,4.850619057142641,-4.850619057142641) q[19];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[19];
rzz(pi/2) q[17],q[19];
rzz(-pi/2) q[17],q[20];
u3(pi/2,5.044769483134489,-5.044769483134489) q[20];
u3(pi/2,1.7065131294299756,-1.7065131294299756) q[20];
rzz(pi/2) q[17],q[20];
u3(pi/2,5.277875658030852,-5.277875658030852) q[18];
u3(pi/2,0.5453804846631881,-0.5453804846631881) q[18];
rzz(pi/2) q[18],q[19];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[19];
u3(pi/2,3.6725218120464684,-3.6725218120464684) q[19];
rzz(pi/2) q[18],q[19];
rzz(-pi/2) q[18],q[20];
u3(pi/2,4.8481057830197685,-4.8481057830197685) q[20];
u3(pi/2,1.3138140477312514,-1.3138140477312514) q[20];
rzz(pi/2) q[18],q[20];
u3(pi/2,2.101725485251572,-2.101725485251572) q[19];
u3(pi/2,3.6624687155549807,-3.6624687155549807) q[19];
rzz(pi/2) q[19],q[20];
u3(pi/2,1.3138140477312514,-1.3138140477312514) q[20];
u3(pi/2,3.670008537923596,-3.670008537923596) q[20];
rzz(pi/2) q[19],q[20];
u3(pi,0,0) q[1];
u3(pi,0,0) q[2];
u3(pi,7*pi/4,-7*pi/4) q[3];
u3(pi,2.1601591086083416,-2.1601591086083416) q[4];
u3(pi,2.061513099285622,-2.061513099285622) q[5];
u3(pi,0.14702653618800232,-0.14702653618800232) q[6];
u3(pi,3.460778467194516,-3.460778467194516) q[7];
u3(pi,5.64481367997014,-5.64481367997014) q[8];
u3(pi,5.387203082375778,-5.387203082375778) q[9];
u3(pi,0.6013008338970863,-0.6013008338970863) q[10];
u3(pi,2.930477627268559,-2.930477627268559) q[11];
u3(pi,2.5126458043411164,-2.5126458043411164) q[12];
u3(pi,1.8246370132049519,-1.8246370132049519) q[13];
u3(pi,3.382866969385489,-3.382866969385489) q[14];
u3(pi,4.006787270388423,-4.006787270388423) q[15];
u3(pi,2.4240528915098842,-2.4240528915098842) q[16];
u3(pi,5.333796007264751,-5.333796007264751) q[17];
u3(pi,5.322486273711827,-5.322486273711827) q[18];
u3(pi,3.7064510127052377,-3.7064510127052377) q[19];
u3(pi/2,3.670008537923596,-3.670008537923596) q[20];
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
measure q[20] -> c[20];
measure q[21] -> c[21];
