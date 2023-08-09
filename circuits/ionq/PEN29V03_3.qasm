OPENQASM 2.0;
include "qelib1.inc";
qreg q[29];
creg c[29];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[27];
u3(pi,1.971035230862236,-1.971035230862236) q[28];
rzz(0.3711685063145751) q[27],q[28];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[26];
rzz(0.7422307325496793) q[26],q[28];
u3(pi/2,0.3518583772020568,-0.3518583772020568) q[25];
rzz(1.4843686153284816) q[25],q[28];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[24];
u3(pi,5.793096853219579,-5.793096853219579) q[28];
rzz(0.17261834264322648) q[24],q[28];
u3(pi/2,2.9757165614802523,-2.9757165614802523) q[23];
rzz(0.3453666058506421) q[23],q[28];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[22];
rzz(0.6906816738738022) q[22],q[28];
u3(pi/2,4.0482562934158075,-4.0482562934158075) q[21];
rzz(1.3811718205514782) q[21],q[28];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[20];
u3(pi,2.1212033597038285,-2.1212033597038285) q[28];
rzz(0.3789977164904761) q[20],q[28];
u3(pi/2,2.0439201804255194,-2.0439201804255194) q[19];
rzz(0.7580790150535008) q[19],q[28];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[18];
rzz(1.5159908659619044) q[18],q[28];
u3(pi/2,0.30913271711323564,-0.30913271711323564) q[17];
u3(pi,2.975088242949534,-2.975088242949534) q[28];
rzz(0.10931240753204062) q[17],q[28];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[16];
rzz(0.2186605192645893) q[16],q[28];
u3(pi/2,5.935725159692555,-5.935725159692555) q[15];
rzz(0.43741539626452913) q[15],q[28];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[14];
rzz(0.8746420770583572) q[14],q[28];
u3(pi/2,3.311238656883642,-3.311238656883642) q[13];
u3(pi,5.28478716186875,-5.28478716186875) q[28];
rzz(1.3920930490488952) q[13],q[28];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[12];
u3(pi,2.3417431639858317,-2.3417431639858317) q[28];
rzz(0.357401717439316) q[12],q[28];
u3(pi/2,5.379663260007161,-5.379663260007161) q[11];
rzz(0.7148460977068678) q[11],q[28];
u3(pi/2,3.669380219392878,-3.669380219392878) q[10];
rzz(1.429606869757264) q[10],q[28];
u3(pi/2,1.0920176063878122,-1.0920176063878122) q[9];
u3(pi,1.7385573744965914,-1.7385573744965914) q[28];
rzz(0.28225829262533036) q[9],q[28];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[8];
rzz(0.5645445925491676) q[8],q[28];
u3(pi/2,1.9879998311916212,-1.9879998311916212) q[7];
rzz(1.1290331705013215) q[7],q[28];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
u3(pi,2.1928316722056755,-2.1928316722056755) q[28];
rzz(0.8835729338221293) q[6],q[28];
u3(pi/2,3.82897312619524,-3.82897312619524) q[5];
u3(pi,6.164433104873892,-6.164433104873892) q[28];
rzz(1.3744046257721234) q[5],q[28];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
u3(pi,0.3933274002294421,-0.3933274002294421) q[28];
rzz(pi/8) q[4],q[28];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
rzz(0.7853830994606742) q[3],q[28];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
u3(pi,9*pi/5,-9*pi/5) q[28];
rzz(pi/2) q[2],q[28];
u3(pi/2,0,0) q[1];
u3(pi/2,pi/4,-pi/4) q[1];
rzz(pi/2) q[0],q[1];
rzz(pi/8) q[0],q[2];
rzz(0.19640131429629326) q[0],q[3];
rzz(pi/32) q[0],q[4];
rzz(0.049100328574073315) q[0],q[5];
rzz(pi/128) q[0],q[6];
rzz(0.012275082143518329) q[0],q[7];
u3(pi/2,0,0) q[0];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[8];
rzz(0) q[0],q[8];
u3(pi,3.568220935947287,-3.568220935947287) q[0];
u3(pi/2,1.0920176063878122,-1.0920176063878122) q[9];
rzz(0) q[0],q[9];
u3(pi,1.2805131656031998,-1.2805131656031998) q[0];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[10];
rzz(0) q[0],q[10];
u3(pi,5.2759907024386985,-5.2759907024386985) q[0];
u3(pi/2,2.2386989249480864,-2.2386989249480864) q[11];
rzz(0) q[0],q[11];
u3(pi,2.988282932094611,-2.988282932094611) q[0];
u3(pi/2,0.5252742916802133,-0.5252742916802133) q[12];
rzz(0) q[0],q[12];
u3(pi,0.7005751617505239,-0.7005751617505239) q[0];
u3(pi/2,3.31186697541436,-3.31186697541436) q[13];
rzz(0) q[0],q[13];
u3(pi,4.696052698586023,-4.696052698586023) q[0];
u3(pi/2,3.666238626739289,-3.666238626739289) q[14];
rzz(0) q[0],q[14];
u3(pi,2.4083449282419354,-2.4083449282419354) q[0];
u3(pi/2,2.79476082463348,-2.79476082463348) q[15];
rzz(0) q[0],q[15];
u3(pi,0.12063715789784804,-0.12063715789784804) q[0];
u3(pi/2,3.666238626739289,-3.666238626739289) q[16];
rzz(0) q[0],q[16];
u3(pi,4.116114694733347,-4.116114694733347) q[0];
u3(pi/2,0.30913271711323564,-0.30913271711323564) q[17];
rzz(0) q[0],q[17];
u3(pi,1.8284069243892596,-1.8284069243892596) q[0];
u3(pi/2,3.666238626739289,-3.666238626739289) q[18];
rzz(0) q[0],q[18];
u3(pi,5.8238844612247584,-5.8238844612247584) q[0];
u3(pi/2,5.185512834015313,-5.185512834015313) q[19];
rzz(0) q[0],q[19];
u3(pi,3.536176690880671,-3.536176690880671) q[0];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[20];
rzz(0) q[0],q[20];
u3(pi,1.2484689205365838,-1.2484689205365838) q[0];
u3(pi/2,0.9066636398260144,-0.9066636398260144) q[21];
rzz(0) q[0],q[21];
u3(pi,5.243946457372083,-5.243946457372083) q[0];
u3(pi/2,3.666238626739289,-3.666238626739289) q[22];
rzz(0) q[0],q[22];
u3(pi/2,2.528982086139784,-2.528982086139784) q[0];
u3(pi/2,6.117309215070045,-6.117309215070045) q[23];
rzz(pi/2) q[0],q[23];
rzz(-pi/2) q[0],q[23];
u3(pi/2,5.6705747397295765,-5.6705747397295765) q[0];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[24];
rzz(0) q[0],q[24];
u3(pi/2,2.528982086139784,-2.528982086139784) q[0];
u3(pi/2,3.4934510307918503,-3.4934510307918503) q[25];
rzz(pi/2) q[0],q[25];
rzz(pi/2) q[0],q[25];
u3(pi/2,3*pi/4,-3*pi/4) q[1];
u3(pi/2,pi,-pi) q[1];
rzz(0.7853830994606742) q[1],q[2];
rzz(pi/8) q[1],q[3];
rzz(0.19640131429629326) q[1],q[4];
rzz(pi/32) q[1],q[5];
rzz(0.049100328574073315) q[1],q[6];
rzz(pi/128) q[1],q[7];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(0) q[1],q[9];
u3(pi,5.1390172627421835,-5.1390172627421835) q[1];
rzz(0) q[1],q[10];
u3(pi,2.851309492398096,-2.851309492398096) q[1];
rzz(0) q[1],q[11];
u3(pi,0.5636017220540089,-0.5636017220540089) q[1];
rzz(0) q[1],q[12];
u3(pi,4.559079258889508,-4.559079258889508) q[1];
rzz(0) q[1],q[13];
u3(pi,2.2713714885454204,-2.2713714885454204) q[1];
rzz(0) q[1],q[14];
u3(pi,6.266849025380919,-6.266849025380919) q[1];
rzz(0) q[1],q[15];
u3(pi,3.979141255036832,-3.979141255036832) q[1];
rzz(0) q[1],q[16];
u3(pi,1.6914334846927446,-1.6914334846927446) q[1];
rzz(0) q[1],q[17];
u3(pi,5.6869110215282435,-5.6869110215282435) q[1];
rzz(0) q[1],q[18];
u3(pi,3.3992032511841566,-3.3992032511841566) q[1];
rzz(0) q[1],q[19];
u3(pi,1.1114954808400688,-1.1114954808400688) q[1];
rzz(0) q[1],q[20];
u3(pi,5.106973017675568,-5.106973017675568) q[1];
rzz(0) q[1],q[21];
u3(pi,2.8192652473314803,-2.8192652473314803) q[1];
rzz(0) q[1],q[22];
u3(pi/2,11*pi/8,-11*pi/8) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,4.230468667324016,-4.230468667324016) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
rzz(0) q[2],q[10];
u3(pi,0.03392920065876977,-0.03392920065876977) q[2];
rzz(0) q[2],q[11];
u3(pi,4.029406737494268,-4.029406737494268) q[2];
rzz(0) q[2],q[12];
u3(pi/2,1.3150706847926874,-1.3150706847926874) q[2];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[13];
rzz(-pi/2) q[2],q[14];
rzz(pi/2) q[2],q[14];
rzz(pi/2) q[2],q[15];
rzz(pi/2) q[2],q[15];
rzz(-pi/2) q[2],q[16];
rzz(pi/2) q[2],q[16];
rzz(pi/2) q[2],q[17];
rzz(pi/2) q[2],q[17];
rzz(pi/2) q[2],q[18];
rzz(-pi/2) q[2],q[18];
rzz(pi/2) q[2],q[19];
rzz(pi/2) q[2],q[19];
rzz(pi/2) q[2],q[20];
rzz(-pi/2) q[2],q[20];
rzz(pi/2) q[2],q[21];
rzz(pi/2) q[2],q[21];
rzz(pi/2) q[2],q[22];
rzz(-pi/2) q[2],q[22];
u3(pi/2,4.456663338382481,-4.456663338382481) q[2];
rzz(0) q[2],q[23];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,3.664981989677853,-3.664981989677853) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[3];
rzz(0) q[3],q[11];
u3(pi,2.8594776332974297,-2.8594776332974297) q[3];
rzz(0) q[3],q[12];
u3(pi,5.044141164603771,-5.044141164603771) q[3];
rzz(0) q[3],q[13];
u3(pi,0.9456193887305276,-0.9456193887305276) q[3];
rzz(0) q[3],q[14];
u3(pi,3.1302829200368696,-3.1302829200368696) q[3];
rzz(0) q[3],q[15];
u3(pi,5.314946451343212,-5.314946451343212) q[3];
rzz(0) q[3],q[16];
u3(pi,1.21579635693925,-1.21579635693925) q[3];
rzz(0) q[3],q[17];
u3(pi,3.4004598882455923,-3.4004598882455923) q[3];
rzz(0) q[3],q[18];
u3(pi,5.585123419551934,-5.585123419551934) q[3];
rzz(0) q[3],q[19];
u3(pi,1.4866016436786902,-1.4866016436786902) q[3];
rzz(0) q[3],q[20];
u3(pi,3.6712651749850327,-3.6712651749850327) q[3];
rzz(0) q[3],q[21];
u3(pi,5.855928706291374,-5.855928706291374) q[3];
rzz(0) q[3],q[22];
u3(pi/2,2.235557332294497,-2.235557332294497) q[3];
rzz(pi/2) q[3],q[23];
rzz(pi/2) q[3],q[23];
u3(pi/2,2.235557332294497,-2.235557332294497) q[3];
rzz(0) q[3],q[24];
u3(pi/2,1.865477717701619,-1.865477717701619) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
u3(pi/2,5.129592484781415,-5.129592484781415) q[7];
rzz(pi/2) q[4],q[7];
u3(pi/2,3.5587961579865177,-3.5587961579865177) q[7];
u3(pi/2,0.22116812281272144,-0.22116812281272144) q[7];
rzz(-pi/2) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,5.373380074699982,-5.373380074699982) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,0.2946813909067226,-0.2946813909067226) q[4];
rzz(0) q[4],q[12];
u3(pi,2.819893565862198,-2.819893565862198) q[4];
rzz(0) q[4],q[13];
u3(pi,4.729353580714075,-4.729353580714075) q[4];
rzz(0) q[4],q[14];
u3(pi/2,0.9720087670206821,-0.9720087670206821) q[4];
rzz(pi/2) q[4],q[15];
rzz(pi/2) q[4],q[15];
rzz(pi/2) q[4],q[16];
rzz(pi/2) q[4],q[16];
rzz(pi/2) q[4],q[17];
rzz(pi/2) q[4],q[17];
rzz(pi/2) q[4],q[18];
rzz(pi/2) q[4],q[18];
rzz(-pi/2) q[4],q[19];
rzz(pi/2) q[4],q[19];
rzz(pi/2) q[4],q[20];
rzz(pi/2) q[4],q[20];
rzz(-pi/2) q[4],q[21];
rzz(pi/2) q[4],q[21];
rzz(pi/2) q[4],q[22];
rzz(pi/2) q[4],q[22];
u3(pi/2,0.9720087670206821,-0.9720087670206821) q[4];
rzz(0) q[4],q[23];
u3(pi/2,4.1136014206104745,-4.1136014206104745) q[4];
rzz(pi/2) q[4],q[24];
rzz(-pi/2) q[4],q[24];
u3(pi/2,0.9720087670206821,-0.9720087670206821) q[4];
rzz(0) q[4],q[25];
u3(pi/2,3.583300580684518,-3.583300580684518) q[5];
rzz(0.7853830994606742) q[5],q[6];
u3(pi/2,1.7919644496076181,-1.7919644496076181) q[7];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,3.660583759962827,-3.660583759962827) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,2.012504253889621,-2.012504253889621) q[5];
rzz(0) q[5],q[13];
u3(pi,5.581353508367626,-5.581353508367626) q[5];
rzz(0) q[5],q[14];
u3(pi,3.293645738023539,-3.293645738023539) q[5];
rzz(0) q[5],q[15];
u3(pi,1.0053096491487339,-1.0053096491487339) q[5];
rzz(0) q[5],q[16];
u3(pi,5.000787185984233,-5.000787185984233) q[5];
rzz(0) q[5],q[17];
u3(pi,2.7130794156401454,-2.7130794156401454) q[5];
rzz(0) q[5],q[18];
u3(pi,0.425371645296058,-0.425371645296058) q[5];
rzz(0) q[5],q[19];
u3(pi,4.4208491821315565,-4.4208491821315565) q[5];
rzz(0) q[5],q[20];
u3(pi,2.1331414117874696,-2.1331414117874696) q[5];
rzz(0) q[5],q[21];
u3(pi,6.128618948622969,-6.128618948622969) q[5];
rzz(0) q[5],q[22];
u3(pi/2,3.4142828959213873,-3.4142828959213873) q[5];
rzz(pi/2) q[5],q[23];
rzz(pi/2) q[5],q[23];
u3(pi/2,3.4142828959213873,-3.4142828959213873) q[5];
rzz(0) q[5],q[24];
u3(pi/2,0.27269024233159406,-0.27269024233159406) q[5];
rzz(pi/2) q[5],q[25];
rzz(-pi/2) q[5],q[25];
u3(pi/2,3.4142828959213873,-3.4142828959213873) q[5];
u3(pi/2,3.666238626739289,-3.666238626739289) q[26];
rzz(0) q[5],q[26];
u3(pi/2,1.2516105131901736,-1.2516105131901736) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,3.3055837901071805,-3.3055837901071805) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,5.963999493574864,-5.963999493574864) q[6];
rzz(0) q[6],q[14];
u3(pi,3.249663440873282,-3.249663440873282) q[6];
rzz(0) q[6],q[15];
u3(pi,0.9619556705291947,-0.9619556705291947) q[6];
rzz(0) q[6],q[16];
u3(pi/2,4.5301766064764815,-4.5301766064764815) q[6];
rzz(pi/2) q[6],q[17];
rzz(pi/2) q[6],q[17];
rzz(pi/2) q[6],q[18];
rzz(-pi/2) q[6],q[18];
rzz(pi/2) q[6],q[19];
rzz(pi/2) q[6],q[19];
rzz(pi/2) q[6],q[20];
rzz(-pi/2) q[6],q[20];
rzz(pi/2) q[6],q[21];
rzz(pi/2) q[6],q[21];
rzz(pi/2) q[6],q[22];
rzz(pi/2) q[6],q[22];
u3(pi/2,4.5301766064764815,-4.5301766064764815) q[6];
rzz(0) q[6],q[23];
u3(pi/2,1.3885839528866886,-1.3885839528866886) q[6];
rzz(-pi/2) q[6],q[24];
rzz(pi/2) q[6],q[24];
u3(pi/2,4.5301766064764815,-4.5301766064764815) q[6];
rzz(0) q[6],q[25];
u3(pi/2,1.3885839528866886,-1.3885839528866886) q[6];
rzz(pi/2) q[6],q[26];
rzz(pi/2) q[6],q[26];
u3(pi/2,1.3885839528866886,-1.3885839528866886) q[6];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[27];
rzz(0) q[6],q[27];
u3(pi/2,1.337690151898534,-1.337690151898534) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,2.90848647869343,-2.90848647869343) q[7];
rzz(0) q[7],q[15];
u3(pi,5.571300411876139,-5.571300411876139) q[7];
rzz(0) q[7],q[16];
u3(pi,1.472778636002895,-1.472778636002895) q[7];
rzz(0) q[7],q[17];
u3(pi,3.6574421673092368,-3.6574421673092368) q[7];
rzz(0) q[7],q[18];
u3(pi,5.842105698615579,-5.842105698615579) q[7];
rzz(0) q[7],q[19];
u3(pi,1.7435839227423353,-1.7435839227423353) q[7];
rzz(0) q[7],q[20];
u3(pi,3.928247454048677,-3.928247454048677) q[7];
rzz(0) q[7],q[21];
u3(pi,6.112910985355019,-6.112910985355019) q[7];
rzz(0) q[7],q[22];
u3(pi/2,2.4925396113581417,-2.4925396113581417) q[7];
rzz(-pi/2) q[7],q[23];
rzz(pi/2) q[7],q[23];
u3(pi/2,5.634132264947936,-5.634132264947936) q[7];
rzz(0) q[7],q[24];
u3(pi/2,2.4925396113581417,-2.4925396113581417) q[7];
rzz(pi/2) q[7],q[25];
rzz(pi/2) q[7],q[25];
u3(pi/2,2.4925396113581417,-2.4925396113581417) q[7];
rzz(0) q[7],q[26];
u3(pi/2,5.634132264947936,-5.634132264947936) q[7];
rzz(pi/2) q[7],q[27];
rzz(-pi/2) q[7],q[27];
u3(pi/2,2.6628139331827088,-2.6628139331827088) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,5.930070292916093,-5.930070292916093) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,1.0920176063878122,-1.0920176063878122) q[8];
rzz(0) q[8],q[16];
u3(pi,4.920990732583052,-4.920990732583052) q[8];
rzz(0) q[8],q[17];
u3(pi,3.15478734273487,-3.15478734273487) q[8];
rzz(0) q[8],q[18];
u3(pi/2,0.7005751617505239,-0.7005751617505239) q[8];
rzz(pi/2) q[8],q[19];
rzz(pi/2) q[8],q[19];
rzz(pi/2) q[8],q[20];
rzz(-pi/2) q[8],q[20];
rzz(pi/2) q[8],q[21];
rzz(pi/2) q[8],q[21];
rzz(pi/2) q[8],q[22];
rzz(pi/2) q[8],q[22];
u3(pi/2,3.8421678153403174,-3.8421678153403174) q[8];
rzz(0) q[8],q[23];
u3(pi/2,0.7005751617505239,-0.7005751617505239) q[8];
rzz(pi/2) q[8],q[24];
rzz(pi/2) q[8],q[24];
u3(pi/2,0.7005751617505239,-0.7005751617505239) q[8];
rzz(0) q[8],q[25];
u3(pi/2,3.8421678153403174,-3.8421678153403174) q[8];
rzz(pi/2) q[8],q[26];
rzz(-pi/2) q[8],q[26];
u3(pi/2,0.7005751617505239,-0.7005751617505239) q[8];
rzz(0) q[8],q[27];
u3(pi/2,6.074583554981224,-6.074583554981224) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,3.660583759962827,-3.660583759962827) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,1.3621945745965343,-1.3621945745965343) q[9];
rzz(0) q[9],q[17];
u3(pi,4.380636796165608,-4.380636796165608) q[9];
rzz(0) q[9],q[18];
u3(pi,0.9921149600036567,-0.9921149600036567) q[9];
rzz(0) q[9],q[19];
u3(pi,3.8867784310212925,-3.8867784310212925) q[9];
rzz(0) q[9],q[20];
u3(pi,0.49825659485934115,-0.49825659485934115) q[9];
rzz(0) q[9],q[21];
u3(pi,3.392920065876977,-3.392920065876977) q[9];
rzz(0) q[9],q[22];
u3(pi/2,0.12817698026646357,-0.12817698026646357) q[9];
rzz(pi/2) q[9],q[23];
rzz(pi/2) q[9],q[23];
u3(pi/2,0.12817698026646357,-0.12817698026646357) q[9];
rzz(0) q[9],q[24];
u3(pi/2,3.2697696338562565,-3.2697696338562565) q[9];
rzz(pi/2) q[9],q[25];
rzz(pi/2) q[9],q[25];
u3(pi/2,3.2697696338562565,-3.2697696338562565) q[9];
rzz(0) q[9],q[26];
u3(pi/2,0.12817698026646357,-0.12817698026646357) q[9];
rzz(pi/2) q[9],q[27];
rzz(pi/2) q[9],q[27];
u3(pi/2,0.6565928646002668,-0.6565928646002668) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,0.30284953180605606,-0.30284953180605606) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi/2,5.3689818449849565,-5.3689818449849565) q[10];
rzz(0) q[10],q[18];
u3(pi,2.5616546497371173,-2.5616546497371173) q[10];
rzz(0) q[10],q[19];
u3(pi,0.08922123136195013,-0.08922123136195013) q[10];
rzz(0) q[10],q[20];
u3(pi/2,3.565707661824415,-3.565707661824415) q[10];
rzz(pi/2) q[10],q[21];
rzz(pi/2) q[10],q[21];
rzz(-pi/2) q[10],q[22];
rzz(pi/2) q[10],q[22];
u3(pi/2,0.4241150082346221,-0.4241150082346221) q[10];
rzz(0) q[10],q[23];
u3(pi/2,3.565707661824415,-3.565707661824415) q[10];
rzz(pi/2) q[10],q[24];
rzz(pi/2) q[10],q[24];
u3(pi/2,3.565707661824415,-3.565707661824415) q[10];
rzz(0) q[10],q[25];
u3(pi/2,0.4241150082346221,-0.4241150082346221) q[10];
rzz(-pi/2) q[10],q[26];
rzz(pi/2) q[10],q[26];
u3(pi/2,3.565707661824415,-3.565707661824415) q[10];
rzz(0) q[10],q[27];
u3(pi/2,3.081274074640869,-3.081274074640869) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,3.660583759962827,-3.660583759962827) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,4.652070401435766,-4.652070401435766) q[11];
rzz(0) q[11],q[19];
u3(pi,1.0316990274388882,-1.0316990274388882) q[11];
rzz(0) q[11],q[20];
u3(pi,3.2163625587452302,-3.2163625587452302) q[11];
rzz(0) q[11],q[21];
u3(pi,5.401026090051572,-5.401026090051572) q[11];
rzz(0) q[11],q[22];
u3(pi/2,1.7806547160546946,-1.7806547160546946) q[11];
rzz(pi/2) q[11],q[23];
rzz(pi/2) q[11],q[23];
u3(pi/2,1.7806547160546946,-1.7806547160546946) q[11];
rzz(0) q[11],q[24];
u3(pi/2,4.922247369644488,-4.922247369644488) q[11];
rzz(pi/2) q[11],q[25];
rzz(-pi/2) q[11],q[25];
u3(pi/2,1.7806547160546946,-1.7806547160546946) q[11];
rzz(0) q[11],q[26];
u3(pi/2,4.922247369644488,-4.922247369644488) q[11];
rzz(pi/2) q[11],q[27];
rzz(pi/2) q[11],q[27];
u3(pi/2,4.867583657472026,-4.867583657472026) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/128) q[12],q[18];
u3(pi/2,5.179229648708133,-5.179229648708133) q[19];
rzz(0.012275082143518329) q[12],q[19];
u3(pi/2,3.2967873306771294,-3.2967873306771294) q[12];
rzz(0) q[12],q[20];
u3(pi,5.899282684910913,-5.899282684910913) q[12];
rzz(0) q[12],q[21];
u3(pi,1.6782387955476674,-1.6782387955476674) q[12];
rzz(0) q[12],q[22];
u3(pi,3.741008531894726,-3.741008531894726) q[12];
rzz(0) q[12],q[23];
u3(pi/2,0.05969026041820607,-0.05969026041820607) q[12];
rzz(pi/2) q[12],q[24];
rzz(-pi/2) q[12],q[24];
u3(pi/2,3.201282914007999,-3.201282914007999) q[12];
rzz(0) q[12],q[25];
u3(pi/2,0.05969026041820607,-0.05969026041820607) q[12];
rzz(pi/2) q[12],q[26];
rzz(pi/2) q[12],q[26];
u3(pi/2,0.05969026041820607,-0.05969026041820607) q[12];
rzz(0) q[12],q[27];
u3(pi/2,6.261822477135176,-6.261822477135176) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/128) q[13],q[19];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi/2,1.5494334967504861,-1.5494334967504861) q[13];
rzz(0) q[13],q[21];
u3(pi,5.272849109785108,-5.272849109785108) q[13];
rzz(0) q[13],q[22];
u3(pi,3.2955306936156927,-3.2955306936156927) q[13];
rzz(0) q[13],q[23];
u3(pi,1.3175839589155591,-1.3175839589155591) q[13];
rzz(0) q[13],q[24];
u3(pi/2,5.040999571950183,-5.040999571950183) q[13];
rzz(pi/2) q[13],q[25];
rzz(pi/2) q[13],q[25];
u3(pi/2,5.040999571950183,-5.040999571950183) q[13];
rzz(0) q[13],q[26];
u3(pi/2,1.899406918360389,-1.899406918360389) q[13];
rzz(-pi/2) q[13],q[27];
rzz(pi/2) q[13],q[27];
u3(pi/2,4.350477506691146,-4.350477506691146) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/128) q[14],q[20];
u3(pi/2,4.041973108108627,-4.041973108108627) q[21];
rzz(0.012275082143518329) q[14],q[21];
u3(pi/2,2.779681179896249,-2.779681179896249) q[14];
rzz(0) q[14],q[22];
u3(pi,6.256167610358714,-6.256167610358714) q[14];
rzz(0) q[14],q[23];
u3(pi,3.7837341919835468,-3.7837341919835468) q[14];
rzz(0) q[14],q[24];
u3(pi,1.3113007736083797,-1.3113007736083797) q[14];
rzz(0) q[14],q[25];
u3(pi/2,4.787158885540127,-4.787158885540127) q[14];
rzz(pi/2) q[14],q[26];
rzz(-pi/2) q[14],q[26];
u3(pi/2,1.6455662319503337,-1.6455662319503337) q[14];
rzz(0) q[14],q[27];
u3(pi/2,3.915681083434318,-3.915681083434318) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/32) q[15],q[19];
rzz(0.049100328574073315) q[15],q[20];
rzz(pi/128) q[15],q[21];
u3(pi/2,3.660583759962827,-3.660583759962827) q[22];
rzz(0.012275082143518329) q[15],q[22];
u3(pi/2,5.486477410229215,-5.486477410229215) q[15];
rzz(0) q[15],q[23];
u3(pi,2.7715130389969156,-2.7715130389969156) q[15];
rzz(0) q[15],q[24];
u3(pi,0.48380526865282814,-0.48380526865282814) q[15];
rzz(0) q[15],q[25];
u3(pi,4.479282805488327,-4.479282805488327) q[15];
rzz(0) q[15],q[26];
u3(pi/2,1.7649467527867457,-1.7649467527867457) q[15];
rzz(pi/2) q[15],q[27];
rzz(pi/2) q[15],q[27];
u3(pi/2,1.865477717701619,-1.865477717701619) q[16];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/8) q[16],q[18];
rzz(0.19640131429629326) q[16],q[19];
rzz(pi/32) q[16],q[20];
rzz(0.049100328574073315) q[16],q[21];
rzz(pi/128) q[16],q[22];
u3(pi/2,2.9700616947037903,-2.9700616947037903) q[23];
rzz(0.012275082143518329) q[16],q[23];
u3(pi/2,0.2946813909067226,-0.2946813909067226) q[16];
rzz(0) q[16],q[24];
u3(pi,4.421477500662275,-4.421477500662275) q[16];
rzz(0) q[16],q[25];
u3(pi,3.2509200779347176,-3.2509200779347176) q[16];
rzz(0) q[16],q[26];
u3(pi,2.0797343366764434,-2.0797343366764434) q[16];
rzz(0) q[16],q[27];
u3(pi/2,1.7580352489488482,-1.7580352489488482) q[17];
rzz(0.7853830994606742) q[17],q[18];
rzz(pi/8) q[17],q[19];
rzz(0.19640131429629326) q[17],q[20];
rzz(pi/32) q[17],q[21];
rzz(0.049100328574073315) q[17],q[22];
rzz(pi/128) q[17],q[23];
u3(pi/2,3.660583759962827,-3.660583759962827) q[24];
rzz(0.012275082143518329) q[17],q[24];
u3(pi/2,3.328831575743745,-3.328831575743745) q[17];
rzz(0) q[17],q[25];
u3(pi,1.247212283475148,-1.247212283475148) q[17];
rzz(0) q[17],q[26];
u3(pi,0.22556635252774715,-0.22556635252774715) q[17];
rzz(0) q[17],q[27];
u3(pi/2,3.600265181013903,-3.600265181013903) q[18];
rzz(0.7853830994606742) q[18],q[19];
rzz(pi/8) q[18],q[20];
rzz(0.19640131429629326) q[18],q[21];
rzz(pi/32) q[18],q[22];
rzz(0.049100328574073315) q[18],q[23];
rzz(pi/128) q[18],q[24];
u3(pi/2,0.34683182895631315,-0.34683182895631315) q[25];
rzz(0.012275082143518329) q[18],q[25];
u3(pi/2,2.0294688542190062,-2.0294688542190062) q[18];
rzz(0) q[18],q[26];
u3(pi,5.197450886098954,-5.197450886098954) q[18];
rzz(0) q[18],q[27];
u3(pi/2,4.3605306031826325,-4.3605306031826325) q[19];
rzz(0.7853830994606742) q[19],q[20];
rzz(pi/8) q[19],q[21];
rzz(0.19640131429629326) q[19],q[22];
rzz(pi/32) q[19],q[23];
rzz(0.049100328574073315) q[19],q[24];
rzz(pi/128) q[19],q[25];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[26];
rzz(0.012275082143518329) q[19],q[26];
u3(pi/2,5.931326929977529,-5.931326929977529) q[19];
rzz(0) q[19],q[27];
u3(pi/2,2.463008640414398,-2.463008640414398) q[20];
rzz(0.7853830994606742) q[20],q[21];
rzz(pi/8) q[20],q[22];
rzz(0.19640131429629326) q[20],q[23];
rzz(pi/32) q[20],q[24];
rzz(0.049100328574073315) q[20],q[25];
rzz(pi/128) q[20],q[26];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[27];
rzz(0.012275082143518329) q[20],q[27];
u3(pi/2,1.0838494654884785,-1.0838494654884785) q[21];
rzz(0.7853830994606742) q[21],q[22];
rzz(pi/8) q[21],q[23];
rzz(0.19640131429629326) q[21],q[24];
rzz(pi/32) q[21],q[25];
rzz(0.049100328574073315) q[21],q[26];
rzz(pi/128) q[21],q[27];
u3(pi/2,1.3936105011324322,-1.3936105011324322) q[22];
rzz(0.7853830994606742) q[22],q[23];
rzz(pi/8) q[22],q[24];
rzz(0.19640131429629326) q[22],q[25];
rzz(pi/32) q[22],q[26];
rzz(0.049100328574073315) q[22],q[27];
u3(pi/2,1.048035309237555,-1.048035309237555) q[23];
rzz(0.7853830994606742) q[23],q[24];
rzz(pi/8) q[23],q[25];
rzz(0.19640131429629326) q[23],q[26];
rzz(pi/32) q[23],q[27];
u3(pi/2,5.052937624033824,-5.052937624033824) q[24];
rzz(0.7853830994606742) q[24],q[25];
rzz(pi/8) q[24],q[26];
rzz(0.19640131429629326) q[24],q[27];
u3(pi/2,0.25446900494077324,-0.25446900494077324) q[25];
rzz(0.7853830994606742) q[25],q[26];
rzz(pi/8) q[25],q[27];
u3(pi/2,5.967769404759171,-5.967769404759171) q[26];
rzz(0.7853830994606742) q[26],q[27];
u3(pi/2,0.10430087609918114,-0.10430087609918114) q[1];
u3(pi/2,1.3150706847926874,-1.3150706847926874) q[2];
u3(pi/2,5.37714998588429,-5.37714998588429) q[3];
u3(pi/2,4.1136014206104745,-4.1136014206104745) q[4];
u3(pi/2,0.27269024233159406,-0.27269024233159406) q[5];
u3(pi/2,4.5301766064764815,-4.5301766064764815) q[6];
u3(pi/2,3.8421678153403174,-3.8421678153403174) q[8];
u3(pi/2,0.4241150082346221,-0.4241150082346221) q[10];
u3(pi/2,3.201282914007999,-3.201282914007999) q[12];
u3(pi/2,4.787158885540127,-4.787158885540127) q[14];
u3(pi/2,6.207158764962713,-6.207158764962713) q[16];
u3(pi/2,4.427132367438737,-4.427132367438737) q[17];
u3(pi/2,2.082875929330033,-2.082875929330033) q[18];
u3(pi/2,2.7897342763877364,-2.7897342763877364) q[19];
u3(pi/2,2.4548404995150643,-2.4548404995150643) q[27];
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
measure q[22] -> c[22];
measure q[23] -> c[23];
measure q[24] -> c[24];
measure q[25] -> c[25];
measure q[26] -> c[26];
measure q[27] -> c[27];
measure q[28] -> c[28];
