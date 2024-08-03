OPENQASM 2.0;
include "qelib1.inc";
qreg q[18];
creg c[18];
u(pi/2,4.0809288570131415,-4.0809288570131415) q[17];
u(pi/2,1.7247343668207962,-1.7247343668207962) q[17];
rzz(-pi/2) q[16],q[17];
u(pi/2,3.295530693615693,-3.295530693615693) q[17];
rzz(pi/8) q[15],q[17];
u(pi/2,3.7617430434084183,-3.7617430434084183) q[16];
rzz(0.7853830994606742) q[15],q[16];
u(pi,1.7479821524573609,-1.7479821524573609) q[17];
rzz(0.19640131429629326) q[14],q[17];
u(pi,1.8013892275683872,-1.8013892275683872) q[16];
rzz(pi/8) q[14],q[16];
u(pi/2,-0.9556724852220151,0.9556724852220151) q[15];
rzz(0.7853830994606742) q[14],q[15];
u(pi,2.042663543364083,-2.042663543364083) q[17];
rzz(pi/32) q[13],q[17];
u(pi,0.9311680625240149,-0.9311680625240149) q[16];
rzz(0.19640131429629326) q[13],q[16];
u(pi,3.673778449107904,-3.673778449107904) q[15];
rzz(pi/8) q[13],q[15];
u(pi/2,4.052026204600115,-4.052026204600115) q[14];
rzz(0.7853830994606742) q[13],q[14];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/32) q[12],q[16];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/8) q[12],q[14];
u(pi/2,0.8896990394966293,-0.8896990394966293) q[13];
rzz(0.7853830994606742) q[12],q[13];
u(pi,0.9820618635121692,-0.9820618635121692) q[17];
rzz(pi/128) q[11],q[17];
u(pi,0.5692565888304708,-0.5692565888304708) q[16];
rzz(0.049100328574073315) q[11],q[16];
u(pi,3.6913713679680074,-3.6913713679680074) q[15];
rzz(pi/32) q[11],q[15];
u(pi,3.1447342462433827,-3.1447342462433827) q[14];
rzz(0.19640131429629326) q[11],q[14];
u(pi,3.5512563356179028,-3.5512563356179028) q[13];
rzz(pi/8) q[11],q[13];
u(pi/2,1.261035291150943,-1.261035291150943) q[12];
rzz(0.7853830994606742) q[11],q[12];
u(pi,4.172035043967245,-4.172035043967245) q[17];
rzz(0.012275082143518329) q[10],q[17];
u(pi,1.1316016738230434,-1.1316016738230434) q[16];
rzz(pi/128) q[10],q[16];
u(pi,1.2955928103404308,-1.2955928103404308) q[15];
rzz(0.049100328574073315) q[10],q[15];
u(pi,3.837769585625291,-3.837769585625291) q[14];
rzz(pi/32) q[10],q[14];
u(pi,0.09927432785343737,-0.09927432785343737) q[13];
rzz(0.19640131429629326) q[10],q[13];
u(pi,3.634822700203391,-3.634822700203391) q[12];
rzz(pi/8) q[10],q[12];
u(pi/2,1.9879998311916212,-1.9879998311916212) q[11];
rzz(0.7853830994606742) q[10],q[11];
u(pi/2,0.09801769079200162,-0.09801769079200162) q[9];
u(pi/2,0.8394335570391926,-0.8394335570391926) q[17];
rzz(0) q[9],q[17];
u(pi/2,3.239610344381794,-3.239610344381794) q[9];
rzz(0.012275082143518329) q[9],q[16];
rzz(pi/128) q[9],q[15];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/32) q[9],q[13];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/8) q[9],q[11];
u(pi/2,3.5713625286008766,-3.5713625286008766) q[10];
rzz(0.7853830994606742) q[9],q[10];
u(pi/2,0.1470265361880021,-0.1470265361880021) q[8];
rzz(0) q[8],q[17];
u(pi,2.020672394788955,-2.020672394788955) q[8];
u(pi/2,0.003769911184307695,-0.003769911184307695) q[16];
rzz(0) q[8],q[16];
u(pi/2,3.8936899348591893,-3.8936899348591893) q[8];
rzz(0.012275082143518329) q[8],q[15];
rzz(pi/128) q[8],q[14];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/32) q[8],q[12];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/8) q[8],q[10];
u(pi/2,0.09801769079200162,-0.09801769079200162) q[9];
rzz(0.7853830994606742) q[8],q[9];
u(pi/2,-3*pi/8,3*pi/8) q[7];
rzz(0) q[7],q[17];
u(pi,3.003362576831842,-3.003362576831842) q[7];
rzz(0) q[7],q[16];
u(pi,1.940875941387774,-1.940875941387774) q[7];
u(pi/2,3.5123005867133887,-3.5123005867133887) q[15];
rzz(0) q[7],q[15];
u(pi/2,-0.16084954386379757,0.16084954386379757) q[7];
u(pi,0.5265309287416491,-0.5265309287416491) q[14];
rzz(0.012275082143518329) q[7],q[14];
u(pi,2.1312564561953153,-2.1312564561953153) q[13];
rzz(pi/128) q[7],q[13];
u(pi,-1.1806105192190444,1.1806105192190444) q[12];
rzz(0.049100328574073315) q[7],q[12];
u(pi,-0.06346017160251383,0.06346017160251383) q[11];
rzz(pi/32) q[7],q[11];
u(pi,2.016902483604647,-2.016902483604647) q[10];
rzz(0.19640131429629326) q[7],q[10];
u(pi,3.8955748904513436,-3.8955748904513436) q[9];
rzz(pi/8) q[7],q[9];
u(pi/2,3.8936899348591893,-3.8936899348591893) q[8];
rzz(0.7853830994606742) q[7],q[8];
u(pi/2,-0.9814335449814514,0.9814335449814514) q[6];
rzz(0) q[6],q[17];
u(pi,1.3420883816135598,-1.3420883816135598) q[6];
rzz(0) q[6],q[16];
u(pi,2.8475395812137885,-2.8475395812137885) q[6];
rzz(0) q[6],q[15];
u(pi,4.353619099344735,-4.353619099344735) q[6];
u(pi/2,1.8975219627682351,-1.8975219627682351) q[14];
rzz(0) q[6],q[14];
u(pi/2,0.3939557187601599,-0.3939557187601599) q[6];
u(pi,3.80446870349724,-3.80446870349724) q[13];
rzz(0.012275082143518329) q[6],q[13];
u(pi,0.23750440461138833,-0.23750440461138833) q[12];
rzz(pi/128) q[6],q[12];
u(pi,-1.0247875236009905,1.0247875236009905) q[11];
rzz(0.049100328574073315) q[6],q[11];
u(pi,2.118690085580957,-2.118690085580957) q[10];
rzz(pi/32) q[6],q[10];
u(pi,2.065911329000648,-2.065911329000648) q[9];
rzz(0.19640131429629326) q[6],q[9];
u(pi,2.4586104106993716,-2.4586104106993716) q[8];
rzz(pi/8) q[6],q[8];
u(pi/2,2.9807431097259958,-2.9807431097259958) q[7];
rzz(0.7853830994606742) q[6],q[7];
u(pi/2,0,0) q[5];
rzz(0) q[5],q[17];
u(pi,2.823035158515788,-2.823035158515788) q[5];
rzz(0) q[5],q[16];
u(pi,-0.9550441666912971,0.9550441666912971) q[5];
rzz(0) q[5],q[15];
u(pi,1.5500618152812038,-1.5500618152812038) q[5];
rzz(0) q[5],q[14];
u(pi,4.055167797253705,-4.055167797253705) q[5];
u(pi/2,0.47186721656918706,-0.47186721656918706) q[13];
rzz(0) q[5],q[13];
u(pi/2,0.595017648589907,-0.595017648589907) q[5];
u(pi,1.9427608969799284,-1.9427608969799284) q[12];
rzz(0.012275082143518329) q[5],q[12];
u(pi,4.381893433227043,-4.381893433227043) q[11];
rzz(pi/128) q[5],q[11];
u(pi,2.5340086343855273,-2.5340086343855273) q[10];
rzz(0.049100328574073315) q[5],q[10];
u(pi,-0.6107256118578558,0.6107256118578558) q[9];
rzz(pi/32) q[5],q[9];
u(pi,2.723760830662351,-2.723760830662351) q[8];
rzz(0.19640131429629326) q[5],q[8];
u(pi,-0.6848671984825749,0.6848671984825749) q[7];
rzz(pi/8) q[5],q[7];
u(pi/2,0.3939557187601599,-0.3939557187601599) q[6];
rzz(0.7853830994606742) q[5],q[6];
u(pi/2,pi/4,-pi/4) q[4];
rzz(0) q[4],q[17];
u(pi,4.033176648678577,-4.033176648678577) q[4];
rzz(0) q[4],q[16];
u(pi,1.1045839770021715,-1.1045839770021715) q[4];
rzz(0) q[4],q[15];
u(pi,4.458548293974634,-4.458548293974634) q[4];
rzz(0) q[4],q[14];
u(pi,1.5299556222982291,-1.5299556222982291) q[4];
rzz(0) q[4],q[13];
u(pi,-1.3992653679088938,1.3992653679088938) q[4];
u(pi/2,4.46608811634325,-4.46608811634325) q[12];
rzz(0) q[4],q[12];
u(pi/2,1.8485131173722342,-1.8485131173722342) q[4];
rzz(0.012275082143518329) q[4],q[11];
rzz(pi/128) q[4],q[10];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/32) q[4],q[8];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/8) q[4],q[6];
u(pi/2,0.595017648589907,-0.595017648589907) q[5];
rzz(0.7853830994606742) q[4],q[5];
u(pi/2,pi/4,-pi/4) q[3];
rzz(0) q[3],q[17];
u(pi,4.033176648678577,-4.033176648678577) q[3];
rzz(0) q[3],q[16];
u(pi,1.1045839770021715,-1.1045839770021715) q[3];
rzz(0) q[3],q[15];
u(pi,4.458548293974634,-4.458548293974634) q[3];
rzz(0) q[3],q[14];
u(pi,1.5299556222982291,-1.5299556222982291) q[3];
rzz(0) q[3],q[13];
u(pi,-1.3992653679088938,1.3992653679088938) q[3];
rzz(0) q[3],q[12];
u(pi,1.9546989490635696,-1.9546989490635696) q[3];
u(pi/2,2.414628113549115,-2.414628113549115) q[11];
rzz(0) q[3],q[11];
u(pi/2,-1.0800795543041708,1.0800795543041708) q[3];
rzz(0.012275082143518329) q[3],q[10];
rzz(pi/128) q[3],q[9];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/32) q[3],q[7];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/8) q[3],q[5];
u(pi/2,1.8485131173722342,-1.8485131173722342) q[4];
rzz(0.7853830994606742) q[3],q[4];
u(pi/2,pi/2,-pi/2) q[2];
rzz(0) q[2],q[17];
u(pi,-1.4017786420317657,1.4017786420317657) q[2];
rzz(0) q[2],q[16];
u(pi,2.077849381084289,-2.077849381084289) q[2];
rzz(0) q[2],q[15];
u(pi,-0.7257079029792421,0.7257079029792421) q[2];
rzz(0) q[2],q[14];
u(pi,2.7539201201368124,-2.7539201201368124) q[2];
rzz(0) q[2],q[13];
u(pi,-0.049637163926718575,0.049637163926718575) q[2];
rzz(0) q[2],q[12];
u(pi,3.4299908591893367,-3.4299908591893367) q[2];
rzz(0) q[2],q[11];
u(pi,0.6264335751258048,-0.6264335751258048) q[2];
u(pi/2,1.293079536217559,-1.293079536217559) q[10];
rzz(0) q[2],q[10];
u(pi/2,3.937043913478729,-3.937043913478729) q[2];
u(pi,2.021929031850391,-2.021929031850391) q[9];
rzz(0.012275082143518329) q[2],q[9];
u(pi,3.184946632209332,-3.184946632209332) q[8];
rzz(pi/128) q[2],q[8];
u(pi,0.8645662982679112,-0.8645662982679112) q[7];
rzz(0.049100328574073315) q[2],q[7];
u(pi,-0.8469733794078083,0.8469733794078083) q[6];
rzz(pi/32) q[2],q[6];
u(pi,3.4186811256364127,-3.4186811256364127) q[5];
rzz(0.19640131429629326) q[2],q[5];
u(pi,0.44107960856400696,-0.44107960856400696) q[4];
rzz(pi/8) q[2],q[4];
u(pi/2,-1.0800795543041708,1.0800795543041708) q[3];
rzz(0.7853830994606742) q[2],q[3];
u(pi/2,-pi/2,pi/2) q[1];
rzz(0) q[1],q[17];
u(pi,2.2267608728644457,-2.2267608728644457) q[1];
rzz(0) q[1],q[16];
u(pi,0.39646899288303183,-0.39646899288303183) q[1];
rzz(0) q[1],q[15];
u(pi,-1.4331945685676637,1.4331945685676637) q[1];
rzz(0) q[1],q[14];
u(pi,3.020327177161227,-3.020327177161227) q[1];
rzz(0) q[1],q[13];
u(pi,1.1900352971798136,-1.1900352971798136) q[1];
rzz(0) q[1],q[12];
u(pi,-0.6396282642708819,0.6396282642708819) q[1];
rzz(0) q[1],q[11];
u(pi,3.8132651629272907,-3.8132651629272907) q[1];
rzz(0) q[1],q[10];
u(pi,1.9836016014765954,-1.9836016014765954) q[1];
u(pi/2,-1.4375927982826893,1.4375927982826893) q[9];
rzz(0) q[1],q[9];
u(pi/2,-0.5020265060436488,0.5020265060436488) q[1];
u(pi,3.2528050335268723,-3.2528050335268723) q[8];
rzz(0.012275082143518329) q[1],q[8];
u(pi,2.19157503514424,-2.19157503514424) q[7];
rzz(pi/128) q[1],q[7];
u(pi,-0.7809999336824225,0.7809999336824225) q[6];
rzz(0.049100328574073315) q[1],q[6];
u(pi,-1.5004246513544852,1.5004246513544852) q[5];
rzz(pi/32) q[1],q[5];
u(pi,0.7338760438785759,-0.7338760438785759) q[4];
rzz(0.19640131429629326) q[1],q[4];
u(pi,3.7956722440671875,-3.7956722440671875) q[3];
rzz(pi/8) q[1],q[3];
u(pi/2,3.937043913478729,-3.937043913478729) q[2];
rzz(0.7853830994606742) q[1],q[2];
rzz(pi/2) q[0],q[17];
rzz(pi/2) q[0],q[16];
rzz(-pi/2) q[0],q[15];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[13];
rzz(-pi/2) q[0],q[12];
rzz(pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[10];
rzz(pi/2) q[0],q[9];
u(pi/2,1.4174866052997146,-1.4174866052997146) q[8];
rzz(-pi/2) q[0],q[8];
u(pi/2,4.586725274241098,-4.586725274241098) q[7];
rzz(pi/2) q[0],q[7];
u(pi/2,3.666238626739289,-3.666238626739289) q[6];
rzz(pi/2) q[0],q[6];
u(pi/2,0.18221237390820821,-0.18221237390820821) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-0.7081149841191393,0.7081149841191393) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,2.3882387352589607,-2.3882387352589607) q[3];
rzz(-pi/2) q[0],q[3];
u(pi/2,0.7954512598889356,-0.7954512598889356) q[2];
rzz(pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u(pi,-pi/2,pi/2) q[0];
u(pi/2,2.6395661475461445,-2.6395661475461445) q[1];
rzz(-pi/2) q[0],q[1];
u(pi/2,2.3662475866838326,-2.3662475866838326) q[2];
u(pi/2,0.010053096491487334,-0.010053096491487334) q[2];
rzz(pi/2) q[0],q[2];
u(pi/2,0.8174424084640641,-0.8174424084640641) q[3];
u(pi/2,4.351734143752581,-4.351734143752581) q[3];
rzz(pi/2) q[0],q[3];
u(pi/2,0.8626813426757569,-0.8626813426757569) q[4];
u(pi/2,4.200309377849553,-4.200309377849553) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-1.3885839528866886,1.3885839528866886) q[5];
u(pi/2,1.8510263914951057,-1.8510263914951057) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-1.046150353645401,1.046150353645401) q[6];
u(pi/2,2.144451145340393,-2.144451145340393) q[6];
rzz(-pi/2) q[0],q[6];
u(pi/2,-pi/25,pi/25) q[7];
u(pi/2,3.0404333701442017,-3.0404333701442017) q[7];
rzz(pi/2) q[0],q[7];
u(pi/2,-0.15330972149518196,0.15330972149518196) q[8];
u(pi/2,3.0002209841782523,-3.0002209841782523) q[8];
rzz(pi/2) q[0],q[8];
rzz(pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[11];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[15];
rzz(-pi/2) q[0],q[16];
rzz(-pi/2) q[0],q[17];
u(pi/2,1.0687698207512475,-1.0687698207512475) q[1];
u(pi/2,-1.5607432303034092,1.5607432303034092) q[2];
rzz(0.7853830994606742) q[1],q[2];
u(pi/2,2.7809378169576844,-2.7809378169576844) q[3];
rzz(pi/8) q[1],q[3];
u(pi/2,2.6295130510546567,-2.6295130510546567) q[4];
rzz(0.19640131429629326) q[1],q[4];
u(pi/2,0.2802300647002094,-0.2802300647002094) q[5];
rzz(pi/32) q[1],q[5];
u(pi/2,3.7152474721352897,-3.7152474721352897) q[6];
rzz(0.049100328574073315) q[1],q[6];
u(pi/2,1.4696370433493051,-1.4696370433493051) q[7];
rzz(pi/128) q[1],q[7];
u(pi/2,1.4294246573833558,-1.4294246573833558) q[8];
rzz(0.012275082143518329) q[1],q[8];
u(pi/2,1.0687698207512475,-1.0687698207512475) q[1];
rzz(0) q[1],q[9];
u(pi,-1.4174866052997146,1.4174866052997146) q[1];
rzz(0) q[1],q[10];
u(pi,3.036035140429176,-3.036035140429176) q[1];
rzz(0) q[1],q[11];
u(pi,1.2063715789784806,-1.2063715789784806) q[1];
rzz(0) q[1],q[12];
u(pi,-0.6239203010029329,0.6239203010029329) q[1];
rzz(0) q[1],q[13];
u(pi,3.8296014447259576,-3.8296014447259576) q[1];
rzz(0) q[1],q[14];
u(pi,1.9993095647445447,-1.9993095647445447) q[1];
rzz(0) q[1],q[15];
u(pi,0.16964600329384893,-0.16964600329384893) q[1];
rzz(0) q[1],q[16];
u(pi,4.623167749022739,-4.623167749022739) q[1];
rzz(0) q[1],q[17];
u(pi/2,0.7954512598889356,-0.7954512598889356) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u(pi/2,-1.4501591688970485,1.4501591688970485) q[9];
rzz(0.012275082143518329) q[2],q[9];
u(pi/2,-0.7753450669059611,0.7753450669059611) q[2];
rzz(0) q[2],q[10];
u(pi,2.047690091609827,-2.047690091609827) q[2];
rzz(0) q[2],q[11];
u(pi,4.5527960735823285,-4.5527960735823285) q[2];
rzz(0) q[2],q[12];
u(pi,0.7747167483752433,-0.7747167483752433) q[2];
rzz(0) q[2],q[13];
u(pi,3.2798227303477443,-3.2798227303477443) q[2];
rzz(0) q[2],q[14];
u(pi,-0.498884913390059,0.498884913390059) q[2];
rzz(0) q[2],q[15];
u(pi,2.006221068582442,-2.006221068582442) q[2];
rzz(0) q[2],q[16];
u(pi,4.511327050554943,-4.511327050554943) q[2];
rzz(0) q[2],q[17];
u(pi/2,2.3882387352589607,-2.3882387352589607) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u(pi/2,4.419592545070121,-4.419592545070121) q[10];
rzz(0.012275082143518329) q[3],q[10];
u(pi/2,3.9590350620538572,-3.9590350620538572) q[3];
rzz(0) q[3],q[11];
u(pi,-0.0006283185307178751,0.0006283185307178751) q[3];
rzz(0) q[3],q[12];
u(pi,1.5054511996002287,-1.5054511996002287) q[3];
rzz(0) q[3],q[13];
u(pi,3.0109023992004573,-3.0109023992004573) q[3];
rzz(0) q[3],q[14];
u(pi,4.516981917331404,-4.516981917331404) q[3];
rzz(0) q[3],q[15];
u(pi,-0.2607521902479528,0.2607521902479528) q[3];
rzz(0) q[3],q[16];
u(pi,1.2453273278829937,-1.2453273278829937) q[3];
rzz(0) q[3],q[17];
u(pi/2,2.432849350939936,-2.432849350939936) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u(pi/2,2.40080510587332,-2.40080510587332) q[11];
rzz(0.012275082143518329) q[4],q[11];
u(pi/2,0.8620530241450393,-0.8620530241450393) q[4];
rzz(0) q[4],q[12];
u(pi,2.785336046672711,-2.785336046672711) q[4];
rzz(0) q[4],q[13];
u(pi,3.4896811196075426,-3.4896811196075426) q[4];
rzz(0) q[4],q[14];
u(pi,4.194026192542374,-4.194026192542374) q[4];
rzz(0) q[4],q[15];
u(pi,-1.3841857231716628,1.3841857231716628) q[4];
rzz(0) q[4],q[16];
u(pi,-0.6798406502368312,0.6798406502368312) q[4];
rzz(0) q[4],q[17];
u(pi/2,0.18221237390820821,-0.18221237390820821) q[5];
rzz(0.7853830994606742) q[5],q[6];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u(pi/2,4.448495197483147,-4.448495197483147) q[12];
rzz(0.012275082143518329) q[5],q[12];
u(pi/2,1.7530087007031048,-1.7530087007031048) q[5];
rzz(0) q[5],q[13];
u(pi,3.6970262347444685,-3.6970262347444685) q[5];
rzz(0) q[5],q[14];
u(pi,4.4428403307066855,-4.4428403307066855) q[5];
rzz(0) q[5],q[15];
u(pi,-1.093902561979966,1.093902561979966) q[5];
rzz(0) q[5],q[16];
u(pi,-0.34746014748703113,0.34746014748703113) q[5];
rzz(0) q[5],q[17];
u(pi/2,0.5246459731494957,-0.5246459731494957) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u(pi/2,0.454902616239802,-0.454902616239802) q[13];
rzz(0.012275082143518329) q[6],q[13];
u(pi/2,-1.046150353645401,1.046150353645401) q[6];
rzz(0) q[6],q[14];
u(pi,1.074424687527709,-1.074424687527709) q[6];
rzz(0) q[6],q[15];
u(pi,2.173353797753419,-2.173353797753419) q[6];
rzz(0) q[6],q[16];
u(pi,3.272282907979129,-3.272282907979129) q[6];
rzz(0) q[6],q[17];
u(pi/2,1.4451326206513047,-1.4451326206513047) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u(pi/2,-1.2622919282123788,1.2622919282123788) q[14];
rzz(0.012275082143518329) q[7],q[14];
u(pi/2,3.015928947446201,-3.015928947446201) q[7];
rzz(0) q[7],q[15];
u(pi,-0.4435928826868787,0.4435928826868787) q[7];
rzz(0) q[7],q[16];
u(pi,2.061513099285622,-2.061513099285622) q[7];
rzz(0) q[7],q[17];
u(pi/2,4.55845094035879,-4.55845094035879) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u(pi/2,0.3524866957327746,-0.3524866957327746) q[15];
rzz(0.012275082143518329) q[8],q[15];
u(pi/2,2.987654613563893,-2.987654613563893) q[8];
rzz(0) q[8],q[16];
u(pi,-0.9902300044115028,0.9902300044115028) q[8];
rzz(0) q[8],q[17];
u(pi/2,1.6794954326091034,-1.6794954326091034) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u(pi/2,3.1271413273832804,-3.1271413273832804) q[16];
rzz(0.012275082143518329) q[9],q[16];
u(pi/2,3.250291759404,-3.250291759404) q[9];
rzz(0) q[9],q[17];
u(pi/2,1.2685751135195584,-1.2685751135195584) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u(pi/2,3.962804973238165,-3.962804973238165) q[17];
rzz(0.012275082143518329) q[10],q[17];
u(pi/2,2.3932652835047046,-2.3932652835047046) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u(pi/2,1.2999910400554562,-1.2999910400554562) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
u(pi/2,3.590212084522416,-3.590212084522416) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
u(pi/2,-1.2685751135195584,1.2685751135195584) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
u(pi/2,3.4877961640153883,-3.4877961640153883) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
u(pi/2,3.1208581420761003,-3.1208581420761003) q[16];
rzz(0.7853830994606742) q[16],q[17];
u(pi,0,0) q[0];
u(pi/2,2.136911322971777,-2.136911322971777) q[1];
u(pi/2,1.0511769018911448,-1.0511769018911448) q[2];
u(pi/2,3.568849254478005,-3.568849254478005) q[3];
u(pi/2,1.2428140537601222,-1.2428140537601222) q[4];
u(pi/2,1.5965573865543328,-1.5965573865543328) q[5];
u(pi/2,-0.8909556765580653,0.8909556765580653) q[6];
u(pi/2,-1.398637049378176,1.398637049378176) q[7];
u(pi/2,1.3144423662619693,-1.3144423662619693) q[8];
u(pi/2,0.10869910581420683,-0.10869910581420683) q[9];
u(pi/2,3.9565217879309857,-3.9565217879309857) q[17];
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
