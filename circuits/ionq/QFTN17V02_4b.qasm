OPENQASM 2.0;
include "qelib1.inc";
qreg q[17];
creg c[17];
u(pi/2,3.3948050214691303,-3.3948050214691303) q[16];
u(pi/2,1.0386105312767855,-1.0386105312767855) q[16];
rzz(pi/2) q[15],q[16];
u(pi/2,-0.5321857955181108,0.5321857955181108) q[16];
rzz(pi/8) q[14],q[16];
u(pi/2,-0.9311680625240146,0.9311680625240146) q[15];
rzz(0.7853830994606742) q[14],q[15];
u(pi,4.271309371820682,-4.271309371820682) q[16];
rzz(0.19640131429629326) q[13],q[16];
u(pi,2.1620440642004954,-2.1620440642004954) q[15];
rzz(pi/8) q[13],q[15];
u(pi/2,2.5836457983122463,-2.5836457983122463) q[14];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/32) q[12],q[16];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/8) q[12],q[14];
u(pi/2,0.9883450488193488,-0.9883450488193488) q[13];
rzz(0.7853830994606742) q[12],q[13];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/32) q[11],q[15];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/8) q[11],q[13];
u(pi/2,-0.6609910943152923,0.6609910943152923) q[12];
rzz(0.7853830994606742) q[11],q[12];
u(pi,4.593636778078995,-4.593636778078995) q[16];
rzz(pi/128) q[10],q[16];
u(pi,4.456663338382481,-4.456663338382481) q[15];
rzz(0.049100328574073315) q[10],q[15];
u(pi,1.1485662741524285,-1.1485662741524285) q[14];
rzz(pi/32) q[10],q[14];
u(pi,3.3841236064469253,-3.3841236064469253) q[13];
rzz(0.19640131429629326) q[10],q[13];
u(pi,4.3222031728088375,-4.3222031728088375) q[12];
rzz(pi/8) q[10],q[12];
u(pi/2,-0.7577521480458582,0.7577521480458582) q[11];
rzz(0.7853830994606742) q[10],q[11];
rzz(0.012275082143518329) q[9],q[16];
rzz(pi/128) q[9],q[15];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/32) q[9],q[13];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/8) q[9],q[11];
u(pi/2,2.0677962845928013,-2.0677962845928013) q[10];
rzz(0.7853830994606742) q[9],q[10];
u(pi/2,2.822406839985071,-2.822406839985071) q[8];
u(pi/2,3.255318307649744,-3.255318307649744) q[16];
rzz(0) q[8],q[16];
u(pi/2,-0.319185813604723,0.319185813604723) q[8];
rzz(0.012275082143518329) q[8],q[15];
rzz(pi/128) q[8],q[14];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/32) q[8],q[12];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/8) q[8],q[10];
u(pi/2,-1.4602122653885359,1.4602122653885359) q[9];
rzz(0.7853830994606742) q[8],q[9];
u(pi/2,-1.1290883997001717,1.1290883997001717) q[7];
rzz(0) q[7],q[16];
u(pi,0.7464424144929351,-0.7464424144929351) q[7];
u(pi/2,0.5164778322501618,-0.5164778322501618) q[15];
rzz(0) q[7],q[15];
u(pi/2,2.6219732286860413,-2.6219732286860413) q[7];
rzz(0.012275082143518329) q[7],q[14];
rzz(pi/128) q[7],q[13];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/32) q[7],q[11];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/8) q[7],q[9];
u(pi/2,-0.319185813604723,0.319185813604723) q[8];
rzz(0.7853830994606742) q[7],q[8];
u(pi/2,0.09801769079200162,-0.09801769079200162) q[6];
rzz(0) q[6],q[16];
u(pi,2.9782298356031234,-2.9782298356031234) q[6];
rzz(0) q[6],q[15];
u(pi,-0.6867521540747288,0.6867521540747288) q[6];
u(pi/2,-0.28588493147667116,0.28588493147667116) q[14];
rzz(0) q[6],q[14];
u(pi/2,2.1928316722056755,-2.1928316722056755) q[6];
rzz(0.012275082143518329) q[6],q[13];
rzz(pi/128) q[6],q[12];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/32) q[6],q[10];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/8) q[6],q[8];
u(pi/2,-0.5196194249037518,0.5196194249037518) q[7];
rzz(0.7853830994606742) q[6],q[7];
u(pi/2,0.19666370011472112,-0.19666370011472112) q[5];
rzz(0) q[5],q[16];
u(pi,2.1400529156253674,-2.1400529156253674) q[5];
rzz(0) q[5],q[15];
u(pi,2.8864953301183025,-2.8864953301183025) q[5];
rzz(0) q[5],q[14];
u(pi,3.6329377446112368,-3.6329377446112368) q[5];
u(pi/2,-0.5032831431050848,0.5032831431050848) q[13];
rzz(0) q[5],q[13];
u(pi/2,-0.7068583470577033,0.7068583470577033) q[5];
u(pi,-1.5255573925832036,1.5255573925832036) q[12];
rzz(0.012275082143518329) q[5],q[12];
u(pi,3.0398050516134836,-3.0398050516134836) q[11];
rzz(pi/128) q[5],q[11];
u(pi,0.5202477434344694,-0.5202477434344694) q[10];
rzz(0.049100328574073315) q[5],q[10];
u(pi,0.6911503837897546,-0.6911503837897546) q[9];
rzz(pi/32) q[5],q[9];
u(pi,1.6248317204366414,-1.6248317204366414) q[8];
rzz(0.19640131429629326) q[5],q[8];
u(pi,2.0037077944595696,-2.0037077944595696) q[7];
rzz(pi/8) q[5],q[7];
u(pi/2,2.1928316722056755,-2.1928316722056755) q[6];
rzz(0.7853830994606742) q[5],q[6];
u(pi/2,-3*pi/8,3*pi/8) q[4];
rzz(0) q[4],q[16];
u(pi,0.7451857774314989,-0.7451857774314989) q[4];
rzz(0) q[4],q[15];
u(pi,1.4495308503663304,-1.4495308503663304) q[4];
rzz(0) q[4],q[14];
u(pi,2.153875923301162,-2.153875923301162) q[4];
rzz(0) q[4],q[13];
u(pi,2.8582209962359935,-2.8582209962359935) q[4];
u(pi/2,3.3514510428495914,-3.3514510428495914) q[12];
rzz(0) q[4],q[12];
u(pi/2,-1.5016812884159212,1.5016812884159212) q[4];
rzz(0.012275082143518329) q[4],q[11];
rzz(pi/128) q[4],q[10];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/32) q[4],q[8];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/8) q[4],q[6];
u(pi/2,-0.7068583470577033,0.7068583470577033) q[5];
rzz(0.7853830994606742) q[4],q[5];
u(pi/2,-pi/4,pi/4) q[3];
rzz(0) q[3],q[16];
u(pi,1.5381237631975626,-1.5381237631975626) q[3];
rzz(0) q[3],q[15];
u(pi,3.0442032813285094,-3.0442032813285094) q[3];
rzz(0) q[3],q[14];
u(pi,4.549654480928738,-4.549654480928738) q[3];
rzz(0) q[3],q[13];
u(pi,-0.227451308119901,0.227451308119901) q[3];
rzz(0) q[3],q[12];
u(pi,1.277999891480328,-1.277999891480328) q[3];
u(pi/2,0.5541769440932396,-0.5541769440932396) q[11];
rzz(0) q[3],q[11];
u(pi/2,3.601521818075339,-3.601521818075339) q[3];
u(pi,-1.4313096129755096,1.4313096129755096) q[10];
rzz(0.012275082143518329) q[3],q[10];
u(pi,-1.382929086110227,1.382929086110227) q[9];
rzz(pi/128) q[3],q[9];
u(pi,-1.275486617357456,1.275486617357456) q[8];
rzz(0.049100328574073315) q[3],q[8];
u(pi,-0.1551946770873358,0.1551946770873358) q[7];
rzz(pi/32) q[3],q[7];
u(pi,-1.139141496191659,1.139141496191659) q[6];
rzz(0.19640131429629326) q[3],q[6];
u(pi,4.134964250654885,-4.134964250654885) q[5];
rzz(pi/8) q[3],q[5];
u(pi/2,1.6399113651738721,-1.6399113651738721) q[4];
rzz(0.7853830994606742) q[3],q[4];
u(pi/2,0,0) q[2];
rzz(0) q[2],q[16];
u(pi,2.823035158515788,-2.823035158515788) q[2];
rzz(0) q[2],q[15];
u(pi,-0.9550441666912971,0.9550441666912971) q[2];
rzz(0) q[2],q[14];
u(pi,1.5500618152812038,-1.5500618152812038) q[2];
rzz(0) q[2],q[13];
u(pi,4.055167797253705,-4.055167797253705) q[2];
rzz(0) q[2],q[12];
u(pi,0.27708847204661957,-0.27708847204661957) q[2];
rzz(0) q[2],q[11];
u(pi,2.7815661354884025,-2.7815661354884025) q[2];
u(pi/2,1.3056459068319177,-1.3056459068319177) q[10];
rzz(0) q[2],q[10];
u(pi/2,-0.6779556946446773,0.6779556946446773) q[2];
u(pi,-1.2968494474018666,1.2968494474018666) q[9];
rzz(0.012275082143518329) q[2],q[9];
u(pi,-1.2717167061731482,1.2717167061731482) q[8];
rzz(pi/128) q[2],q[8];
u(pi,-0.034557519189487795,0.034557519189487795) q[7];
rzz(0.049100328574073315) q[2],q[7];
u(pi,0.011309733552923307,-0.011309733552923307) q[6];
rzz(pi/32) q[2],q[6];
u(pi,4.354875736406171,-4.354875736406171) q[5];
rzz(0.19640131429629326) q[2],q[5];
u(pi,3.5581678394558,-3.5581678394558) q[4];
rzz(pi/8) q[2],q[4];
u(pi/2,3.601521818075339,-3.601521818075339) q[3];
rzz(0.7853830994606742) q[2],q[3];
u(pi/2,pi/2,-pi/2) q[1];
rzz(0) q[1],q[16];
u(pi,-1.4017786420317657,1.4017786420317657) q[1];
rzz(0) q[1],q[15];
u(pi,2.077849381084289,-2.077849381084289) q[1];
rzz(0) q[1],q[14];
u(pi,-0.7257079029792421,0.7257079029792421) q[1];
rzz(0) q[1],q[13];
u(pi,2.7539201201368124,-2.7539201201368124) q[1];
rzz(0) q[1],q[12];
u(pi,-0.049637163926718575,0.049637163926718575) q[1];
rzz(0) q[1],q[11];
u(pi,3.4299908591893367,-3.4299908591893367) q[1];
rzz(0) q[1],q[10];
u(pi,0.6264335751258048,-0.6264335751258048) q[1];
u(pi/2,3.015300628915483,-3.015300628915483) q[9];
rzz(0) q[1],q[9];
u(pi/2,3.937043913478729,-3.937043913478729) q[1];
rzz(0.012275082143518329) q[1],q[8];
rzz(pi/128) q[1],q[7];
rzz(0.049100328574073315) q[1],q[6];
rzz(pi/32) q[1],q[5];
rzz(0.19640131429629326) q[1],q[4];
rzz(pi/8) q[1],q[3];
u(pi/2,-0.6779556946446773,0.6779556946446773) q[2];
rzz(0.7853830994606742) q[1],q[2];
rzz(-pi/2) q[0],q[16];
rzz(pi/2) q[0],q[15];
rzz(-pi/2) q[0],q[14];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[11];
rzz(pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[9];
u(pi/2,3.5763890768466204,-3.5763890768466204) q[8];
rzz(pi/2) q[0],q[8];
u(pi/2,-1.5148759775609983,1.5148759775609983) q[7];
rzz(pi/2) q[0],q[7];
u(pi/2,1.351513159574329,-1.351513159574329) q[6];
rzz(-pi/2) q[0],q[6];
u(pi/2,2.8745572780346604,-2.8745572780346604) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-0.8067609934418589,0.8067609934418589) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,0.45992916448554544,-0.45992916448554544) q[3];
rzz(-pi/2) q[0],q[3];
u(pi/2,2.4636369589451155,-2.4636369589451155) q[2];
rzz(pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u(pi,-pi/2,pi/2) q[0];
u(pi/2,0.7954512598889356,-0.7954512598889356) q[1];
rzz(pi/2) q[0],q[1];
u(pi/2,4.034433285740012,-4.034433285740012) q[2];
u(pi/2,1.6782387955476676,-1.6782387955476676) q[2];
rzz(-pi/2) q[0],q[2];
u(pi/2,-1.110867162309351,1.110867162309351) q[3];
u(pi/2,2.4234245729791666,-2.4234245729791666) q[3];
rzz(-pi/2) q[0],q[3];
u(pi/2,3.9056279869428314,-3.9056279869428314) q[4];
u(pi/2,0.9606990334677588,-0.9606990334677588) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,4.445353604829557,-4.445353604829557) q[5];
u(pi/2,1.4017786420317657,-1.4017786420317657) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-0.2192831672205675,0.2192831672205675) q[6];
u(pi/2,2.971318331765226,-2.971318331765226) q[6];
rzz(pi/2) q[0],q[6];
u(pi/2,0.055920349233898436,-0.055920349233898436) q[7];
u(pi/2,3.222017425521692,-3.222017425521692) q[7];
rzz(-pi/2) q[0],q[7];
u(pi/2,-1.1359999035380692,1.1359999035380692) q[8];
u(pi/2,2.0175308021353655,-2.0175308021353655) q[8];
rzz(pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[12];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[15];
rzz(pi/2) q[0],q[16];
u(pi/2,2.3662475866838326,-2.3662475866838326) q[1];
u(pi/2,3.2490351223425638,-3.2490351223425638) q[2];
rzz(0.7853830994606742) q[1],q[2];
u(pi/2,3.9942208997740636,-3.9942208997740636) q[3];
rzz(pi/8) q[1],q[3];
u(pi/2,-0.6100972933271378,0.6100972933271378) q[4];
rzz(0.19640131429629326) q[1],q[4];
u(pi/2,-0.16901768476313084,0.16901768476313084) q[5];
rzz(pi/32) q[1],q[5];
u(pi/2,1.4005220049703295,-1.4005220049703295) q[6];
rzz(0.049100328574073315) q[1],q[6];
u(pi/2,-1.4903715548629979,1.4903715548629979) q[7];
rzz(pi/128) q[1],q[7];
u(pi/2,0.4467344753404685,-0.4467344753404685) q[8];
rzz(0.012275082143518329) q[1],q[8];
u(pi/2,-0.7753450669059611,0.7753450669059611) q[1];
rzz(0) q[1],q[9];
u(pi,3.0222121327533813,-3.0222121327533813) q[1];
rzz(0) q[1],q[10];
u(pi,1.1919202527719674,-1.1919202527719674) q[1];
rzz(0) q[1],q[11];
u(pi,-0.6377433086787281,0.6377433086787281) q[1];
rzz(0) q[1],q[12];
u(pi,3.815150118519444,-3.815150118519444) q[1];
rzz(0) q[1],q[13];
u(pi,1.9854865570687488,-1.9854865570687488) q[1];
rzz(0) q[1],q[14];
u(pi,0.15582299561805368,-0.15582299561805368) q[1];
rzz(0) q[1],q[15];
u(pi,4.608716422816227,-4.608716422816227) q[1];
rzz(0) q[1],q[16];
u(pi/2,2.4636369589451155,-2.4636369589451155) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u(pi/2,-0.13885839528866883,0.13885839528866883) q[9];
rzz(0.012275082143518329) q[2],q[9];
u(pi/2,4.034433285740012,-4.034433285740012) q[2];
rzz(0) q[2],q[10];
u(pi,0.5742831370762143,-0.5742831370762143) q[2];
rzz(0) q[2],q[11];
u(pi,3.079389119048715,-3.079389119048715) q[2];
rzz(0) q[2],q[12];
u(pi,-0.69869020615837,0.69869020615837) q[2];
rzz(0) q[2],q[13];
u(pi,1.806415775814131,-1.806415775814131) q[2];
rzz(0) q[2],q[14];
u(pi,4.310893439255914,-4.310893439255914) q[2];
rzz(0) q[2],q[15];
u(pi,0.5328141140488287,-0.5328141140488287) q[2];
rzz(0) q[2],q[16];
u(pi/2,0.45992916448554544,-0.45992916448554544) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u(pi/2,1.2899379435639688,-1.2899379435639688) q[10];
rzz(0.012275082143518329) q[3],q[10];
u(pi/2,-1.110867162309351,1.110867162309351) q[3];
rzz(0) q[3],q[11];
u(pi,1.2132830828163779,-1.2132830828163779) q[3];
rzz(0) q[3],q[12];
u(pi,2.718734282416607,-2.718734282416607) q[3];
rzz(0) q[3],q[13];
u(pi,4.224813800547554,-4.224813800547554) q[3];
rzz(0) q[3],q[14];
u(pi,-0.5529203070318036,0.5529203070318036) q[3];
rzz(0) q[3],q[15];
u(pi,0.9531592110991434,-0.9531592110991434) q[3];
rzz(0) q[3],q[16];
u(pi/2,-0.8061326749111409,0.8061326749111409) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u(pi/2,3.6819465900072377,-3.6819465900072377) q[11];
rzz(0.012275082143518329) q[4],q[11];
u(pi/2,0.7646636518837555,-0.7646636518837555) q[4];
rzz(0) q[4],q[12];
u(pi,2.9204245307770718,-2.9204245307770718) q[4];
rzz(0) q[4],q[13];
u(pi,4.089725316443193,-4.089725316443193) q[4];
rzz(0) q[4],q[14];
u(pi,-1.0241592050702726,1.0241592050702726) q[4];
rzz(0) q[4],q[15];
u(pi,0.14576989912656635,-0.14576989912656635) q[4];
rzz(0) q[4],q[16];
u(pi/2,-0.26703537555513246,0.26703537555513246) q[5];
rzz(0.7853830994606742) q[5],q[6];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u(pi/2,3.333858123989488,-3.333858123989488) q[12];
rzz(0.012275082143518329) q[5],q[12];
u(pi/2,4.445353604829557,-4.445353604829557) q[5];
rzz(0) q[5],q[13];
u(pi,1.838460020880747,-1.838460020880747) q[5];
rzz(0) q[5],q[14];
u(pi,-0.23373449342708064,0.23373449342708064) q[5];
rzz(0) q[5],q[15];
u(pi,3.9766279809139604,-3.9766279809139604) q[5];
rzz(0) q[5],q[16];
u(pi/2,4.493105813164122,-4.493105813164122) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u(pi/2,-0.5208760619651878,0.5208760619651878) q[13];
rzz(0.012275082143518329) q[6],q[13];
u(pi/2,-0.2192831672205675,0.2192831672205675) q[6];
rzz(0) q[6],q[14];
u(pi,2.7319289715616835,-2.7319289715616835) q[6];
rzz(0) q[6],q[15];
u(pi,-0.7910530301739099,0.7910530301739099) q[6];
rzz(0) q[6],q[16];
u(pi/2,1.6267166760287952,-1.6267166760287952) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u(pi/2,2.8374864847223016,-2.8374864847223016) q[14];
rzz(0.012275082143518329) q[7],q[14];
u(pi/2,0.055920349233898436,-0.055920349233898436) q[7];
rzz(0) q[7],q[15];
u(pi,1.9792033717615691,-1.9792033717615691) q[7];
rzz(0) q[7],q[16];
u(pi/2,3.5757607583159032,-3.5757607583159032) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u(pi/2,0.49825659485934093,-0.49825659485934093) q[15];
rzz(0.012275082143518329) q[8],q[15];
u(pi/2,-1.136628222068787,1.136628222068787) q[8];
rzz(0) q[8],q[16];
u(pi/2,2.9907962062174827,-2.9907962062174827) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u(pi/2,0.0948760981384118,-0.0948760981384118) q[16];
rzz(0.012275082143518329) q[9],q[16];
u(pi/2,4.422734137723711,-4.422734137723711) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u(pi/2,3.673778449107904,-3.673778449107904) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
u(pi/2,3.326946620151591,-3.326946620151591) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
u(pi/2,-0.5271592472723674,0.5271592472723674) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
u(pi/2,-0.3103893541746716,0.3103893541746716) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
u(pi/2,0.4919734095521613,-0.4919734095521613) q[15];
rzz(0.7853830994606742) q[15],q[16];
u(pi,pi/2,-pi/2) q[0];
u(pi/2,2.123088315295982,-2.123088315295982) q[1];
u(pi/2,3.3564775910953353,-3.3564775910953353) q[2];
u(pi/2,3.276681137694154,-3.276681137694154) q[3];
u(pi/2,2.3015307780198824,-2.3015307780198824) q[4];
u(pi/2,1.3697343969651499,-1.3697343969651499) q[5];
u(pi/2,2.1601591086083416,-2.1601591086083416) q[6];
u(pi/2,3.901858075758523,-3.901858075758523) q[7];
u(pi/2,2.0049644315210062,-2.0049644315210062) q[8];
u(pi/2,0.08859291283123216,-0.08859291283123216) q[16];
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
