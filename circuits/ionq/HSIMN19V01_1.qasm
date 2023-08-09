OPENQASM 2.0;
include "qelib1.inc";
qreg q[19];
creg c[19];
u3(pi/2,5.085610187631157,-5.085610187631157) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,0.3989822670059037,-0.3989822670059037) q[0];
u3(pi/2,0.1445132620651305,-0.1445132620651305) q[1];
u3(pi/2,3.275424500632718,-3.275424500632718) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,4.846220827427615,-4.846220827427615) q[1];
u3(pi/2,4.85690224244982,-4.85690224244982) q[1];
u3(pi/2,3.5512563356179023,-3.5512563356179023) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,4.731866854836946,-4.731866854836946) q[10];
rzz(pi/2) q[11],q[10];
u3(pi/2,3.2050528251923067,-3.2050528251923067) q[10];
u3(pi/2,3.0297519551219967,-3.0297519551219967) q[11];
u3(pi/2,6.160663193689585,-6.160663193689585) q[11];
rzz(pi/2) q[10],q[11];
u3(pi/2,1.4482742133048947,-1.4482742133048947) q[11];
u3(pi/2,1.4589556283271,-1.4589556283271) q[11];
u3(pi/2,0.07351326809400116,-0.07351326809400116) q[10];
rzz(-pi/2) q[11],q[10];
u3(pi/2,3*pi/2,-3*pi/2) q[13];
u3(pi/2,1.5362388076054088,-1.5362388076054088) q[13];
u3(pi/2,4.770194285210742,-4.770194285210742) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,3.229557247890307,-3.229557247890307) q[12];
u3(pi/2,3.0844156672944587,-3.0844156672944587) q[13];
u3(pi/2,6.215326905862047,-6.215326905862047) q[13];
rzz(pi/2) q[12],q[13];
u3(pi/2,1.502937925477357,-1.502937925477357) q[13];
u3(pi/2,1.5136193404995624,-1.5136193404995624) q[13];
u3(pi/2,0.0986460093227195,-0.0986460093227195) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,4.81103498970741,-4.81103498970741) q[12];
u3(pi/2,1.4589556283271,-1.4589556283271) q[11];
rzz(-pi/2) q[12],q[11];
u3(pi/2,3.0297519551219967,-3.0297519551219967) q[11];
u3(pi/2,4.81103498970741,-4.81103498970741) q[12];
u3(pi/2,1.6587609210954108,-1.6587609210954108) q[12];
rzz(pi/2) q[11],q[12];
u3(pi/2,3.229557247890307,-3.229557247890307) q[12];
u3(pi/2,3.240238662912513,-3.240238662912513) q[12];
u3(pi/2,6.182026023733995,-6.182026023733995) q[11];
rzz(pi/2) q[12],q[11];
u3(pi/2,pi/2,-pi/2) q[15];
u3(pi/2,4.6847429650331,-4.6847429650331) q[15];
u3(pi/2,4.655840312620073,-4.655840312620073) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,3.1145749567689207,-3.1145749567689207) q[14];
u3(pi/2,3.146619201835537,-3.146619201835537) q[15];
u3(pi/2,6.277530440403124,-6.277530440403124) q[15];
rzz(pi/2) q[14],q[15];
u3(pi/2,1.5651414600184348,-1.5651414600184348) q[15];
u3(pi/2,1.5758228750406404,-1.5758228750406404) q[15];
u3(pi/2,6.266849025380919,-6.266849025380919) q[14];
rzz(-pi/2) q[15],q[14];
u3(pi/2,1.5544600449962296,-1.5544600449962296) q[14];
u3(pi/2,4.655211994089355,-4.655211994089355) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,3.0844156672944587,-3.0844156672944587) q[13];
u3(pi/2,4.696052698586023,-4.696052698586023) q[14];
u3(pi/2,1.5437786299740244,-1.5437786299740244) q[14];
rzz(pi/2) q[13],q[14];
u3(pi/2,3.1145749567689207,-3.1145749567689207) q[14];
u3(pi/2,3.125256371791126,-3.125256371791126) q[14];
u3(pi/2,6.23606141737574,-6.23606141737574) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,3*pi/2,-3*pi/2) q[17];
u3(pi/2,1.5155042960917162,-1.5155042960917162) q[17];
u3(pi/2,4.654583675558637,-4.654583675558637) q[16];
rzz(-pi/2) q[17],q[16];
u3(pi/2,6.166318060466046,-6.166318060466046) q[16];
u3(pi/2,pi/5000,-pi/5000) q[17];
u3(pi/2,3.132167875629024,-3.132167875629024) q[17];
rzz(pi/2) q[16],q[17];
u3(pi/2,4.702964202423921,-4.702964202423921) q[17];
u3(pi/2,4.713017298915408,-4.713017298915408) q[17];
u3(pi/2,3.035406821898458,-3.035406821898458) q[16];
rzz(pi/2) q[17],q[16];
u3(pi/2,1.4646104951035617,-1.4646104951035617) q[16];
u3(pi/2,1.5758228750406404,-1.5758228750406404) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,pi/625,-pi/625) q[15];
u3(pi/2,4.606203148693354,-4.606203148693354) q[16];
u3(pi/2,1.4539290800813562,-1.4539290800813562) q[16];
rzz(pi/2) q[15],q[16];
u3(pi/2,3.024725406876253,-3.024725406876253) q[16];
u3(pi/2,3.035406821898458,-3.035406821898458) q[16];
u3(pi/2,3.1573006168577415,-3.1573006168577415) q[15];
rzz(-pi/2) q[16],q[15];
u3(pi/2,1.063114953974786,-1.063114953974786) q[18];
u3(pi/2,4.195282829603809,-4.195282829603809) q[18];
u3(pi/2,1.5714246453256144,-1.5714246453256144) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,pi/5000,-pi/5000) q[17];
u3(pi/2,2.5823891612508096,-2.5823891612508096) q[18];
u3(pi/2,2.5930705762730155,-2.5930705762730155) q[18];
rzz(pi/2) q[17],q[18];
u3(pi/2,4.163866903067912,-4.163866903067912) q[18];
u3(pi/2,4.174548318090117,-4.174548318090117) q[18];
u3(pi/2,6.273760529218817,-6.273760529218817) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,4.795955344970178,-4.795955344970178) q[3];
u3(pi/2,1.6160352610065896,-1.6160352610065896) q[3];
u3(pi/2,4.743176588389869,-4.743176588389869) q[2];
rzz(-pi/2) q[3],q[2];
u3(pi/2,6.276273803341689,-6.276273803341689) q[2];
u3(pi/2,6.271247255095945,-6.271247255095945) q[3];
u3(pi/2,3.1189731864839465,-3.1189731864839465) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,4.689769513278843,-4.689769513278843) q[3];
u3(pi/2,4.7004509283010485,-4.7004509283010485) q[3];
u3(pi/2,3.1453625647741013,-3.1453625647741013) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,1.5745662379792043,-1.5745662379792043) q[2];
u3(pi/2,1.7153095888600272,-1.7153095888600272) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,0.1445132620651305,-0.1445132620651305) q[1];
u3(pi/2,4.716158891568997,-4.716158891568997) q[2];
u3(pi/2,1.563884822956999,-1.563884822956999) q[2];
rzz(-pi/2) q[1],q[2];
u3(pi/2,6.276273803341689,-6.276273803341689) q[2];
u3(pi/2,0.0037699111843077513,-0.0037699111843077513) q[2];
u3(pi/2,0.15456635855661782,-0.15456635855661782) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,1.487229962209408,-1.487229962209408) q[5];
u3(pi/2,4.586725274241098,-4.586725274241098) q[5];
u3(pi/2,4.715530573038279,-4.715530573038279) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,3.0938404452552284,-3.0938404452552284) q[4];
u3(pi/2,3.01027408066974,-3.01027408066974) q[5];
u3(pi/2,6.141185319237328,-6.141185319237328) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,1.428796338852638,-1.428796338852638) q[5];
u3(pi/2,1.4394777538748431,-1.4394777538748431) q[5];
u3(pi/2,6.2461145138672265,-6.2461145138672265) q[4];
rzz(-pi/2) q[5],q[4];
u3(pi/2,1.533725533482537,-1.533725533482537) q[4];
u3(pi/2,1.5588582747112552,-1.5588582747112552) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,6.271247255095945,-6.271247255095945) q[3];
u3(pi/2,4.67531818707233,-4.67531818707233) q[4];
u3(pi/2,1.5230441184603318,-1.5230441184603318) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.0938404452552284,-3.0938404452552284) q[4];
u3(pi/2,3.1045218602774334,-3.1045218602774334) q[4];
u3(pi/2,3.1397076979976393,-3.1397076979976393) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,4.670919957357304,-4.670919957357304) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,3.129654601506152,-3.129654601506152) q[6];
u3(pi/2,3.042318325736356,-3.042318325736356) q[7];
u3(pi/2,6.173229564303944,-6.173229564303944) q[7];
rzz(pi/2) q[6],q[7];
u3(pi/2,1.460840583919254,-1.460840583919254) q[7];
u3(pi/2,1.471521998941459,-1.471521998941459) q[7];
u3(pi/2,6.281300351587433,-6.281300351587433) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,4.710504024792536,-4.710504024792536) q[6];
u3(pi/2,1.4394777538748431,-1.4394777538748431) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,6.151866734259532,-6.151866734259532) q[5];
u3(pi/2,1.5689113712027427,-1.5689113712027427) q[6];
u3(pi/2,4.7004509283010485,-4.7004509283010485) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,6.271247255095945,-6.271247255095945) q[6];
u3(pi/2,6.281300351587433,-6.281300351587433) q[6];
u3(pi/2,3.020955495691945,-3.020955495691945) q[5];
rzz(-pi/2) q[6],q[5];
u3(pi/2,1.487229962209408,-1.487229962209408) q[9];
u3(pi/2,4.587981911302534,-4.587981911302534) q[9];
u3(pi/2,4.695424380055305,-4.695424380055305) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,3.0932121267245103,-3.0932121267245103) q[8];
u3(pi/2,2.959380279681585,-2.959380279681585) q[9];
u3(pi/2,6.090291518249173,-6.090291518249173) q[9];
rzz(pi/2) q[8],q[9];
u3(pi/2,1.3779025378644831,-1.3779025378644831) q[9];
u3(pi/2,1.3885839528866886,-1.3885839528866886) q[9];
u3(pi/2,6.245486195336508,-6.245486195336508) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,1.6443095948888977,-1.6443095948888977) q[10];
u3(pi/2,4.5301766064764815,-4.5301766064764815) q[9];
rzz(-pi/2) q[10],q[9];
u3(pi/2,4.674689868541612,-4.674689868541612) q[8];
u3(pi/2,4.613114652531252,-4.613114652531252) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,3.042318325736356,-3.042318325736356) q[7];
u3(pi/2,1.533097214951819,-1.533097214951819) q[8];
u3(pi/2,4.664008453519407,-4.664008453519407) q[8];
rzz(pi/2) q[7],q[8];
u3(pi/2,6.234804780314303,-6.234804780314303) q[8];
u3(pi/2,6.245486195336508,-6.245486195336508) q[8];
u3(pi/2,6.193964075817636,-6.193964075817636) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,6.100972933271378,-6.100972933271378) q[9];
u3(pi/2,1.6443095948888977,-1.6443095948888977) q[10];
u3(pi/2,4.775849151987203,-4.775849151987203) q[10];
rzz(pi/2) q[9],q[10];
u3(pi/2,0.06346017160251381,-0.06346017160251381) q[10];
u3(pi/2,0.07351326809400116,-0.07351326809400116) q[10];
u3(pi/2,2.9700616947037903,-2.9700616947037903) q[9];
rzz(-pi/2) q[10],q[9];
u3(pi/2,4.866955338941307,-4.866955338941307) q[1];
u3(pi/2,5.121424343882081,-5.121424343882081) q[0];
u3(pi/2,0.41846014145816046,-0.41846014145816046) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,2.015017528012493,-2.015017528012493) q[0];
u3(pi/2,1.7806547160546946,-1.7806547160546946) q[1];
u3(pi/2,4.911565954622282,-4.911565954622282) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,0.19917697423759287,-0.19917697423759287) q[1];
u3(pi/2,3.3307165313358986,-3.3307165313358986) q[1];
u3(pi/2,2.004336112990288,-2.004336112990288) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,4.611229696939098,-4.611229696939098) q[11];
u3(pi/2,0.09361946107697583,-0.09361946107697583) q[10];
rzz(-pi/2) q[11],q[10];
u3(pi/2,1.7077697664914115,-1.7077697664914115) q[10];
u3(pi/2,4.58295536305679,-4.58295536305679) q[11];
u3(pi/2,1.430681294444792,-1.430681294444792) q[11];
rzz(pi/2) q[10],q[11];
u3(pi/2,3.0014776212396885,-3.0014776212396885) q[11];
u3(pi/2,3.0121590362618935,-3.0121590362618935) q[11];
u3(pi/2,4.86004383510341,-4.86004383510341) q[10];
rzz(pi/2) q[11],q[10];
u3(pi/2,3.094468763785946,-3.094468763785946) q[13];
u3(pi/2,4.630707571391355,-4.630707571391355) q[13];
u3(pi/2,0.15582299561805374,-0.15582299561805374) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,4.898999584007923,-4.898999584007923) q[12];
u3(pi/2,6.178884431080405,-6.178884431080405) q[13];
u3(pi/2,3.0266103624684066,-3.0266103624684066) q[13];
rzz(-pi/2) q[12],q[13];
u3(pi/2,1.45581403567351,-1.45581403567351) q[13];
u3(pi/2,1.4664954506957153,-1.4664954506957153) q[13];
u3(pi/2,4.909680999030129,-4.909680999030129) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,3.338884672235232,-3.338884672235232) q[12];
u3(pi/2,6.153751689851687,-6.153751689851687) q[11];
rzz(pi/2) q[12],q[11];
u3(pi/2,4.58295536305679,-4.58295536305679) q[11];
u3(pi/2,0.197292018645439,-0.197292018645439) q[12];
u3(pi/2,3.3282032572130267,-3.3282032572130267) q[12];
rzz(pi/2) q[11],q[12];
u3(pi/2,4.898999584007923,-4.898999584007923) q[12];
u3(pi/2,4.909680999030129,-4.909680999030129) q[12];
u3(pi/2,1.4520441244892024,-1.4520441244892024) q[11];
rzz(-pi/2) q[12],q[11];
u3(pi/2,3.1573006168577415,-3.1573006168577415) q[15];
u3(pi/2,4.755114640473511,-4.755114640473511) q[15];
u3(pi/2,6.231034869129996,-6.231034869129996) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,4.689769513278843,-4.689769513278843) q[14];
u3(pi/2,0.07539822368615504,-0.07539822368615504) q[15];
u3(pi/2,3.2069377807844606,-3.2069377807844606) q[15];
rzz(pi/2) q[14],q[15];
u3(pi/2,4.777734107579358,-4.777734107579358) q[15];
u3(pi/2,4.7877872040708445,-4.7877872040708445) q[15];
u3(pi/2,1.5588582747112552,-1.5588582747112552) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,6.271247255095945,-6.271247255095945) q[14];
u3(pi/2,4.608088104285509,-4.608088104285509) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,3.037291777490612,-3.037291777490612) q[13];
u3(pi/2,3.129654601506152,-3.129654601506152) q[14];
u3(pi/2,6.260565840073739,-6.260565840073739) q[14];
rzz(-pi/2) q[13],q[14];
u3(pi/2,4.689769513278843,-4.689769513278843) q[14];
u3(pi/2,4.7004509283010485,-4.7004509283010485) q[14];
u3(pi/2,3.0473448739820994,-3.0473448739820994) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,3.132167875629024,-3.132167875629024) q[17];
u3(pi/2,4.757627914596383,-4.757627914596383) q[17];
u3(pi/2,2.9983360285860985,-2.9983360285860985) q[16];
rzz(pi/2) q[17],q[16];
u3(pi/2,1.368477759903714,-1.368477759903714) q[16];
u3(pi/2,6.27250389215738,-6.27250389215738) q[17];
u3(pi/2,0,0) q[17];
rzz(pi/2) q[16],q[17];
u3(pi/2,pi/2,-pi/2) q[17];
u3(pi/2,1.5814777418171018,-1.5814777418171018) q[17];
u3(pi/2,1.3577963448815085,-1.3577963448815085) q[16];
rzz(-pi/2) q[17],q[16];
u3(pi/2,2.9285926716764052,-2.9285926716764052) q[16];
u3(pi/2,1.6461945504810516,-1.6461945504810516) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,0.07539822368615504,-0.07539822368615504) q[15];
u3(pi/2,6.070185325266198,-6.070185325266198) q[16];
u3(pi/2,6.080866740288403,-6.080866740288403) q[16];
rzz(pi/2) q[15],q[16];
u3(pi/2,1.368477759903714,-1.368477759903714) q[16];
u3(pi/2,1.379159174925919,-1.379159174925919) q[16];
u3(pi/2,0.06534512719466769,-0.06534512719466769) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,1.032955664500324,-1.032955664500324) q[18];
u3(pi/2,4.165123540129348,-4.165123540129348) q[18];
u3(pi/2,1.5814777418171018,-1.5814777418171018) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,0.010681415022205296,-0.010681415022205296) q[17];
u3(pi/2,2.552229871776348,-2.552229871776348) q[18];
u3(pi/2,2.562911286798553,-2.562911286798553) q[18];
rzz(-pi/2) q[17],q[18];
u3(pi/2,0.9921149600036567,-0.9921149600036567) q[18];
u3(pi/2,1.002168056495144,-1.002168056495144) q[18];
u3(pi/2,pi,-pi) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,6.281300351587433,-6.281300351587433) q[3];
u3(pi/2,1.5305839408289472,-1.5305839408289472) q[3];
u3(pi/2,3.1968846842929737,-3.1968846842929737) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,1.5883892456549995,-1.5883892456549995) q[2];
u3(pi/2,3.0442032813285094,-3.0442032813285094) q[3];
u3(pi/2,6.1751145198960975,-6.1751145198960975) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,1.4627255395114078,-1.4627255395114078) q[3];
u3(pi/2,1.4734069545336128,-1.4734069545336128) q[3];
u3(pi/2,4.740663314266998,-4.740663314266998) q[2];
rzz(-pi/2) q[3],q[2];
u3(pi/2,0.028274333882308135,-0.028274333882308135) q[2];
u3(pi/2,0.18912387774610553,-0.18912387774610553) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,4.901512858130795,-4.901512858130795) q[1];
u3(pi/2,3.169866987472101,-3.169866987472101) q[2];
u3(pi/2,0.017592918860102842,-0.017592918860102842) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,1.5883892456549995,-1.5883892456549995) q[2];
u3(pi/2,4.719928802753305,-4.719928802753305) q[2];
u3(pi/2,4.89083144310859,-4.89083144310859) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3.020955495691945,-3.020955495691945) q[5];
u3(pi/2,4.633220845514227,-4.633220845514227) q[5];
u3(pi/2,6.269990618034509,-6.269990618034509) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,4.648928808782176,-4.648928808782176) q[4];
u3(pi/2,6.198362305532662,-6.198362305532662) q[5];
u3(pi/2,3.0460882369206637,-3.0460882369206637) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,4.61688456371556,-4.61688456371556) q[5];
u3(pi/2,4.627565978737765,-4.627565978737765) q[5];
u3(pi/2,1.51738925168387,-1.51738925168387) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,6.22977823206856,-6.22977823206856) q[4];
u3(pi/2,1.4734069545336128,-1.4734069545336128) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,6.1857959349183025,-6.1857959349183025) q[3];
u3(pi/2,3.0881855784787664,-3.0881855784787664) q[4];
u3(pi/2,6.219725135577073,-6.219725135577073) q[4];
rzz(-pi/2) q[3],q[4];
u3(pi/2,4.648928808782176,-4.648928808782176) q[4];
u3(pi/2,4.658981905273664,-4.658981905273664) q[4];
u3(pi/2,6.1964773499405075,-6.1964773499405075) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,4.623167749022739,-4.623167749022739) q[7];
u3(pi/2,6.239831328560047,-6.239831328560047) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,4.698565972708895,-4.698565972708895) q[6];
u3(pi/2,1.4664954506957153,-1.4664954506957153) q[7];
u3(pi/2,4.597406689263304,-4.597406689263304) q[7];
rzz(pi/2) q[6],q[7];
u3(pi/2,6.1682030160582,-6.1682030160582) q[7];
u3(pi/2,6.178884431080405,-6.178884431080405) q[7];
u3(pi/2,1.5676547341413067,-1.5676547341413067) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,6.280043714525997,-6.280043714525997) q[6];
u3(pi/2,1.485973325147972,-1.485973325147972) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,6.198362305532662,-6.198362305532662) q[5];
u3(pi/2,3.1384510609362035,-3.1384510609362035) q[6];
u3(pi/2,6.269362299503792,-6.269362299503792) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,1.5569733191191015,-1.5569733191191015) q[6];
u3(pi/2,1.5676547341413067,-1.5676547341413067) q[6];
u3(pi/2,3.066822748434356,-3.066822748434356) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,2.9700616947037903,-2.9700616947037903) q[9];
u3(pi/2,4.581070407464636,-4.581070407464636) q[9];
u3(pi/2,3.0869289414173307,-3.0869289414173307) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,1.4847166880865363,-1.4847166880865363) q[8];
u3(pi/2,6.094061429433481,-6.094061429433481) q[9];
u3(pi/2,2.9424156793522003,-2.9424156793522003) q[9];
rzz(-pi/2) q[8],q[9];
u3(pi/2,1.3716193525573037,-1.3716193525573037) q[9];
u3(pi/2,1.381672449048791,-1.381672449048791) q[9];
u3(pi/2,1.4953981031087415,-1.4953981031087415) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,3.289247508308513,-3.289247508308513) q[10];
u3(pi/2,4.523265102638584,-4.523265102638584) q[9];
rzz(pi/2) q[10],q[9];
u3(pi/2,6.2077870834934314,-6.2077870834934314) q[8];
u3(pi/2,3.037291777490612,-3.037291777490612) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,1.4664954506957153,-1.4664954506957153) q[7];
u3(pi/2,3.066194429903638,-3.066194429903638) q[8];
u3(pi/2,6.197105668471226,-6.197105668471226) q[8];
rzz(-pi/2) q[7],q[8];
u3(pi/2,4.626309341676329,-4.626309341676329) q[8];
u3(pi/2,4.636990756698535,-4.636990756698535) q[8];
u3(pi/2,1.4771768657179207,-1.4771768657179207) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,2.9524687758436876,-2.9524687758436876) q[9];
u3(pi/2,0.14765485471872028,-0.14765485471872028) q[10];
u3(pi/2,3.2785660932863085,-3.2785660932863085) q[10];
rzz(pi/2) q[9],q[10];
u3(pi/2,4.849362420081205,-4.849362420081205) q[10];
u3(pi/2,4.86004383510341,-4.86004383510341) q[10];
u3(pi/2,6.104742844455686,-6.104742844455686) q[9];
rzz(pi/2) q[10],q[9];
u3(pi/2,3.320035116313693,-3.320035116313693) q[1];
u3(pi/2,3.5751324397851842,-3.5751324397851842) q[0];
u3(pi/2,5.154725226010132,-5.154725226010132) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,0.46809730538487915,-0.46809730538487915) q[0];
u3(pi/2,0.23373449342708058,-0.23373449342708058) q[1];
u3(pi/2,3.3646457319946683,-3.3646457319946683) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,1.7938494051997718,-1.7938494051997718) q[1];
u3(pi/2,4.9253889622980775,-4.9253889622980775) q[1];
u3(pi/2,3.599008543952467,-3.599008543952467) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.022840451284099,-3.022840451284099) q[11];
u3(pi/2,1.7379290559658735,-1.7379290559658735) q[10];
rzz(pi/2) q[11],q[10];
u3(pi/2,0.21048670779051615,-0.21048670779051615) q[10];
u3(pi/2,6.136158770991584,-6.136158770991584) q[11];
u3(pi/2,2.9838847023795854,-2.9838847023795854) q[11];
rzz(pi/2) q[10],q[11];
u3(pi/2,4.554681029174482,-4.554681029174482) q[11];
u3(pi/2,4.565362444196688,-4.565362444196688) q[11];
u3(pi/2,3.3627607764025145,-3.3627607764025145) q[10];
rzz(-pi/2) q[11],q[10];
u3(pi/2,6.188937527571892,-6.188937527571892) q[13];
u3(pi/2,1.441991027997715,-1.441991027997715) q[13];
u3(pi/2,4.966857985325463,-4.966857985325463) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,3.426849266535746,-3.426849266535746) q[12];
u3(pi/2,2.990167887686765,-2.990167887686765) q[13];
u3(pi/2,6.121079126254353,-6.121079126254353) q[13];
rzz(pi/2) q[12],q[13];
u3(pi/2,1.4086901458696632,-1.4086901458696632) q[13];
u3(pi/2,1.4193715608918684,-1.4193715608918684) q[13];
u3(pi/2,0.29593802796815855,-0.29593802796815855) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,5.008327008352849,-5.008327008352849) q[12];
u3(pi/2,4.565362444196688,-4.565362444196688) q[11];
rzz(pi/2) q[12],q[11];
u3(pi/2,2.994566117401791,-2.994566117401791) q[11];
u3(pi/2,1.8667343547630548,-1.8667343547630548) q[12];
u3(pi/2,4.997645593330643,-4.997645593330643) q[12];
rzz(-pi/2) q[11],q[12];
u3(pi/2,3.426849266535746,-3.426849266535746) q[12];
u3(pi/2,3.4375306815579516,-3.4375306815579516) q[12];
u3(pi/2,3.0052475324239962,-3.0052475324239962) q[11];
rzz(pi/2) q[12],q[11];
u3(pi/2,3.2069377807844606,-3.2069377807844606) q[15];
u3(pi/2,4.7500880922277675,-4.7500880922277675) q[15];
u3(pi/2,1.5230441184603318,-1.5230441184603318) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,6.264964069788765,-6.264964069788765) q[14];
u3(pi/2,pi/625,-pi/625) q[15];
u3(pi/2,pi/200,-pi/200) q[15];
rzz(-pi/2) q[14],q[15];
u3(pi/2,4.7280969436526386,-4.7280969436526386) q[15];
u3(pi/2,4.738150040144126,-4.738150040144126) q[15];
u3(pi/2,3.112690001176767,-3.112690001176767) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,1.5418936743818705,-1.5418936743818705) q[14];
u3(pi/2,4.560964214481662,-4.560964214481662) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,2.990167887686765,-2.990167887686765) q[13];
u3(pi/2,4.683486327971663,-4.683486327971663) q[14];
u3(pi/2,4.694167742993868,-4.694167742993868) q[14];
rzz(pi/2) q[13],q[14];
u3(pi/2,6.264964069788765,-6.264964069788765) q[14];
u3(pi/2,6.275017166280253,-6.275017166280253) q[14];
u3(pi/2,2.97948647266456,-2.97948647266456) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,0,0) q[17];
u3(pi/2,1.6260883574980767,-1.6260883574980767) q[17];
u3(pi/2,4.483681035203353,-4.483681035203353) q[16];
rzz(-pi/2) q[17],q[16];
u3(pi/2,5.995415420110762,-5.995415420110762) q[16];
u3(pi/2,6.282556988648868,-6.282556988648868) q[17];
u3(pi/2,0.00942477796076938,-0.00942477796076938) q[17];
rzz(pi/2) q[16],q[17];
u3(pi/2,1.580221104755666,-1.580221104755666) q[17];
u3(pi/2,1.590902519777871,-1.590902519777871) q[17];
u3(pi/2,5.984734005088556,-5.984734005088556) q[16];
rzz(pi/2) q[17],q[16];
u3(pi/2,4.413937678293659,-4.413937678293659) q[16];
u3(pi/2,1.5965573865543328,-1.5965573865543328) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,0.025761059759436305,-0.025761059759436305) q[15];
u3(pi/2,1.2723450247038663,-1.2723450247038663) q[16];
u3(pi/2,1.2830264397260716,-1.2830264397260716) q[16];
rzz(-pi/2) q[15],q[16];
u3(pi/2,5.995415420110762,-5.995415420110762) q[16];
u3(pi/2,6.006096835132967,-6.006096835132967) q[16];
u3(pi/2,3.1573006168577415,-3.1573006168577415) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,4.143760710084937,-4.143760710084937) q[18];
u3(pi/2,0.9933715970650925,-0.9933715970650925) q[18];
u3(pi/2,4.732495173367664,-4.732495173367664) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,3.161698846572768,-3.161698846572768) q[17];
u3(pi/2,5.663663235891679,-5.663663235891679) q[18];
u3(pi/2,5.674344650913884,-5.674344650913884) q[18];
rzz(pi/2) q[17],q[18];
u3(pi/2,0.9619556705291947,-0.9619556705291947) q[18];
u3(pi/2,0.9720087670206821,-0.9720087670206821) q[18];
u3(pi/2,3.1510174315505624,-3.1510174315505624) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,3.054884696350715,-3.054884696350715) q[3];
u3(pi/2,4.586725274241098,-4.586725274241098) q[3];
u3(pi/2,1.5261857111139214,-1.5261857111139214) q[2];
rzz(-pi/2) q[3],q[2];
u3(pi/2,3.0630528372500483,-3.0630528372500483) q[2];
u3(pi/2,2.9587519611508672,-2.9587519611508672) q[3];
u3(pi/2,6.089663199718455,-6.089663199718455) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,1.3772742193337653,-1.3772742193337653) q[3];
u3(pi/2,4.5088137764320715,-4.5088137764320715) q[3];
u3(pi/2,3.052999740758561,-3.052999740758561) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,1.4822034139636644,-1.4822034139636644) q[2];
u3(pi/2,1.7837963087082844,-1.7837963087082844) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,0.21299998191338798,-0.21299998191338798) q[1];
u3(pi/2,4.6237960675534575,-4.6237960675534575) q[2];
u3(pi/2,1.471521998941459,-1.471521998941459) q[2];
rzz(-pi/2) q[1],q[2];
u3(pi/2,6.183910979326148,-6.183910979326148) q[2];
u3(pi/2,3.03163691071415,-3.03163691071415) q[2];
u3(pi/2,3.343911220480976,-3.343911220480976) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,6.208415402024149,-6.208415402024149) q[5];
u3(pi/2,1.5374954446668447,-1.5374954446668447) q[5];
u3(pi/2,1.5412653558511524,-1.5412653558511524) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,6.203388853778405,-6.203388853778405) q[4];
u3(pi/2,3.10263690468528,-3.10263690468528) q[5];
u3(pi/2,6.233548143252867,-6.233548143252867) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,1.521159162868178,-1.521159162868178) q[5];
u3(pi/2,1.531840577890383,-1.531840577890383) q[5];
u3(pi/2,3.0718492966801,-3.0718492966801) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,1.501052969885203,-1.501052969885203) q[4];
u3(pi/2,1.367221122842278,-1.367221122842278) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,6.079610103226968,-6.079610103226968) q[3];
u3(pi/2,4.642645623474996,-4.642645623474996) q[4];
u3(pi/2,1.490999873393716,-1.490999873393716) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.0617962001886125,-3.0617962001886125) q[4];
u3(pi/2,6.192707438756201,-6.192707438756201) q[4];
u3(pi/2,6.068928688204762,-6.068928688204762) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,6.18956584610261,-6.18956584610261) q[7];
u3(pi/2,4.6671500461729964,-4.6671500461729964) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,3.1258846903218442,-3.1258846903218442) q[6];
u3(pi/2,3.032265229244868,-3.032265229244868) q[7];
u3(pi/2,6.163804786343174,-6.163804786343174) q[7];
rzz(-pi/2) q[6],q[7];
u3(pi/2,4.593008459548278,-4.593008459548278) q[7];
u3(pi/2,4.603061556039765,-4.603061556039765) q[7];
u3(pi/2,3.1365661053440492,-3.1365661053440492) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,1.565769778549153,-1.565769778549153) q[6];
u3(pi/2,4.6734332314801765,-4.6734332314801765) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,3.10263690468528,-3.10263690468528) q[5];
u3(pi/2,4.707362432138946,-4.707362432138946) q[6];
u3(pi/2,1.5550883635269477,-1.5550883635269477) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,3.1258846903218442,-3.1258846903218442) q[6];
u3(pi/2,3.1365661053440492,-3.1365661053440492) q[6];
u3(pi/2,6.254910973297278,-6.254910973297278) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,2.963150190865893,-2.963150190865893) q[9];
u3(pi/2,4.574787222157457,-4.574787222157457) q[9];
u3(pi/2,1.4784335027793567,-1.4784335027793567) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,6.160034875158867,-6.160034875158867) q[8];
u3(pi/2,6.087778244126301,-6.087778244126301) q[9];
u3(pi/2,2.9355041755143025,-2.9355041755143025) q[9];
rzz(pi/2) q[8],q[9];
u3(pi/2,4.506300502309199,-4.506300502309199) q[9];
u3(pi/2,4.516981917331404,-4.516981917331404) q[9];
u3(pi/2,3.0284953180605605,-3.0284953180605605) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,4.933557103197411,-4.933557103197411) q[10];
u3(pi/2,1.3753892637416114,-1.3753892637416114) q[9];
rzz(-pi/2) q[10],q[9];
u3(pi/2,1.4576989912656642,-1.4576989912656642) q[8];
u3(pi/2,1.4614689024499719,-1.4614689024499719) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,6.173857882834661,-6.173857882834661) q[7];
u3(pi/2,4.599291644855457,-4.599291644855457) q[8];
u3(pi/2,1.4476458947741766,-1.4476458947741766) q[8];
rzz(pi/2) q[7],q[8];
u3(pi/2,3.018442221569073,-3.018442221569073) q[8];
u3(pi/2,3.0284953180605605,-3.0284953180605605) q[8];
u3(pi/2,3.0429466442670736,-3.0429466442670736) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,2.946185590536508,-2.946185590536508) q[9];
u3(pi/2,4.933557103197411,-4.933557103197411) q[10];
u3(pi/2,1.7812830345854125,-1.7812830345854125) q[10];
rzz(-pi/2) q[9],q[10];
u3(pi/2,0.21048670779051615,-0.21048670779051615) q[10];
u3(pi/2,0.22116812281272144,-0.22116812281272144) q[10];
u3(pi/2,2.9562386870279953,-2.9562386870279953) q[9];
rzz(pi/2) q[10],q[9];
u3(pi/2,2.0282122171575705,-2.0282122171575705) q[0];
u3(pi/2,1.7731148936860792,-1.7731148936860792) q[1];
u3(pi/2,4.498132361409866,-4.498132361409866) q[3];
u3(pi/2,4.6841146465023815,-4.6841146465023815) q[5];
u3(pi,5.10948629179844,-5.10948629179844) q[6];
u3(pi/2,1.472150317472177,-1.472150317472177) q[7];
u3(pi,5.355158837309161,-5.355158837309161) q[8];
u3(pi/2,1.3854423602330987,-1.3854423602330987) q[9];
u3(pi,1.7932210866690539,-1.7932210866690539) q[10];
u3(pi/2,1.4344512056290994,-1.4344512056290994) q[11];
u3(pi,5.10948629179844,-5.10948629179844) q[12];
u3(pi/2,4.550282799459456,-4.550282799459456) q[13];
u3(pi,0.4800353574685204,-0.4800353574685204) q[14];
u3(pi/2,4.7280969436526386,-4.7280969436526386) q[15];
u3(pi,1.318212277446277,-1.318212277446277) q[16];
u3(pi/2,4.7218137583454585,-4.7218137583454585) q[17];
u3(pi,2.5139024414025526,-2.5139024414025526) q[18];
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
