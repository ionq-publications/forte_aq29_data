OPENQASM 2.0;
include "qelib1.inc";
qreg q[23];
creg c[23];
u3(pi/2,2.111150263212341,-2.111150263212341) q[22];
u3(pi/2,4.467344753404686,-4.467344753404686) q[22];
rzz(pi/2) q[21],q[22];
u3(pi/2,2.8965484266097894,-2.8965484266097894) q[22];
rzz(pi/8) q[20],q[22];
u3(pi/2,5.326256184896136,-5.326256184896136) q[21];
rzz(0.7853830994606742) q[20],q[21];
u3(pi,0.18158405537749003,-0.18158405537749003) q[22];
rzz(0.19640131429629326) q[19],q[22];
u3(pi,3.9181943575571903,-3.9181943575571903) q[21];
rzz(pi/8) q[19],q[21];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[20];
rzz(0.7853830994606742) q[19],q[20];
u3(pi,2.7683714463433255,-2.7683714463433255) q[22];
rzz(pi/32) q[18],q[22];
u3(pi,7*pi/5,-7*pi/5) q[21];
rzz(0.19640131429629326) q[18],q[21];
u3(pi,0.0043982297150257105,-0.0043982297150257105) q[20];
rzz(pi/8) q[18],q[20];
u3(pi/2,2.4560971365765005,-2.4560971365765005) q[19];
rzz(0.7853830994606742) q[18],q[19];
u3(pi,5.925672063201068,-5.925672063201068) q[22];
rzz(0.049100328574073315) q[17],q[22];
u3(pi,5.329397777549724,-5.329397777549724) q[21];
rzz(pi/32) q[17],q[21];
u3(pi,5.191167700791774,-5.191167700791774) q[20];
rzz(0.19640131429629326) q[17],q[20];
u3(pi,1.048035309237555,-1.048035309237555) q[19];
rzz(pi/8) q[17],q[19];
u3(pi/2,2.3580794457844987,-2.3580794457844987) q[18];
rzz(0.7853830994606742) q[17],q[18];
u3(pi,3.353335998441745,-3.353335998441745) q[22];
rzz(pi/128) q[16],q[22];
u3(pi,6.0488224952217875,-6.0488224952217875) q[21];
rzz(0.049100328574073315) q[16],q[21];
u3(pi,3.1359377868133316,-3.1359377868133316) q[20];
rzz(pi/32) q[16],q[20];
u3(pi,1.5280706667060753,-1.5280706667060753) q[19];
rzz(0.19640131429629326) q[16],q[19];
u3(pi,6.187052571979739,-6.187052571979739) q[18];
rzz(pi/8) q[16],q[18];
u3(pi/2,5.77424729729804,-5.77424729729804) q[17];
rzz(0.7853830994606742) q[16],q[17];
u3(pi,5.965884449167017,-5.965884449167017) q[22];
rzz(0.012275082143518329) q[15],q[22];
u3(pi,0.485690224244982,-0.485690224244982) q[21];
rzz(pi/128) q[15],q[21];
u3(pi,0.44799111240190453,-0.44799111240190453) q[20];
rzz(0.049100328574073315) q[15],q[20];
u3(pi,2.041406906302648,-2.041406906302648) q[19];
rzz(pi/32) q[15],q[19];
u3(pi,3.6096899589746725,-3.6096899589746725) q[18];
rzz(0.19640131429629326) q[15],q[18];
u3(pi,1.8302918799814134,-1.8302918799814134) q[17];
rzz(pi/8) q[15],q[17];
u3(pi/2,3.1485041574276904,-3.1485041574276904) q[16];
rzz(0.7853830994606742) q[15],q[16];
u3(pi/2,0.028274333882308135,-0.028274333882308135) q[14];
u3(pi/2,2.7237608306623504,-2.7237608306623504) q[22];
rzz(0) q[14],q[22];
u3(pi/2,3.169866987472101,-3.169866987472101) q[14];
rzz(0.012275082143518329) q[14],q[21];
rzz(pi/128) q[14],q[20];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/32) q[14],q[18];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/8) q[14],q[16];
u3(pi/2,4.224185482016836,-4.224185482016836) q[15];
rzz(0.7853830994606742) q[14],q[15];
u3(pi/2,0.10807078728348889,-0.10807078728348889) q[13];
rzz(0) q[13],q[22];
u3(pi,4.477397849896173,-4.477397849896173) q[13];
u3(pi/2,5.812574727671835,-5.812574727671835) q[21];
rzz(0) q[13],q[21];
u3(pi/2,2.5641679238599893,-2.5641679238599893) q[13];
rzz(0.012275082143518329) q[13],q[20];
rzz(pi/128) q[13],q[19];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/32) q[13],q[17];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/8) q[13],q[15];
u3(pi/2,0.028274333882308135,-0.028274333882308135) q[14];
rzz(0.7853830994606742) q[13],q[14];
u3(pi/2,0.11184069846779664,-0.11184069846779664) q[12];
rzz(0) q[12],q[22];
u3(pi,3.851592593301086,-3.851592593301086) q[12];
rzz(0) q[12],q[21];
u3(pi,1.9056901036675686,-1.9056901036675686) q[12];
u3(pi/2,4.016212048349192,-4.016212048349192) q[20];
rzz(0) q[12],q[20];
u3(pi/2,5.645441998500858,-5.645441998500858) q[12];
u3(pi,2.549716597653476,-2.549716597653476) q[19];
rzz(0.012275082143518329) q[12],q[19];
u3(pi,0.7715751557216532,-0.7715751557216532) q[18];
rzz(pi/128) q[12],q[18];
u3(pi,4.503787228186328,-4.503787228186328) q[17];
rzz(0.049100328574073315) q[12],q[17];
u3(pi,5.674344650913884,-5.674344650913884) q[16];
rzz(pi/32) q[12],q[16];
u3(pi,0.19666370011472106,-0.19666370011472106) q[15];
rzz(0.19640131429629326) q[12],q[15];
u3(pi,2.630141369585375,-2.630141369585375) q[14];
rzz(pi/8) q[12],q[14];
u3(pi/2,2.5641679238599893,-2.5641679238599893) q[13];
rzz(0.7853830994606742) q[12],q[13];
u3(pi/2,0.43228314913395555,-0.43228314913395555) q[11];
rzz(0) q[11],q[22];
u3(pi,4.001132403611961,-4.001132403611961) q[11];
rzz(0) q[11],q[21];
u3(pi,1.713424633267873,-1.713424633267873) q[11];
rzz(0) q[11],q[20];
u3(pi,5.7089021701033715,-5.7089021701033715) q[11];
u3(pi/2,1.2911945806254048,-1.2911945806254048) q[19];
rzz(0) q[11],q[19];
u3(pi/2,2.9939377988710727,-2.9939377988710727) q[11];
rzz(0.012275082143518329) q[11],q[18];
rzz(pi/128) q[11],q[17];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/32) q[11],q[15];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/8) q[11],q[13];
u3(pi/2,5.645441998500858,-5.645441998500858) q[12];
rzz(0.7853830994606742) q[11],q[12];
u3(pi/2,0.44799111240190453,-0.44799111240190453) q[10];
rzz(0) q[10],q[22];
u3(pi,2.7865926837341464,-2.7865926837341464) q[10];
rzz(0) q[10],q[21];
u3(pi,4.322831491339555,-4.322831491339555) q[10];
rzz(0) q[10],q[20];
u3(pi,5.859070298944964,-5.859070298944964) q[10];
rzz(0) q[10],q[19];
u3(pi,1.1121237993707866,-1.1121237993707866) q[10];
u3(pi/2,4.33979609166894,-4.33979609166894) q[18];
rzz(0) q[10],q[18];
u3(pi/2,3.4513536892337466,-3.4513536892337466) q[10];
rzz(0.012275082143518329) q[10],q[17];
rzz(pi/128) q[10],q[16];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/32) q[10],q[14];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/8) q[10],q[12];
u3(pi/2,2.9939377988710727,-2.9939377988710727) q[11];
rzz(0.7853830994606742) q[10],q[11];
u3(pi/2,1.730389233597258,-1.730389233597258) q[9];
rzz(0) q[9],q[22];
u3(pi,3.910654535188574,-3.910654535188574) q[9];
rzz(0) q[9],q[21];
u3(pi,5.130220803312132,-5.130220803312132) q[9];
rzz(0) q[9],q[20];
u3(pi,0.06660176425610362,-0.06660176425610362) q[9];
rzz(0) q[9],q[19];
u3(pi,1.2861680323796612,-1.2861680323796612) q[9];
rzz(0) q[9],q[18];
u3(pi,2.505734300503219,-2.505734300503219) q[9];
u3(pi/2,1.6964600329384885,-1.6964600329384885) q[17];
rzz(0) q[9],q[17];
u3(pi/2,4.686627920625253,-4.686627920625253) q[9];
rzz(0.012275082143518329) q[9],q[16];
rzz(pi/128) q[9],q[15];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/32) q[9],q[13];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/8) q[9],q[11];
u3(pi/2,3.4513536892337466,-3.4513536892337466) q[10];
rzz(0.7853830994606742) q[9],q[10];
u3(pi/2,1.7919644496076181,-1.7919644496076181) q[8];
rzz(0) q[8],q[22];
u3(pi,4.032548330147859,-4.032548330147859) q[8];
rzz(0) q[8],q[21];
u3(pi,5.372123437638546,-5.372123437638546) q[8];
rzz(0) q[8],q[20];
u3(pi,0.42851323794964774,-0.42851323794964774) q[8];
rzz(0) q[8],q[19];
u3(pi,1.7680883454403356,-1.7680883454403356) q[8];
rzz(0) q[8],q[18];
u3(pi,3.1076634529310234,-3.1076634529310234) q[8];
rzz(0) q[8],q[17];
u3(pi,4.447238560421711,-4.447238560421711) q[8];
u3(pi/2,1.9163715186897738,-1.9163715186897738) q[16];
rzz(0) q[8],q[16];
u3(pi/2,0.4046371337823653,-0.4046371337823653) q[8];
rzz(0.012275082143518329) q[8],q[15];
rzz(pi/128) q[8],q[14];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/32) q[8],q[12];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/8) q[8],q[10];
u3(pi/2,4.686627920625253,-4.686627920625253) q[9];
rzz(0.7853830994606742) q[8],q[9];
u3(pi/2,1.03107070890817,-1.03107070890817) q[7];
rzz(0) q[7],q[22];
u3(pi,4.599291644855457,-4.599291644855457) q[7];
rzz(0) q[7],q[21];
u3(pi,2.3115838745113697,-2.3115838745113697) q[7];
rzz(0) q[7],q[20];
u3(pi,0.023876104167282426,-0.023876104167282426) q[7];
rzz(0) q[7],q[19];
u3(pi,4.019353641002781,-4.019353641002781) q[7];
rzz(0) q[7],q[18];
u3(pi,1.7316458706586941,-1.7316458706586941) q[7];
rzz(0) q[7],q[17];
u3(pi,5.727123407494193,-5.727123407494193) q[7];
rzz(0) q[7],q[16];
u3(pi,3.4394156371501055,-3.4394156371501055) q[7];
u3(pi/2,2.4523272253921924,-2.4523272253921924) q[15];
rzz(0) q[7],q[15];
u3(pi/2,0.7244512659178063,-0.7244512659178063) q[7];
rzz(0.012275082143518329) q[7],q[14];
rzz(pi/128) q[7],q[13];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/32) q[7],q[11];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/8) q[7],q[9];
u3(pi/2,0.4046371337823653,-0.4046371337823653) q[8];
rzz(0.7853830994606742) q[7],q[8];
u3(pi/2,4.025008507779242,-4.025008507779242) q[6];
rzz(0) q[6],q[22];
u3(pi,0.19038051480754148,-0.19038051480754148) q[6];
rzz(0) q[6],q[21];
u3(pi,1.944645852572082,-1.944645852572082) q[6];
rzz(0) q[6],q[20];
u3(pi,3.6995395088673404,-3.6995395088673404) q[6];
rzz(0) q[6],q[19];
u3(pi,5.453804846631881,-5.453804846631881) q[6];
rzz(0) q[6],q[18];
u3(pi,0.925513195747553,-0.925513195747553) q[6];
rzz(0) q[6],q[17];
u3(pi,2.6804068520428115,-2.6804068520428115) q[6];
rzz(0) q[6],q[16];
u3(pi,4.434672189807352,-4.434672189807352) q[6];
rzz(0) q[6],q[15];
u3(pi,6.18956584610261,-6.18956584610261) q[6];
u3(pi/2,5.2320084052884415,-5.2320084052884415) q[14];
rzz(0) q[6],q[14];
u3(pi/2,2.354309534600191,-2.354309534600191) q[6];
u3(pi,0.0043982297150257105,-0.0043982297150257105) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi,4.688512876217407,-4.688512876217407) q[12];
rzz(pi/128) q[6],q[12];
u3(pi,0.9123185066024759,-0.9123185066024759) q[11];
rzz(0.049100328574073315) q[6],q[11];
u3(pi,2.1978582204514194,-2.1978582204514194) q[10];
rzz(pi/32) q[6],q[10];
u3(pi,1.9716635493929544,-1.9716635493929544) q[9];
rzz(0.19640131429629326) q[6],q[9];
u3(pi,5.340707511102648,-5.340707511102648) q[8];
rzz(pi/8) q[6],q[8];
u3(pi/2,3.8660439195075993,-3.8660439195075993) q[7];
rzz(0.7853830994606742) q[6],q[7];
u3(pi/2,2.552229871776348,-2.552229871776348) q[5];
rzz(0) q[5],q[22];
u3(pi,5.5939198789819855,-5.5939198789819855) q[5];
rzz(0) q[5],q[21];
u3(pi,2.2525219326238815,-2.2525219326238815) q[5];
rzz(0) q[5],q[20];
u3(pi,5.194309293445364,-5.194309293445364) q[5];
rzz(0) q[5],q[19];
u3(pi,1.8522830285565421,-1.8522830285565421) q[5];
rzz(0) q[5],q[18];
u3(pi,4.7940703893780245,-4.7940703893780245) q[5];
rzz(0) q[5],q[17];
u3(pi,1.4520441244892024,-1.4520441244892024) q[5];
rzz(0) q[5],q[16];
u3(pi,4.393831485310685,-4.393831485310685) q[5];
rzz(0) q[5],q[15];
u3(pi,1.0524335389525807,-1.0524335389525807) q[5];
rzz(0) q[5],q[14];
u3(pi,3.9935925812433455,-3.9935925812433455) q[5];
u3(pi/2,3.727813842749649,-3.727813842749649) q[13];
rzz(0) q[5],q[13];
u3(pi/2,0.7520972812693965,-0.7520972812693965) q[5];
u3(pi,0.7621503777608838,-0.7621503777608838) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi,5.54051280387096,-5.54051280387096) q[11];
rzz(pi/128) q[5],q[11];
u3(pi,2.569822790636451,-2.569822790636451) q[10];
rzz(0.049100328574073315) q[5],q[10];
u3(pi,5.5669021821611135,-5.5669021821611135) q[9];
rzz(pi/32) q[5],q[9];
u3(pi,5.9394950708768635,-5.9394950708768635) q[8];
rzz(0.19640131429629326) q[5],q[8];
u3(pi,1.151707866806018,-1.151707866806018) q[7];
rzz(pi/8) q[5],q[7];
u3(pi/2,2.354309534600191,-2.354309534600191) q[6];
rzz(0.7853830994606742) q[5],q[6];
u3(pi/2,9*pi/8,-9*pi/8) q[4];
rzz(0) q[4],q[22];
u3(pi,0.26954864967800424,-0.26954864967800424) q[4];
rzz(0) q[4],q[21];
u3(pi,3.1642121206956397,-3.1642121206956397) q[4];
rzz(0) q[4],q[20];
u3(pi,6.058875591713275,-6.058875591713275) q[4];
rzz(0) q[4],q[19];
u3(pi,2.670353755551324,-2.670353755551324) q[4];
rzz(0) q[4],q[18];
u3(pi,5.56501722656896,-5.56501722656896) q[4];
rzz(0) q[4],q[17];
u3(pi,2.1771237089377267,-2.1771237089377267) q[4];
rzz(0) q[4],q[16];
u3(pi,5.071787179955362,-5.071787179955362) q[4];
rzz(0) q[4],q[15];
u3(pi,1.6832653437934113,-1.6832653437934113) q[4];
rzz(0) q[4],q[14];
u3(pi,4.577928814811047,-4.577928814811047) q[4];
rzz(0) q[4],q[13];
u3(pi,1.1894069786490957,-1.1894069786490957) q[4];
u3(pi/2,0.9336813366468866,-0.9336813366468866) q[12];
rzz(0) q[4],q[12];
u3(pi/2,4.207849200218169,-4.207849200218169) q[4];
rzz(0.012275082143518329) q[4],q[11];
rzz(pi/128) q[4],q[10];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/32) q[4],q[8];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/8) q[4],q[6];
u3(pi/2,0.7520972812693965,-0.7520972812693965) q[5];
rzz(0.7853830994606742) q[4],q[5];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(0) q[3],q[22];
u3(pi,1.2120264457549421,-1.2120264457549421) q[3];
rzz(0) q[3],q[21];
u3(pi,5.207503982590441,-5.207503982590441) q[3];
rzz(0) q[3],q[20];
u3(pi,2.9197962122463537,-2.9197962122463537) q[3];
rzz(0) q[3],q[19];
u3(pi,0.6320884419022663,-0.6320884419022663) q[3];
rzz(0) q[3],q[18];
u3(pi,4.627565978737765,-4.627565978737765) q[3];
rzz(0) q[3],q[17];
u3(pi,2.339858208393678,-2.339858208393678) q[3];
rzz(0) q[3],q[16];
u3(pi,0.05215043804959057,-0.05215043804959057) q[3];
rzz(0) q[3],q[15];
u3(pi,4.047627974885089,-4.047627974885089) q[3];
rzz(0) q[3],q[14];
u3(pi,1.7599202045410023,-1.7599202045410023) q[3];
rzz(0) q[3],q[13];
u3(pi,5.755397741376501,-5.755397741376501) q[3];
rzz(0) q[3],q[12];
u3(pi,3.467689971032413,-3.467689971032413) q[3];
u3(pi/2,2.8261767511693776,-2.8261767511693776) q[11];
rzz(0) q[3],q[11];
u3(pi/2,0.7527255998001144,-0.7527255998001144) q[3];
rzz(0.012275082143518329) q[3],q[10];
rzz(pi/128) q[3],q[9];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/32) q[3],q[7];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/8) q[3],q[5];
u3(pi/2,4.207849200218169,-4.207849200218169) q[4];
rzz(0.7853830994606742) q[3],q[4];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(0) q[2],q[22];
u3(pi,5.1390172627421835,-5.1390172627421835) q[2];
rzz(0) q[2],q[21];
u3(pi,2.851309492398096,-2.851309492398096) q[2];
rzz(0) q[2],q[20];
u3(pi,0.5636017220540089,-0.5636017220540089) q[2];
rzz(0) q[2],q[19];
u3(pi,4.559079258889508,-4.559079258889508) q[2];
rzz(0) q[2],q[18];
u3(pi,2.2713714885454204,-2.2713714885454204) q[2];
rzz(0) q[2],q[17];
u3(pi,6.266849025380919,-6.266849025380919) q[2];
rzz(0) q[2],q[16];
u3(pi,3.979141255036832,-3.979141255036832) q[2];
rzz(0) q[2],q[15];
u3(pi,1.6914334846927446,-1.6914334846927446) q[2];
rzz(0) q[2],q[14];
u3(pi,5.6869110215282435,-5.6869110215282435) q[2];
rzz(0) q[2],q[13];
u3(pi,3.3992032511841566,-3.3992032511841566) q[2];
rzz(0) q[2],q[12];
u3(pi,1.1114954808400688,-1.1114954808400688) q[2];
rzz(0) q[2],q[11];
u3(pi,5.106973017675568,-5.106973017675568) q[2];
u3(pi/2,1.0543184945447346,-1.0543184945447346) q[10];
rzz(0) q[2],q[10];
u3(pi/2,2.3920086464432684,-2.3920086464432684) q[2];
u3(pi,2.4787166036823467,-2.4787166036823467) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi,0.22682298958918307,-0.22682298958918307) q[8];
rzz(pi/128) q[2],q[8];
u3(pi,5.147185403641517,-5.147185403641517) q[7];
rzz(0.049100328574073315) q[2],q[7];
u3(pi,1.39738041231674,-1.39738041231674) q[6];
rzz(pi/32) q[2],q[6];
u3(pi,3.2000262769465633,-3.2000262769465633) q[5];
rzz(0.19640131429629326) q[2],q[5];
u3(pi,2.8180086102700446,-2.8180086102700446) q[4];
rzz(pi/8) q[2],q[4];
u3(pi/2,3.8943182533899074,-3.8943182533899074) q[3];
rzz(0.7853830994606742) q[2],q[3];
u3(pi/2,pi,-pi) q[1];
rzz(0) q[1],q[22];
u3(pi,5.74345968929286,-5.74345968929286) q[1];
rzz(0) q[1],q[21];
u3(pi,1.5230441184603318,-1.5230441184603318) q[1];
rzz(0) q[1],q[20];
u3(pi,3.5858138548073897,-3.5858138548073897) q[1];
rzz(0) q[1],q[19];
u3(pi,5.648583591154448,-5.648583591154448) q[1];
rzz(0) q[1],q[18];
u3(pi,1.4275397017912022,-1.4275397017912022) q[1];
rzz(0) q[1],q[17];
u3(pi,3.4903094381382602,-3.4903094381382602) q[1];
rzz(0) q[1],q[16];
u3(pi,5.553079174485318,-5.553079174485318) q[1];
rzz(0) q[1],q[15];
u3(pi,1.3326636036527904,-1.3326636036527904) q[1];
rzz(0) q[1],q[14];
u3(pi,3.3948050214691303,-3.3948050214691303) q[1];
rzz(0) q[1],q[13];
u3(pi,5.457574757816189,-5.457574757816189) q[1];
rzz(0) q[1],q[12];
u3(pi,1.2371591869836605,-1.2371591869836605) q[1];
rzz(0) q[1],q[11];
u3(pi,3.2999289233307185,-3.2999289233307185) q[1];
rzz(0) q[1],q[10];
u3(pi,5.362070341147059,-5.362070341147059) q[1];
u3(pi/2,5.647326954093012,-5.647326954093012) q[9];
rzz(0) q[1],q[9];
u3(pi/2,1.6813803882012572,-1.6813803882012572) q[1];
rzz(0.012275082143518329) q[1],q[8];
rzz(pi/128) q[1],q[7];
rzz(0.049100328574073315) q[1],q[6];
rzz(pi/32) q[1],q[5];
rzz(0.19640131429629326) q[1],q[4];
rzz(pi/8) q[1],q[3];
u3(pi/2,5.5336013000330615,-5.5336013000330615) q[2];
rzz(0.7853830994606742) q[1],q[2];
rzz(pi/2) q[0],q[22];
rzz(pi/2) q[0],q[21];
rzz(-pi/2) q[0],q[20];
rzz(pi/2) q[0],q[19];
rzz(pi/2) q[0],q[18];
rzz(pi/2) q[0],q[17];
rzz(pi/2) q[0],q[16];
rzz(-pi/2) q[0],q[15];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[11];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[9];
u3(pi/2,5.135875670088594,-5.135875670088594) q[8];
rzz(pi/2) q[0],q[8];
u3(pi/2,2.4322210324092177,-2.4322210324092177) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,0.4398229715025711,-0.4398229715025711) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,5.648583591154448,-5.648583591154448) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,1.4275397017912022,-1.4275397017912022) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,0.7527255998001144,-0.7527255998001144) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,2.3920086464432684,-2.3920086464432684) q[2];
rzz(-pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi,pi,-pi) q[0];
u3(pi/2,1.6813803882012572,-1.6813803882012572) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,3.9628049732381654,-3.9628049732381654) q[2];
u3(pi/2,0.03581415625092364,-0.03581415625092364) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,5.465114580184804,-5.465114580184804) q[3];
u3(pi/2,1.930822844896287,-1.930822844896287) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,6.139928682175891,-6.139928682175891) q[4];
u3(pi/2,2.8023006470020957,-2.8023006470020957) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,4.077787264359552,-4.077787264359552) q[5];
u3(pi/2,0.8375486014470389,-0.8375486014470389) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,2.0106192982974678,-2.0106192982974678) q[6];
u3(pi/2,5.10320310649126,-5.10320310649126) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,0.8614247056143213,-0.8614247056143213) q[7];
u3(pi/2,3.978512936506114,-3.978512936506114) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,3.5650793432936974,-3.5650793432936974) q[8];
u3(pi/2,0.41092031908954496,-0.41092031908954496) q[8];
rzz(pi/2) q[0],q[8];
rzz(pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[11];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[15];
rzz(-pi/2) q[0],q[16];
rzz(pi/2) q[0],q[17];
rzz(pi/2) q[0],q[18];
rzz(pi/2) q[0],q[19];
rzz(-pi/2) q[0],q[20];
rzz(pi/2) q[0],q[21];
rzz(pi/2) q[0],q[22];
u3(pi/2,3.2521767149961534,-3.2521767149961534) q[1];
u3(pi/2,4.748203136635613,-4.748203136635613) q[2];
rzz(0.7853830994606742) q[1],q[2];
u3(pi/2,3.5016191716911833,-3.5016191716911833) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,1.231504320207199,-1.231504320207199) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi/2,5.5499375818317285,-5.5499375818317285) q[5];
rzz(pi/32) q[1],q[5];
u3(pi/2,3.5324067796963634,-3.5324067796963634) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi/2,5.54930926330101,-5.54930926330101) q[7];
rzz(pi/128) q[1],q[7];
u3(pi/2,5.123309299474235,-5.123309299474235) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi/2,3.2521767149961534,-3.2521767149961534) q[1];
rzz(0) q[1],q[9];
u3(pi,0.5372123437638546,-0.5372123437638546) q[1];
rzz(0) q[1],q[10];
u3(pi,4.532689880599354,-4.532689880599354) q[1];
rzz(0) q[1],q[11];
u3(pi,2.244982110255266,-2.244982110255266) q[1];
rzz(0) q[1],q[12];
u3(pi,6.2404596470907645,-6.2404596470907645) q[1];
rzz(0) q[1],q[13];
u3(pi,3.9527518767466776,-3.9527518767466776) q[1];
rzz(0) q[1],q[14];
u3(pi,1.6650441064025905,-1.6650441064025905) q[1];
rzz(0) q[1],q[15];
u3(pi,5.66052164323809,-5.66052164323809) q[1];
rzz(0) q[1],q[16];
u3(pi,3.3728138728940023,-3.3728138728940023) q[1];
rzz(0) q[1],q[17];
u3(pi,1.0851061025499145,-1.0851061025499145) q[1];
rzz(0) q[1],q[18];
u3(pi,5.080583639385413,-5.080583639385413) q[1];
rzz(0) q[1],q[19];
u3(pi,2.792875869041326,-2.792875869041326) q[1];
rzz(0) q[1],q[20];
u3(pi,0.5051680986972388,-0.5051680986972388) q[1];
rzz(0) q[1],q[21];
u3(pi,4.500645635532738,-4.500645635532738) q[1];
rzz(0) q[1],q[22];
u3(pi/2,0.8212123196483719,-0.8212123196483719) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,5.634760583478653,-5.634760583478653) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi/2,2.3920086464432684,-2.3920086464432684) q[2];
rzz(0) q[2],q[10];
u3(pi,5.960857900921273,-5.960857900921273) q[2];
rzz(0) q[2],q[11];
u3(pi,3.673150130577186,-3.673150130577186) q[2];
rzz(0) q[2],q[12];
u3(pi,1.3854423602330987,-1.3854423602330987) q[2];
rzz(0) q[2],q[13];
u3(pi,5.380919897068598,-5.380919897068598) q[2];
rzz(0) q[2],q[14];
u3(pi,3.0932121267245103,-3.0932121267245103) q[2];
rzz(0) q[2],q[15];
u3(pi,0.8055043563804231,-0.8055043563804231) q[2];
rzz(0) q[2],q[16];
u3(pi,4.800981893215922,-4.800981893215922) q[2];
rzz(0) q[2],q[17];
u3(pi,4*pi/5,-4*pi/5) q[2];
rzz(0) q[2],q[18];
u3(pi,0.2249380339970292,-0.2249380339970292) q[2];
rzz(0) q[2],q[19];
u3(pi,4.220415570832528,-4.220415570832528) q[2];
rzz(0) q[2],q[20];
u3(pi,1.9327078004884406,-1.9327078004884406) q[2];
rzz(0) q[2],q[21];
u3(pi,5.9281853373239395,-5.9281853373239395) q[2];
rzz(0) q[2],q[22];
u3(pi/2,3.108920089992459,-3.108920089992459) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,4.186486370173759,-4.186486370173759) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,1.5381237631975626,-1.5381237631975626) q[3];
rzz(0) q[3],q[11];
u3(pi,4.706734113608228,-4.706734113608228) q[3];
rzz(0) q[3],q[12];
u3(pi,1.6185485351294613,-1.6185485351294613) q[3];
rzz(0) q[3],q[13];
u3(pi,4.813548263830281,-4.813548263830281) q[3];
rzz(0) q[3],q[14];
u3(pi,1.7253626853515145,-1.7253626853515145) q[3];
rzz(0) q[3],q[15];
u3(pi,4.920362414052334,-4.920362414052334) q[3];
rzz(0) q[3],q[16];
u3(pi,1.8321768355735675,-1.8321768355735675) q[3];
rzz(0) q[3],q[17];
u3(pi,5.027176564274387,-5.027176564274387) q[3];
rzz(0) q[3],q[18];
u3(pi,1.9389909857956202,-1.9389909857956202) q[3];
rzz(0) q[3],q[19];
u3(pi,5.1339907144964405,-5.1339907144964405) q[3];
rzz(0) q[3],q[20];
u3(pi,2.045805136017673,-2.045805136017673) q[3];
rzz(0) q[3],q[21];
u3(pi,5.240804864718492,-5.240804864718492) q[3];
rzz(0) q[3],q[22];
u3(pi/2,1.034840620092478,-1.034840620092478) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,2.809212150839993,-2.809212150839993) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,2.6056369468873744,-2.6056369468873744) q[4];
rzz(0) q[4],q[12];
u3(pi,0.4498760679940584,-0.4498760679940584) q[4];
rzz(0) q[4],q[13];
u3(pi,5.56187563391537,-5.56187563391537) q[4];
rzz(0) q[4],q[14];
u3(pi,4.391318211187812,-4.391318211187812) q[4];
rzz(0) q[4],q[15];
u3(pi,3.2207607884602556,-3.2207607884602556) q[4];
rzz(0) q[4],q[16];
u3(pi,2.050203365732699,-2.050203365732699) q[4];
rzz(0) q[4],q[17];
u3(pi,0.8790176244744241,-0.8790176244744241) q[4];
rzz(0) q[4],q[18];
u3(pi,5.991645508926453,-5.991645508926453) q[4];
rzz(0) q[4],q[19];
u3(pi,4.8210880861988965,-4.8210880861988965) q[4];
rzz(0) q[4],q[20];
u3(pi,3.6499023449406214,-3.6499023449406214) q[4];
rzz(0) q[4],q[21];
u3(pi,2.479344922213065,-2.479344922213065) q[4];
rzz(0) q[4],q[22];
u3(pi/2,5.4519198910397275,-5.4519198910397275) q[5];
rzz(0.7853830994606742) q[5],q[6];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,4.059566026968731,-4.059566026968731) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,3.8811235642448305,-3.8811235642448305) q[5];
rzz(0) q[5],q[13];
u3(pi,0.12377875055143785,-0.12377875055143785) q[5];
rzz(0) q[5],q[14];
u3(pi,2.033238765403314,-2.033238765403314) q[5];
rzz(0) q[5],q[15];
u3(pi,3.94269878025519,-3.94269878025519) q[5];
rzz(0) q[5],q[16];
u3(pi,5.852158795107067,-5.852158795107067) q[5];
rzz(0) q[5],q[17];
u3(pi,1.4784335027793567,-1.4784335027793567) q[5];
rzz(0) q[5],q[18];
u3(pi,3.387893517631233,-3.387893517631233) q[5];
rzz(0) q[5],q[19];
u3(pi,5.2973535324831085,-5.2973535324831085) q[5];
rzz(0) q[5],q[20];
u3(pi,0.9236282401553991,-0.9236282401553991) q[5];
rzz(0) q[5],q[21];
u3(pi,2.833088255007276,-2.833088255007276) q[5];
rzz(0) q[5],q[22];
u3(pi/2,0.3418052807105695,-0.3418052807105695) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,3.7095926053588277,-3.7095926053588277) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,1.9126016075054662,-1.9126016075054662) q[6];
rzz(0) q[6],q[14];
u3(pi,5.389088037967931,-5.389088037967931) q[6];
rzz(0) q[6],q[15];
u3(pi,2.916654619592764,-2.916654619592764) q[6];
rzz(0) q[6],q[16];
u3(pi,0.44359288268687874,-0.44359288268687874) q[6];
rzz(0) q[6],q[17];
u3(pi,4.254344771491298,-4.254344771491298) q[6];
rzz(0) q[6],q[18];
u3(pi,1.7819113531161308,-1.7819113531161308) q[6];
rzz(0) q[6],q[19];
u3(pi,5.592663241920549,-5.592663241920549) q[6];
rzz(0) q[6],q[20];
u3(pi,3.1202298235453823,-3.1202298235453823) q[6];
rzz(0) q[6],q[21];
u3(pi,0.6477964051702153,-0.6477964051702153) q[6];
rzz(0) q[6],q[22];
u3(pi/2,5.52480484060301,-5.52480484060301) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,2.0728228328385456,-2.0728228328385456) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,3.9540085138081134,-3.9540085138081134) q[7];
rzz(0) q[7],q[15];
u3(pi,0.8394335570391926,-0.8394335570391926) q[7];
rzz(0) q[7],q[16];
u3(pi,4.034433285740012,-4.034433285740012) q[7];
rzz(0) q[7],q[17];
u3(pi,0.9462477072612457,-0.9462477072612457) q[7];
rzz(0) q[7],q[18];
u3(pi,4.1412474359620655,-4.1412474359620655) q[7];
rzz(0) q[7],q[19];
u3(pi,1.0530618574832986,-1.0530618574832986) q[7];
rzz(0) q[7],q[20];
u3(pi,4.248061586184119,-4.248061586184119) q[7];
rzz(0) q[7],q[21];
u3(pi,1.1598760077053516,-1.1598760077053516) q[7];
rzz(0) q[7],q[22];
u3(pi/2,5.1113712473905935,-5.1113712473905935) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,2.434105988001372,-2.434105988001372) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,0.3989822670059037,-0.3989822670059037) q[8];
rzz(0) q[8],q[16];
u3(pi,4.122397880040527,-4.122397880040527) q[8];
rzz(0) q[8],q[17];
u3(pi,2.1444511453403927,-2.1444511453403927) q[8];
rzz(0) q[8],q[18];
u3(pi,0.167132729170977,-0.167132729170977) q[8];
rzz(0) q[8],q[19];
u3(pi,4.47237130165043,-4.47237130165043) q[8];
rzz(0) q[8],q[20];
u3(pi,2.4950528854810137,-2.4950528854810137) q[8];
rzz(0) q[8],q[21];
u3(pi,0.51710615078088,-0.51710615078088) q[8];
rzz(0) q[8],q[22];
u3(pi/2,5.622822531395012,-5.622822531395012) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,1.8981502812989528,-1.8981502812989528) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,4.052026204600115,-4.052026204600115) q[9];
rzz(0) q[9],q[17];
u3(pi,1.337061833367816,-1.337061833367816) q[9];
rzz(0) q[9],q[18];
u3(pi,5.332539370203315,-5.332539370203315) q[9];
rzz(0) q[9],q[19];
u3(pi,3.0448315998592275,-3.0448315998592275) q[9];
rzz(0) q[9],q[20];
u3(pi,0.7571238295151401,-0.7571238295151401) q[9];
rzz(0) q[9],q[21];
u3(pi,4.752601366350639,-4.752601366350639) q[9];
rzz(0) q[9],q[22];
u3(pi/2,4.177061592212989,-4.177061592212989) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,4.81983144913746,-4.81983144913746) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi/2,5.750999511661475,-5.750999511661475) q[10];
rzz(0) q[10],q[18];
u3(pi,3.036663458959894,-3.036663458959894) q[10];
rzz(0) q[10],q[19];
u3(pi,0.7489556886158066,-0.7489556886158066) q[10];
rzz(0) q[10],q[20];
u3(pi,4.7444332254513055,-4.7444332254513055) q[10];
rzz(0) q[10],q[21];
u3(pi,2.456725455107218,-2.456725455107218) q[10];
rzz(0) q[10],q[22];
u3(pi/2,2.8016723284713776,-2.8016723284713776) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,1.1799822006883263,-1.1799822006883263) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,1.2308760016764808,-1.2308760016764808) q[11];
rzz(0) q[11],q[19];
u3(pi,4.799096937623768,-4.799096937623768) q[11];
rzz(0) q[11],q[20];
u3(pi,2.5113891672796806,-2.5113891672796806) q[11];
rzz(0) q[11],q[21];
u3(pi,0.22368139693559327,-0.22368139693559327) q[11];
rzz(0) q[11],q[22];
u3(pi/2,4.052654523130833,-4.052654523130833) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/128) q[12],q[18];
u3(pi/2,4.414565996824377,-4.414565996824377) q[19];
rzz(0.012275082143518329) q[12],q[19];
u3(pi/2,5.624079168456448,-5.624079168456448) q[12];
rzz(0) q[12],q[20];
u3(pi,2.3593360828459344,-2.3593360828459344) q[12];
rzz(0) q[12],q[21];
u3(pi,5.25399955386357,-5.25399955386357) q[12];
rzz(0) q[12],q[22];
u3(pi/2,0.561716766461855,-0.561716766461855) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/128) q[13],q[19];
u3(pi/2,3.9979908109583704,-3.9979908109583704) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi/2,5.274105746846545,-5.274105746846545) q[13];
rzz(0) q[13],q[21];
u3(pi,1.6543626913803848,-1.6543626913803848) q[13];
rzz(0) q[13],q[22];
u3(pi/2,5.208132301121159,-5.208132301121159) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/128) q[14],q[20];
u3(pi/2,5.794353490281015,-5.794353490281015) q[21];
rzz(0.012275082143518329) q[14],q[21];
u3(pi/2,0.4957433207364693,-0.4957433207364693) q[14];
rzz(0) q[14],q[22];
u3(pi/2,5.569415456283985,-5.569415456283985) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/32) q[15],q[19];
rzz(0.049100328574073315) q[15],q[20];
rzz(pi/128) q[15],q[21];
u3(pi/2,2.7055395932715296,-2.7055395932715296) q[22];
rzz(0.012275082143518329) q[15],q[22];
u3(pi/2,1.8918670959917734,-1.8918670959917734) q[16];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/8) q[16],q[18];
rzz(0.19640131429629326) q[16],q[19];
rzz(pi/32) q[16],q[20];
rzz(0.049100328574073315) q[16],q[21];
rzz(pi/128) q[16],q[22];
u3(pi/2,1.671955610240488,-1.671955610240488) q[17];
rzz(0.7853830994606742) q[17],q[18];
rzz(pi/8) q[17],q[19];
rzz(0.19640131429629326) q[17],q[20];
rzz(pi/32) q[17],q[21];
rzz(0.049100328574073315) q[17],q[22];
u3(pi/2,1.1736990153811466,-1.1736990153811466) q[18];
rzz(0.7853830994606742) q[18],q[19];
rzz(pi/8) q[18],q[20];
rzz(0.19640131429629326) q[18],q[21];
rzz(pi/32) q[18],q[22];
u3(pi/2,1.2666901579274046,-1.2666901579274046) q[19];
rzz(0.7853830994606742) q[19],q[20];
rzz(pi/8) q[19],q[21];
rzz(0.19640131429629326) q[19],q[22];
u3(pi/2,3.991707625651191,-3.991707625651191) q[20];
rzz(0.7853830994606742) q[20],q[21];
rzz(pi/8) q[20],q[22];
u3(pi/2,5.788070304973835,-5.788070304973835) q[21];
rzz(0.7853830994606742) q[21],q[22];
u3(pi/2,1.7856812643004385,-1.7856812643004385) q[1];
u3(pi/2,3.213849284622358,-3.213849284622358) q[2];
u3(pi/2,2.126229907949572,-2.126229907949572) q[3];
u3(pi/2,0.32295572478903073,-0.32295572478903073) q[4];
u3(pi/2,5.3583004299627515,-5.3583004299627515) q[5];
u3(pi/2,4.12428283563268,-4.12428283563268) q[6];
u3(pi/2,4.327858039585299,-4.327858039585299) q[7];
u3(pi/2,4.2405217638155035,-4.2405217638155035) q[8];
u3(pi/2,2.0376369951183397,-2.0376369951183397) q[9];
u3(pi/2,6.024946391054505,-6.024946391054505) q[10];
u3(pi/2,3.7925306514135984,-3.7925306514135984) q[11];
u3(pi/2,1.988628149722339,-1.988628149722339) q[12];
u3(pi/2,4.317176624563094,-4.317176624563094) q[13];
u3(pi/2,3.637335974326262,-3.637335974326262) q[14];
u3(pi/2,2.69925640796435,-2.69925640796435) q[22];
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
