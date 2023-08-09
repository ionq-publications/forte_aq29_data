OPENQASM 2.0;
include "qelib1.inc";
qreg q[23];
creg c[23];
u3(pi/2,2.111150263212341,-2.111150263212341) q[22];
u3(pi/2,4.467344753404686,-4.467344753404686) q[22];
rzz(-pi/2) q[21],q[22];
u3(pi/2,6.0381410801995825,-6.0381410801995825) q[22];
rzz(pi/8) q[20],q[22];
u3(pi/2,2.184663531306342,-2.184663531306342) q[21];
rzz(0.7853830994606742) q[20],q[21];
u3(pi,3.3231767089672832,-3.3231767089672832) q[22];
rzz(0.19640131429629326) q[19],q[22];
u3(pi,0.7766017039673969,-0.7766017039673969) q[21];
rzz(pi/8) q[19],q[21];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[20];
rzz(0.7853830994606742) q[19],q[20];
u3(pi,5.909964099933119,-5.909964099933119) q[22];
rzz(pi/32) q[18],q[22];
u3(pi,2*pi/5,-2*pi/5) q[21];
rzz(0.19640131429629326) q[18],q[21];
u3(pi,0.0043982297150257105,-0.0043982297150257105) q[20];
rzz(pi/8) q[18],q[20];
u3(pi/2,2.4560971365765005,-2.4560971365765005) q[19];
rzz(0.7853830994606742) q[18],q[19];
u3(pi,2.7840794096112744,-2.7840794096112744) q[22];
rzz(0.049100328574073315) q[17],q[22];
u3(pi,2.187805123959932,-2.187805123959932) q[21];
rzz(pi/32) q[17],q[21];
u3(pi,5.191167700791774,-5.191167700791774) q[20];
rzz(0.19640131429629326) q[17],q[20];
u3(pi,1.048035309237555,-1.048035309237555) q[19];
rzz(pi/8) q[17],q[19];
u3(pi/2,2.3580794457844987,-2.3580794457844987) q[18];
rzz(0.7853830994606742) q[17],q[18];
u3(pi,0.21174334485195206,-0.21174334485195206) q[22];
rzz(pi/128) q[16],q[22];
u3(pi,2.9072298416319944,-2.9072298416319944) q[21];
rzz(0.049100328574073315) q[16],q[21];
u3(pi,3.1359377868133316,-3.1359377868133316) q[20];
rzz(pi/32) q[16],q[20];
u3(pi,1.5280706667060753,-1.5280706667060753) q[19];
rzz(0.19640131429629326) q[16],q[19];
u3(pi,6.187052571979739,-6.187052571979739) q[18];
rzz(pi/8) q[16],q[18];
u3(pi/2,5.77424729729804,-5.77424729729804) q[17];
rzz(0.7853830994606742) q[16],q[17];
u3(pi,2.824291795577224,-2.824291795577224) q[22];
rzz(0.012275082143518329) q[15],q[22];
u3(pi,3.627282877834775,-3.627282877834775) q[21];
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
u3(pi/2,5.865353484252144,-5.865353484252144) q[22];
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
u3(pi/2,2.670982074082042,-2.670982074082042) q[21];
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
rzz(pi/2) q[0],q[20];
rzz(pi/2) q[0],q[19];
rzz(-pi/2) q[0],q[18];
rzz(pi/2) q[0],q[17];
rzz(pi/2) q[0],q[16];
rzz(pi/2) q[0],q[15];
rzz(-pi/2) q[0],q[14];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[11];
rzz(pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[9];
u3(pi/2,5.135875670088594,-5.135875670088594) q[8];
rzz(pi/2) q[0],q[8];
u3(pi/2,2.4322210324092177,-2.4322210324092177) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,0.4398229715025711,-0.4398229715025711) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,5.648583591154448,-5.648583591154448) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,1.4275397017912022,-1.4275397017912022) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,0.7527255998001144,-0.7527255998001144) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,2.3920086464432684,-2.3920086464432684) q[2];
rzz(pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi,0,0) q[0];
u3(pi/2,1.6813803882012572,-1.6813803882012572) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,0.8212123196483719,-0.8212123196483719) q[2];
u3(pi/2,3.177406809840717,-3.177406809840717) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,5.465114580184804,-5.465114580184804) q[3];
u3(pi/2,1.930822844896287,-1.930822844896287) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,6.139928682175891,-6.139928682175891) q[4];
u3(pi/2,2.8023006470020957,-2.8023006470020957) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,0.9361946107697583,-0.9361946107697583) q[5];
u3(pi/2,3.979141255036832,-3.979141255036832) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,5.15221195188726,-5.15221195188726) q[6];
u3(pi/2,1.9616104529014666,-1.9616104529014666) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,0.8614247056143213,-0.8614247056143213) q[7];
u3(pi/2,3.978512936506114,-3.978512936506114) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,3.5650793432936974,-3.5650793432936974) q[8];
u3(pi/2,0.41092031908954496,-0.41092031908954496) q[8];
rzz(pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[11];
rzz(pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[13];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[15];
rzz(pi/2) q[0],q[16];
rzz(pi/2) q[0],q[17];
rzz(pi/2) q[0],q[18];
rzz(pi/2) q[0],q[19];
rzz(pi/2) q[0],q[20];
rzz(pi/2) q[0],q[21];
rzz(-pi/2) q[0],q[22];
u3(pi/2,0.11058406140636072,-0.11058406140636072) q[1];
u3(pi/2,1.60661048304582,-1.60661048304582) q[2];
rzz(0.7853830994606742) q[1],q[2];
u3(pi/2,0.36002651810139025,-0.36002651810139025) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,1.231504320207199,-1.231504320207199) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi/2,5.5499375818317285,-5.5499375818317285) q[5];
rzz(pi/32) q[1],q[5];
u3(pi/2,0.39081412610657024,-0.39081412610657024) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi/2,2.4077166097112173,-2.4077166097112173) q[7];
rzz(pi/128) q[1],q[7];
u3(pi/2,5.123309299474235,-5.123309299474235) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi/2,0.11058406140636072,-0.11058406140636072) q[1];
rzz(0) q[1],q[9];
u3(pi,3.678804997353648,-3.678804997353648) q[1];
rzz(0) q[1],q[10];
u3(pi,1.3910972270095605,-1.3910972270095605) q[1];
rzz(0) q[1],q[11];
u3(pi,5.386574763845059,-5.386574763845059) q[1];
rzz(0) q[1],q[12];
u3(pi,3.0988669935009723,-3.0988669935009723) q[1];
rzz(0) q[1],q[13];
u3(pi,0.8111592231568845,-0.8111592231568845) q[1];
rzz(0) q[1],q[14];
u3(pi,4.806636759992384,-4.806636759992384) q[1];
rzz(0) q[1],q[15];
u3(pi,2.518928989648296,-2.518928989648296) q[1];
rzz(0) q[1],q[16];
u3(pi,0.23122121930420877,-0.23122121930420877) q[1];
rzz(0) q[1],q[17];
u3(pi,4.226698756139707,-4.226698756139707) q[1];
rzz(0) q[1],q[18];
u3(pi,1.9389909857956202,-1.9389909857956202) q[1];
rzz(0) q[1],q[19];
u3(pi,5.93446852263112,-5.93446852263112) q[1];
rzz(0) q[1],q[20];
u3(pi,3.6467607522870322,-3.6467607522870322) q[1];
rzz(0) q[1],q[21];
u3(pi,1.3590529819429444,-1.3590529819429444) q[1];
rzz(0) q[1],q[22];
u3(pi/2,3.9628049732381654,-3.9628049732381654) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,5.634760583478653,-5.634760583478653) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi/2,5.5336013000330615,-5.5336013000330615) q[2];
rzz(0) q[2],q[10];
u3(pi,2.8192652473314803,-2.8192652473314803) q[2];
rzz(0) q[2],q[11];
u3(pi,0.531557476987393,-0.531557476987393) q[2];
rzz(0) q[2],q[12];
u3(pi,4.527035013822892,-4.527035013822892) q[2];
rzz(0) q[2],q[13];
u3(pi,2.2393272434788045,-2.2393272434788045) q[2];
rzz(0) q[2],q[14];
u3(pi,6.234804780314303,-6.234804780314303) q[2];
rzz(0) q[2],q[15];
u3(pi,3.947097009970216,-3.947097009970216) q[2];
rzz(0) q[2],q[16];
u3(pi,1.6593892396261287,-1.6593892396261287) q[2];
rzz(0) q[2],q[17];
u3(pi,9*pi/5,-9*pi/5) q[2];
rzz(0) q[2],q[18];
u3(pi,3.3665306875868226,-3.3665306875868226) q[2];
rzz(0) q[2],q[19];
u3(pi,1.0788229172427348,-1.0788229172427348) q[2];
rzz(0) q[2],q[20];
u3(pi,5.074300454078234,-5.074300454078234) q[2];
rzz(0) q[2],q[21];
u3(pi,2.7865926837341464,-2.7865926837341464) q[2];
rzz(0) q[2],q[22];
u3(pi/2,6.250512743582252,-6.250512743582252) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,4.186486370173759,-4.186486370173759) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,4.679716416787356,-4.679716416787356) q[3];
rzz(0) q[3],q[11];
u3(pi,1.5651414600184348,-1.5651414600184348) q[3];
rzz(0) q[3],q[12];
u3(pi,4.760141188719255,-4.760141188719255) q[3];
rzz(0) q[3],q[13];
u3(pi,1.671955610240488,-1.671955610240488) q[3];
rzz(0) q[3],q[14];
u3(pi,4.866955338941307,-4.866955338941307) q[3];
rzz(0) q[3],q[15];
u3(pi,1.778769760462541,-1.778769760462541) q[3];
rzz(0) q[3],q[16];
u3(pi,4.97376948916336,-4.97376948916336) q[3];
rzz(0) q[3],q[17];
u3(pi,1.8855839106845937,-1.8855839106845937) q[3];
rzz(0) q[3],q[18];
u3(pi,5.080583639385413,-5.080583639385413) q[3];
rzz(0) q[3],q[19];
u3(pi,1.9923980609066467,-1.9923980609066467) q[3];
rzz(0) q[3],q[20];
u3(pi,5.187397789607466,-5.187397789607466) q[3];
rzz(0) q[3],q[21];
u3(pi,2.0992122111287,-2.0992122111287) q[3];
rzz(0) q[3],q[22];
u3(pi/2,1.034840620092478,-1.034840620092478) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,5.9508048044297865,-5.9508048044297865) q[11];
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
u3(pi/2,3.4833979343003625,-3.4833979343003625) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,0.5679999517690345,-0.5679999517690345) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,5.054194261095259,-5.054194261095259) q[6];
rzz(0) q[6],q[14];
u3(pi,2.247495384378138,-2.247495384378138) q[6];
rzz(0) q[6],q[15];
u3(pi,6.058247273182556,-6.058247273182556) q[6];
rzz(0) q[6],q[16];
u3(pi,3.585185536276672,-3.585185536276672) q[6];
rzz(0) q[6],q[17];
u3(pi,1.1127521179015047,-1.1127521179015047) q[6];
rzz(0) q[6],q[18];
u3(pi,4.923504006705923,-4.923504006705923) q[6];
rzz(0) q[6],q[19];
u3(pi,2.4510705883307566,-2.4510705883307566) q[6];
rzz(0) q[6],q[20];
u3(pi,6.261822477135176,-6.261822477135176) q[6];
rzz(0) q[6],q[21];
u3(pi,3.7893890587600083,-3.7893890587600083) q[6];
rzz(0) q[6],q[22];
u3(pi/2,2.3832121870132172,-2.3832121870132172) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,5.214415486428338,-5.214415486428338) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,0.8124158602183205,-0.8124158602183205) q[7];
rzz(0) q[7],q[15];
u3(pi,3.981026210628986,-3.981026210628986) q[7];
rzz(0) q[7],q[16];
u3(pi,0.8928406321502192,-0.8928406321502192) q[7];
rzz(0) q[7],q[17];
u3(pi,4.087840360851039,-4.087840360851039) q[7];
rzz(0) q[7],q[18];
u3(pi,0.9996547823722721,-0.9996547823722721) q[7];
rzz(0) q[7],q[19];
u3(pi,4.194654511073091,-4.194654511073091) q[7];
rzz(0) q[7],q[20];
u3(pi,1.106468932594325,-1.106468932594325) q[7];
rzz(0) q[7],q[21];
u3(pi,4.301468661295145,-4.301468661295145) q[7];
rzz(0) q[7],q[22];
u3(pi/2,5.1113712473905935,-5.1113712473905935) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,5.575698641591164,-5.575698641591164) q[15];
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
u3(pi/2,5.039742934888746,-5.039742934888746) q[16];
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
u3(pi/2,5.94326498206117,-5.94326498206117) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,4.321574854278119,-4.321574854278119) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,4.372468655266274,-4.372468655266274) q[11];
rzz(0) q[11],q[19];
u3(pi,1.6575042840339747,-1.6575042840339747) q[11];
rzz(0) q[11],q[20];
u3(pi,5.652981820869474,-5.652981820869474) q[11];
rzz(0) q[11],q[21];
u3(pi,3.365274050525386,-3.365274050525386) q[11];
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
u3(pi/2,3.7033094200516485,-3.7033094200516485) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/128) q[13],q[19];
u3(pi/2,3.9979908109583704,-3.9979908109583704) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi/2,2.1325130932567515,-2.1325130932567515) q[13];
rzz(0) q[13],q[21];
u3(pi,4.795955344970178,-4.795955344970178) q[13];
rzz(0) q[13],q[22];
u3(pi/2,2.066539647531366,-2.066539647531366) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/128) q[14],q[20];
u3(pi/2,2.6527608366912214,-2.6527608366912214) q[21];
rzz(0.012275082143518329) q[14],q[21];
u3(pi/2,3.637335974326262,-3.637335974326262) q[14];
rzz(0) q[14],q[22];
u3(pi/2,2.4278228026941924,-2.4278228026941924) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/32) q[15],q[19];
rzz(0.049100328574073315) q[15],q[20];
rzz(pi/128) q[15],q[21];
u3(pi/2,2.7055395932715296,-2.7055395932715296) q[22];
rzz(0.012275082143518329) q[15],q[22];
u3(pi/2,5.033459749581567,-5.033459749581567) q[16];
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
u3(pi/2,4.315291668970939,-4.315291668970939) q[18];
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
u3(pi/2,2.646477651384042,-2.646477651384042) q[21];
rzz(0.7853830994606742) q[21],q[22];
u3(pi/2,4.927273917890232,-4.927273917890232) q[1];
u3(pi/2,0.07225663103256524,-0.07225663103256524) q[2];
u3(pi/2,5.267822561539365,-5.267822561539365) q[3];
u3(pi/2,0.32295572478903073,-0.32295572478903073) q[4];
u3(pi/2,5.3583004299627515,-5.3583004299627515) q[5];
u3(pi/2,0.9826901820428874,-0.9826901820428874) q[6];
u3(pi/2,1.186265385995506,-1.186265385995506) q[7];
u3(pi/2,4.2405217638155035,-4.2405217638155035) q[8];
u3(pi/2,2.0376369951183397,-2.0376369951183397) q[9];
u3(pi/2,6.024946391054505,-6.024946391054505) q[10];
u3(pi/2,0.6509379978238051,-0.6509379978238051) q[11];
u3(pi/2,1.988628149722339,-1.988628149722339) q[12];
u3(pi/2,1.1755839709733005,-1.1755839709733005) q[13];
u3(pi/2,0.4957433207364693,-0.4957433207364693) q[14];
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
