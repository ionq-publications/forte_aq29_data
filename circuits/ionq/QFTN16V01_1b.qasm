OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u(pi/2,3.739123576302571,-3.739123576302571) q[15];
u(pi/2,1.3829290861102272,-1.3829290861102272) q[15];
rzz(-pi/2) q[14],q[15];
u(pi/2,-0.18786724068466976,0.18786724068466976) q[15];
rzz(pi/8) q[13],q[15];
u(pi/2,0.9167167363175013,-0.9167167363175013) q[14];
rzz(0.7853830994606742) q[13],q[14];
u(pi,3.6096899589746725,-3.6096899589746725) q[15];
rzz(0.19640131429629326) q[12],q[15];
u(pi,-0.5636017220540088,0.5636017220540088) q[14];
rzz(pi/8) q[12],q[14];
u(pi/2,0.819327364056218,-0.819327364056218) q[13];
rzz(0.7853830994606742) q[12],q[13];
u(pi,2.5635396053292716,-2.5635396053292716) q[15];
rzz(pi/32) q[11],q[15];
u(pi,-0.20169024836046456,0.20169024836046456) q[14];
rzz(0.19640131429629326) q[11],q[14];
u(pi,3.31186697541436,-3.31186697541436) q[13];
rzz(pi/8) q[11],q[13];
u(pi/2,2.0954422999443922,-2.0954422999443922) q[12];
rzz(0.7853830994606742) q[11],q[12];
u(pi,0.3367787324648257,-0.3367787324648257) q[15];
rzz(0.049100328574073315) q[10],q[15];
u(pi,-0.3304955471576463,0.3304955471576463) q[14];
rzz(pi/32) q[10],q[14];
u(pi,-1.00342469355658,1.00342469355658) q[13];
rzz(0.19640131429629326) q[10],q[13];
u(pi,2.2110529095964964,-2.2110529095964964) q[12];
rzz(pi/8) q[10],q[12];
u(pi/2,1.7058848108992573,-1.7058848108992573) q[11];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/128) q[9],q[15];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/32) q[9],q[13];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/8) q[9],q[11];
u(pi/2,2.9084864786934306,-2.9084864786934306) q[10];
rzz(0.7853830994606742) q[9],q[10];
rzz(0.012275082143518329) q[8],q[15];
rzz(pi/128) q[8],q[14];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/32) q[8],q[12];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/8) q[8],q[10];
u(pi/2,3.656813848778519,-3.656813848778519) q[9];
rzz(0.7853830994606742) q[8],q[9];
u(pi/2,0.589362781813445,-0.589362781813445) q[7];
u(pi/2,2.954353731435841,-2.954353731435841) q[15];
rzz(0) q[7],q[15];
u(pi/2,3.7309554354032386,-3.7309554354032386) q[7];
u(pi,-0.9952565526572464,0.9952565526572464) q[14];
rzz(0.012275082143518329) q[7],q[14];
u(pi,0.46495571273128933,-0.46495571273128933) q[13];
rzz(pi/128) q[7],q[13];
u(pi,4.227955393201144,-4.227955393201144) q[12];
rzz(0.049100328574073315) q[7],q[12];
u(pi,-0.3455751918948773,0.3455751918948773) q[11];
rzz(pi/32) q[7],q[11];
u(pi,0.8570264758992958,-0.8570264758992958) q[10];
rzz(0.19640131429629326) q[7],q[10];
u(pi,2.322265289533575,-2.322265289533575) q[9];
rzz(pi/8) q[7],q[9];
u(pi/2,0.5397256178867265,-0.5397256178867265) q[8];
rzz(0.7853830994606742) q[7],q[8];
u(pi/2,0.5887344632827274,-0.5887344632827274) q[6];
rzz(0) q[6],q[15];
u(pi,3.388521836161951,-3.388521836161951) q[6];
u(pi/2,3.452610326295183,-3.452610326295183) q[14];
rzz(0) q[6],q[14];
u(pi/2,-0.09550441666912968,0.09550441666912968) q[6];
rzz(0.012275082143518329) q[6],q[13];
rzz(pi/128) q[6],q[12];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/32) q[6],q[10];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/8) q[6],q[8];
u(pi/2,3.7309554354032386,-3.7309554354032386) q[7];
rzz(0.7853830994606742) q[6],q[7];
u(pi/2,pi/4,-pi/4) q[5];
rzz(0) q[5],q[15];
u(pi,4.033176648678577,-4.033176648678577) q[5];
rzz(0) q[5],q[14];
u(pi,1.1045839770021715,-1.1045839770021715) q[5];
u(pi/2,2.456725455107218,-2.456725455107218) q[13];
rzz(0) q[5],q[13];
u(pi/2,4.352362462283299,-4.352362462283299) q[5];
u(pi,2.796645780225634,-2.796645780225634) q[12];
rzz(0.012275082143518329) q[5],q[12];
u(pi,-0.5906194188748811,0.5906194188748811) q[11];
rzz(pi/128) q[5],q[11];
u(pi,0.5120796025351364,-0.5120796025351364) q[10];
rzz(0.049100328574073315) q[5],q[10];
u(pi,0.16901768476313084,-0.16901768476313084) q[9];
rzz(pi/32) q[5],q[9];
u(pi,3.490937756668978,-3.490937756668978) q[8];
rzz(0.19640131429629326) q[5],q[8];
u(pi,2.4900263372352702,-2.4900263372352702) q[7];
rzz(pi/8) q[5],q[7];
u(pi/2,3.0460882369206637,-3.0460882369206637) q[6];
rzz(0.7853830994606742) q[5],q[6];
u(pi/2,5*pi/8,-5*pi/8) q[4];
rzz(0) q[4],q[15];
u(pi,4.7004509283010485,-4.7004509283010485) q[4];
rzz(0) q[4],q[14];
u(pi,0.7495840071465247,-0.7495840071465247) q[4];
rzz(0) q[4],q[13];
u(pi,3.0825307117023053,-3.0825307117023053) q[4];
u(pi/2,-0.5353273881717007,0.5353273881717007) q[12];
rzz(0) q[4],q[12];
u(pi/2,-0.4636990756698536,0.4636990756698536) q[4];
u(pi,-1.4514158059584845,1.4514158059584845) q[11];
rzz(0.012275082143518329) q[4],q[11];
u(pi,0.8124158602183207,-0.8124158602183207) q[10];
rzz(pi/128) q[4],q[10];
u(pi,1.5205308443374599,-1.5205308443374599) q[9];
rzz(0.049100328574073315) q[4],q[9];
u(pi,-0.03204424506661585,0.03204424506661585) q[8];
rzz(pi/32) q[4],q[8];
u(pi,3.0517431036971248,-3.0517431036971248) q[7];
rzz(0.19640131429629326) q[4],q[7];
u(pi,0.11184069846779643,-0.11184069846779643) q[6];
rzz(pi/8) q[4],q[6];
u(pi/2,4.352362462283299,-4.352362462283299) q[5];
rzz(0.7853830994606742) q[4],q[5];
u(pi/2,pi/2,-pi/2) q[3];
rzz(0) q[3],q[15];
u(pi,-1.4017786420317657,1.4017786420317657) q[3];
rzz(0) q[3],q[14];
u(pi,2.077849381084289,-2.077849381084289) q[3];
rzz(0) q[3],q[13];
u(pi,-0.7257079029792421,0.7257079029792421) q[3];
rzz(0) q[3],q[12];
u(pi,2.7539201201368124,-2.7539201201368124) q[3];
u(pi/2,2.1639290197926497,-2.1639290197926497) q[11];
rzz(0) q[3],q[11];
u(pi/2,-0.21865484868984963,0.21865484868984963) q[3];
u(pi,1.0675131836898117,-1.0675131836898117) q[10];
rzz(0.012275082143518329) q[3],q[10];
u(pi,2.4121148394262435,-2.4121148394262435) q[9];
rzz(pi/128) q[3],q[9];
u(pi,1.7888228569540279,-1.7888228569540279) q[8];
rzz(0.049100328574073315) q[3],q[8];
u(pi,3.1654687577570755,-3.1654687577570755) q[7];
rzz(pi/32) q[3],q[7];
u(pi,4.550911117990174,-4.550911117990174) q[6];
rzz(0.19640131429629326) q[3],q[6];
u(pi,2.944928953475072,-2.944928953475072) q[5];
rzz(pi/8) q[3],q[5];
u(pi/2,-0.4636990756698536,0.4636990756698536) q[4];
rzz(0.7853830994606742) q[3],q[4];
u(pi/2,pi/2,-pi/2) q[2];
rzz(0) q[2],q[15];
u(pi,-1.4017786420317657,1.4017786420317657) q[2];
rzz(0) q[2],q[14];
u(pi,2.077849381084289,-2.077849381084289) q[2];
rzz(0) q[2],q[13];
u(pi,-0.7257079029792421,0.7257079029792421) q[2];
rzz(0) q[2],q[12];
u(pi,2.7539201201368124,-2.7539201201368124) q[2];
rzz(0) q[2],q[11];
u(pi,-0.049637163926718575,0.049637163926718575) q[2];
u(pi/2,-0.4121769561509807,0.4121769561509807) q[10];
rzz(0) q[2],q[10];
u(pi/2,3.260973174426205,-3.260973174426205) q[2];
u(pi,3.1038935417467153,-3.1038935417467153) q[9];
rzz(0.012275082143518329) q[2],q[9];
u(pi,3.7516899469169305,-3.7516899469169305) q[8];
rzz(pi/128) q[2],q[8];
u(pi,3.071220978149382,-3.071220978149382) q[7];
rzz(0.049100328574073315) q[2],q[7];
u(pi,1.2478406020058657,-1.2478406020058657) q[6];
rzz(pi/32) q[2],q[6];
u(pi,3.2377253887896407,-3.2377253887896407) q[5];
rzz(0.19640131429629326) q[2],q[5];
u(pi,4.412681041232224,-4.412681041232224) q[4];
rzz(pi/8) q[2],q[4];
u(pi/2,-0.21865484868984963,0.21865484868984963) q[3];
rzz(0.7853830994606742) q[2],q[3];
u(pi/2,-pi/2,pi/2) q[1];
rzz(0) q[1],q[15];
u(pi,2.2267608728644457,-2.2267608728644457) q[1];
rzz(0) q[1],q[14];
u(pi,0.39646899288303183,-0.39646899288303183) q[1];
rzz(0) q[1],q[13];
u(pi,-1.4331945685676637,1.4331945685676637) q[1];
rzz(0) q[1],q[12];
u(pi,3.020327177161227,-3.020327177161227) q[1];
rzz(0) q[1],q[11];
u(pi,1.1900352971798136,-1.1900352971798136) q[1];
rzz(0) q[1],q[10];
u(pi,-0.6396282642708819,0.6396282642708819) q[1];
u(pi/2,-1.2088848531013525,1.2088848531013525) q[9];
rzz(0) q[1],q[9];
u(pi/2,3.15792893538846,-3.15792893538846) q[1];
u(pi,0.2431592713878501,-0.2431592713878501) q[8];
rzz(0.012275082143518329) q[1],q[8];
u(pi,3.184946632209332,-3.184946632209332) q[7];
rzz(pi/128) q[1],q[7];
u(pi,2.9468139090672256,-2.9468139090672256) q[6];
rzz(0.049100328574073315) q[1],q[6];
u(pi,3.3897784732233864,-3.3897784732233864) q[5];
rzz(pi/32) q[1],q[5];
u(pi,4.471742983119712,-4.471742983119712) q[4];
rzz(0.19640131429629326) q[1],q[4];
u(pi,-1.4979113772316133,1.4979113772316133) q[3];
rzz(pi/8) q[1],q[3];
u(pi/2,0.11938052083641226,-0.11938052083641226) q[2];
rzz(0.7853830994606742) q[1],q[2];
rzz(-pi/2) q[0],q[15];
rzz(pi/2) q[0],q[14];
rzz(-pi/2) q[0],q[13];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[10];
rzz(pi/2) q[0],q[9];
u(pi/2,3.0668227484343555,-3.0668227484343555) q[8];
rzz(pi/2) q[0],q[8];
u(pi/2,1.7052564923685396,-1.7052564923685396) q[7];
rzz(pi/2) q[0],q[7];
u(pi/2,-0.3851592593301085,0.3851592593301085) q[6];
rzz(pi/2) q[0],q[6];
u(pi/2,1.8428582505957727,-1.8428582505957727) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,2.796017461694915,-2.796017461694915) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3.5060174014062095,-3.5060174014062095) q[3];
rzz(pi/2) q[0],q[3];
u(pi/2,3.260973174426205,-3.260973174426205) q[2];
rzz(pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u(pi,-pi/2,pi/2) q[0];
u(pi/2,3.15792893538846,-3.15792893538846) q[1];
rzz(-pi/2) q[0],q[1];
u(pi/2,-1.4514158059584845,1.4514158059584845) q[2];
u(pi/2,2.4755750110287567,-2.4755750110287567) q[2];
rzz(pi/2) q[0],q[2];
u(pi/2,-1.2063715789784806,1.2063715789784806) q[3];
u(pi/2,2.327920156310037,-2.327920156310037) q[3];
rzz(pi/2) q[0],q[3];
u(pi/2,4.366813788489812,-4.366813788489812) q[4];
u(pi/2,1.4212565164840223,-1.4212565164840223) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,0.2720619238008761,-0.2720619238008761) q[5];
u(pi/2,3.5123005867133887,-3.5123005867133887) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,1.1856370674647878,-1.1856370674647878) q[6];
u(pi/2,4.376238566450582,-4.376238566450582) q[6];
rzz(-pi/2) q[0],q[6];
u(pi/2,3.2760528191634357,-3.2760528191634357) q[7];
u(pi/2,0.1589645882716435,-0.1589645882716435) q[7];
rzz(pi/2) q[0],q[7];
u(pi/2,4.637619075229252,-4.637619075229252) q[8];
u(pi/2,1.5085927922538187,-1.5085927922538187) q[8];
rzz(pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[11];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[15];
u(pi/2,1.587132608593564,-1.587132608593564) q[1];
u(pi/2,4.046371337823653,-4.046371337823653) q[2];
rzz(0.7853830994606742) q[1],q[2];
u(pi/2,3.8987164831049332,-3.8987164831049332) q[3];
rzz(pi/8) q[1],q[3];
u(pi/2,2.992052843278919,-2.992052843278919) q[4];
rzz(0.19640131429629326) q[1],q[4];
u(pi/2,-1.200088393671301,1.200088393671301) q[5];
rzz(pi/32) q[1],q[5];
u(pi/2,2.8054422396556857,-2.8054422396556857) q[6];
rzz(0.049100328574073315) q[1],q[6];
u(pi/2,1.72976091506654,-1.72976091506654) q[7];
rzz(pi/128) q[1],q[7];
u(pi/2,3.079389119048715,-3.079389119048715) q[8];
rzz(0.012275082143518329) q[1],q[8];
u(pi/2,-1.5544600449962296,1.5544600449962296) q[1];
rzz(0) q[1],q[9];
u(pi,2.242468836132394,-2.242468836132394) q[1];
rzz(0) q[1],q[10];
u(pi,0.4128052746816986,-0.4128052746816986) q[1];
rzz(0) q[1],q[11];
u(pi,-1.4174866052997146,1.4174866052997146) q[1];
rzz(0) q[1],q[12];
u(pi,3.036035140429176,-3.036035140429176) q[1];
rzz(0) q[1],q[13];
u(pi,1.2063715789784806,-1.2063715789784806) q[1];
rzz(0) q[1],q[14];
u(pi,-0.6239203010029329,0.6239203010029329) q[1];
rzz(0) q[1],q[15];
u(pi/2,3.260973174426205,-3.260973174426205) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u(pi/2,-1.2214512237157116,1.2214512237157116) q[9];
rzz(0.012275082143518329) q[2],q[9];
u(pi/2,1.6901768476313088,-1.6901768476313088) q[2];
rzz(0) q[2],q[10];
u(pi,4.513212006147097,-4.513212006147097) q[2];
rzz(0) q[2],q[11];
u(pi,0.7351326809400116,-0.7351326809400116) q[2];
rzz(0) q[2],q[12];
u(pi,3.240238662912513,-3.240238662912513) q[2];
rzz(0) q[2],q[13];
u(pi,-0.5378406622945726,0.5378406622945726) q[2];
rzz(0) q[2],q[14];
u(pi,1.9666370011472103,-1.9666370011472103) q[2];
rzz(0) q[2],q[15];
u(pi/2,0.364424747816416,-0.364424747816416) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u(pi/2,2.7199909194780423,-2.7199909194780423) q[10];
rzz(0.012275082143518329) q[3],q[10];
u(pi/2,1.935221074611313,-1.935221074611313) q[3];
rzz(0) q[3],q[11];
u(pi,4.67217659441874,-4.67217659441874) q[3];
rzz(0) q[3],q[12];
u(pi,0.7213096732642166,-0.7213096732642166) q[3];
rzz(0) q[3],q[13];
u(pi,3.054256377819997,-3.054256377819997) q[3];
rzz(0) q[3],q[14];
u(pi,-0.8966105433345269,0.8966105433345269) q[3];
rzz(0) q[3],q[15];
u(pi/2,2.795389143164198,-2.795389143164198) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u(pi/2,2.1475927379939823,-2.1475927379939823) q[11];
rzz(0.012275082143518329) q[4],q[11];
u(pi/2,1.2245928163693014,-1.2245928163693014) q[4];
rzz(0) q[4],q[12];
u(pi,3.3803536952626176,-3.3803536952626176) q[4];
rzz(0) q[4],q[13];
u(pi,4.549654480928738,-4.549654480928738) q[4];
rzz(0) q[4],q[14];
u(pi,-0.5636017220540088,0.5636017220540088) q[4];
rzz(0) q[4],q[15];
u(pi/2,1.8434865691264908,-1.8434865691264908) q[5];
rzz(0.7853830994606742) q[5],q[6];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u(pi/2,-0.5516636699703676,0.5516636699703676) q[12];
rzz(0.012275082143518329) q[5],q[12];
u(pi/2,3.414282895921387,-3.414282895921387) q[5];
rzz(0) q[5],q[13];
u(pi,-0.04586725274241088,0.04586725274241088) q[5];
rzz(0) q[5],q[14];
u(pi,2.4592387292300897,-2.4592387292300897) q[5];
rzz(0) q[5],q[15];
u(pi/2,2.756433394259684,-2.756433394259684) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u(pi/2,2.439132536247115,-2.439132536247115) q[13];
rzz(0.012275082143518329) q[6],q[13];
u(pi/2,1.1856370674647878,-1.1856370674647878) q[6];
rzz(0) q[6],q[14];
u(pi,3.5807873065616462,-3.5807873065616462) q[6];
rzz(0) q[6],q[15];
u(pi/2,1.7052564923685396,-1.7052564923685396) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u(pi/2,0.29279643531456867,-0.29279643531456867) q[14];
rzz(0.012275082143518329) q[7],q[14];
u(pi/2,3.2760528191634357,-3.2760528191634357) q[7];
rzz(0) q[7],q[15];
u(pi/2,-0.07414158662471926,0.07414158662471926) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u(pi/2,-0.20546015954477248,0.20546015954477248) q[15];
rzz(0.012275082143518329) q[8],q[15];
u(pi/2,1.907575059259722,-1.907575059259722) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u(pi/2,-0.43102651207251963,0.43102651207251963) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
u(pi/2,2.1394245970946493,-2.1394245970946493) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
u(pi/2,2.583017479781528,-2.583017479781528) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
u(pi/2,2.4322210324092177,-2.4322210324092177) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
u(pi/2,0.28651325000738903,-0.28651325000738903) q[14];
rzz(0.7853830994606742) q[14],q[15];
u(pi/2,3.173636898656409,-3.173636898656409) q[1];
u(pi/2,-1.4928848289858696,1.4928848289858696) q[2];
u(pi/2,1.8403449764729012,-1.8403449764729012) q[3];
u(pi/2,1.5915308383085889,-1.5915308383085889) q[4];
u(pi/2,-1.00028310090299,1.00028310090299) q[5];
u(pi/2,-0.3066194429903639,0.3066194429903639) q[6];
u(pi/2,0.13446016557364304,-0.13446016557364304) q[7];
u(pi/2,2.9298493087378414,-2.9298493087378414) q[15];
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