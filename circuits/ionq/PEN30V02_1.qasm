OPENQASM 2.0;
include "qelib1.inc";
qreg q[30];
creg c[30];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[28];
rzz(0.36160994363065263) q[28],q[29];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[27];
rzz(0.7233485040645432) q[27],q[29];
u3(pi/2,0.7728317927830891,-0.7728317927830891) q[26];
rzz(1.446741257461612) q[26],q[29];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[25];
u3(pi,2.235557332294497,-2.235557332294497) q[29];
rzz(0.24819563711063614) q[25],q[29];
u3(pi/2,4.658981905273664,-4.658981905273664) q[24];
rzz(0.4964797728863239) q[24],q[29];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[23];
rzz(0.9927825484425445) q[23],q[29];
u3(pi/2,1.3540264336972008,-1.3540264336972008) q[22];
u3(pi,5.809433135018245,-5.809433135018245) q[29];
rzz(1.1560941741769235) q[22],q[29];
u3(pi/2,3.669380219392878,-3.669380219392878) q[21];
u3(pi,2.5151590784639883,-2.5151590784639883) q[29];
rzz(0.8295263718185014) q[21],q[29];
u3(pi/2,0.7012034802812419,-0.7012034802812419) q[20];
u3(pi,5.164778322501619,-5.164778322501619) q[29];
rzz(1.4825368154840266) q[20],q[29];
u3(pi/2,3.669380219392878,-3.669380219392878) q[19];
u3(pi,4.167008495721501,-4.167008495721501) q[29];
rzz(0.17659848920791243) q[19],q[29];
u3(pi/2,4.372468655266274,-4.372468655266274) q[18];
rzz(0.3532181213343836) q[18],q[29];
u3(pi/2,3.669380219392878,-3.669380219392878) q[17];
rzz(0.7063939568316497) q[17],q[29];
u3(pi/2,0.20860175219836227,-0.20860175219836227) q[16];
rzz(1.4128724853375343) q[16],q[29];
u3(pi/2,3.669380219392878,-3.669380219392878) q[15];
u3(pi,1.2503538761287378,-1.2503538761287378) q[29];
rzz(0.3159692154265287) q[15],q[29];
u3(pi/2,2.4020617429347557,-2.4020617429347557) q[14];
rzz(0.6319886806275516) q[14],q[29];
u3(pi/2,3.669380219392878,-3.669380219392878) q[13];
rzz(1.263876861706115) q[13],q[29];
u3(pi/2,4.892716398700744,-4.892716398700744) q[12];
u3(pi,1.9126016075054662,-1.9126016075054662) q[29];
rzz(0.6135923151542565) q[12],q[29];
u3(pi/2,3.669380219392878,-3.669380219392878) q[11];
rzz(1.2272305132692187) q[11],q[29];
u3(pi/2,5.433070335118188,-5.433070335118188) q[10];
u3(pi,5.568158819222549,-5.568158819222549) q[29];
rzz(0.6871447903245744) q[10],q[29];
u3(pi/2,3.6756634047000576,-3.6756634047000576) q[9];
rzz(1.3744046257721234) q[9],q[29];
u3(pi/2,1.3131857292005336,-1.3131857292005336) q[8];
u3(pi,0.49197340955216157,-0.49197340955216157) q[29];
rzz(pi/8) q[8],q[29];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[7];
rzz(0.7853830994606742) q[7],q[29];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
u3(pi,0.925513195747553,-0.925513195747553) q[29];
rzz(-pi/2) q[6],q[29];
u3(pi/2,0,0) q[1];
u3(pi/2,3*pi/4,-3*pi/4) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(pi/8) q[0],q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(0.19640131429629326) q[0],q[3];
u3(pi/2,3.730955435403238,-3.730955435403238) q[4];
rzz(pi/32) q[0],q[4];
u3(pi/2,3.82897312619524,-3.82897312619524) q[5];
rzz(0.049100328574073315) q[0],q[5];
rzz(pi/128) q[0],q[6];
rzz(0.012275082143518329) q[0],q[7];
u3(pi/2,0,0) q[0];
u3(pi/2,1.3131857292005336,-1.3131857292005336) q[8];
rzz(0) q[0],q[8];
u3(pi,3.568220935947287,-3.568220935947287) q[0];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
rzz(0) q[0],q[9];
u3(pi,1.2805131656031998,-1.2805131656031998) q[0];
u3(pi/2,5.433070335118188,-5.433070335118188) q[10];
rzz(0) q[0],q[10];
u3(pi,5.2759907024386985,-5.2759907024386985) q[0];
u3(pi/2,0.5284158843338032,-0.5284158843338032) q[11];
rzz(0) q[0],q[11];
u3(pi,2.988282932094611,-2.988282932094611) q[0];
u3(pi/2,4.892716398700744,-4.892716398700744) q[12];
rzz(0) q[0],q[12];
u3(pi,0.7005751617505239,-0.7005751617505239) q[0];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[13];
rzz(0) q[0],q[13];
u3(pi,4.696052698586023,-4.696052698586023) q[0];
u3(pi/2,5.543654396524548,-5.543654396524548) q[14];
rzz(0) q[0],q[14];
u3(pi,2.4083449282419354,-2.4083449282419354) q[0];
u3(pi/2,3.669380219392878,-3.669380219392878) q[15];
rzz(0) q[0],q[15];
u3(pi,0.12063715789784804,-0.12063715789784804) q[0];
u3(pi/2,3.3501944057881556,-3.3501944057881556) q[16];
rzz(0) q[0],q[16];
u3(pi,4.116114694733347,-4.116114694733347) q[0];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[17];
rzz(0) q[0],q[17];
u3(pi,1.8284069243892596,-1.8284069243892596) q[0];
u3(pi/2,1.2308760016764808,-1.2308760016764808) q[18];
rzz(0) q[0],q[18];
u3(pi,5.8238844612247584,-5.8238844612247584) q[0];
u3(pi/2,3.669380219392878,-3.669380219392878) q[19];
rzz(0) q[0],q[19];
u3(pi,3.536176690880671,-3.536176690880671) q[0];
u3(pi/2,0.7012034802812419,-0.7012034802812419) q[20];
rzz(0) q[0],q[20];
u3(pi,1.2484689205365838,-1.2484689205365838) q[0];
u3(pi/2,3.669380219392878,-3.669380219392878) q[21];
rzz(0) q[0],q[21];
u3(pi,5.243946457372083,-5.243946457372083) q[0];
u3(pi/2,1.3540264336972008,-1.3540264336972008) q[22];
rzz(0) q[0],q[22];
u3(pi/2,2.528982086139784,-2.528982086139784) q[0];
u3(pi/2,3.666238626739289,-3.666238626739289) q[23];
rzz(pi/2) q[0],q[23];
rzz(pi/2) q[0],q[23];
u3(pi/2,2.528982086139784,-2.528982086139784) q[0];
u3(pi/2,1.51738925168387,-1.51738925168387) q[24];
rzz(0) q[0],q[24];
u3(pi/2,5.6705747397295765,-5.6705747397295765) q[0];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[25];
rzz(pi/2) q[0],q[25];
rzz(-pi/2) q[0],q[25];
u3(pi/2,5*pi/4,-5*pi/4) q[1];
u3(pi/2,0,0) q[1];
rzz(0.7853830994606742) q[1],q[2];
rzz(pi/8) q[1],q[3];
rzz(0.19640131429629326) q[1],q[4];
rzz(pi/32) q[1],q[5];
rzz(0.049100328574073315) q[1],q[6];
rzz(pi/128) q[1],q[7];
u3(pi/2,4.454778382790327,-4.454778382790327) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(0) q[1],q[9];
u3(pi,4.172663362497963,-4.172663362497963) q[1];
rzz(0) q[1],q[10];
u3(pi,6.235433098845021,-6.235433098845021) q[1];
rzz(0) q[1],q[11];
u3(pi,2.015017528012493,-2.015017528012493) q[1];
rzz(0) q[1],q[12];
u3(pi,4.077787264359552,-4.077787264359552) q[1];
rzz(0) q[1],q[13];
u3(pi,6.139928682175891,-6.139928682175891) q[1];
rzz(0) q[1],q[14];
u3(pi,1.9195131113433637,-1.9195131113433637) q[1];
rzz(0) q[1],q[15];
u3(pi,3.982282847690422,-3.982282847690422) q[1];
rzz(0) q[1],q[16];
u3(pi,6.04505258403748,-6.04505258403748) q[1];
rzz(0) q[1],q[17];
u3(pi,1.824008694674234,-1.824008694674234) q[1];
rzz(0) q[1],q[18];
u3(pi,3.8867784310212925,-3.8867784310212925) q[1];
rzz(0) q[1],q[19];
u3(pi,5.94954816736835,-5.94954816736835) q[1];
rzz(0) q[1],q[20];
u3(pi,1.7291325965358222,-1.7291325965358222) q[1];
rzz(0) q[1],q[21];
u3(pi,3.7912740143521626,-3.7912740143521626) q[1];
rzz(0) q[1],q[22];
u3(pi/2,5*pi/8,-5*pi/8) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,3.6725218120464684,-3.6725218120464684) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi/2,pi/8,-pi/8) q[2];
rzz(0) q[2],q[10];
u3(pi,3.960920017646011,-3.960920017646011) q[2];
rzz(0) q[2],q[11];
u3(pi,1.6732122473019237,-1.6732122473019237) q[2];
rzz(0) q[2],q[12];
u3(pi/2,5.242061501779929,-5.242061501779929) q[2];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[14];
rzz(-pi/2) q[2],q[14];
rzz(pi/2) q[2],q[15];
rzz(pi/2) q[2],q[15];
rzz(pi/2) q[2],q[16];
rzz(pi/2) q[2],q[16];
rzz(-pi/2) q[2],q[17];
rzz(pi/2) q[2],q[17];
rzz(pi/2) q[2],q[18];
rzz(pi/2) q[2],q[18];
rzz(-pi/2) q[2],q[19];
rzz(pi/2) q[2],q[19];
rzz(pi/2) q[2],q[20];
rzz(pi/2) q[2],q[20];
rzz(pi/2) q[2],q[21];
rzz(-pi/2) q[2],q[21];
rzz(pi/2) q[2],q[22];
rzz(pi/2) q[2],q[22];
u3(pi/2,5.242061501779929,-5.242061501779929) q[2];
rzz(0) q[2],q[23];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,2.2870794518133692,-2.2870794518133692) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[3];
rzz(0) q[3],q[11];
u3(pi,4.1575837177607315,-4.1575837177607315) q[3];
rzz(0) q[3],q[12];
u3(pi,1.8698759474166446,-1.8698759474166446) q[3];
rzz(0) q[3],q[13];
u3(pi,5.865353484252144,-5.865353484252144) q[3];
rzz(0) q[3],q[14];
u3(pi,3.5776457139080566,-3.5776457139080566) q[3];
rzz(0) q[3],q[15];
u3(pi,1.289937943563969,-1.289937943563969) q[3];
rzz(0) q[3],q[16];
u3(pi,5.285415480399467,-5.285415480399467) q[3];
rzz(0) q[3],q[17];
u3(pi,2.997707710055381,-2.997707710055381) q[3];
rzz(0) q[3],q[18];
u3(pi,0.7099999397112933,-0.7099999397112933) q[3];
rzz(0) q[3],q[19];
u3(pi,4.705477476546792,-4.705477476546792) q[3];
rzz(0) q[3],q[20];
u3(pi,2.4177697062027046,-2.4177697062027046) q[3];
rzz(0) q[3],q[21];
u3(pi,0.1294336173278995,-0.1294336173278995) q[3];
rzz(0) q[3],q[22];
u3(pi/2,3.6982828718059046,-3.6982828718059046) q[3];
rzz(pi/2) q[3],q[23];
rzz(-pi/2) q[3],q[23];
u3(pi/2,0.5566902182161113,-0.5566902182161113) q[3];
rzz(0) q[3],q[24];
u3(pi/2,2.061513099285622,-2.061513099285622) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
rzz(pi/2) q[4],q[7];
u3(pi/2,4.466716434873968,-4.466716434873968) q[7];
u3(pi/2,1.1290883997001717,-1.1290883997001717) q[7];
rzz(pi/2) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,3.6630970340856988,-3.6630970340856988) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,3.6323094260805187,-3.6323094260805187) q[4];
rzz(0) q[4],q[12];
u3(pi,0.9179733733789376,-0.9179733733789376) q[4];
rzz(0) q[4],q[13];
u3(pi,4.913450910214436,-4.913450910214436) q[4];
rzz(0) q[4],q[14];
u3(pi/2,2.198486538982137,-2.198486538982137) q[4];
rzz(pi/2) q[4],q[15];
rzz(-pi/2) q[4],q[15];
rzz(pi/2) q[4],q[16];
rzz(pi/2) q[4],q[16];
rzz(pi/2) q[4],q[17];
rzz(-pi/2) q[4],q[17];
rzz(pi/2) q[4],q[18];
rzz(pi/2) q[4],q[18];
rzz(pi/2) q[4],q[19];
rzz(pi/2) q[4],q[19];
rzz(pi/2) q[4],q[20];
rzz(pi/2) q[4],q[20];
rzz(pi/2) q[4],q[21];
rzz(pi/2) q[4],q[21];
rzz(-pi/2) q[4],q[22];
rzz(pi/2) q[4],q[22];
u3(pi/2,5.34007919257193,-5.34007919257193) q[4];
rzz(0) q[4],q[23];
u3(pi/2,2.198486538982137,-2.198486538982137) q[4];
rzz(pi/2) q[4],q[24];
rzz(pi/2) q[4],q[24];
u3(pi/2,2.198486538982137,-2.198486538982137) q[4];
rzz(0) q[4],q[25];
u3(pi/2,2.2091679540043425,-2.2091679540043425) q[5];
rzz(0.7853830994606742) q[5],q[6];
u3(pi/2,5.841477380084861,-5.841477380084861) q[7];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,1.745468878334489,-1.745468878334489) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,0.638371627209446,-0.638371627209446) q[5];
rzz(0) q[5],q[13];
u3(pi,4.765167736964998,-4.765167736964998) q[5];
rzz(0) q[5],q[14];
u3(pi,3.5946103142374417,-3.5946103142374417) q[5];
rzz(0) q[5],q[15];
u3(pi,2.423424572979166,-2.423424572979166) q[5];
rzz(0) q[5],q[16];
u3(pi,1.2528671502516096,-1.2528671502516096) q[5];
rzz(0) q[5],q[17];
u3(pi,0.08230972752405258,-0.08230972752405258) q[5];
rzz(0) q[5],q[18];
u3(pi,5.194309293445364,-5.194309293445364) q[5];
rzz(0) q[5],q[19];
u3(pi,4.023751870717807,-4.023751870717807) q[5];
rzz(0) q[5],q[20];
u3(pi,2.85319444799025,-2.85319444799025) q[5];
rzz(0) q[5],q[21];
u3(pi,1.6820087067319751,-1.6820087067319751) q[5];
rzz(0) q[5],q[22];
u3(pi/2,5.809433135018245,-5.809433135018245) q[5];
rzz(-pi/2) q[5],q[23];
rzz(pi/2) q[5],q[23];
u3(pi/2,2.667840481428452,-2.667840481428452) q[5];
rzz(0) q[5],q[24];
u3(pi/2,5.809433135018245,-5.809433135018245) q[5];
rzz(pi/2) q[5],q[25];
rzz(pi/2) q[5],q[25];
u3(pi/2,5.809433135018245,-5.809433135018245) q[5];
u3(pi/2,3.914424446372882,-3.914424446372882) q[26];
rzz(0) q[5],q[26];
u3(pi/2,3.70582269417452,-3.70582269417452) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,0.5215043804959056,-0.5215043804959056) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[6];
rzz(0) q[6],q[14];
u3(pi,5.7038756218576285,-5.7038756218576285) q[6];
rzz(0) q[6],q[15];
u3(pi,3.4161678515135407,-3.4161678515135407) q[6];
rzz(0) q[6],q[16];
u3(pi/2,0.7012034802812419,-0.7012034802812419) q[6];
rzz(pi/2) q[6],q[17];
rzz(-pi/2) q[6],q[17];
rzz(pi/2) q[6],q[18];
rzz(pi/2) q[6],q[18];
rzz(pi/2) q[6],q[19];
rzz(-pi/2) q[6],q[19];
rzz(pi/2) q[6],q[20];
rzz(pi/2) q[6],q[20];
rzz(pi/2) q[6],q[21];
rzz(-pi/2) q[6],q[21];
rzz(pi/2) q[6],q[22];
rzz(pi/2) q[6],q[22];
u3(pi/2,3.842796133871035,-3.842796133871035) q[6];
rzz(0) q[6],q[23];
u3(pi/2,0.7012034802812419,-0.7012034802812419) q[6];
rzz(pi/2) q[6],q[24];
rzz(-pi/2) q[6],q[24];
u3(pi/2,3.842796133871035,-3.842796133871035) q[6];
rzz(0) q[6],q[25];
u3(pi/2,0.7012034802812419,-0.7012034802812419) q[6];
rzz(pi/2) q[6],q[26];
rzz(pi/2) q[6],q[26];
u3(pi/2,0.7012034802812419,-0.7012034802812419) q[6];
u3(pi/2,3.666238626739289,-3.666238626739289) q[27];
rzz(0) q[6],q[27];
u3(pi/2,5.043512846073054,-5.043512846073054) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,2.396406876158294,-2.396406876158294) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,3.472716519278157,-3.472716519278157) q[7];
rzz(0) q[7],q[15];
u3(pi,0.7583804665765761,-0.7583804665765761) q[7];
rzz(0) q[7],q[16];
u3(pi,4.753858003412075,-4.753858003412075) q[7];
rzz(0) q[7],q[17];
u3(pi,2.466150233067988,-2.466150233067988) q[7];
rzz(0) q[7],q[18];
u3(pi,0.17844246272390027,-0.17844246272390027) q[7];
rzz(0) q[7],q[19];
u3(pi,4.173919999559399,-4.173919999559399) q[7];
rzz(0) q[7],q[20];
u3(pi,1.8862122292153118,-1.8862122292153118) q[7];
rzz(0) q[7],q[21];
u3(pi,5.881689766050811,-5.881689766050811) q[7];
rzz(0) q[7],q[22];
u3(pi/2,3.1667253948185117,-3.1667253948185117) q[7];
rzz(pi/2) q[7],q[23];
rzz(pi/2) q[7],q[23];
u3(pi/2,3.1667253948185117,-3.1667253948185117) q[7];
rzz(0) q[7],q[24];
u3(pi/2,pi/125,-pi/125) q[7];
rzz(-pi/2) q[7],q[25];
rzz(pi/2) q[7],q[25];
u3(pi/2,3.1667253948185117,-3.1667253948185117) q[7];
rzz(0) q[7],q[26];
u3(pi/2,pi/125,-pi/125) q[7];
rzz(pi/2) q[7],q[27];
rzz(pi/2) q[7],q[27];
u3(pi/2,pi/125,-pi/125) q[7];
u3(pi/2,3.666238626739289,-3.666238626739289) q[28];
rzz(0) q[7],q[28];
u3(pi/2,0.12252211349000193,-0.12252211349000193) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,3.6630970340856988,-3.6630970340856988) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,1.6933184402848986,-1.6933184402848986) q[8];
rzz(0) q[8],q[16];
u3(pi,5.433070335118188,-5.433070335118188) q[8];
rzz(0) q[8],q[17];
u3(pi,3.4871678454846706,-3.4871678454846706) q[8];
rzz(0) q[8],q[18];
u3(pi/2,0.9437344331383739,-0.9437344331383739) q[8];
rzz(-pi/2) q[8],q[19];
rzz(pi/2) q[8],q[19];
rzz(pi/2) q[8],q[20];
rzz(pi/2) q[8],q[20];
rzz(pi/2) q[8],q[21];
rzz(-pi/2) q[8],q[21];
rzz(pi/2) q[8],q[22];
rzz(pi/2) q[8],q[22];
u3(pi/2,0.9437344331383739,-0.9437344331383739) q[8];
rzz(0) q[8],q[23];
u3(pi/2,4.085327086728167,-4.085327086728167) q[8];
rzz(pi/2) q[8],q[24];
rzz(-pi/2) q[8],q[24];
u3(pi/2,0.9437344331383739,-0.9437344331383739) q[8];
rzz(0) q[8],q[25];
u3(pi/2,4.085327086728167,-4.085327086728167) q[8];
rzz(pi/2) q[8],q[26];
rzz(pi/2) q[8],q[26];
u3(pi/2,4.085327086728167,-4.085327086728167) q[8];
rzz(0) q[8],q[27];
u3(pi/2,0.9437344331383739,-0.9437344331383739) q[8];
rzz(pi/2) q[8],q[28];
rzz(pi/2) q[8],q[28];
u3(pi/2,0.7181680806106266,-0.7181680806106266) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,0.20294688542190065,-0.20294688542190065) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,5.430557060995316,-5.430557060995316) q[9];
rzz(0) q[9],q[17];
u3(pi,2.715592689763017,-2.715592689763017) q[9];
rzz(0) q[9],q[18];
u3(pi,0.4278849194189298,-0.4278849194189298) q[9];
rzz(0) q[9],q[19];
u3(pi,4.423362456254428,-4.423362456254428) q[9];
rzz(0) q[9],q[20];
u3(pi,2.135654685910341,-2.135654685910341) q[9];
rzz(0) q[9],q[21];
u3(pi,6.1311322227458405,-6.1311322227458405) q[9];
rzz(0) q[9],q[22];
u3(pi/2,3.4161678515135407,-3.4161678515135407) q[9];
rzz(pi/2) q[9],q[23];
rzz(pi/2) q[9],q[23];
u3(pi/2,3.4161678515135407,-3.4161678515135407) q[9];
rzz(0) q[9],q[24];
u3(pi/2,0.27457519792374796,-0.27457519792374796) q[9];
rzz(pi/2) q[9],q[25];
rzz(-pi/2) q[9],q[25];
u3(pi/2,3.4161678515135407,-3.4161678515135407) q[9];
rzz(0) q[9],q[26];
u3(pi/2,0.27457519792374796,-0.27457519792374796) q[9];
rzz(pi/2) q[9],q[27];
rzz(pi/2) q[9],q[27];
u3(pi/2,0.27457519792374796,-0.27457519792374796) q[9];
rzz(0) q[9],q[28];
u3(pi/2,3.1629554836342035,-3.1629554836342035) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,3.6630970340856988,-3.6630970340856988) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi/2,4.7337518104291,-4.7337518104291) q[10];
rzz(0) q[10],q[18];
u3(pi,1.4690087248185872,-1.4690087248185872) q[10];
rzz(0) q[10],q[19];
u3(pi,4.363672195836223,-4.363672195836223) q[10];
rzz(0) q[10],q[20];
u3(pi/2,1.0989291102257097,-1.0989291102257097) q[10];
rzz(pi/2) q[10],q[21];
rzz(pi/2) q[10],q[21];
rzz(pi/2) q[10],q[22];
rzz(pi/2) q[10],q[22];
u3(pi/2,1.0989291102257097,-1.0989291102257097) q[10];
rzz(0) q[10],q[23];
u3(pi/2,4.2405217638155035,-4.2405217638155035) q[10];
rzz(pi/2) q[10],q[24];
rzz(pi/2) q[10],q[24];
u3(pi/2,4.2405217638155035,-4.2405217638155035) q[10];
rzz(0) q[10],q[25];
u3(pi/2,1.0989291102257097,-1.0989291102257097) q[10];
rzz(-pi/2) q[10],q[26];
rzz(pi/2) q[10],q[26];
u3(pi/2,4.2405217638155035,-4.2405217638155035) q[10];
rzz(0) q[10],q[27];
u3(pi/2,1.0989291102257097,-1.0989291102257097) q[10];
rzz(pi/2) q[10],q[28];
rzz(pi/2) q[10],q[28];
u3(pi/2,3.312495293945078,-3.312495293945078) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,1.2258494534307371,-1.2258494534307371) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,1.7416989671501812,-1.7416989671501812) q[11];
rzz(0) q[11],q[19];
u3(pi,4.910309317560847,-4.910309317560847) q[11];
rzz(0) q[11],q[20];
u3(pi,1.82212373908208,-1.82212373908208) q[11];
rzz(0) q[11],q[21];
u3(pi,5.017123467782899,-5.017123467782899) q[11];
rzz(0) q[11],q[22];
u3(pi/2,1.901920192483261,-1.901920192483261) q[11];
rzz(-pi/2) q[11],q[23];
rzz(pi/2) q[11],q[23];
u3(pi/2,5.043512846073054,-5.043512846073054) q[11];
rzz(0) q[11],q[24];
u3(pi/2,1.901920192483261,-1.901920192483261) q[11];
rzz(pi/2) q[11],q[25];
rzz(pi/2) q[11],q[25];
u3(pi/2,1.901920192483261,-1.901920192483261) q[11];
rzz(0) q[11],q[26];
u3(pi/2,5.043512846073054,-5.043512846073054) q[11];
rzz(pi/2) q[11],q[27];
rzz(-pi/2) q[11],q[27];
u3(pi/2,1.901920192483261,-1.901920192483261) q[11];
rzz(0) q[11],q[28];
u3(pi/2,3.923849224333652,-3.923849224333652) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/128) q[12],q[18];
u3(pi/2,3.6630970340856988,-3.6630970340856988) q[19];
rzz(0.012275082143518329) q[12],q[19];
u3(pi/2,5.494645551128548,-5.494645551128548) q[12];
rzz(0) q[12],q[20];
u3(pi,2.229274146987317,-2.229274146987317) q[12];
rzz(0) q[12],q[21];
u3(pi,5.123937618004953,-5.123937618004953) q[12];
rzz(0) q[12],q[22];
u3(pi,1.7360441003737197,-1.7360441003737197) q[12];
rzz(0) q[12],q[23];
u3(pi/2,4.753858003412075,-4.753858003412075) q[12];
rzz(pi/2) q[12],q[24];
rzz(pi/2) q[12],q[24];
u3(pi/2,4.753858003412075,-4.753858003412075) q[12];
rzz(0) q[12],q[25];
u3(pi/2,1.6122653498222819,-1.6122653498222819) q[12];
rzz(pi/2) q[12],q[26];
rzz(-pi/2) q[12],q[26];
u3(pi/2,4.753858003412075,-4.753858003412075) q[12];
rzz(0) q[12],q[27];
u3(pi/2,1.6122653498222819,-1.6122653498222819) q[12];
rzz(pi/2) q[12],q[28];
rzz(pi/2) q[12],q[28];
u3(pi/2,3.963433291768883,-3.963433291768883) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/128) q[13],q[19];
u3(pi/2,3.837769585625291,-3.837769585625291) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi/2,2.3926369649739865,-2.3926369649739865) q[13];
rzz(0) q[13],q[21];
u3(pi,4.648928808782176,-4.648928808782176) q[13];
rzz(0) q[13],q[22];
u3(pi,6.018663205747325,-6.018663205747325) q[13];
rzz(0) q[13],q[23];
u3(pi,1.1052122955328891,-1.1052122955328891) q[13];
rzz(0) q[13],q[24];
u3(pi/2,3.3615041393410787,-3.3615041393410787) q[13];
rzz(pi/2) q[13],q[25];
rzz(pi/2) q[13],q[25];
u3(pi/2,3.3615041393410787,-3.3615041393410787) q[13];
rzz(0) q[13],q[26];
u3(pi/2,0.21991148575128555,-0.21991148575128555) q[13];
rzz(-pi/2) q[13],q[27];
rzz(pi/2) q[13],q[27];
u3(pi/2,3.3615041393410787,-3.3615041393410787) q[13];
rzz(0) q[13],q[28];
u3(pi/2,0.18786724068466962,-0.18786724068466962) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/128) q[14],q[20];
u3(pi/2,0.5215043804959056,-0.5215043804959056) q[21];
rzz(0.012275082143518329) q[14],q[21];
u3(pi/2,1.758663567479566,-1.758663567479566) q[14];
rzz(0) q[14],q[22];
u3(pi,4.926645599359514,-4.926645599359514) q[14];
rzz(0) q[14],q[23];
u3(pi,1.8384600208807471,-1.8384600208807471) q[14];
rzz(0) q[14],q[24];
u3(pi,5.033459749581567,-5.033459749581567) q[14];
rzz(0) q[14],q[25];
u3(pi/2,1.9188847928126456,-1.9188847928126456) q[14];
rzz(pi/2) q[14],q[26];
rzz(-pi/2) q[14],q[26];
u3(pi/2,5.060477446402439,-5.060477446402439) q[14];
rzz(0) q[14],q[27];
u3(pi/2,1.9188847928126456,-1.9188847928126456) q[14];
rzz(pi/2) q[14],q[28];
rzz(pi/2) q[14],q[28];
u3(pi/2,4.911565954622282,-4.911565954622282) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/32) q[15],q[19];
rzz(0.049100328574073315) q[15],q[20];
rzz(pi/128) q[15],q[21];
u3(pi/2,4.49059253904125,-4.49059253904125) q[22];
rzz(0.012275082143518329) q[15],q[22];
u3(pi/2,3.3407696278273855,-3.3407696278273855) q[15];
rzz(0) q[15],q[23];
u3(pi,0.7816282522131405,-0.7816282522131405) q[15];
rzz(0) q[15],q[24];
u3(pi,5.086866824692593,-5.086866824692593) q[15];
rzz(0) q[15],q[25];
u3(pi,3.108920089992459,-3.108920089992459) q[15];
rzz(0) q[15],q[26];
u3(pi/2,0.5497787143782138,-0.5497787143782138) q[15];
rzz(pi/2) q[15],q[27];
rzz(pi/2) q[15],q[27];
u3(pi/2,0.5497787143782138,-0.5497787143782138) q[15];
rzz(0) q[15],q[28];
u3(pi/2,0.038955748904513435,-0.038955748904513435) q[16];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/8) q[16],q[18];
rzz(0.19640131429629326) q[16],q[19];
rzz(pi/32) q[16],q[20];
rzz(0.049100328574073315) q[16],q[21];
rzz(pi/128) q[16],q[22];
u3(pi/2,3.660583759962827,-3.660583759962827) q[23];
rzz(0.012275082143518329) q[16],q[23];
u3(pi/2,1.60975207569941,-1.60975207569941) q[16];
rzz(0) q[16],q[24];
u3(pi,5.178601330177416,-5.178601330177416) q[16];
rzz(0) q[16],q[25];
u3(pi,2.890893559833328,-2.890893559833328) q[16];
rzz(0) q[16],q[26];
u3(pi,0.6031857894892403,-0.6031857894892403) q[16];
rzz(0) q[16],q[27];
u3(pi/2,4.171406725436528,-4.171406725436528) q[16];
rzz(-pi/2) q[16],q[28];
rzz(pi/2) q[16],q[28];
u3(pi/2,2.792875869041326,-2.792875869041326) q[17];
rzz(0.7853830994606742) q[17],q[18];
rzz(pi/8) q[17],q[19];
rzz(0.19640131429629326) q[17],q[20];
rzz(pi/32) q[17],q[21];
rzz(0.049100328574073315) q[17],q[22];
rzz(pi/128) q[17],q[23];
u3(pi/2,1.5123627034381264,-1.5123627034381264) q[24];
rzz(0.012275082143518329) q[17],q[24];
u3(pi/2,1.2220795422464295,-1.2220795422464295) q[17];
rzz(0) q[17],q[25];
u3(pi,3.477743067523901,-3.477743067523901) q[17];
rzz(0) q[17],q[26];
u3(pi,4.84747746448905,-4.84747746448905) q[17];
rzz(0) q[17],q[27];
u3(pi,6.217840179984918,-6.217840179984918) q[17];
rzz(0) q[17],q[28];
u3(pi/2,0.0018849555921538756,-0.0018849555921538756) q[18];
rzz(0.7853830994606742) q[18],q[19];
rzz(pi/8) q[18],q[20];
rzz(0.19640131429629326) q[18],q[21];
rzz(pi/32) q[18],q[22];
rzz(0.049100328574073315) q[18],q[23];
rzz(pi/128) q[18],q[24];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[25];
rzz(0.012275082143518329) q[18],q[25];
u3(pi/2,1.5726812823870506,-1.5726812823870506) q[18];
rzz(0) q[18],q[26];
u3(pi,5.77424729729804,-5.77424729729804) q[18];
rzz(0) q[18],q[27];
u3(pi,4.752601366350639,-4.752601366350639) q[18];
rzz(0) q[18],q[28];
u3(pi/2,5.40479600123588,-5.40479600123588) q[19];
rzz(0.7853830994606742) q[19],q[20];
rzz(pi/8) q[19],q[21];
rzz(0.19640131429629326) q[19],q[22];
rzz(pi/32) q[19],q[23];
rzz(0.049100328574073315) q[19],q[24];
rzz(pi/128) q[19],q[25];
u3(pi/2,0.7678052445373454,-0.7678052445373454) q[26];
rzz(0.012275082143518329) q[19],q[26];
u3(pi/2,3.8339996744409834,-3.8339996744409834) q[19];
rzz(0) q[19],q[27];
u3(pi,0.7187963991413446,-0.7187963991413446) q[19];
rzz(0) q[19],q[28];
u3(pi/2,3.919450994618626,-3.919450994618626) q[20];
rzz(0.7853830994606742) q[20],q[21];
rzz(pi/8) q[20],q[22];
rzz(0.19640131429629326) q[20],q[23];
rzz(pi/32) q[20],q[24];
rzz(0.049100328574073315) q[20],q[25];
rzz(pi/128) q[20],q[26];
u3(pi/2,3.660583759962827,-3.660583759962827) q[27];
rzz(0.012275082143518329) q[20],q[27];
u3(pi/2,5.490247321413523,-5.490247321413523) q[20];
rzz(0) q[20],q[28];
u3(pi/2,2.916026301062046,-2.916026301062046) q[21];
rzz(0.7853830994606742) q[21],q[22];
rzz(pi/8) q[21],q[23];
rzz(0.19640131429629326) q[21],q[24];
rzz(pi/32) q[21],q[25];
rzz(0.049100328574073315) q[21],q[26];
rzz(pi/128) q[21],q[27];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[28];
rzz(0.012275082143518329) q[21],q[28];
u3(pi/2,4.898999584007923,-4.898999584007923) q[22];
rzz(0.7853830994606742) q[22],q[23];
rzz(pi/8) q[22],q[24];
rzz(0.19640131429629326) q[22],q[25];
rzz(pi/32) q[22],q[26];
rzz(0.049100328574073315) q[22],q[27];
rzz(pi/128) q[22],q[28];
u3(pi/2,3.0768758449258433,-3.0768758449258433) q[23];
rzz(0.7853830994606742) q[23],q[24];
rzz(pi/8) q[23],q[25];
rzz(0.19640131429629326) q[23],q[26];
rzz(pi/32) q[23],q[27];
rzz(0.049100328574073315) q[23],q[28];
u3(pi/2,0.43165483060323756,-0.43165483060323756) q[24];
rzz(0.7853830994606742) q[24],q[25];
rzz(pi/8) q[24],q[26];
rzz(0.19640131429629326) q[24],q[27];
rzz(pi/32) q[24],q[28];
u3(pi/2,2.332318386025062,-2.332318386025062) q[25];
rzz(0.7853830994606742) q[25],q[26];
rzz(pi/8) q[25],q[27];
rzz(0.19640131429629326) q[25],q[28];
u3(pi/2,4.026893463371397,-4.026893463371397) q[26];
rzz(0.7853830994606742) q[26],q[27];
rzz(pi/8) q[26],q[28];
u3(pi/2,1.3603096190043804,-1.3603096190043804) q[27];
rzz(0.7853830994606742) q[27],q[28];
u3(pi/2,0.11058406140636072,-0.11058406140636072) q[1];
u3(pi/2,2.1004688481901357,-2.1004688481901357) q[2];
u3(pi/2,3.6982828718059046,-3.6982828718059046) q[3];
u3(pi/2,5.34007919257193,-5.34007919257193) q[4];
u3(pi/2,2.667840481428452,-2.667840481428452) q[5];
u3(pi/2,3.842796133871035,-3.842796133871035) q[6];
u3(pi/2,3.1667253948185117,-3.1667253948185117) q[7];
u3(pi/2,3.4161678515135407,-3.4161678515135407) q[9];
u3(pi/2,5.043512846073054,-5.043512846073054) q[11];
u3(pi/2,0.21991148575128555,-0.21991148575128555) q[13];
u3(pi/2,3.691371367968007,-3.691371367968007) q[15];
u3(pi/2,2.190318398082804,-2.190318398082804) q[17];
u3(pi/2,2.670353755551324,-2.670353755551324) q[18];
u3(pi/2,3.88740674955201,-3.88740674955201) q[19];
u3(pi/2,2.3486546678237294,-2.3486546678237294) q[20];
u3(pi/2,4.863813746287717,-4.863813746287717) q[28];
u3(pi,1.00342469355658,-1.00342469355658) q[29];
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
measure q[29] -> c[29];
