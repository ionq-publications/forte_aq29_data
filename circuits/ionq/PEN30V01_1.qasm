OPENQASM 2.0;
include "qelib1.inc";
qreg q[30];
creg c[30];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[28];
rzz(1.3811718205514782) q[28],q[29];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[27];
u3(pi,1.8466281617800804,-1.8466281617800804) q[29];
rzz(0.3791520472295837) q[27],q[29];
u3(pi/2,5.1829995598924405,-5.1829995598924405) q[26];
rzz(0.7585155707686435) q[26],q[29];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[25];
rzz(1.5169168975204395) q[25],q[29];
u3(pi/2,3.450725370703029,-3.450725370703029) q[24];
u3(pi,5.252742916802134,-5.252742916802134) q[29];
rzz(0.1078257430565089) q[24],q[29];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[23];
rzz(0.21562360447821716) q[23],q[29];
u3(pi/2,2.8035572840635314,-2.8035572840635314) q[22];
rzz(0.4311383684789507) q[22],q[29];
u3(pi/2,3.669380219392878,-3.669380219392878) q[21];
rzz(0.8622616416052009) q[21],q[29];
u3(pi/2,0.21676989309769573,-0.21676989309769573) q[20];
u3(pi,0.9745220411435538,-0.9745220411435538) q[29];
rzz(1.4167484723977906) q[20],q[29];
u3(pi/2,3.669380219392878,-3.669380219392878) q[19];
u3(pi,0.09801769079200154,-0.09801769079200154) q[29];
rzz(0.30784983204517896) q[19],q[29];
u3(pi/2,2.4347343065320897,-2.4347343065320897) q[18];
rzz(0.6155573813590769) q[18],q[29];
u3(pi/2,3.669380219392878,-3.669380219392878) q[17];
rzz(1.2311211872751302) q[17],q[29];
u3(pi/2,5.024034971620797,-5.024034971620797) q[16];
u3(pi,4.461061568097506,-4.461061568097506) q[29];
rzz(0.6790650695504762) q[16],q[29];
u3(pi/2,3.669380219392878,-3.669380219392878) q[15];
rzz(1.3582895906360852) q[15],q[29];
u3(pi/2,2.8154953361471726,-2.8154953361471726) q[14];
u3(pi,6.030601257830967,-6.030601257830967) q[29];
rzz(0.4252335565909694) q[14],q[29];
u3(pi/2,3.669380219392878,-3.669380219392878) q[13];
rzz(0.8506628815281473) q[13],q[29];
u3(pi/2,0.26326546437082465,-0.26326546437082465) q[12];
u3(pi,0.9500176184455534,-0.9500176184455534) q[29];
rzz(1.4404167130871481) q[12],q[29];
u3(pi/2,3.669380219392878,-3.669380219392878) q[11];
u3(pi,5.350760607594136,-5.350760607594136) q[29];
rzz(0.260840547541335) q[11],q[29];
u3(pi/2,5.764822519337271,-5.764822519337271) q[10];
rzz(0.5215144964242503) q[10],q[29];
u3(pi/2,3.6756634047000576,-3.6756634047000576) q[9];
rzz(1.0431061660720358) q[9],q[29];
u3(pi/2,2.6383095104847083,-2.6383095104847083) q[8];
u3(pi,3.069336022557228,-3.069336022557228) q[29];
rzz(1.0554307597157746) q[8],q[29];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[7];
u3(pi,0.7288494956328321,-0.7288494956328321) q[29];
rzz(1.0308539547230358) q[7],q[29];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
u3(pi,6.077096829104096,-6.077096829104096) q[29];
rzz(1.0800365773166698) q[6],q[29];
u3(pi/2,3.82897312619524,-3.82897312619524) q[5];
u3(pi,4.917849139929462,-4.917849139929462) q[29];
rzz(0.9818236679571741) q[5],q[29];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
u3(pi,3.25343335205759,-3.25343335205759) q[29];
rzz(1.1780281300577933) q[4],q[29];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
u3(pi,1.472778636002895,-1.472778636002895) q[29];
rzz(0.7853830994606742) q[3],q[29];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
u3(pi,5.167291596624492,-5.167291596624492) q[29];
rzz(-pi/2) q[2],q[29];
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
u3(pi/2,2.6383095104847083,-2.6383095104847083) q[8];
rzz(0) q[0],q[8];
u3(pi,3.568220935947287,-3.568220935947287) q[0];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
rzz(0) q[0],q[9];
u3(pi,1.2805131656031998,-1.2805131656031998) q[0];
u3(pi/2,2.623229865747477,-2.623229865747477) q[10];
rzz(0) q[0],q[10];
u3(pi,5.2759907024386985,-5.2759907024386985) q[0];
u3(pi/2,3.670008537923596,-3.670008537923596) q[11];
rzz(0) q[0],q[11];
u3(pi,2.988282932094611,-2.988282932094611) q[0];
u3(pi/2,0.26326546437082465,-0.26326546437082465) q[12];
rzz(0) q[0],q[12];
u3(pi,0.7005751617505239,-0.7005751617505239) q[0];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[13];
rzz(0) q[0],q[13];
u3(pi,4.696052698586023,-4.696052698586023) q[0];
u3(pi/2,2.8154953361471726,-2.8154953361471726) q[14];
rzz(0) q[0],q[14];
u3(pi,2.4083449282419354,-2.4083449282419354) q[0];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[15];
rzz(0) q[0],q[15];
u3(pi,0.12063715789784804,-0.12063715789784804) q[0];
u3(pi/2,5.024034971620797,-5.024034971620797) q[16];
rzz(0) q[0],q[16];
u3(pi,4.116114694733347,-4.116114694733347) q[0];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[17];
rzz(0) q[0],q[17];
u3(pi,1.8284069243892596,-1.8284069243892596) q[0];
u3(pi/2,5.576326960121882,-5.576326960121882) q[18];
rzz(0) q[0],q[18];
u3(pi,5.8238844612247584,-5.8238844612247584) q[0];
u3(pi/2,3.669380219392878,-3.669380219392878) q[19];
rzz(0) q[0],q[19];
u3(pi,3.536176690880671,-3.536176690880671) q[0];
u3(pi/2,0.21676989309769573,-0.21676989309769573) q[20];
rzz(0) q[0],q[20];
u3(pi,1.2484689205365838,-1.2484689205365838) q[0];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[21];
rzz(0) q[0],q[21];
u3(pi,5.243946457372083,-5.243946457372083) q[0];
u3(pi/2,5.945149937653325,-5.945149937653325) q[22];
rzz(0) q[0],q[22];
u3(pi/2,2.528982086139784,-2.528982086139784) q[0];
u3(pi/2,3.666238626739289,-3.666238626739289) q[23];
rzz(pi/2) q[0],q[23];
rzz(pi/2) q[0],q[23];
u3(pi/2,2.528982086139784,-2.528982086139784) q[0];
u3(pi/2,3.450725370703029,-3.450725370703029) q[24];
rzz(0) q[0],q[24];
u3(pi/2,5.6705747397295765,-5.6705747397295765) q[0];
u3(pi/2,3.666238626739289,-3.666238626739289) q[25];
rzz(pi/2) q[0],q[25];
rzz(-pi/2) q[0],q[25];
u3(pi/2,3*pi/4,-3*pi/4) q[1];
u3(pi/2,pi,-pi) q[1];
rzz(0.7853830994606742) q[1],q[2];
rzz(pi/8) q[1],q[3];
rzz(0.19640131429629326) q[1],q[4];
rzz(pi/32) q[1],q[5];
rzz(0.049100328574073315) q[1],q[6];
rzz(pi/128) q[1],q[7];
u3(pi/2,5.779902164074501,-5.779902164074501) q[8];
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
u3(pi/2,3*pi/8,-3*pi/8) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,3.6725218120464684,-3.6725218120464684) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi/2,15*pi/8,-15*pi/8) q[2];
rzz(0) q[2],q[10];
u3(pi,3.175521854248563,-3.175521854248563) q[2];
rzz(0) q[2],q[11];
u3(pi,0.8878140839044756,-0.8878140839044756) q[2];
rzz(0) q[2],q[12];
u3(pi/2,4.456663338382481,-4.456663338382481) q[2];
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
u3(pi/2,4.456663338382481,-4.456663338382481) q[2];
rzz(0) q[2],q[23];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,5.759795971091527,-5.759795971091527) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[3];
rzz(0) q[3],q[11];
u3(pi,4.323459809870274,-4.323459809870274) q[3];
rzz(0) q[3],q[12];
u3(pi,3.1529023871427166,-3.1529023871427166) q[3];
rzz(0) q[3],q[13];
u3(pi,1.9817166458844415,-1.9817166458844415) q[3];
rzz(0) q[3],q[14];
u3(pi,0.8111592231568845,-0.8111592231568845) q[3];
rzz(0) q[3],q[15];
u3(pi,5.923787107608914,-5.923787107608914) q[3];
rzz(0) q[3],q[16];
u3(pi,4.752601366350639,-4.752601366350639) q[3];
rzz(0) q[3],q[17];
u3(pi,3.5820439436230824,-3.5820439436230824) q[3];
rzz(0) q[3],q[18];
u3(pi,2.411486520895525,-2.411486520895525) q[3];
rzz(0) q[3],q[19];
u3(pi,1.2403007796372503,-1.2403007796372503) q[3];
rzz(0) q[3],q[20];
u3(pi,0.06974335690969341,-0.06974335690969341) q[3];
rzz(0) q[3],q[21];
u3(pi,5.1823712413617224,-5.1823712413617224) q[3];
rzz(0) q[3],q[22];
u3(pi/2,3.0259820439376885,-3.0259820439376885) q[3];
rzz(pi/2) q[3],q[23];
rzz(-pi/2) q[3],q[23];
u3(pi/2,6.167574697527482,-6.167574697527482) q[3];
rzz(0) q[3],q[24];
u3(pi/2,3.436274044496516,-3.436274044496516) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[7];
rzz(pi/2) q[4],q[7];
u3(pi/2,1.3251237812841747,-1.3251237812841747) q[7];
u3(pi/2,4.270681053289964,-4.270681053289964) q[7];
rzz(pi/2) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,0.5215043804959056,-0.5215043804959056) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,5.007070371291412,-5.007070371291412) q[4];
rzz(0) q[4],q[12];
u3(pi,0.8834158541894498,-0.8834158541894498) q[4];
rzz(0) q[4],q[13];
u3(pi,2.0608847807549044,-2.0608847807549044) q[4];
rzz(0) q[4],q[14];
u3(pi/2,4.221043889363246,-4.221043889363246) q[4];
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
u3(pi/2,1.079451235773453,-1.079451235773453) q[4];
rzz(0) q[4],q[23];
u3(pi/2,4.221043889363246,-4.221043889363246) q[4];
rzz(pi/2) q[4],q[24];
rzz(pi/2) q[4],q[24];
u3(pi/2,4.221043889363246,-4.221043889363246) q[4];
rzz(0) q[4],q[25];
u3(pi/2,1.2271060904921731,-1.2271060904921731) q[5];
rzz(0.7853830994606742) q[5],q[6];
u3(pi/2,2.6998847264950685,-2.6998847264950685) q[7];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,3.3992032511841566,-3.3992032511841566) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,5.9394950708768635,-5.9394950708768635) q[5];
rzz(0) q[5],q[13];
u3(pi,2.368760860806704,-2.368760860806704) q[5];
rzz(0) q[5],q[14];
u3(pi,4.652070401435766,-4.652070401435766) q[5];
rzz(0) q[5],q[15];
u3(pi,0.6521946348852411,-0.6521946348852411) q[5];
rzz(0) q[5],q[16];
u3(pi,2.9361324940450206,-2.9361324940450206) q[5];
rzz(0) q[5],q[17];
u3(pi,5.219442034674082,-5.219442034674082) q[5];
rzz(0) q[5],q[18];
u3(pi,1.2195662681235577,-1.2195662681235577) q[5];
rzz(0) q[5],q[19];
u3(pi,3.5028758087526195,-3.5028758087526195) q[5];
rzz(0) q[5],q[20];
u3(pi,5.7861853493816815,-5.7861853493816815) q[5];
rzz(0) q[5],q[21];
u3(pi,1.7863095828311564,-1.7863095828311564) q[5];
rzz(0) q[5],q[22];
u3(pi/2,4.498760679940584,-4.498760679940584) q[5];
rzz(-pi/2) q[5],q[23];
rzz(pi/2) q[5],q[23];
u3(pi/2,1.3571680263507906,-1.3571680263507906) q[5];
rzz(0) q[5],q[24];
u3(pi/2,4.498760679940584,-4.498760679940584) q[5];
rzz(pi/2) q[5],q[25];
rzz(pi/2) q[5],q[25];
u3(pi/2,4.498760679940584,-4.498760679940584) q[5];
u3(pi/2,2.041406906302648,-2.041406906302648) q[26];
rzz(0) q[5],q[26];
u3(pi/2,3.2151059216837945,-3.2151059216837945) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,0.5215043804959056,-0.5215043804959056) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,1.6443095948888977,-1.6443095948888977) q[6];
rzz(0) q[6],q[14];
u3(pi,4.356760691998325,-4.356760691998325) q[6];
rzz(0) q[6],q[15];
u3(pi,0.35688492544780054,-0.35688492544780054) q[6];
rzz(0) q[6],q[16];
u3(pi/2,3.069336022557228,-3.069336022557228) q[6];
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
u3(pi/2,6.210928676147021,-6.210928676147021) q[6];
rzz(0) q[6],q[23];
u3(pi/2,3.069336022557228,-3.069336022557228) q[6];
rzz(pi/2) q[6],q[24];
rzz(-pi/2) q[6],q[24];
u3(pi/2,6.210928676147021,-6.210928676147021) q[6];
rzz(0) q[6],q[25];
u3(pi/2,3.069336022557228,-3.069336022557228) q[6];
rzz(pi/2) q[6],q[26];
rzz(pi/2) q[6],q[26];
u3(pi/2,3.069336022557228,-3.069336022557228) q[6];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[27];
rzz(0) q[6],q[27];
u3(pi/2,3.2276722922981538,-3.2276722922981538) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,5.951433122960505,-5.951433122960505) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,1.6568759655032568,-1.6568759655032568) q[7];
rzz(0) q[7],q[15];
u3(pi,5.3802915785378795,-5.3802915785378795) q[7];
rzz(0) q[7],q[16];
u3(pi,3.4023448438377457,-3.4023448438377457) q[7];
rzz(0) q[7],q[17];
u3(pi,1.4250264276683302,-1.4250264276683302) q[7];
rzz(0) q[7],q[18];
u3(pi,5.730265000147783,-5.730265000147783) q[7];
rzz(0) q[7],q[19];
u3(pi,3.7523182654476486,-3.7523182654476486) q[7];
rzz(0) q[7],q[20];
u3(pi,1.7749998492782328,-1.7749998492782328) q[7];
rzz(0) q[7],q[21];
u3(pi,6.080238421757685,-6.080238421757685) q[7];
rzz(0) q[7],q[22];
u3(pi/2,3.5204687276127222,-3.5204687276127222) q[7];
rzz(pi/2) q[7],q[23];
rzz(pi/2) q[7],q[23];
u3(pi/2,3.5204687276127222,-3.5204687276127222) q[7];
rzz(0) q[7],q[24];
u3(pi/2,0.37887607402292905,-0.37887607402292905) q[7];
rzz(-pi/2) q[7],q[25];
rzz(pi/2) q[7],q[25];
u3(pi/2,3.5204687276127222,-3.5204687276127222) q[7];
rzz(0) q[7],q[26];
u3(pi/2,0.37887607402292905,-0.37887607402292905) q[7];
rzz(pi/2) q[7],q[27];
rzz(pi/2) q[7],q[27];
u3(pi/2,0.37887607402292905,-0.37887607402292905) q[7];
u3(pi/2,3.666238626739289,-3.666238626739289) q[28];
rzz(0) q[7],q[28];
u3(pi/2,2.110521944681623,-2.110521944681623) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,0.5215043804959056,-0.5215043804959056) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[8];
rzz(0) q[8],q[16];
u3(pi,0.41657518586600656,-0.41657518586600656) q[8];
rzz(0) q[8],q[17];
u3(pi,3.311238656883642,-3.311238656883642) q[8];
rzz(0) q[8],q[18];
u3(pi/2,0.04649557127312894,-0.04649557127312894) q[8];
rzz(-pi/2) q[8],q[19];
rzz(pi/2) q[8],q[19];
rzz(pi/2) q[8],q[20];
rzz(pi/2) q[8],q[20];
rzz(pi/2) q[8],q[21];
rzz(-pi/2) q[8],q[21];
rzz(pi/2) q[8],q[22];
rzz(pi/2) q[8],q[22];
u3(pi/2,0.04649557127312894,-0.04649557127312894) q[8];
rzz(0) q[8],q[23];
u3(pi/2,3.1880882248629216,-3.1880882248629216) q[8];
rzz(pi/2) q[8],q[24];
rzz(-pi/2) q[8],q[24];
u3(pi/2,0.04649557127312894,-0.04649557127312894) q[8];
rzz(0) q[8],q[25];
u3(pi/2,3.1880882248629216,-3.1880882248629216) q[8];
rzz(pi/2) q[8],q[26];
rzz(pi/2) q[8],q[26];
u3(pi/2,3.1880882248629216,-3.1880882248629216) q[8];
rzz(0) q[8],q[27];
u3(pi/2,0.04649557127312894,-0.04649557127312894) q[8];
rzz(pi/2) q[8],q[28];
rzz(pi/2) q[8],q[28];
u3(pi/2,1.049291946298991,-1.049291946298991) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,1.8774157697852605,-1.8774157697852605) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,5.761680926683681,-5.761680926683681) q[9];
rzz(0) q[9],q[17];
u3(pi,3.2019112325387176,-3.2019112325387176) q[9];
rzz(0) q[9],q[18];
u3(pi,1.2245928163693014,-1.2245928163693014) q[9];
rzz(0) q[9],q[19];
u3(pi,5.529831388848754,-5.529831388848754) q[9];
rzz(0) q[9],q[20];
u3(pi,3.5518846541486204,-3.5518846541486204) q[9];
rzz(0) q[9],q[21];
u3(pi,1.5745662379792043,-1.5745662379792043) q[9];
rzz(0) q[9],q[22];
u3(pi/2,5.297981851013827,-5.297981851013827) q[9];
rzz(pi/2) q[9],q[23];
rzz(pi/2) q[9],q[23];
u3(pi/2,5.297981851013827,-5.297981851013827) q[9];
rzz(0) q[9],q[24];
u3(pi/2,2.156389197424034,-2.156389197424034) q[9];
rzz(pi/2) q[9],q[25];
rzz(-pi/2) q[9],q[25];
u3(pi/2,5.297981851013827,-5.297981851013827) q[9];
rzz(0) q[9],q[26];
u3(pi/2,2.156389197424034,-2.156389197424034) q[9];
rzz(pi/2) q[9],q[27];
rzz(pi/2) q[9],q[27];
u3(pi/2,2.156389197424034,-2.156389197424034) q[9];
rzz(0) q[9],q[28];
u3(pi/2,3.659955441432109,-3.659955441432109) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,3.6630970340856988,-3.6630970340856988) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi/2,5.230751768227005,-5.230751768227005) q[10];
rzz(0) q[10],q[18];
u3(pi,2.1161768114580846,-2.1161768114580846) q[10];
rzz(0) q[10],q[19];
u3(pi,5.311176540158905,-5.311176540158905) q[10];
rzz(0) q[10],q[20];
u3(pi/2,2.195973264859265,-2.195973264859265) q[10];
rzz(pi/2) q[10],q[21];
rzz(pi/2) q[10],q[21];
rzz(pi/2) q[10],q[22];
rzz(pi/2) q[10],q[22];
u3(pi/2,2.195973264859265,-2.195973264859265) q[10];
rzz(0) q[10],q[23];
u3(pi/2,5.337565918449059,-5.337565918449059) q[10];
rzz(pi/2) q[10],q[24];
rzz(pi/2) q[10],q[24];
u3(pi/2,5.337565918449059,-5.337565918449059) q[10];
rzz(0) q[10],q[25];
u3(pi/2,2.195973264859265,-2.195973264859265) q[10];
rzz(-pi/2) q[10],q[26];
rzz(pi/2) q[10],q[26];
u3(pi/2,5.337565918449059,-5.337565918449059) q[10];
rzz(0) q[10],q[27];
u3(pi/2,2.195973264859265,-2.195973264859265) q[10];
rzz(pi/2) q[10],q[28];
rzz(pi/2) q[10],q[28];
u3(pi/2,1.8246370132049519,-1.8246370132049519) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,5.571300411876139,-5.571300411876139) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,0.25384068641005525,-0.25384068641005525) q[11];
rzz(0) q[11],q[19];
u3(pi,3.822689940888061,-3.822689940888061) q[11];
rzz(0) q[11],q[20];
u3(pi,1.534353852013255,-1.534353852013255) q[11];
rzz(0) q[11],q[21];
u3(pi,5.529831388848754,-5.529831388848754) q[11];
rzz(0) q[11],q[22];
u3(pi/2,2.8154953361471726,-2.8154953361471726) q[11];
rzz(-pi/2) q[11],q[23];
rzz(pi/2) q[11],q[23];
u3(pi/2,5.957087989736966,-5.957087989736966) q[11];
rzz(0) q[11],q[24];
u3(pi/2,2.8154953361471726,-2.8154953361471726) q[11];
rzz(pi/2) q[11],q[25];
rzz(pi/2) q[11],q[25];
u3(pi/2,2.8154953361471726,-2.8154953361471726) q[11];
rzz(0) q[11],q[26];
u3(pi/2,5.957087989736966,-5.957087989736966) q[11];
rzz(pi/2) q[11],q[27];
rzz(-pi/2) q[11],q[27];
u3(pi/2,2.8154953361471726,-2.8154953361471726) q[11];
rzz(0) q[11],q[28];
u3(pi/2,0.12063715789784804,-0.12063715789784804) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/128) q[12],q[18];
u3(pi/2,3.6630970340856988,-3.6630970340856988) q[19];
rzz(0.012275082143518329) q[12],q[19];
u3(pi/2,1.6914334846927446,-1.6914334846927446) q[12];
rzz(0) q[12],q[20];
u3(pi,5.26028273917075,-5.26028273917075) q[12];
rzz(0) q[12],q[21];
u3(pi,2.9725749688266623,-2.9725749688266623) q[12];
rzz(0) q[12],q[22];
u3(pi,0.6848671984825749,-0.6848671984825749) q[12];
rzz(0) q[12],q[23];
u3(pi/2,4.253088134429862,-4.253088134429862) q[12];
rzz(pi/2) q[12],q[24];
rzz(pi/2) q[12],q[24];
u3(pi/2,4.253088134429862,-4.253088134429862) q[12];
rzz(0) q[12],q[25];
u3(pi/2,1.1114954808400688,-1.1114954808400688) q[12];
rzz(pi/2) q[12],q[26];
rzz(-pi/2) q[12],q[26];
u3(pi/2,4.253088134429862,-4.253088134429862) q[12];
rzz(0) q[12],q[27];
u3(pi/2,1.1114954808400688,-1.1114954808400688) q[12];
rzz(pi/2) q[12],q[28];
rzz(pi/2) q[12],q[28];
u3(pi/2,4.3768668849812995,-4.3768668849812995) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/128) q[13],q[19];
u3(pi/2,3.353335998441745,-3.353335998441745) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi/2,2.8060705581864034,-2.8060705581864034) q[13];
rzz(0) q[13],q[21];
u3(pi,5.145300448049363,-5.145300448049363) q[13];
rzz(0) q[13],q[22];
u3(pi,0.39835394847518574,-0.39835394847518574) q[13];
rzz(0) q[13],q[23];
u3(pi,1.9345927560805947,-1.9345927560805947) q[13];
rzz(0) q[13],q[24];
u3(pi/2,4.273194327412837,-4.273194327412837) q[13];
rzz(pi/2) q[13],q[25];
rzz(pi/2) q[13],q[25];
u3(pi/2,4.273194327412837,-4.273194327412837) q[13];
rzz(0) q[13],q[26];
u3(pi/2,1.1316016738230434,-1.1316016738230434) q[13];
rzz(-pi/2) q[13],q[27];
rzz(pi/2) q[13],q[27];
u3(pi/2,4.273194327412837,-4.273194327412837) q[13];
rzz(0) q[13],q[28];
u3(pi/2,0.8080176305032948,-0.8080176305032948) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/128) q[14],q[20];
u3(pi/2,3.6630970340856988,-3.6630970340856988) q[21];
rzz(0.012275082143518329) q[14],q[21];
u3(pi/2,2.3788139572981915,-2.3788139572981915) q[14];
rzz(0) q[14],q[22];
u3(pi,6.117937533600763,-6.117937533600763) q[14];
rzz(0) q[14],q[23];
u3(pi,4.172663362497963,-4.172663362497963) q[14];
rzz(0) q[14],q[24];
u3(pi,2.2267608728644452,-2.2267608728644452) q[14];
rzz(0) q[14],q[25];
u3(pi/2,5.966512767697735,-5.966512767697735) q[14];
rzz(pi/2) q[14],q[26];
rzz(-pi/2) q[14],q[26];
u3(pi/2,2.824920114107942,-2.824920114107942) q[14];
rzz(0) q[14],q[27];
u3(pi/2,5.966512767697735,-5.966512767697735) q[14];
rzz(pi/2) q[14],q[28];
rzz(pi/2) q[14],q[28];
u3(pi/2,0.30284953180605606,-0.30284953180605606) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/32) q[15],q[19];
rzz(0.049100328574073315) q[15],q[20];
rzz(pi/128) q[15],q[21];
u3(pi/2,2.798530735817788,-2.798530735817788) q[22];
rzz(0.012275082143518329) q[15],q[22];
u3(pi/2,5.015238512190746,-5.015238512190746) q[15];
rzz(0) q[15],q[23];
u3(pi,1.0706547763434016,-1.0706547763434016) q[15];
rzz(0) q[15],q[24];
u3(pi,2.60689358394881,-2.60689358394881) q[15];
rzz(0) q[15],q[25];
u3(pi,4.143132391554219,-4.143132391554219) q[15];
rzz(0) q[15],q[26];
u3(pi/2,0.19917697423759287,-0.19917697423759287) q[15];
rzz(pi/2) q[15],q[27];
rzz(pi/2) q[15],q[27];
u3(pi/2,0.19917697423759287,-0.19917697423759287) q[15];
rzz(0) q[15],q[28];
u3(pi/2,4.121141242979091,-4.121141242979091) q[16];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/8) q[16],q[18];
rzz(0.19640131429629326) q[16],q[19];
rzz(pi/32) q[16],q[20];
rzz(0.049100328574073315) q[16],q[21];
rzz(pi/128) q[16],q[22];
u3(pi/2,3.660583759962827,-3.660583759962827) q[23];
rzz(0.012275082143518329) q[16],q[23];
u3(pi/2,5.691937569773987,-5.691937569773987) q[16];
rzz(0) q[16],q[24];
u3(pi,1.747353833926643,-1.747353833926643) q[16];
rzz(0) q[16],q[25];
u3(pi,3.2835926415320515,-3.2835926415320515) q[16];
rzz(0) q[16],q[26];
u3(pi,4.81983144913746,-4.81983144913746) q[16];
rzz(0) q[16],q[27];
u3(pi/2,0.8758760318208343,-0.8758760318208343) q[16];
rzz(-pi/2) q[16],q[28];
rzz(pi/2) q[16],q[28];
u3(pi/2,0.8551415203071416,-0.8551415203071416) q[17];
rzz(0.7853830994606742) q[17],q[18];
rzz(pi/8) q[17],q[19];
rzz(0.19640131429629326) q[17],q[20];
rzz(pi/32) q[17],q[21];
rzz(0.049100328574073315) q[17],q[22];
rzz(pi/128) q[17],q[23];
u3(pi/2,3.445070503926567,-3.445070503926567) q[24];
rzz(0.012275082143518329) q[17],q[24];
u3(pi/2,5.567530500691832,-5.567530500691832) q[17];
rzz(0) q[17],q[25];
u3(pi,3.4111413032677977,-3.4111413032677977) q[17];
rzz(0) q[17],q[26];
u3(pi,2.2405838805402403,-2.2405838805402403) q[17];
rzz(0) q[17],q[27];
u3(pi,1.0700264578126837,-1.0700264578126837) q[17];
rzz(0) q[17],q[28];
u3(pi/2,3.3784687396704633,-3.3784687396704633) q[18];
rzz(0.7853830994606742) q[18],q[19];
rzz(pi/8) q[18],q[20];
rzz(0.19640131429629326) q[18],q[21];
rzz(pi/32) q[18],q[22];
rzz(0.049100328574073315) q[18],q[23];
rzz(pi/128) q[18],q[24];
u3(pi/2,3.660583759962827,-3.660583759962827) q[25];
rzz(0.012275082143518329) q[18],q[25];
u3(pi/2,4.94926506646536,-4.94926506646536) q[18];
rzz(0) q[18],q[26];
u3(pi,0.8469733794078083,-0.8469733794078083) q[18];
rzz(0) q[18],q[27];
u3(pi,2.066539647531366,-2.066539647531366) q[18];
rzz(0) q[18],q[28];
u3(pi/2,4.920362414052334,-4.920362414052334) q[19];
rzz(0.7853830994606742) q[19],q[20];
rzz(pi/8) q[19],q[21];
rzz(0.19640131429629326) q[19],q[22];
rzz(pi/32) q[19],q[23];
rzz(0.049100328574073315) q[19],q[24];
rzz(pi/128) q[19],q[25];
u3(pi/2,5.1779730116466975,-5.1779730116466975) q[26];
rzz(0.012275082143518329) q[19],q[26];
u3(pi/2,3.3495660872574375,-3.3495660872574375) q[19];
rzz(0) q[19],q[27];
u3(pi,0.8055043563804231,-0.8055043563804231) q[19];
rzz(0) q[19],q[28];
u3(pi/2,0.05152211951887261,-0.05152211951887261) q[20];
rzz(0.7853830994606742) q[20],q[21];
rzz(pi/8) q[20],q[22];
rzz(0.19640131429629326) q[20],q[23];
rzz(pi/32) q[20],q[24];
rzz(0.049100328574073315) q[20],q[25];
rzz(pi/128) q[20],q[26];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[27];
rzz(0.012275082143518329) q[20],q[27];
u3(pi/2,1.622318446313769,-1.622318446313769) q[20];
rzz(0) q[20],q[28];
u3(pi/2,1.2239644978385833,-1.2239644978385833) q[21];
rzz(0.7853830994606742) q[21],q[22];
rzz(pi/8) q[21],q[23];
rzz(0.19640131429629326) q[21],q[24];
rzz(pi/32) q[21],q[25];
rzz(0.049100328574073315) q[21],q[26];
rzz(pi/128) q[21],q[27];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[28];
rzz(0.012275082143518329) q[21],q[28];
u3(pi/2,0.790424711643192,-0.790424711643192) q[22];
rzz(0.7853830994606742) q[22],q[23];
rzz(pi/8) q[22],q[24];
rzz(0.19640131429629326) q[22],q[25];
rzz(pi/32) q[22],q[26];
rzz(0.049100328574073315) q[22],q[27];
rzz(pi/128) q[22],q[28];
u3(pi/2,1.868619310355209,-1.868619310355209) q[23];
rzz(0.7853830994606742) q[23],q[24];
rzz(pi/8) q[23],q[25];
rzz(0.19640131429629326) q[23],q[26];
rzz(pi/32) q[23],q[27];
rzz(0.049100328574073315) q[23],q[28];
u3(pi/2,4.902141176661513,-4.902141176661513) q[24];
rzz(0.7853830994606742) q[24],q[25];
rzz(pi/8) q[24],q[26];
rzz(0.19640131429629326) q[24],q[27];
rzz(pi/32) q[24],q[28];
u3(pi/2,3.6008934995446213,-3.6008934995446213) q[25];
rzz(0.7853830994606742) q[25],q[26];
rzz(pi/8) q[25],q[27];
rzz(0.19640131429629326) q[25],q[28];
u3(pi/2,4.359273966121196,-4.359273966121196) q[26];
rzz(0.7853830994606742) q[26],q[27];
rzz(pi/8) q[26],q[28];
u3(pi/2,2.463008640414398,-2.463008640414398) q[27];
rzz(0.7853830994606742) q[27],q[28];
u3(pi/2,0.10430087609918114,-0.10430087609918114) q[1];
u3(pi/2,1.3150706847926874,-1.3150706847926874) q[2];
u3(pi/2,3.0259820439376885,-3.0259820439376885) q[3];
u3(pi/2,1.079451235773453,-1.079451235773453) q[4];
u3(pi/2,1.3571680263507906,-1.3571680263507906) q[5];
u3(pi/2,6.210928676147021,-6.210928676147021) q[6];
u3(pi/2,3.5204687276127222,-3.5204687276127222) q[7];
u3(pi/2,5.297981851013827,-5.297981851013827) q[9];
u3(pi/2,5.957087989736966,-5.957087989736966) q[11];
u3(pi/2,1.1316016738230434,-1.1316016738230434) q[13];
u3(pi/2,3.3407696278273855,-3.3407696278273855) q[15];
u3(pi/2,5.196822567568235,-5.196822567568235) q[17];
u3(pi/2,4.246804949122682,-4.246804949122682) q[18];
u3(pi/2,4.545256251213713,-4.545256251213713) q[19];
u3(pi/2,4.763911099903562,-4.763911099903562) q[20];
u3(pi/2,3.8440527709324708,-3.8440527709324708) q[28];
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
