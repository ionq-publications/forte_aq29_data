OPENQASM 2.0;
include "qelib1.inc";
qreg q[23];
creg c[23];
u3(pi/2,3.669380219392878,-3.669380219392878) q[21];
u3(pi,3.7787076437378033,-3.7787076437378033) q[22];
rzz(0.7710102345306846) q[21],q[22];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[20];
rzz(1.5420266265829703) q[20],q[22];
u3(pi/2,0.4128052746816988,-0.4128052746816988) q[19];
u3(pi,5.842105698615579,-5.842105698615579) q[22];
rzz(0.05754448980395116) q[19],q[22];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[18];
rzz(0.11510558292507657) q[18],q[22];
u3(pi/2,0.06723008278682156,-0.06723008278682156) q[17];
rzz(0.2302982350905474) q[17],q[22];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[16];
rzz(0.4605924332345348) q[16],q[22];
u3(pi/2,4.9687429409176165,-4.9687429409176165) q[15];
rzz(0.9211929403621896) q[15],q[22];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[14];
u3(pi,0.4128052746816988,-0.4128052746816988) q[22];
rzz(1.2993474062605521) q[14],q[22];
u3(pi/2,5.724610133371321,-5.724610133371321) q[13];
u3(pi,6.096574703556353,-6.096574703556353) q[22];
rzz(0.5430934052113746) q[13],q[22];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[12];
rzz(1.0860364066744588) q[12],q[22];
u3(pi/2,2.4648935960065517,-2.4648935960065517) q[11];
u3(pi,0.6697875537453439,-0.6697875537453439) q[22];
rzz(0.9694442221057035) q[11],q[22];
u3(pi/2,3.669380219392878,-3.669380219392878) q[10];
u3(pi,4.203450970503144,-4.203450970503144) q[22];
rzz(1.2026409377023426) q[10],q[22];
u3(pi/2,2.0005662018059804,-2.0005662018059804) q[9];
u3(pi,4.909052680499411,-4.909052680499411) q[22];
rzz(0.736411450521692) q[9],q[22];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[8];
rzz(1.472482132488249) q[8],q[22];
u3(pi/2,5.620309257272139,-5.620309257272139) q[7];
u3(pi,2.6477342884454775,-2.6477342884454775) q[22];
rzz(0.19640131429629326) q[7],q[22];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
rzz(pi/8) q[6],q[22];
u3(pi/2,5.301751762198135,-5.301751762198135) q[5];
rzz(0.7853830994606742) q[5],q[22];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
u3(pi,1.1278317626387357,-1.1278317626387357) q[22];
rzz(-pi/2) q[4],q[22];
u3(pi/2,pi/2,-pi/2) q[1];
u3(pi/2,5*pi/4,-5*pi/4) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(pi/8) q[0],q[2];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(0.19640131429629326) q[0],q[3];
rzz(pi/32) q[0],q[4];
rzz(0.049100328574073315) q[0],q[5];
rzz(pi/128) q[0],q[6];
rzz(0.012275082143518329) q[0],q[7];
u3(pi/2,0,0) q[0];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[8];
rzz(0) q[0],q[8];
u3(pi,3.568220935947287,-3.568220935947287) q[0];
u3(pi/2,2.0005662018059804,-2.0005662018059804) q[9];
rzz(0) q[0],q[9];
u3(pi,1.2805131656031998,-1.2805131656031998) q[0];
u3(pi/2,3.669380219392878,-3.669380219392878) q[10];
rzz(0) q[0],q[10];
u3(pi,5.2759907024386985,-5.2759907024386985) q[0];
u3(pi/2,2.466150233067988,-2.466150233067988) q[11];
rzz(0) q[0],q[11];
u3(pi,2.988282932094611,-2.988282932094611) q[0];
u3(pi/2,3.6668669452700065,-3.6668669452700065) q[12];
rzz(0) q[0],q[12];
u3(pi,0.7005751617505239,-0.7005751617505239) q[0];
u3(pi/2,5.724610133371321,-5.724610133371321) q[13];
rzz(0) q[0],q[13];
u3(pi,4.696052698586023,-4.696052698586023) q[0];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[14];
rzz(0) q[0],q[14];
u3(pi,2.4083449282419354,-2.4083449282419354) q[0];
u3(pi/2,1.8271502873278236,-1.8271502873278236) q[15];
rzz(0) q[0],q[15];
u3(pi,0.12063715789784804,-0.12063715789784804) q[0];
u3(pi/2,3.666238626739289,-3.666238626739289) q[16];
rzz(0) q[0],q[16];
u3(pi,4.116114694733347,-4.116114694733347) q[0];
u3(pi/2,3.208822736376615,-3.208822736376615) q[17];
rzz(0) q[0],q[17];
u3(pi,1.8284069243892596,-1.8284069243892596) q[0];
u3(pi/2,3.666238626739289,-3.666238626739289) q[18];
rzz(0) q[0],q[18];
u3(pi,5.8238844612247584,-5.8238844612247584) q[0];
u3(pi/2,0.4128052746816988,-0.4128052746816988) q[19];
rzz(0) q[0],q[19];
u3(pi,3.536176690880671,-3.536176690880671) q[0];
u3(pi/2,3.666238626739289,-3.666238626739289) q[20];
rzz(0) q[0],q[20];
u3(pi,1.2484689205365838,-1.2484689205365838) q[0];
u3(pi/2,3.669380219392878,-3.669380219392878) q[21];
rzz(0) q[0],q[21];
u3(pi/2,7*pi/4,-7*pi/4) q[1];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(0.7853830994606742) q[1],q[2];
rzz(pi/8) q[1],q[3];
rzz(0.19640131429629326) q[1],q[4];
rzz(pi/32) q[1],q[5];
rzz(0.049100328574073315) q[1],q[6];
rzz(pi/128) q[1],q[7];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi/2,0,0) q[1];
rzz(0) q[1],q[9];
u3(pi,3.568220935947287,-3.568220935947287) q[1];
rzz(0) q[1],q[10];
u3(pi,1.2805131656031998,-1.2805131656031998) q[1];
rzz(0) q[1],q[11];
u3(pi,5.2759907024386985,-5.2759907024386985) q[1];
rzz(0) q[1],q[12];
u3(pi,2.988282932094611,-2.988282932094611) q[1];
rzz(0) q[1],q[13];
u3(pi,0.7005751617505239,-0.7005751617505239) q[1];
rzz(0) q[1],q[14];
u3(pi,4.696052698586023,-4.696052698586023) q[1];
rzz(0) q[1],q[15];
u3(pi,2.4083449282419354,-2.4083449282419354) q[1];
rzz(0) q[1],q[16];
u3(pi,0.12063715789784804,-0.12063715789784804) q[1];
rzz(0) q[1],q[17];
u3(pi,4.116114694733347,-4.116114694733347) q[1];
rzz(0) q[1],q[18];
u3(pi,1.8284069243892596,-1.8284069243892596) q[1];
rzz(0) q[1],q[19];
u3(pi,5.8238844612247584,-5.8238844612247584) q[1];
rzz(0) q[1],q[20];
u3(pi,3.536176690880671,-3.536176690880671) q[1];
rzz(0) q[1],q[21];
u3(pi/2,5*pi/8,-5*pi/8) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,5.1390172627421835,-5.1390172627421835) q[9];
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
rzz(-pi/2) q[2],q[16];
rzz(pi/2) q[2],q[17];
rzz(pi/2) q[2],q[17];
rzz(pi/2) q[2],q[18];
rzz(pi/2) q[2],q[18];
rzz(-pi/2) q[2],q[19];
rzz(pi/2) q[2],q[19];
rzz(pi/2) q[2],q[20];
rzz(pi/2) q[2],q[20];
rzz(-pi/2) q[2],q[21];
rzz(pi/2) q[2],q[21];
u3(pi/2,5.301751762198135,-5.301751762198135) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,0.5233893360880595,-0.5233893360880595) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,0.5887344632827273,-0.5887344632827273) q[3];
rzz(0) q[3],q[11];
u3(pi,2.927964353145687,-2.927964353145687) q[3];
rzz(0) q[3],q[12];
u3(pi,4.464203160751096,-4.464203160751096) q[3];
rzz(0) q[3],q[13];
u3(pi,6.000441968356505,-6.000441968356505) q[3];
rzz(0) q[3],q[14];
u3(pi,1.2534954687823274,-1.2534954687823274) q[3];
rzz(0) q[3],q[15];
u3(pi,2.7897342763877364,-2.7897342763877364) q[3];
rzz(0) q[3],q[16];
u3(pi,4.325973083993145,-4.325973083993145) q[3];
rzz(0) q[3],q[17];
u3(pi,5.862211891598554,-5.862211891598554) q[3];
rzz(0) q[3],q[18];
u3(pi,1.1152653920243765,-1.1152653920243765) q[3];
rzz(0) q[3],q[19];
u3(pi,2.6515041996297852,-2.6515041996297852) q[3];
rzz(0) q[3],q[20];
u3(pi,4.187743007235194,-4.187743007235194) q[3];
rzz(0) q[3],q[21];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
u3(pi/2,5.620309257272139,-5.620309257272139) q[7];
rzz(pi/2) q[4],q[7];
u3(pi/2,4.049512930477243,-4.049512930477243) q[7];
u3(pi/2,0.7118848953034471,-0.7118848953034471) q[7];
rzz(pi/2) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,5.600831382819883,-5.600831382819883) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[4];
rzz(0) q[4],q[12];
u3(pi,4.783388974355819,-4.783388974355819) q[4];
rzz(0) q[4],q[13];
u3(pi,0.409663682028109,-0.409663682028109) q[4];
rzz(0) q[4],q[14];
u3(pi/2,2.9355041755143025,-2.9355041755143025) q[4];
rzz(pi/2) q[4],q[15];
rzz(-pi/2) q[4],q[15];
rzz(pi/2) q[4],q[16];
rzz(pi/2) q[4],q[16];
rzz(pi/2) q[4],q[17];
rzz(-pi/2) q[4],q[17];
rzz(pi/2) q[4],q[18];
rzz(pi/2) q[4],q[18];
rzz(pi/2) q[4],q[19];
rzz(-pi/2) q[4],q[19];
rzz(pi/2) q[4],q[20];
rzz(pi/2) q[4],q[20];
rzz(pi/2) q[4],q[21];
rzz(-pi/2) q[4],q[21];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[5];
rzz(0.7853830994606742) q[5],q[6];
u3(pi/2,5.424273875688137,-5.424273875688137) q[7];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,4.466716434873968,-4.466716434873968) q[5];
rzz(0) q[5],q[13];
u3(pi,0.6320884419022663,-0.6320884419022663) q[5];
rzz(0) q[5],q[14];
u3(pi,2.386353779666807,-2.386353779666807) q[5];
rzz(0) q[5],q[15];
u3(pi,4.1412474359620655,-4.1412474359620655) q[5];
rzz(0) q[5],q[16];
u3(pi,5.896141092257324,-5.896141092257324) q[5];
rzz(0) q[5],q[17];
u3(pi,1.367221122842278,-1.367221122842278) q[5];
rzz(0) q[5],q[18];
u3(pi,3.1221147791375365,-3.1221147791375365) q[5];
rzz(0) q[5],q[19];
u3(pi,4.876380116902077,-4.876380116902077) q[5];
rzz(0) q[5],q[20];
u3(pi,0.34808846601774907,-0.34808846601774907) q[5];
rzz(0) q[5],q[21];
u3(pi/2,1.7423272856808991,-1.7423272856808991) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,5.718326948064141,-5.718326948064141) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,0.1715309588860027,-0.1715309588860027) q[6];
rzz(0) q[6],q[14];
u3(pi,2.7740263131197875,-2.7740263131197875) q[6];
rzz(0) q[6],q[15];
u3(pi,4.836796049466845,-4.836796049466845) q[6];
rzz(0) q[6],q[16];
u3(pi/2,1.155477777990326,-1.155477777990326) q[6];
rzz(pi/2) q[6],q[17];
rzz(pi/2) q[6],q[17];
rzz(pi/2) q[6],q[18];
rzz(pi/2) q[6],q[18];
rzz(pi/2) q[6],q[19];
rzz(pi/2) q[6],q[19];
rzz(pi/2) q[6],q[20];
rzz(pi/2) q[6],q[20];
rzz(-pi/2) q[6],q[21];
rzz(pi/2) q[6],q[21];
u3(pi/2,0.5032831431050849,-0.5032831431050849) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,3.660583759962827,-3.660583759962827) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,2.0740794698999814,-2.0740794698999814) q[7];
rzz(0) q[7],q[15];
u3(pi,4.599291644855457,-4.599291644855457) q[7];
rzz(0) q[7],q[16];
u3(pi,0.22556635252774715,-0.22556635252774715) q[7];
rzz(0) q[7],q[17];
u3(pi,2.1350263673796235,-2.1350263673796235) q[7];
rzz(0) q[7],q[18];
u3(pi,4.0444863822315,-4.0444863822315) q[7];
rzz(0) q[7],q[19];
u3(pi,5.953946397083376,-5.953946397083376) q[7];
rzz(0) q[7],q[20];
u3(pi,1.580221104755666,-1.580221104755666) q[7];
rzz(0) q[7],q[21];
u3(pi/2,3.571362528600877,-3.571362528600877) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,1.8208671020206442,-1.8208671020206442) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,2.0005662018059804,-2.0005662018059804) q[8];
rzz(0) q[8],q[16];
u3(pi,4.3391677731382226,-4.3391677731382226) q[8];
rzz(0) q[8],q[17];
u3(pi,5.875406580743632,-5.875406580743632) q[8];
rzz(0) q[8],q[18];
u3(pi/2,1.9314511634270048,-1.9314511634270048) q[8];
rzz(pi/2) q[8],q[19];
rzz(pi/2) q[8],q[19];
rzz(-pi/2) q[8],q[20];
rzz(pi/2) q[8],q[20];
rzz(pi/2) q[8],q[21];
rzz(pi/2) q[8],q[21];
u3(pi/2,1.1535928223981722,-1.1535928223981722) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,3.660583759962827,-3.660583759962827) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,2.7243891491930685,-2.7243891491930685) q[9];
rzz(0) q[9],q[17];
u3(pi,0.16461945504810516,-0.16461945504810516) q[9];
rzz(0) q[9],q[18];
u3(pi,4.4704863460582756,-4.4704863460582756) q[9];
rzz(0) q[9],q[19];
u3(pi,2.4925396113581417,-2.4925396113581417) q[9];
rzz(0) q[9],q[20];
u3(pi,0.5145928766580081,-0.5145928766580081) q[9];
rzz(0) q[9],q[21];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,0.060946897479641986,-0.060946897479641986) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi/2,5.595804834574139,-5.595804834574139) q[10];
rzz(0) q[10],q[18];
u3(pi,2.881468781872558,-2.881468781872558) q[10];
rzz(0) q[10],q[19];
u3(pi,0.5937610115284709,-0.5937610115284709) q[10];
rzz(0) q[10],q[20];
u3(pi/2,4.161981947475758,-4.161981947475758) q[10];
rzz(pi/2) q[10],q[21];
rzz(-pi/2) q[10],q[21];
u3(pi/2,1.851026391495106,-1.851026391495106) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,3.660583759962827,-3.660583759962827) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,3.4218227182900023,-3.4218227182900023) q[11];
rzz(0) q[11],q[19];
u3(pi,0.7068583470577035,-0.7068583470577035) q[11];
rzz(0) q[11],q[20];
u3(pi,4.702335883893202,-4.702335883893202) q[11];
rzz(0) q[11],q[21];
u3(pi/2,4.139362480369912,-4.139362480369912) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/128) q[12],q[18];
u3(pi/2,3.5481147429643123,-3.5481147429643123) q[19];
rzz(0.012275082143518329) q[12],q[19];
u3(pi/2,2.5685661535750146,-2.5685661535750146) q[12];
rzz(0) q[12],q[20];
u3(pi,5.7365481854549625,-5.7365481854549625) q[12];
rzz(0) q[12],q[21];
u3(pi/2,0.4567875718319559,-0.4567875718319559) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/128) q[13],q[19];
u3(pi/2,3.660583759962827,-3.660583759962827) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi/2,2.0275838986268524,-2.0275838986268524) q[13];
rzz(0) q[13],q[21];
u3(pi/2,0.24127431579569608,-0.24127431579569608) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/128) q[14],q[20];
u3(pi/2,3.6630970340856988,-3.6630970340856988) q[21];
rzz(0.012275082143518329) q[14],q[21];
u3(pi/2,5.605857931065627,-5.605857931065627) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/32) q[15],q[19];
rzz(0.049100328574073315) q[15],q[20];
rzz(pi/128) q[15],q[21];
u3(pi/2,1.6235750833752052,-1.6235750833752052) q[16];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/8) q[16],q[18];
rzz(0.19640131429629326) q[16],q[19];
rzz(pi/32) q[16],q[20];
rzz(0.049100328574073315) q[16],q[21];
u3(pi/2,4.537088110314379,-4.537088110314379) q[17];
rzz(0.7853830994606742) q[17],q[18];
rzz(pi/8) q[17],q[19];
rzz(0.19640131429629326) q[17],q[20];
rzz(pi/32) q[17],q[21];
u3(pi/2,1.9685219567393646,-1.9685219567393646) q[18];
rzz(0.7853830994606742) q[18],q[19];
rzz(pi/8) q[18],q[20];
rzz(0.19640131429629326) q[18],q[21];
u3(pi/2,5.055450898156695,-5.055450898156695) q[19];
rzz(0.7853830994606742) q[19],q[20];
rzz(pi/8) q[19],q[21];
u3(pi/2,3.626026240773339,-3.626026240773339) q[20];
rzz(0.7853830994606742) q[20],q[21];
u3(pi/2,4.816689856483871,-4.816689856483871) q[0];
u3(pi/2,0.8212123196483719,-0.8212123196483719) q[1];
u3(pi/2,0.24315927138784998,-0.24315927138784998) q[3];
u3(pi/2,2.796017461694916,-2.796017461694916) q[5];
u3(pi/2,4.106061598241859,-4.106061598241859) q[7];
u3(pi/2,4.238636808223348,-4.238636808223348) q[9];
u3(pi/2,1.9873715126609033,-1.9873715126609033) q[11];
u3(pi/2,2.6219732286860413,-2.6219732286860413) q[12];
u3(pi/2,5.169176552216645,-5.169176552216645) q[13];
u3(pi/2,5.9991853312950685,-5.9991853312950685) q[21];
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
