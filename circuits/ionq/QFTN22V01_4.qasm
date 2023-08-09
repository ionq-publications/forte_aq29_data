OPENQASM 2.0;
include "qelib1.inc";
qreg q[22];
creg c[22];
u3(pi/2,0.4511327050554943,-0.4511327050554943) q[21];
u3(pi/2,2.807327195247839,-2.807327195247839) q[21];
rzz(-pi/2) q[20],q[21];
u3(pi/2,4.378123522042736,-4.378123522042736) q[21];
rzz(pi/8) q[19],q[21];
u3(pi/2,1.7077697664914115,-1.7077697664914115) q[20];
rzz(0.7853830994606742) q[19],q[20];
rzz(0.19640131429629326) q[18],q[21];
rzz(pi/8) q[18],q[20];
u3(pi/2,4.94612347381177,-4.94612347381177) q[19];
rzz(0.7853830994606742) q[18],q[19];
u3(pi,3.003362576831842,-3.003362576831842) q[21];
rzz(pi/32) q[17],q[21];
u3(pi,5.909335781402401,-5.909335781402401) q[20];
rzz(0.19640131429629326) q[17],q[20];
u3(pi,3.749804991324777,-3.749804991324777) q[19];
rzz(pi/8) q[17],q[19];
u3(pi/2,3.6894864123758535,-3.6894864123758535) q[18];
rzz(0.7853830994606742) q[17],q[18];
u3(pi,3.28422096006277,-3.28422096006277) q[21];
rzz(0.049100328574073315) q[16],q[21];
u3(pi,4.254344771491298,-4.254344771491298) q[20];
rzz(pi/32) q[16],q[20];
u3(pi,4.363043877305505,-4.363043877305505) q[19];
rzz(0.19640131429629326) q[16],q[19];
u3(pi,0.44799111240190453,-0.44799111240190453) q[18];
rzz(pi/8) q[16],q[18];
u3(pi/2,4.075902308767398,-4.075902308767398) q[17];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/128) q[15],q[21];
rzz(0.049100328574073315) q[15],q[20];
rzz(pi/32) q[15],q[19];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/8) q[15],q[17];
u3(pi/2,2.1909467166135217,-2.1909467166135217) q[16];
rzz(0.7853830994606742) q[15],q[16];
rzz(0.012275082143518329) q[14],q[21];
rzz(pi/128) q[14],q[20];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/32) q[14],q[18];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/8) q[14],q[16];
u3(pi/2,0.5956459671206248,-0.5956459671206248) q[15];
rzz(0.7853830994606742) q[14],q[15];
u3(pi/2,2.382583868482499,-2.382583868482499) q[13];
u3(pi/2,1.7982476349147976,-1.7982476349147976) q[21];
rzz(0) q[13],q[21];
u3(pi/2,5.524176522072292,-5.524176522072292) q[13];
u3(pi,1.4162299682382786,-1.4162299682382786) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi,5.2156721234897745,-5.2156721234897745) q[19];
rzz(pi/128) q[13],q[19];
u3(pi,4.087212042320321,-4.087212042320321) q[18];
rzz(0.049100328574073315) q[13],q[18];
u3(pi,3.1189731864839465,-3.1189731864839465) q[17];
rzz(pi/32) q[13],q[17];
u3(pi,0.035185837720205684,-0.035185837720205684) q[16];
rzz(0.19640131429629326) q[13],q[16];
u3(pi,5.504070329089317,-5.504070329089317) q[15];
rzz(pi/8) q[13],q[15];
u3(pi/2,5.622822531395012,-5.622822531395012) q[14];
rzz(0.7853830994606742) q[13],q[14];
u3(pi/2,3.6417342040412883,-3.6417342040412883) q[12];
rzz(0) q[12],q[21];
u3(pi,0.5265309287416493,-0.5265309287416493) q[12];
u3(pi/2,4.434672189807352,-4.434672189807352) q[20];
rzz(0) q[12],q[20];
u3(pi/2,3.6951412791523146,-3.6951412791523146) q[12];
rzz(0.012275082143518329) q[12],q[19];
rzz(pi/128) q[12],q[18];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/32) q[12],q[16];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/8) q[12],q[14];
u3(pi/2,2.382583868482499,-2.382583868482499) q[13];
rzz(0.7853830994606742) q[12],q[13];
u3(pi/2,3.2458935296889737,-3.2458935296889737) q[11];
rzz(0) q[11],q[21];
u3(pi,0.5309291584566751,-0.5309291584566751) q[11];
rzz(0) q[11],q[20];
u3(pi,4.526406695292174,-4.526406695292174) q[11];
u3(pi/2,4.258743001206323,-4.258743001206323) q[19];
rzz(0) q[11],q[19];
u3(pi/2,1.8120706425905926,-1.8120706425905926) q[11];
u3(pi,3.3106103383529244,-3.3106103383529244) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi,3.1472475203662547,-3.1472475203662547) q[17];
rzz(pi/128) q[11],q[17];
u3(pi,5.8955127737266055,-5.8955127737266055) q[16];
rzz(0.049100328574073315) q[11],q[16];
u3(pi,5.115141158574901,-5.115141158574901) q[15];
rzz(pi/32) q[11],q[15];
u3(pi,4.290787246272939,-4.290787246272939) q[14];
rzz(0.19640131429629326) q[11],q[14];
u3(pi,4.984450904185565,-4.984450904185565) q[13];
rzz(pi/8) q[11],q[13];
u3(pi/2,3.6951412791523146,-3.6951412791523146) q[12];
rzz(0.7853830994606742) q[11],q[12];
u3(pi/2,2.0005662018059804,-2.0005662018059804) q[10];
rzz(0) q[10],q[21];
u3(pi,4.3391677731382226,-4.3391677731382226) q[10];
rzz(0) q[10],q[20];
u3(pi,5.875406580743632,-5.875406580743632) q[10];
rzz(0) q[10],q[19];
u3(pi,1.1284600811694538,-1.1284600811694538) q[10];
u3(pi/2,1.9358493931420304,-1.9358493931420304) q[18];
rzz(0) q[10],q[18];
u3(pi/2,3.467689971032413,-3.467689971032413) q[10];
u3(pi,1.5890175641857174,-1.5890175641857174) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi,1.51738925168387,-1.51738925168387) q[16];
rzz(pi/128) q[10],q[16];
u3(pi,2.4196546617948584,-2.4196546617948584) q[15];
rzz(0.049100328574073315) q[10],q[15];
u3(pi,4.787158885540127,-4.787158885540127) q[14];
rzz(pi/32) q[10],q[14];
u3(pi,0.7640353333530376,-0.7640353333530376) q[13];
rzz(0.19640131429629326) q[10],q[13];
u3(pi,2.261318392053933,-2.261318392053933) q[12];
rzz(pi/8) q[10],q[12];
u3(pi/2,4.953663296180386,-4.953663296180386) q[11];
rzz(0.7853830994606742) q[10],q[11];
u3(pi/2,0.41720350439672454,-0.41720350439672454) q[9];
rzz(0) q[9],q[21];
u3(pi,3.607176684851801,-3.607176684851801) q[9];
rzz(0) q[9],q[20];
u3(pi,0.561716766461855,-0.561716766461855) q[9];
rzz(0) q[9],q[19];
u3(pi,3.8000704737822137,-3.8000704737822137) q[9];
rzz(0) q[9],q[18];
u3(pi,0.7546105553922683,-0.7546105553922683) q[9];
u3(pi/2,5.328141140488289,-5.328141140488289) q[17];
rzz(0) q[9],q[17];
u3(pi/2,3.9445837358473446,-3.9445837358473446) q[9];
rzz(0.012275082143518329) q[9],q[16];
rzz(pi/128) q[9],q[15];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/32) q[9],q[13];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/8) q[9],q[11];
u3(pi/2,0.32609731744262055,-0.32609731744262055) q[10];
rzz(0.7853830994606742) q[9],q[10];
u3(pi/2,1.7178228629828987,-1.7178228629828987) q[8];
rzz(0) q[8],q[21];
u3(pi,3.958406743523139,-3.958406743523139) q[8];
rzz(0) q[8],q[20];
u3(pi,5.297981851013827,-5.297981851013827) q[8];
rzz(0) q[8],q[19];
u3(pi,0.35499996985564664,-0.35499996985564664) q[8];
rzz(0) q[8],q[18];
u3(pi,1.6945750773463344,-1.6945750773463344) q[8];
rzz(0) q[8],q[17];
u3(pi,3.034150184837022,-3.034150184837022) q[8];
u3(pi/2,1.6889202105698726,-1.6889202105698726) q[16];
rzz(0) q[8],q[16];
u3(pi/2,5.274734065377263,-5.274734065377263) q[8];
rzz(0.012275082143518329) q[8],q[15];
rzz(pi/128) q[8],q[14];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/32) q[8],q[12];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/8) q[8],q[10];
u3(pi/2,0.8029910822575511,-0.8029910822575511) q[9];
rzz(0.7853830994606742) q[8],q[9];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[7];
rzz(0) q[7],q[21];
u3(pi,5.392229630621521,-5.392229630621521) q[7];
rzz(0) q[7],q[20];
u3(pi,3.414911214452105,-3.414911214452105) q[7];
rzz(0) q[7],q[19];
u3(pi,1.4369644797519714,-1.4369644797519714) q[7];
rzz(0) q[7],q[18];
u3(pi,5.742831370762142,-5.742831370762142) q[7];
rzz(0) q[7],q[17];
u3(pi,3.764884636062008,-3.764884636062008) q[7];
rzz(0) q[7],q[16];
u3(pi,1.7869379013618742,-1.7869379013618742) q[7];
u3(pi/2,5.021521697497925,-5.021521697497925) q[15];
rzz(0) q[7],q[15];
u3(pi/2,5.510981832927215,-5.510981832927215) q[7];
u3(pi,5.284158843338032,-5.284158843338032) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi,3.2427519370353846,-3.2427519370353846) q[13];
rzz(pi/128) q[7],q[13];
u3(pi,2.57924756859722,-2.57924756859722) q[12];
rzz(0.049100328574073315) q[7],q[12];
u3(pi,0.8300087790784233,-0.8300087790784233) q[11];
rzz(pi/32) q[7],q[11];
u3(pi,5.350132289063418,-5.350132289063418) q[10];
rzz(0.19640131429629326) q[7],q[10];
u3(pi,3.844681089463189,-3.844681089463189) q[9];
rzz(pi/8) q[7],q[9];
u3(pi/2,5.274734065377263,-5.274734065377263) q[8];
rzz(0.7853830994606742) q[7],q[8];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[6];
rzz(0) q[6],q[21];
u3(pi,3.5072740384676453,-3.5072740384676453) q[6];
rzz(0) q[6],q[20];
u3(pi,5.416734053319521,-5.416734053319521) q[6];
rzz(0) q[6],q[19];
u3(pi,1.0430087609918113,-1.0430087609918113) q[6];
rzz(0) q[6],q[18];
u3(pi,2.9524687758436876,-2.9524687758436876) q[6];
rzz(0) q[6],q[17];
u3(pi,4.861928790695564,-4.861928790695564) q[6];
rzz(0) q[6],q[16];
u3(pi,0.48820349836785387,-0.48820349836785387) q[6];
rzz(0) q[6],q[15];
u3(pi,2.39766351321973,-2.39766351321973) q[6];
u3(pi/2,3.95212355821596,-3.95212355821596) q[14];
rzz(0) q[6],q[14];
u3(pi/2,4.923504006705923,-4.923504006705923) q[6];
rzz(0.012275082143518329) q[6],q[13];
rzz(pi/128) q[6],q[12];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/32) q[6],q[10];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/8) q[6],q[8];
u3(pi/2,5.510981832927215,-5.510981832927215) q[7];
rzz(0.7853830994606742) q[6],q[7];
u3(pi/2,3*pi/8,-3*pi/8) q[5];
rzz(0) q[5],q[21];
u3(pi,5.304893354851725,-5.304893354851725) q[5];
rzz(0) q[5],q[20];
u3(pi,4.134335932124168,-4.134335932124168) q[5];
rzz(0) q[5],q[19];
u3(pi,2.9637785093966107,-2.9637785093966107) q[5];
rzz(0) q[5],q[18];
u3(pi,1.792592768138336,-1.792592768138336) q[5];
rzz(0) q[5],q[17];
u3(pi,0.6220353454107791,-0.6220353454107791) q[5];
rzz(0) q[5],q[16];
u3(pi,5.734663229862808,-5.734663229862808) q[5];
rzz(0) q[5],q[15];
u3(pi,4.5634774886045335,-4.5634774886045335) q[5];
rzz(0) q[5],q[14];
u3(pi,3.392920065876977,-3.392920065876977) q[5];
u3(pi/2,6.260565840073739,-6.260565840073739) q[13];
rzz(0) q[5],q[13];
u3(pi/2,1.2365308684529426,-1.2365308684529426) q[5];
u3(pi,1.6160352610065896,-1.6160352610065896) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi,4.768937648149306,-4.768937648149306) q[11];
rzz(pi/128) q[5],q[11];
u3(pi,3.968459840014627,-3.968459840014627) q[10];
rzz(0.049100328574073315) q[5],q[10];
u3(pi,2.787849320795582,-2.787849320795582) q[9];
rzz(pi/32) q[5],q[9];
u3(pi,2.7149643712322993,-2.7149643712322993) q[8];
rzz(0.19640131429629326) q[5],q[8];
u3(pi,4.553424392113047,-4.553424392113047) q[7];
rzz(pi/8) q[5],q[7];
u3(pi/2,1.7819113531161308,-1.7819113531161308) q[6];
rzz(0.7853830994606742) q[5],q[6];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(0) q[4],q[21];
u3(pi,5.924415426139632,-5.924415426139632) q[4];
rzz(0) q[4],q[20];
u3(pi,3.6367076557955444,-3.6367076557955444) q[4];
rzz(0) q[4],q[19];
u3(pi,1.348999885451457,-1.348999885451457) q[4];
rzz(0) q[4],q[18];
u3(pi,5.344477422286956,-5.344477422286956) q[4];
rzz(0) q[4],q[17];
u3(pi,3.0567696519428686,-3.0567696519428686) q[4];
rzz(0) q[4],q[16];
u3(pi,0.7690618815987813,-0.7690618815987813) q[4];
rzz(0) q[4],q[15];
u3(pi,4.76453941843428,-4.76453941843428) q[4];
rzz(0) q[4],q[14];
u3(pi,2.476831648090193,-2.476831648090193) q[4];
rzz(0) q[4],q[13];
u3(pi,0.18912387774610553,-0.18912387774610553) q[4];
u3(pi/2,5.184256196953876,-5.184256196953876) q[12];
rzz(0) q[4],q[12];
u3(pi/2,3.75797313222411,-3.75797313222411) q[4];
u3(pi,5.352017244655571,-5.352017244655571) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi,0.27457519792374796,-0.27457519792374796) q[10];
rzz(pi/128) q[4],q[10];
u3(pi,3.4570085560102086,-3.4570085560102086) q[9];
rzz(0.049100328574073315) q[4],q[9];
u3(pi,0.5824512779755476,-0.5824512779755476) q[8];
rzz(pi/32) q[4],q[8];
u3(pi,5.406052638297316,-5.406052638297316) q[7];
rzz(0.19640131429629326) q[4],q[7];
u3(pi,4.494362450225558,-4.494362450225558) q[6];
rzz(pi/8) q[4],q[6];
u3(pi/2,1.2365308684529426,-1.2365308684529426) q[5];
rzz(0.7853830994606742) q[4],q[5];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(0) q[3],q[21];
u3(pi,5.1390172627421835,-5.1390172627421835) q[3];
rzz(0) q[3],q[20];
u3(pi,2.851309492398096,-2.851309492398096) q[3];
rzz(0) q[3],q[19];
u3(pi,0.5636017220540089,-0.5636017220540089) q[3];
rzz(0) q[3],q[18];
u3(pi,4.559079258889508,-4.559079258889508) q[3];
rzz(0) q[3],q[17];
u3(pi,2.2713714885454204,-2.2713714885454204) q[3];
rzz(0) q[3],q[16];
u3(pi,6.266849025380919,-6.266849025380919) q[3];
rzz(0) q[3],q[15];
u3(pi,3.979141255036832,-3.979141255036832) q[3];
rzz(0) q[3],q[14];
u3(pi,1.6914334846927446,-1.6914334846927446) q[3];
rzz(0) q[3],q[13];
u3(pi,5.6869110215282435,-5.6869110215282435) q[3];
rzz(0) q[3],q[12];
u3(pi,3.3992032511841566,-3.3992032511841566) q[3];
u3(pi/2,4.155698762168578,-4.155698762168578) q[11];
rzz(0) q[3],q[11];
u3(pi/2,0.684238879951857,-0.684238879951857) q[3];
rzz(0.012275082143518329) q[3],q[10];
rzz(pi/128) q[3],q[9];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/32) q[3],q[7];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/8) q[3],q[5];
u3(pi/2,3.75797313222411,-3.75797313222411) q[4];
rzz(0.7853830994606742) q[3],q[4];
u3(pi/2,pi,-pi) q[2];
rzz(0) q[2],q[21];
u3(pi,5.74345968929286,-5.74345968929286) q[2];
rzz(0) q[2],q[20];
u3(pi,1.5230441184603318,-1.5230441184603318) q[2];
rzz(0) q[2],q[19];
u3(pi,3.5858138548073897,-3.5858138548073897) q[2];
rzz(0) q[2],q[18];
u3(pi,5.648583591154448,-5.648583591154448) q[2];
rzz(0) q[2],q[17];
u3(pi,1.4275397017912022,-1.4275397017912022) q[2];
rzz(0) q[2],q[16];
u3(pi,3.4903094381382602,-3.4903094381382602) q[2];
rzz(0) q[2],q[15];
u3(pi,5.553079174485318,-5.553079174485318) q[2];
rzz(0) q[2],q[14];
u3(pi,1.3326636036527904,-1.3326636036527904) q[2];
rzz(0) q[2],q[13];
u3(pi,3.3948050214691303,-3.3948050214691303) q[2];
rzz(0) q[2],q[12];
u3(pi,5.457574757816189,-5.457574757816189) q[2];
rzz(0) q[2],q[11];
u3(pi,1.2371591869836605,-1.2371591869836605) q[2];
u3(pi/2,2.9870262950331754,-2.9870262950331754) q[10];
rzz(0) q[2],q[10];
u3(pi/2,3.8390262226867273,-3.8390262226867273) q[2];
rzz(0.012275082143518329) q[2],q[9];
rzz(pi/128) q[2],q[8];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/32) q[2],q[6];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/8) q[2],q[4];
u3(pi/2,0.684238879951857,-0.684238879951857) q[3];
rzz(0.7853830994606742) q[2],q[3];
u3(pi/2,0,0) q[1];
rzz(0) q[1],q[21];
u3(pi,3.568220935947287,-3.568220935947287) q[1];
rzz(0) q[1],q[20];
u3(pi,1.2805131656031998,-1.2805131656031998) q[1];
rzz(0) q[1],q[19];
u3(pi,5.2759907024386985,-5.2759907024386985) q[1];
rzz(0) q[1],q[18];
u3(pi,2.988282932094611,-2.988282932094611) q[1];
rzz(0) q[1],q[17];
u3(pi,0.7005751617505239,-0.7005751617505239) q[1];
rzz(0) q[1],q[16];
u3(pi,4.696052698586023,-4.696052698586023) q[1];
rzz(0) q[1],q[15];
u3(pi,2.4083449282419354,-2.4083449282419354) q[1];
rzz(0) q[1],q[14];
u3(pi,0.12063715789784804,-0.12063715789784804) q[1];
rzz(0) q[1],q[13];
u3(pi,4.116114694733347,-4.116114694733347) q[1];
rzz(0) q[1],q[12];
u3(pi,1.8284069243892596,-1.8284069243892596) q[1];
rzz(0) q[1],q[11];
u3(pi,5.8238844612247584,-5.8238844612247584) q[1];
rzz(0) q[1],q[10];
u3(pi,3.536176690880671,-3.536176690880671) q[1];
u3(pi/2,1.9415042599184922,-1.9415042599184922) q[9];
rzz(0) q[1],q[9];
u3(pi/2,0.8212123196483719,-0.8212123196483719) q[1];
u3(pi,4.027521781902115,-4.027521781902115) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi,6.258680884481586,-6.258680884481586) q[7];
rzz(pi/128) q[1],q[7];
u3(pi,1.521159162868178,-1.521159162868178) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi,0.04084070449666731,-0.04084070449666731) q[5];
rzz(pi/32) q[1],q[5];
u3(pi,1.0430087609918113,-1.0430087609918113) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi,5.559362359792498,-5.559362359792498) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,0.6974335690969341,-0.6974335690969341) q[2];
rzz(0.7853830994606742) q[1],q[2];
rzz(pi/2) q[0],q[21];
rzz(pi/2) q[0],q[20];
rzz(pi/2) q[0],q[19];
rzz(pi/2) q[0],q[18];
rzz(-pi/2) q[0],q[17];
rzz(pi/2) q[0],q[16];
rzz(pi/2) q[0],q[15];
rzz(pi/2) q[0],q[14];
rzz(-pi/2) q[0],q[13];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[11];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[9];
u3(pi/2,0.7621503777608838,-0.7621503777608838) q[8];
rzz(-pi/2) q[0],q[8];
u3(pi/2,5.302380080728852,-5.302380080728852) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,5.260911057701468,-5.260911057701468) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,5.12770752918926,-5.12770752918926) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,4.611229696939098,-4.611229696939098) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,4.151300532453552,-4.151300532453552) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3.8390262226867273,-3.8390262226867273) q[2];
rzz(pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi,0,0) q[0];
u3(pi/2,3.9628049732381654,-3.9628049732381654) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,5.409822549481624,-5.409822549481624) q[2];
u3(pi/2,1.4828317324943823,-1.4828317324943823) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,5.722096859248449,-5.722096859248449) q[3];
u3(pi/2,2.187805123959932,-2.187805123959932) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3.0404333701442017,-3.0404333701442017) q[4];
u3(pi/2,5.985990642149992,-5.985990642149992) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,0.4153185488045707,-0.4153185488045707) q[5];
u3(pi/2,3.458893511602362,-3.458893511602362) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,0.5485220773167779,-0.5485220773167779) q[6];
u3(pi/2,3.6411058855105702,-3.6411058855105702) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,0.5899911003441631,-0.5899911003441631) q[7];
u3(pi/2,3.707079331235956,-3.707079331235956) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,5.474539358145574,-5.474539358145574) q[8];
u3(pi/2,2.3210086524721394,-2.3210086524721394) q[8];
rzz(pi/2) q[0],q[8];
rzz(pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[11];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[13];
rzz(-pi/2) q[0],q[14];
rzz(pi/2) q[0],q[15];
rzz(pi/2) q[0],q[16];
rzz(pi/2) q[0],q[17];
rzz(-pi/2) q[0],q[18];
rzz(pi/2) q[0],q[19];
rzz(pi/2) q[0],q[20];
rzz(pi/2) q[0],q[21];
u3(pi/2,5.5336013000330615,-5.5336013000330615) q[1];
u3(pi/2,3.053628059289279,-3.053628059289279) q[2];
rzz(0.7853830994606742) q[1],q[2];
u3(pi/2,3.758601450754828,-3.758601450754828) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,1.273601661765302,-1.273601661765302) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi/2,1.8880971848074657,-1.8880971848074657) q[5];
rzz(pi/32) q[1],q[5];
u3(pi/2,5.211902212305467,-5.211902212305467) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi/2,5.277875658030852,-5.277875658030852) q[7];
rzz(pi/128) q[1],q[7];
u3(pi/2,3.8918049792670355,-3.8918049792670355) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi/2,2.3920086464432684,-2.3920086464432684) q[1];
rzz(0) q[1],q[9];
u3(pi,4.994504000677053,-4.994504000677053) q[1];
rzz(0) q[1],q[10];
u3(pi,0.7740884298445251,-0.7740884298445251) q[1];
rzz(0) q[1],q[11];
u3(pi,2.8362298476608654,-2.8362298476608654) q[1];
rzz(0) q[1],q[12];
u3(pi,4.898999584007923,-4.898999584007923) q[1];
rzz(0) q[1],q[13];
u3(pi,0.6785840131753953,-0.6785840131753953) q[1];
rzz(0) q[1],q[14];
u3(pi,2.7407254309917355,-2.7407254309917355) q[1];
rzz(0) q[1],q[15];
u3(pi,4.803495167338793,-4.803495167338793) q[1];
rzz(0) q[1],q[16];
u3(pi,0.5830795965062655,-0.5830795965062655) q[1];
rzz(0) q[1],q[17];
u3(pi,2.6458493328533237,-2.6458493328533237) q[1];
rzz(0) q[1],q[18];
u3(pi,4.707990750669664,-4.707990750669664) q[1];
rzz(0) q[1],q[19];
u3(pi,0.4875751798371359,-0.4875751798371359) q[1];
rzz(0) q[1],q[20];
u3(pi,2.550344916184194,-2.550344916184194) q[1];
rzz(0) q[1],q[21];
u3(pi/2,2.2682298958918308,-2.2682298958918308) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,5.070530542893926,-5.070530542893926) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi/2,0.6974335690969341,-0.6974335690969341) q[2];
rzz(0) q[2],q[10];
u3(pi,3.2999289233307185,-3.2999289233307185) q[2];
rzz(0) q[2],q[11];
u3(pi,5.362070341147059,-5.362070341147059) q[2];
rzz(0) q[2],q[12];
u3(pi,1.1416547703145308,-1.1416547703145308) q[2];
rzz(0) q[2],q[13];
u3(pi,3.204424506661589,-3.204424506661589) q[2];
rzz(0) q[2],q[14];
u3(pi,5.267194243008648,-5.267194243008648) q[2];
rzz(0) q[2],q[15];
u3(pi,1.046150353645401,-1.046150353645401) q[2];
rzz(0) q[2],q[16];
u3(pi,3.108920089992459,-3.108920089992459) q[2];
rzz(0) q[2],q[17];
u3(pi,5.1716898263395175,-5.1716898263395175) q[2];
rzz(0) q[2],q[18];
u3(pi,0.9512742555069894,-0.9512742555069894) q[2];
rzz(0) q[2],q[19];
u3(pi,3.0134156733233297,-3.0134156733233297) q[2];
rzz(0) q[2],q[20];
u3(pi,5.076185409670387,-5.076185409670387) q[2];
rzz(0) q[2],q[21];
u3(pi/2,3.365902369056104,-3.365902369056104) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,6.112910985355019,-6.112910985355019) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,4.936698695851001,-4.936698695851001) q[3];
rzz(0) q[3],q[11];
u3(pi,2.22236264314942,-2.22236264314942) q[3];
rzz(0) q[3],q[12];
u3(pi,6.217840179984918,-6.217840179984918) q[3];
rzz(0) q[3],q[13];
u3(pi,3.930132409640831,-3.930132409640831) q[3];
rzz(0) q[3],q[14];
u3(pi,1.6424246392967439,-1.6424246392967439) q[3];
rzz(0) q[3],q[15];
u3(pi,5.6379021761322425,-5.6379021761322425) q[3];
rzz(0) q[3],q[16];
u3(pi,3.3495660872574375,-3.3495660872574375) q[3];
rzz(0) q[3],q[17];
u3(pi,1.0618583169133502,-1.0618583169133502) q[3];
rzz(0) q[3],q[18];
u3(pi,5.057335853748849,-5.057335853748849) q[3];
rzz(0) q[3],q[19];
u3(pi,2.7696280834047617,-2.7696280834047617) q[3];
rzz(0) q[3],q[20];
u3(pi,0.4819203130606743,-0.4819203130606743) q[3];
rzz(0) q[3],q[21];
u3(pi/2,1.076937961650581,-1.076937961650581) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,1.000911419433708,-1.000911419433708) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,5.789326942035271,-5.789326942035271) q[4];
rzz(0) q[4],q[12];
u3(pi,3.0749908893336895,-3.0749908893336895) q[4];
rzz(0) q[4],q[13];
u3(pi,0.7872831189896021,-0.7872831189896021) q[4];
rzz(0) q[4],q[14];
u3(pi,4.782760655825101,-4.782760655825101) q[4];
rzz(0) q[4],q[15];
u3(pi,2.4950528854810137,-2.4950528854810137) q[4];
rzz(0) q[4],q[16];
u3(pi,0.20734511513692636,-0.20734511513692636) q[4];
rzz(0) q[4],q[17];
u3(pi,4.202822651972426,-4.202822651972426) q[4];
rzz(0) q[4],q[18];
u3(pi,1.9151148816283379,-1.9151148816283379) q[4];
rzz(0) q[4],q[19];
u3(pi,5.910592418463836,-5.910592418463836) q[4];
rzz(0) q[4],q[20];
u3(pi,3.6228846481197494,-3.6228846481197494) q[4];
rzz(0) q[4],q[21];
u3(pi/2,4.9310438290745395,-4.9310438290745395) q[5];
rzz(0.7853830994606742) q[5],q[6];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,5.166663278093774,-5.166663278093774) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,0.21865484868984958,-0.21865484868984958) q[5];
rzz(0) q[5],q[13];
u3(pi,3.7875041031678545,-3.7875041031678545) q[5];
rzz(0) q[5],q[14];
u3(pi,1.4997963328237671,-1.4997963328237671) q[5];
rzz(0) q[5],q[15];
u3(pi,5.4952738696592665,-5.4952738696592665) q[5];
rzz(0) q[5],q[16];
u3(pi,3.2075660993151787,-3.2075660993151787) q[5];
rzz(0) q[5],q[17];
u3(pi,0.9198583289710914,-0.9198583289710914) q[5];
rzz(0) q[5],q[18];
u3(pi,4.915335865806591,-4.915335865806591) q[5];
rzz(0) q[5],q[19];
u3(pi,2.627628095462503,-2.627628095462503) q[5];
rzz(0) q[5],q[20];
u3(pi,0.33992032511841563,-0.33992032511841563) q[5];
rzz(0) q[5],q[21];
u3(pi/2,5.162893366909466,-5.162893366909466) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,6.243601239744355,-6.243601239744355) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,3.5920970401145693,-3.5920970401145693) q[6];
rzz(0) q[6],q[14];
u3(pi,0.8771326688822703,-0.8771326688822703) q[6];
rzz(0) q[6],q[15];
u3(pi,4.872610205717769,-4.872610205717769) q[6];
rzz(0) q[6],q[16];
u3(pi,2.5849024353736816,-2.5849024353736816) q[6];
rzz(0) q[6],q[17];
u3(pi,0.29719466502959446,-0.29719466502959446) q[6];
rzz(0) q[6],q[18];
u3(pi,4.292672201865093,-4.292672201865093) q[6];
rzz(0) q[6],q[19];
u3(pi,2.004964431521006,-2.004964431521006) q[6];
rzz(0) q[6],q[20];
u3(pi,6.000441968356505,-6.000441968356505) q[6];
rzz(0) q[6],q[21];
u3(pi/2,5.253371235332851,-5.253371235332851) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,3.9339023208251387,-3.9339023208251387) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,0.5409822549481623,-0.5409822549481623) q[7];
rzz(0) q[7],q[15];
u3(pi,4.10920319089545,-4.10920319089545) q[7];
rzz(0) q[7],q[16];
u3(pi,1.821495420551362,-1.821495420551362) q[7];
rzz(0) q[7],q[17];
u3(pi,5.81697295738686,-5.81697295738686) q[7];
rzz(0) q[7],q[18];
u3(pi,3.5292651870427734,-3.5292651870427734) q[7];
rzz(0) q[7],q[19];
u3(pi,1.2415574166986862,-1.2415574166986862) q[7];
rzz(0) q[7],q[20];
u3(pi,5.237034953534185,-5.237034953534185) q[7];
rzz(0) q[7],q[21];
u3(pi/2,0.7376459550628834,-0.7376459550628834) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,5.003300460107105,-5.003300460107105) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,5.450034935447572,-5.450034935447572) q[8];
rzz(0) q[8],q[16];
u3(pi,2.2085396354736244,-2.2085396354736244) q[8];
rzz(0) q[8],q[17];
u3(pi,5.150326996295107,-5.150326996295107) q[8];
rzz(0) q[8],q[18];
u3(pi,1.808300731406285,-1.808300731406285) q[8];
rzz(0) q[8],q[19];
u3(pi,4.7500880922277675,-4.7500880922277675) q[8];
rzz(0) q[8],q[20];
u3(pi,1.4086901458696632,-1.4086901458696632) q[8];
rzz(0) q[8],q[21];
u3(pi/2,1.9169998372204917,-1.9169998372204917) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,1.670698973179052,-1.670698973179052) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,3.4877961640153887,-3.4877961640153887) q[9];
rzz(0) q[9],q[17];
u3(pi,5.826397735347631,-5.826397735347631) q[9];
rzz(0) q[9],q[18];
u3(pi,1.079451235773453,-1.079451235773453) q[9];
rzz(0) q[9],q[19];
u3(pi,2.6156900433788617,-2.6156900433788617) q[9];
rzz(0) q[9],q[20];
u3(pi,4.151928850984271,-4.151928850984271) q[9];
rzz(0) q[9],q[21];
u3(pi/2,2.962521872335175,-2.962521872335175) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,2.1683272495076755,-2.1683272495076755) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi/2,1.3948671381938682,-1.3948671381938682) q[10];
rzz(0) q[10],q[18];
u3(pi,4.107318235303295,-4.107318235303295) q[10];
rzz(0) q[10],q[19];
u3(pi,0.10744246875277093,-0.10744246875277093) q[10];
rzz(0) q[10],q[20];
u3(pi,2.3907520093818326,-2.3907520093818326) q[10];
rzz(0) q[10],q[21];
u3(pi/2,0.9927432785343746,-0.9927432785343746) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,1.9176281557512098,-1.9176281557512098) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,2.565424560921425,-2.565424560921425) q[11];
rzz(0) q[11],q[19];
u3(pi,5.58323846395978,-5.58323846395978) q[11];
rzz(0) q[11],q[20];
u3(pi,2.1953449463285475,-2.1953449463285475) q[11];
rzz(0) q[11],q[21];
u3(pi/2,2.018159120666083,-2.018159120666083) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/128) q[12],q[18];
u3(pi/2,4.2398934452847845,-4.2398934452847845) q[19];
rzz(0.012275082143518329) q[12],q[19];
u3(pi/2,0.44799111240190453,-0.44799111240190453) q[12];
rzz(0) q[12],q[20];
u3(pi,4.01684036687991,-4.01684036687991) q[12];
rzz(0) q[12],q[21];
u3(pi/2,3.095097082316664,-3.095097082316664) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/128) q[13],q[19];
u3(pi/2,4.415822633885813,-4.415822633885813) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi/2,4.665893409111561,-4.665893409111561) q[13];
rzz(0) q[13],q[21];
u3(pi/2,3.9276191355179595,-3.9276191355179595) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/128) q[14],q[20];
u3(pi/2,4.92161905111377,-4.92161905111377) q[21];
rzz(0.012275082143518329) q[14],q[21];
u3(pi/2,4.997017274799925,-4.997017274799925) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/32) q[15],q[19];
rzz(0.049100328574073315) q[15],q[20];
rzz(pi/128) q[15],q[21];
u3(pi/2,1.6644157878718726,-1.6644157878718726) q[16];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/8) q[16],q[18];
rzz(0.19640131429629326) q[16],q[19];
rzz(pi/32) q[16],q[20];
rzz(0.049100328574073315) q[16],q[21];
u3(pi/2,5.303636717790289,-5.303636717790289) q[17];
rzz(0.7853830994606742) q[17],q[18];
rzz(pi/8) q[17],q[19];
rzz(0.19640131429629326) q[17],q[20];
rzz(pi/32) q[17],q[21];
u3(pi/2,1.9113449704440304,-1.9113449704440304) q[18];
rzz(0.7853830994606742) q[18],q[19];
rzz(pi/8) q[18],q[20];
rzz(0.19640131429629326) q[18],q[21];
u3(pi/2,4.233610259977605,-4.233610259977605) q[19];
rzz(0.7853830994606742) q[19],q[20];
rzz(pi/8) q[19],q[21];
u3(pi/2,4.4101677671093515,-4.4101677671093515) q[20];
rzz(0.7853830994606742) q[20],q[21];
u3(pi,0,0) q[0];
u3(pi/2,5.15221195188726,-5.15221195188726) q[1];
u3(pi/2,1.3948671381938682,-1.3948671381938682) q[2];
u3(pi/2,4.05076956753868,-4.05076956753868) q[3];
u3(pi/2,0.9079202768874501,-0.9079202768874501) q[4];
u3(pi/2,3.9081412610657025,-3.9081412610657025) q[5];
u3(pi/2,3.285477597124206,-3.285477597124206) q[6];
u3(pi/2,2.522698900832604,-2.522698900832604) q[7];
u3(pi/2,4.449751834544584,-4.449751834544584) q[8];
u3(pi/2,0.2079734336676443,-0.2079734336676443) q[9];
u3(pi/2,5.10320310649126,-5.10320310649126) q[10];
u3(pi/2,5.213158849366903,-5.213158849366903) q[11];
u3(pi/2,1.3018759956476103,-1.3018759956476103) q[12];
u3(pi/2,1.5243007555217676,-1.5243007555217676) q[13];
u3(pi/2,1.773743212216797,-1.773743212216797) q[21];
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
