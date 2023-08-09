OPENQASM 2.0;
include "qelib1.inc";
qreg q[21];
creg c[21];
u3(pi/2,3.254061670588308,-3.254061670588308) q[20];
u3(pi/2,5.6102561607806525,-5.6102561607806525) q[20];
rzz(-pi/2) q[19],q[20];
u3(pi/2,4.039459833985756,-4.039459833985756) q[20];
rzz(pi/8) q[18],q[20];
u3(pi/2,0.13634512116579703,-0.13634512116579703) q[19];
rzz(0.7853830994606742) q[18],q[19];
rzz(0.19640131429629326) q[17],q[20];
rzz(pi/8) q[17],q[19];
u3(pi/2,3.5914687215838517,-3.5914687215838517) q[18];
rzz(0.7853830994606742) q[17],q[18];
rzz(pi/32) q[16],q[20];
rzz(0.19640131429629326) q[16],q[19];
rzz(pi/8) q[16],q[18];
u3(pi/2,0.546637121724624,-0.546637121724624) q[17];
rzz(0.7853830994606742) q[16],q[17];
rzz(0.049100328574073315) q[15],q[20];
rzz(pi/32) q[15],q[19];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/8) q[15],q[17];
u3(pi/2,5.259026102109313,-5.259026102109313) q[16];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/128) q[14],q[20];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/32) q[14],q[18];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/8) q[14],q[16];
u3(pi/2,2.185920168367778,-2.185920168367778) q[15];
rzz(0.7853830994606742) q[14],q[15];
u3(pi,0.2814867017616455,-0.2814867017616455) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi,5.166034959563056,-5.166034959563056) q[19];
rzz(pi/128) q[13],q[19];
u3(pi,0.8765043503515524,-0.8765043503515524) q[18];
rzz(0.049100328574073315) q[13],q[18];
u3(pi,0.7181680806106266,-0.7181680806106266) q[17];
rzz(pi/32) q[13],q[17];
u3(pi,4.631964208452791,-4.631964208452791) q[16];
rzz(0.19640131429629326) q[13],q[16];
u3(pi,0.7778583410288328,-0.7778583410288328) q[15];
rzz(pi/8) q[13],q[15];
u3(pi/2,0.9098052324796042,-0.9098052324796042) q[14];
rzz(0.7853830994606742) q[13],q[14];
u3(pi/2,5.971539315943479,-5.971539315943479) q[12];
u3(pi/2,2.807327195247839,-2.807327195247839) q[20];
rzz(0) q[12],q[20];
u3(pi/2,2.8299466623536857,-2.8299466623536857) q[12];
u3(pi,5.646070317031576,-5.646070317031576) q[19];
rzz(0.012275082143518329) q[12],q[19];
u3(pi,5.042884527542336,-5.042884527542336) q[18];
rzz(pi/128) q[12],q[18];
u3(pi,3.0724776152108175,-3.0724776152108175) q[17];
rzz(0.049100328574073315) q[12],q[17];
u3(pi,5.892999499603734,-5.892999499603734) q[16];
rzz(pi/32) q[12],q[16];
u3(pi,1.1642742374203774,-1.1642742374203774) q[15];
rzz(0.19640131429629326) q[12],q[15];
u3(pi,3.928247454048677,-3.928247454048677) q[14];
rzz(pi/8) q[12],q[14];
u3(pi/2,2.461123684822244,-2.461123684822244) q[13];
rzz(0.7853830994606742) q[12],q[13];
u3(pi/2,3.5619377506401073,-3.5619377506401073) q[11];
rzz(0) q[11],q[20];
u3(pi,0.44673447534046856,-0.44673447534046856) q[11];
u3(pi/2,4.238008489692631,-4.238008489692631) q[19];
rzz(0) q[11],q[19];
u3(pi/2,3.615344825751134,-3.615344825751134) q[11];
rzz(0.012275082143518329) q[11],q[18];
rzz(pi/128) q[11],q[17];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/32) q[11],q[15];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/8) q[11],q[13];
u3(pi/2,5.971539315943479,-5.971539315943479) q[12];
rzz(0.7853830994606742) q[11],q[12];
u3(pi/2,1.9942830164988008,-1.9942830164988008) q[10];
rzz(0) q[10],q[20];
u3(pi,4.2348668970390415,-4.2348668970390415) q[10];
rzz(0) q[10],q[19];
u3(pi,5.574442004529729,-5.574442004529729) q[10];
u3(pi/2,2.4994511151960395,-2.4994511151960395) q[18];
rzz(0) q[10],q[18];
u3(pi/2,1.531840577890383,-1.531840577890383) q[10];
u3(pi,2.5402918196927065,-2.5402918196927065) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi,4.160725310414322,-4.160725310414322) q[16];
rzz(pi/128) q[10],q[16];
u3(pi,0.24441590844928587,-0.24441590844928587) q[15];
rzz(0.049100328574073315) q[10],q[15];
u3(pi,2.5459466864691684,-2.5459466864691684) q[14];
rzz(pi/32) q[10],q[14];
u3(pi,0.31792917654328706,-0.31792917654328706) q[13];
rzz(0.19640131429629326) q[10],q[13];
u3(pi,1.928937889304133,-1.928937889304133) q[12];
rzz(pi/8) q[10],q[12];
u3(pi/2,0.47375217216134075,-0.47375217216134075) q[11];
rzz(0.7853830994606742) q[10],q[11];
u3(pi/2,1.6813803882012572,-1.6813803882012572) q[9];
rzz(0) q[9],q[20];
u3(pi,5.249601324148545,-5.249601324148545) q[9];
rzz(0) q[9],q[19];
u3(pi,2.961893553804457,-2.961893553804457) q[9];
rzz(0) q[9],q[18];
u3(pi,0.6741857834603696,-0.6741857834603696) q[9];
u3(pi/2,6.109141074170712,-6.109141074170712) q[17];
rzz(0) q[9],q[17];
u3(pi/2,4.243035037938375,-4.243035037938375) q[9];
rzz(0.012275082143518329) q[9],q[16];
rzz(pi/128) q[9],q[15];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/32) q[9],q[13];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/8) q[9],q[11];
u3(pi/2,4.6734332314801765,-4.6734332314801765) q[10];
rzz(0.7853830994606742) q[9],q[10];
u3(pi/2,1.6933184402848986,-1.6933184402848986) q[8];
rzz(0) q[8],q[20];
u3(pi,5.262167694762904,-5.262167694762904) q[8];
rzz(0) q[8],q[19];
u3(pi,2.974459924418816,-2.974459924418816) q[8];
rzz(0) q[8],q[18];
u3(pi,0.6867521540747288,-0.6867521540747288) q[8];
rzz(0) q[8],q[17];
u3(pi,4.682229690910227,-4.682229690910227) q[8];
u3(pi/2,0.5409822549481623,-0.5409822549481623) q[16];
rzz(0) q[8],q[16];
u3(pi/2,1.9672653196779284,-1.9672653196779284) q[8];
rzz(0.012275082143518329) q[8],q[15];
rzz(pi/128) q[8],q[14];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/32) q[8],q[12];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/8) q[8],q[10];
u3(pi/2,1.1014423843485814,-1.1014423843485814) q[9];
rzz(0.7853830994606742) q[8],q[9];
u3(pi/2,0.4417079270947249,-0.4417079270947249) q[7];
rzz(0) q[7],q[20];
u3(pi,3.0435749627977917,-3.0435749627977917) q[7];
rzz(0) q[7],q[19];
u3(pi,5.10634469914485,-5.10634469914485) q[7];
rzz(0) q[7],q[18];
u3(pi,0.8859291283123216,-0.8859291283123216) q[7];
rzz(0) q[7],q[17];
u3(pi,2.94869886465938,-2.94869886465938) q[7];
rzz(0) q[7],q[16];
u3(pi,5.01084028247572,-5.01084028247572) q[7];
u3(pi/2,3.8132651629272907,-3.8132651629272907) q[15];
rzz(0) q[7],q[15];
u3(pi/2,1.3301503295299184,-1.3301503295299184) q[7];
u3(pi,4.6011766004476105,-4.6011766004476105) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi,4.88454825780141,-4.88454825780141) q[13];
rzz(pi/128) q[7],q[13];
u3(pi,3.741008531894726,-3.741008531894726) q[12];
rzz(0.049100328574073315) q[7],q[12];
u3(pi,4.041973108108627,-4.041973108108627) q[11];
rzz(pi/32) q[7],q[11];
u3(pi,3.235840433197487,-3.235840433197487) q[10];
rzz(0.19640131429629326) q[7],q[10];
u3(pi,4.577300496280329,-4.577300496280329) q[9];
rzz(pi/8) q[7],q[9];
u3(pi/2,1.9672653196779284,-1.9672653196779284) q[8];
rzz(0.7853830994606742) q[7],q[8];
u3(pi/2,0.4907167724907257,-0.4907167724907257) q[6];
rzz(0) q[6],q[20];
u3(pi,3.659327122901391,-3.659327122901391) q[6];
rzz(0) q[6],q[19];
u3(pi,0.5711415444226243,-0.5711415444226243) q[6];
rzz(0) q[6],q[18];
u3(pi,3.766141273123444,-3.766141273123444) q[6];
rzz(0) q[6],q[17];
u3(pi,0.6779556946446773,-0.6779556946446773) q[6];
rzz(0) q[6],q[16];
u3(pi,3.8729554233454966,-3.8729554233454966) q[6];
rzz(0) q[6],q[15];
u3(pi,0.7847698448667303,-0.7847698448667303) q[6];
u3(pi/2,4.772707559333614,-4.772707559333614) q[14];
rzz(0) q[6],q[14];
u3(pi/2,3.9527518767466776,-3.9527518767466776) q[6];
u3(pi,2.596212168926605,-2.596212168926605) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi,0.5045397801665208,-0.5045397801665208) q[12];
rzz(pi/128) q[6],q[12];
u3(pi,1.7542653377645405,-1.7542653377645405) q[11];
rzz(0.049100328574073315) q[6],q[11];
u3(pi,3.6266545593040576,-3.6266545593040576) q[10];
rzz(pi/32) q[6],q[10];
u3(pi,2.1972299019207013,-2.1972299019207013) q[9];
rzz(0.19640131429629326) q[6],q[9];
u3(pi,0.7137698508956011,-0.7137698508956011) q[8];
rzz(pi/8) q[6],q[8];
u3(pi/2,4.471742983119712,-4.471742983119712) q[7];
rzz(0.7853830994606742) q[6],q[7];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[5];
rzz(0) q[5],q[20];
u3(pi,4.292672201865093,-4.292672201865093) q[5];
rzz(0) q[5],q[19];
u3(pi,6.2021322167169695,-6.2021322167169695) q[5];
rzz(0) q[5],q[18];
u3(pi,1.8284069243892596,-1.8284069243892596) q[5];
rzz(0) q[5],q[17];
u3(pi,3.737866939241136,-3.737866939241136) q[5];
rzz(0) q[5],q[16];
u3(pi,5.647326954093012,-5.647326954093012) q[5];
rzz(0) q[5],q[15];
u3(pi,1.273601661765302,-1.273601661765302) q[5];
rzz(0) q[5],q[14];
u3(pi,3.1830616766171786,-3.1830616766171786) q[5];
u3(pi/2,6.16506142340461,-6.16506142340461) q[13];
rzz(0) q[5],q[13];
u3(pi/2,5.7089021701033715,-5.7089021701033715) q[5];
u3(pi,3.5518846541486204,-3.5518846541486204) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi,4.844335871835461,-4.844335871835461) q[11];
rzz(pi/128) q[5],q[11];
u3(pi,3.963433291768883,-3.963433291768883) q[10];
rzz(0.049100328574073315) q[5],q[10];
u3(pi,5.1496986777643885,-5.1496986777643885) q[9];
rzz(pi/32) q[5],q[9];
u3(pi,1.643052957827462,-1.643052957827462) q[8];
rzz(0.19640131429629326) q[5],q[8];
u3(pi,0.790424711643192,-0.790424711643192) q[7];
rzz(pi/8) q[5],q[7];
u3(pi/2,3.9527518767466776,-3.9527518767466776) q[6];
rzz(0.7853830994606742) q[5],q[6];
u3(pi/2,5*pi/8,-5*pi/8) q[4];
rzz(0) q[4],q[20];
u3(pi,4.302096979825863,-4.302096979825863) q[4];
rzz(0) q[4],q[19];
u3(pi,5.838335787431272,-5.838335787431272) q[4];
rzz(0) q[4],q[18];
u3(pi,1.091389287857094,-1.091389287857094) q[4];
rzz(0) q[4],q[17];
u3(pi,2.627628095462503,-2.627628095462503) q[4];
rzz(0) q[4],q[16];
u3(pi,4.163866903067912,-4.163866903067912) q[4];
rzz(0) q[4],q[15];
u3(pi,5.700105710673321,-5.700105710673321) q[4];
rzz(0) q[4],q[14];
u3(pi,0.9531592110991433,-0.9531592110991433) q[4];
rzz(0) q[4],q[13];
u3(pi,2.489398018704552,-2.489398018704552) q[4];
u3(pi/2,6.264335751258048,-6.264335751258048) q[12];
rzz(0) q[4],q[12];
u3(pi/2,4.828627908567512,-4.828627908567512) q[4];
rzz(0.012275082143518329) q[4],q[11];
rzz(pi/128) q[4],q[10];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/32) q[4],q[8];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/8) q[4],q[6];
u3(pi/2,5.7089021701033715,-5.7089021701033715) q[5];
rzz(0.7853830994606742) q[4],q[5];
u3(pi/2,pi,-pi) q[3];
rzz(0) q[3],q[20];
u3(pi,5.74345968929286,-5.74345968929286) q[3];
rzz(0) q[3],q[19];
u3(pi,1.5230441184603318,-1.5230441184603318) q[3];
rzz(0) q[3],q[18];
u3(pi,3.5858138548073897,-3.5858138548073897) q[3];
rzz(0) q[3],q[17];
u3(pi,5.648583591154448,-5.648583591154448) q[3];
rzz(0) q[3],q[16];
u3(pi,1.4275397017912022,-1.4275397017912022) q[3];
rzz(0) q[3],q[15];
u3(pi,3.4903094381382602,-3.4903094381382602) q[3];
rzz(0) q[3],q[14];
u3(pi,5.553079174485318,-5.553079174485318) q[3];
rzz(0) q[3],q[13];
u3(pi,1.3326636036527904,-1.3326636036527904) q[3];
rzz(0) q[3],q[12];
u3(pi,3.3948050214691303,-3.3948050214691303) q[3];
u3(pi/2,1.2245928163693014,-1.2245928163693014) q[11];
rzz(0) q[3],q[11];
u3(pi/2,5.997300375702915,-5.997300375702915) q[3];
rzz(0.012275082143518329) q[3],q[10];
rzz(pi/128) q[3],q[9];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/32) q[3],q[7];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/8) q[3],q[5];
u3(pi/2,4.828627908567512,-4.828627908567512) q[4];
rzz(0.7853830994606742) q[3],q[4];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(0) q[2],q[20];
u3(pi,5.1390172627421835,-5.1390172627421835) q[2];
rzz(0) q[2],q[19];
u3(pi,2.851309492398096,-2.851309492398096) q[2];
rzz(0) q[2],q[18];
u3(pi,0.5636017220540089,-0.5636017220540089) q[2];
rzz(0) q[2],q[17];
u3(pi,4.559079258889508,-4.559079258889508) q[2];
rzz(0) q[2],q[16];
u3(pi,2.2713714885454204,-2.2713714885454204) q[2];
rzz(0) q[2],q[15];
u3(pi,6.266849025380919,-6.266849025380919) q[2];
rzz(0) q[2],q[14];
u3(pi,3.979141255036832,-3.979141255036832) q[2];
rzz(0) q[2],q[13];
u3(pi,1.6914334846927446,-1.6914334846927446) q[2];
rzz(0) q[2],q[12];
u3(pi,5.6869110215282435,-5.6869110215282435) q[2];
rzz(0) q[2],q[11];
u3(pi,3.3992032511841566,-3.3992032511841566) q[2];
u3(pi/2,2.4711767813137313,-2.4711767813137313) q[10];
rzz(0) q[2],q[10];
u3(pi/2,0.684238879951857,-0.684238879951857) q[2];
u3(pi,0.9135751436639118,-0.9135751436639118) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi,2.463008640414398,-2.463008640414398) q[8];
rzz(pi/128) q[2],q[8];
u3(pi,4.378123522042736,-4.378123522042736) q[7];
rzz(0.049100328574073315) q[2],q[7];
u3(pi,2.9964510729939446,-2.9964510729939446) q[6];
rzz(pi/32) q[2],q[6];
u3(pi,1.7643184342560279,-1.7643184342560279) q[5];
rzz(0.19640131429629326) q[2],q[5];
u3(pi,3.481512978708209,-3.481512978708209) q[4];
rzz(pi/8) q[2],q[4];
u3(pi/2,5.997300375702915,-5.997300375702915) q[3];
rzz(0.7853830994606742) q[2],q[3];
u3(pi/2,pi,-pi) q[1];
rzz(0) q[1],q[20];
u3(pi,5.74345968929286,-5.74345968929286) q[1];
rzz(0) q[1],q[19];
u3(pi,1.5230441184603318,-1.5230441184603318) q[1];
rzz(0) q[1],q[18];
u3(pi,3.5858138548073897,-3.5858138548073897) q[1];
rzz(0) q[1],q[17];
u3(pi,5.648583591154448,-5.648583591154448) q[1];
rzz(0) q[1],q[16];
u3(pi,1.4275397017912022,-1.4275397017912022) q[1];
rzz(0) q[1],q[15];
u3(pi,3.4903094381382602,-3.4903094381382602) q[1];
rzz(0) q[1],q[14];
u3(pi,5.553079174485318,-5.553079174485318) q[1];
rzz(0) q[1],q[13];
u3(pi,1.3326636036527904,-1.3326636036527904) q[1];
rzz(0) q[1],q[12];
u3(pi,3.3948050214691303,-3.3948050214691303) q[1];
rzz(0) q[1],q[11];
u3(pi,5.457574757816189,-5.457574757816189) q[1];
rzz(0) q[1],q[10];
u3(pi,1.2371591869836605,-1.2371591869836605) q[1];
u3(pi/2,3.576389076846621,-3.576389076846621) q[9];
rzz(0) q[1],q[9];
u3(pi/2,3.8390262226867273,-3.8390262226867273) q[1];
u3(pi,3.2835926415320515,-3.2835926415320515) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi,1.792592768138336,-1.792592768138336) q[7];
rzz(pi/128) q[1],q[7];
u3(pi,3.773052776961342,-3.773052776961342) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi,5.3312827331418795,-5.3312827331418795) q[5];
rzz(pi/32) q[1],q[5];
u3(pi,4.023123552187089,-4.023123552187089) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi,4.6841146465023815,-4.6841146465023815) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,0.684238879951857,-0.684238879951857) q[2];
rzz(0.7853830994606742) q[1],q[2];
rzz(pi/2) q[0],q[20];
rzz(pi/2) q[0],q[19];
rzz(pi/2) q[0],q[18];
rzz(pi/2) q[0],q[17];
rzz(-pi/2) q[0],q[16];
rzz(pi/2) q[0],q[15];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[13];
rzz(-pi/2) q[0],q[12];
rzz(pi/2) q[0],q[11];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[9];
u3(pi/2,2.324778563656447,-2.324778563656447) q[8];
rzz(pi/2) q[0],q[8];
u3(pi/2,4.505043865247763,-4.505043865247763) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,2.3649909496223964,-2.3649909496223964) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,3.4174244885749774,-3.4174244885749774) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,2.7696280834047617,-2.7696280834047617) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,3.370928917301848,-3.370928917301848) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,3.82583153354165,-3.82583153354165) q[2];
rzz(pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi,0,0) q[0];
u3(pi/2,0.6974335690969341,-0.6974335690969341) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,2.2550352067467534,-2.2550352067467534) q[2];
u3(pi/2,4.611229696939098,-4.611229696939098) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,4.9417252440967445,-4.9417252440967445) q[3];
u3(pi/2,1.4074335088082273,-1.4074335088082273) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,1.198831756609865,-1.198831756609865) q[4];
u3(pi/2,4.143760710084937,-4.143760710084937) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,1.8466281617800804,-1.8466281617800804) q[5];
u3(pi/2,4.890203124577872,-4.890203124577872) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,0.7941946228274998,-0.7941946228274998) q[6];
u3(pi/2,3.8867784310212925,-3.8867784310212925) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,6.0758401920426595,-6.0758401920426595) q[7];
u3(pi/2,2.9097431157548663,-2.9097431157548663) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,0.7539822368615503,-0.7539822368615503) q[8];
u3(pi/2,3.8830085198369844,-3.8830085198369844) q[8];
rzz(pi/2) q[0],q[8];
rzz(pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[11];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[14];
rzz(-pi/2) q[0],q[15];
rzz(pi/2) q[0],q[16];
rzz(pi/2) q[0],q[17];
rzz(pi/2) q[0],q[18];
rzz(pi/2) q[0],q[19];
rzz(pi/2) q[0],q[20];
u3(pi/2,2.2682298958918308,-2.2682298958918308) q[1];
u3(pi/2,6.182026023733995,-6.182026023733995) q[2];
rzz(0.7853830994606742) q[1],q[2];
u3(pi/2,6.1198224891929165,-6.1198224891929165) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,5.7145570368798335,-5.7145570368798335) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi/2,0.17781414419318228,-0.17781414419318228) q[5];
rzz(pi/32) q[1],q[5];
u3(pi/2,5.457574757816189,-5.457574757816189) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi/2,1.33894678895997,-1.33894678895997) q[7];
rzz(pi/128) q[1],q[7];
u3(pi/2,5.453804846631881,-5.453804846631881) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi/2,5.409822549481624,-5.409822549481624) q[1];
rzz(0) q[1],q[9];
u3(pi,2.6954864967800423,-2.6954864967800423) q[1];
rzz(0) q[1],q[10];
u3(pi,0.40777872643595514,-0.40777872643595514) q[1];
rzz(0) q[1],q[11];
u3(pi,4.403256263271454,-4.403256263271454) q[1];
rzz(0) q[1],q[12];
u3(pi,2.1155484929273665,-2.1155484929273665) q[1];
rzz(0) q[1],q[13];
u3(pi,6.110397711232148,-6.110397711232148) q[1];
rzz(0) q[1],q[14];
u3(pi,3.822689940888061,-3.822689940888061) q[1];
rzz(0) q[1],q[15];
u3(pi,1.5349821705439728,-1.5349821705439728) q[1];
rzz(0) q[1],q[16];
u3(pi,5.530459707379472,-5.530459707379472) q[1];
rzz(0) q[1],q[17];
u3(pi,3.2427519370353846,-3.2427519370353846) q[1];
rzz(0) q[1],q[18];
u3(pi,0.9550441666912971,-0.9550441666912971) q[1];
rzz(0) q[1],q[19];
u3(pi,4.950521703526796,-4.950521703526796) q[1];
rzz(0) q[1],q[20];
u3(pi/2,2.2550352067467534,-2.2550352067467534) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,3.5644510247629793,-3.5644510247629793) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi/2,3.82583153354165,-3.82583153354165) q[2];
rzz(0) q[2],q[10];
u3(pi,1.1114954808400688,-1.1114954808400688) q[2];
rzz(0) q[2],q[11];
u3(pi,5.106973017675568,-5.106973017675568) q[2];
rzz(0) q[2],q[12];
u3(pi,2.8192652473314803,-2.8192652473314803) q[2];
rzz(0) q[2],q[13];
u3(pi,0.531557476987393,-0.531557476987393) q[2];
rzz(0) q[2],q[14];
u3(pi,4.527035013822892,-4.527035013822892) q[2];
rzz(0) q[2],q[15];
u3(pi,2.2393272434788045,-2.2393272434788045) q[2];
rzz(0) q[2],q[16];
u3(pi,6.234804780314303,-6.234804780314303) q[2];
rzz(0) q[2],q[17];
u3(pi,3.947097009970216,-3.947097009970216) q[2];
rzz(0) q[2],q[18];
u3(pi,1.6593892396261287,-1.6593892396261287) q[2];
rzz(0) q[2],q[19];
u3(pi,9*pi/5,-9*pi/5) q[2];
rzz(0) q[2],q[20];
u3(pi/2,5.727123407494193,-5.727123407494193) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,5.597061471635576,-5.597061471635576) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,4.156327080699296,-4.156327080699296) q[3];
rzz(0) q[3],q[11];
u3(pi,1.441991027997715,-1.441991027997715) q[3];
rzz(0) q[3],q[12];
u3(pi,5.437468564833214,-5.437468564833214) q[3];
rzz(0) q[3],q[13];
u3(pi,3.149760794489126,-3.149760794489126) q[3];
rzz(0) q[3],q[14];
u3(pi,0.8614247056143213,-0.8614247056143213) q[3];
rzz(0) q[3],q[15];
u3(pi,4.85690224244982,-4.85690224244982) q[3];
rzz(0) q[3],q[16];
u3(pi,2.5691944721057327,-2.5691944721057327) q[3];
rzz(0) q[3],q[17];
u3(pi,0.2814867017616455,-0.2814867017616455) q[3];
rzz(0) q[3],q[18];
u3(pi,4.276964238597144,-4.276964238597144) q[3];
rzz(0) q[3],q[19];
u3(pi,1.989256468253057,-1.989256468253057) q[3];
rzz(0) q[3],q[20];
u3(pi/2,5.518521655295831,-5.518521655295831) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,1.2107698086935064,-1.2107698086935064) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,0.8061326749111409,-0.8061326749111409) q[4];
rzz(0) q[4],q[12];
u3(pi,3.1447342462433827,-3.1447342462433827) q[4];
rzz(0) q[4],q[13];
u3(pi,4.680973053848792,-4.680973053848792) q[4];
rzz(0) q[4],q[14];
u3(pi,6.217211861454201,-6.217211861454201) q[4];
rzz(0) q[4],q[15];
u3(pi,1.4702653618800232,-1.4702653618800232) q[4];
rzz(0) q[4],q[16];
u3(pi,3.006504169485432,-3.006504169485432) q[4];
rzz(0) q[4],q[17];
u3(pi,4.542742977090841,-4.542742977090841) q[4];
rzz(0) q[4],q[18];
u3(pi,6.07898178469625,-6.07898178469625) q[4];
rzz(0) q[4],q[19];
u3(pi,1.3320352851220723,-1.3320352851220723) q[4];
rzz(0) q[4],q[20];
u3(pi/2,0.07916813487046279,-0.07916813487046279) q[5];
rzz(0.7853830994606742) q[5],q[6];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,6.246742832397945,-6.246742832397945) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,4.791557115255152,-4.791557115255152) q[5];
rzz(0) q[5],q[13];
u3(pi,1.03421230156176,-1.03421230156176) q[5];
rzz(0) q[5],q[14];
u3(pi,2.9436723164136365,-2.9436723164136365) q[5];
rzz(0) q[5],q[15];
u3(pi,4.853132331265512,-4.853132331265512) q[5];
rzz(0) q[5],q[16];
u3(pi,0.47940703893780245,-0.47940703893780245) q[5];
rzz(0) q[5],q[17];
u3(pi,2.388867053789679,-2.388867053789679) q[5];
rzz(0) q[5],q[18];
u3(pi,4.298327068641555,-4.298327068641555) q[5];
rzz(0) q[5],q[19];
u3(pi,6.2077870834934314,-6.2077870834934314) q[5];
rzz(0) q[5],q[20];
u3(pi/2,2.2669732588303946,-2.2669732588303946) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,6.147468504544507,-6.147468504544507) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,3.837769585625291,-3.837769585625291) q[6];
rzz(0) q[6],q[14];
u3(pi,0.26703537555513246,-0.26703537555513246) q[6];
rzz(0) q[6],q[15];
u3(pi,2.550344916184194,-2.550344916184194) q[6];
rzz(0) q[6],q[16];
u3(pi,4.833654456813256,-4.833654456813256) q[6];
rzz(0) q[6],q[17];
u3(pi,0.8337786902627312,-0.8337786902627312) q[6];
rzz(0) q[6],q[18];
u3(pi,3.1177165494225103,-3.1177165494225103) q[6];
rzz(0) q[6],q[19];
u3(pi,5.401026090051572,-5.401026090051572) q[6];
rzz(0) q[6],q[20];
u3(pi/2,1.3144423662619695,-1.3144423662619695) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,1.6128936683529997,-1.6128936683529997) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,6.026831346646659,-6.026831346646659) q[7];
rzz(0) q[7],q[15];
u3(pi,3.883636838367702,-3.883636838367702) q[7];
rzz(0) q[7],q[16];
u3(pi,2.7375838383381454,-2.7375838383381454) q[7];
rzz(0) q[7],q[17];
u3(pi,1.5921591568393072,-1.5921591568393072) q[7];
rzz(0) q[7],q[18];
u3(pi,0.44673447534046856,-0.44673447534046856) q[7];
rzz(0) q[7],q[19];
u3(pi,5.584495101021217,-5.584495101021217) q[7];
rzz(0) q[7],q[20];
u3(pi/2,2.300274140958446,-2.300274140958446) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,0.653451271946677,-0.653451271946677) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,3.8710704677533427,-3.8710704677533427) q[8];
rzz(0) q[8],q[16];
u3(pi,0.2506990937564655,-0.2506990937564655) q[8];
rzz(0) q[8],q[17];
u3(pi,2.4353626250628078,-2.4353626250628078) q[8];
rzz(0) q[8],q[18];
u3(pi,4.62002615636915,-4.62002615636915) q[8];
rzz(0) q[8],q[19];
u3(pi,0.5215043804959056,-0.5215043804959056) q[8];
rzz(0) q[8],q[20];
u3(pi/2,0.410292000558827,-0.410292000558827) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,3.664353671147135,-3.664353671147135) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,5.128964166250697,-5.128964166250697) q[9];
rzz(0) q[9],q[17];
u3(pi,2.414628113549115,-2.414628113549115) q[9];
rzz(0) q[9],q[18];
u3(pi,0.12692034320502762,-0.12692034320502762) q[9];
rzz(0) q[9],q[19];
u3(pi,4.121769561509809,-4.121769561509809) q[9];
rzz(0) q[9],q[20];
u3(pi/2,2.446672358615731,-2.446672358615731) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,6.090291518249173,-6.090291518249173) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi/2,4.017468685410627,-4.017468685410627) q[10];
rzz(0) q[10],q[18];
u3(pi,1.910716651913312,-1.910716651913312) q[10];
rzz(0) q[10],q[19];
u3(pi,0.8394335570391926,-0.8394335570391926) q[10];
rzz(0) q[10],q[20];
u3(pi/2,1.2032299863248908,-1.2032299863248908) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,5.622822531395012,-5.622822531395012) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,5.91561896670958,-5.91561896670958) q[11];
rzz(0) q[11],q[19];
u3(pi,2.234300695233061,-2.234300695233061) q[11];
rzz(0) q[11],q[20];
u3(pi/2,6.239831328560047,-6.239831328560047) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/128) q[12],q[18];
u3(pi/2,1.078194598712017,-1.078194598712017) q[19];
rzz(0.012275082143518329) q[12],q[19];
u3(pi/2,1.5280706667060753,-1.5280706667060753) q[12];
rzz(0) q[12],q[20];
u3(pi/2,2.9995926656475342,-2.9995926656475342) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/128) q[13],q[19];
u3(pi/2,2.7891059578570183,-2.7891059578570183) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi/2,1.60661048304582,-1.60661048304582) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/128) q[14],q[20];
u3(pi/2,0.6471680866394973,-0.6471680866394973) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/32) q[15],q[19];
rzz(0.049100328574073315) q[15],q[20];
u3(pi/2,3.6580704858399553,-3.6580704858399553) q[16];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/8) q[16],q[18];
rzz(0.19640131429629326) q[16],q[19];
rzz(pi/32) q[16],q[20];
u3(pi/2,6.084636651472712,-6.084636651472712) q[17];
rzz(0.7853830994606742) q[17],q[18];
rzz(pi/8) q[17],q[19];
rzz(0.19640131429629326) q[17],q[20];
u3(pi/2,5.616539346087833,-5.616539346087833) q[18];
rzz(0.7853830994606742) q[18],q[19];
rzz(pi/8) q[18],q[20];
u3(pi/2,4.213504066994631,-4.213504066994631) q[19];
rzz(0.7853830994606742) q[19],q[20];
u3(pi,0,0) q[0];
u3(pi/2,2.236185650825215,-2.236185650825215) q[1];
u3(pi/2,2.9399024052293283,-2.9399024052293283) q[2];
u3(pi/2,5.558105722731062,-5.558105722731062) q[3];
u3(pi/2,3.6712651749850327,-3.6712651749850327) q[4];
u3(pi/2,2.4504422698000385,-2.4504422698000385) q[5];
u3(pi/2,1.8302918799814134,-1.8302918799814134) q[6];
u3(pi/2,3.440672274211541,-3.440672274211541) q[7];
u3(pi/2,3.1843183136786144,-3.1843183136786144) q[8];
u3(pi/2,1.4074335088082273,-1.4074335088082273) q[9];
u3(pi/2,5.016495149252181,-5.016495149252181) q[10];
u3(pi/2,4.836167730936128,-4.836167730936128) q[11];
u3(pi/2,4.669663320295868,-4.669663320295868) q[12];
u3(pi/2,2.7828227725498387,-2.7828227725498387) q[20];
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
