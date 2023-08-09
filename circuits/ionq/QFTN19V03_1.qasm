OPENQASM 2.0;
include "qelib1.inc";
qreg q[19];
creg c[19];
u3(pi/2,5.659265006176653,-5.659265006176653) q[18];
u3(pi/2,1.732274189189412,-1.732274189189412) q[18];
rzz(-pi/2) q[17],q[18];
u3(pi/2,3.303070515984308,-3.303070515984308) q[18];
rzz(pi/8) q[16],q[18];
u3(pi/2,5.364583615269931,-5.364583615269931) q[17];
rzz(0.7853830994606742) q[16],q[17];
u3(pi,2.1067520334973153,-2.1067520334973153) q[18];
rzz(0.19640131429629326) q[15],q[18];
u3(pi,1.744840559803771,-1.744840559803771) q[17];
rzz(pi/8) q[15],q[17];
u3(pi/2,0.6465397681087794,-0.6465397681087794) q[16];
rzz(0.7853830994606742) q[15],q[16];
u3(pi,2.85633604064384,-2.85633604064384) q[18];
rzz(pi/32) q[14],q[18];
u3(pi,3.8685571936304712,-3.8685571936304712) q[17];
rzz(0.19640131429629326) q[14],q[17];
u3(pi,5.521034929418702,-5.521034929418702) q[16];
rzz(pi/8) q[14],q[16];
u3(pi/2,5.750999511661475,-5.750999511661475) q[15];
rzz(0.7853830994606742) q[14],q[15];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/32) q[13],q[17];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/8) q[13],q[15];
u3(pi/2,2.5849024353736816,-2.5849024353736816) q[14];
rzz(0.7853830994606742) q[13],q[14];
u3(pi,1.5601149117726911,-1.5601149117726911) q[18];
rzz(pi/128) q[12],q[18];
u3(pi,1.8126989611213105,-1.8126989611213105) q[17];
rzz(0.049100328574073315) q[12],q[17];
u3(pi,4.540229702967969,-4.540229702967969) q[16];
rzz(pi/32) q[12],q[16];
u3(pi,4.492477494633404,-4.492477494633404) q[15];
rzz(0.19640131429629326) q[12],q[15];
u3(pi,6.061388865836147,-6.061388865836147) q[14];
rzz(pi/8) q[12],q[14];
u3(pi/2,4.15507044363786,-4.15507044363786) q[13];
rzz(0.7853830994606742) q[12],q[13];
u3(pi,4.062079301091602,-4.062079301091602) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi,2.064026373408494,-2.064026373408494) q[17];
rzz(pi/128) q[11],q[17];
u3(pi,2.3291767933714724,-2.3291767933714724) q[16];
rzz(0.049100328574073315) q[11],q[16];
u3(pi,5.419247327442394,-5.419247327442394) q[15];
rzz(pi/32) q[11],q[15];
u3(pi,3.6813182714765196,-3.6813182714765196) q[14];
rzz(0.19640131429629326) q[11],q[14];
u3(pi,2.746380297768197,-2.746380297768197) q[13];
rzz(pi/8) q[11],q[13];
u3(pi/2,4.057681071376577,-4.057681071376577) q[12];
rzz(0.7853830994606742) q[11],q[12];
u3(pi/2,3.6630970340856988,-3.6630970340856988) q[10];
u3(pi/2,0.3813893481458009,-0.3813893481458009) q[18];
rzz(0) q[10],q[18];
u3(pi/2,0.5215043804959056,-0.5215043804959056) q[10];
rzz(0.012275082143518329) q[10],q[17];
rzz(pi/128) q[10],q[16];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/32) q[10],q[14];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/8) q[10],q[12];
u3(pi/2,0.91106186954104,-0.91106186954104) q[11];
rzz(0.7853830994606742) q[10],q[11];
u3(pi/2,3.6448757966948775,-3.6448757966948775) q[9];
rzz(0) q[9],q[18];
u3(pi,0.9299114254625788,-0.9299114254625788) q[9];
u3(pi/2,0.6892654281976006,-0.6892654281976006) q[17];
rzz(0) q[9],q[17];
u3(pi/2,4.498760679940584,-4.498760679940584) q[9];
rzz(0.012275082143518329) q[9],q[16];
rzz(pi/128) q[9],q[15];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/32) q[9],q[13];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/8) q[9],q[11];
u3(pi/2,3.6630970340856988,-3.6630970340856988) q[10];
rzz(0.7853830994606742) q[9],q[10];
u3(pi/2,2.0860175219836226,-2.0860175219836226) q[8];
rzz(0) q[8],q[18];
u3(pi,6.226008320884252,-6.226008320884252) q[8];
rzz(0) q[8],q[17];
u3(pi,5.080583639385413,-5.080583639385413) q[8];
u3(pi/2,5.974680908597068,-5.974680908597068) q[16];
rzz(0) q[8],q[16];
u3(pi/2,2.9367608125757383,-2.9367608125757383) q[8];
rzz(0.012275082143518329) q[8],q[15];
rzz(pi/128) q[8],q[14];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/32) q[8],q[12];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/8) q[8],q[10];
u3(pi/2,1.3571680263507906,-1.3571680263507906) q[9];
rzz(0.7853830994606742) q[8],q[9];
u3(pi/2,2.012504253889621,-2.012504253889621) q[7];
rzz(0) q[7],q[18];
u3(pi,4.351734143752581,-4.351734143752581) q[7];
rzz(0) q[7],q[17];
u3(pi,5.88797295135799,-5.88797295135799) q[7];
rzz(0) q[7],q[16];
u3(pi,1.141026451783813,-1.141026451783813) q[7];
u3(pi/2,4.462318205158942,-4.462318205158942) q[15];
rzz(0) q[7],q[15];
u3(pi/2,3.479628023116055,-3.479628023116055) q[7];
rzz(0.012275082143518329) q[7],q[14];
rzz(pi/128) q[7],q[13];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/32) q[7],q[11];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/8) q[7],q[9];
u3(pi/2,6.078353466165532,-6.078353466165532) q[8];
rzz(0.7853830994606742) q[7],q[8];
u3(pi/2,2.061513099285622,-2.061513099285622) q[6];
rzz(0) q[6],q[18];
u3(pi,5.537999529748087,-5.537999529748087) q[6];
rzz(0) q[6],q[17];
u3(pi,3.0655661113729202,-3.0655661113729202) q[6];
rzz(0) q[6],q[16];
u3(pi,0.593132692997753,-0.593132692997753) q[6];
rzz(0) q[6],q[15];
u3(pi,4.4038845818021715,-4.4038845818021715) q[6];
u3(pi/2,0.9663539002442203,-0.9663539002442203) q[14];
rzz(0) q[6],q[14];
u3(pi/2,1.5965573865543328,-1.5965573865543328) q[6];
rzz(0.012275082143518329) q[6],q[13];
rzz(pi/128) q[6],q[12];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/32) q[6],q[10];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/8) q[6],q[8];
u3(pi/2,0.33803536952626173,-0.33803536952626173) q[7];
rzz(0.7853830994606742) q[6],q[7];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[5];
rzz(0) q[5],q[18];
u3(pi,4.292672201865093,-4.292672201865093) q[5];
rzz(0) q[5],q[17];
u3(pi,6.2021322167169695,-6.2021322167169695) q[5];
rzz(0) q[5],q[16];
u3(pi,1.8284069243892596,-1.8284069243892596) q[5];
rzz(0) q[5],q[15];
u3(pi,3.737866939241136,-3.737866939241136) q[5];
rzz(0) q[5],q[14];
u3(pi,5.647326954093012,-5.647326954093012) q[5];
u3(pi/2,1.338318470429252,-1.338318470429252) q[13];
rzz(0) q[5],q[13];
u3(pi/2,1.8899821403996195,-1.8899821403996195) q[5];
u3(pi,2.7256457862545047,-2.7256457862545047) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi,4.387548300003505,-4.387548300003505) q[11];
rzz(pi/128) q[5],q[11];
u3(pi,2.4133714764876792,-2.4133714764876792) q[10];
rzz(0.049100328574073315) q[5],q[10];
u3(pi,3.8823802013062663,-3.8823802013062663) q[9];
rzz(pi/32) q[5],q[9];
u3(pi,5.119539388289927,-5.119539388289927) q[8];
rzz(0.19640131429629326) q[5],q[8];
u3(pi,3.050486466635689,-3.050486466635689) q[7];
rzz(pi/8) q[5],q[7];
u3(pi/2,1.5965573865543328,-1.5965573865543328) q[6];
rzz(0.7853830994606742) q[5],q[6];
u3(pi/2,5*pi/8,-5*pi/8) q[4];
rzz(0) q[4],q[18];
u3(pi,4.302096979825863,-4.302096979825863) q[4];
rzz(0) q[4],q[17];
u3(pi,5.838335787431272,-5.838335787431272) q[4];
rzz(0) q[4],q[16];
u3(pi,1.091389287857094,-1.091389287857094) q[4];
rzz(0) q[4],q[15];
u3(pi,2.627628095462503,-2.627628095462503) q[4];
rzz(0) q[4],q[14];
u3(pi,4.163866903067912,-4.163866903067912) q[4];
rzz(0) q[4],q[13];
u3(pi,5.700105710673321,-5.700105710673321) q[4];
u3(pi/2,1.3936105011324322,-1.3936105011324322) q[12];
rzz(0) q[4],q[12];
u3(pi/2,1.7561502933566946,-1.7561502933566946) q[4];
u3(pi,2.0074777056438777,-2.0074777056438777) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi,3.3476811316652837,-3.3476811316652837) q[10];
rzz(pi/128) q[4],q[10];
u3(pi,5.868495076905734,-5.868495076905734) q[9];
rzz(0.049100328574073315) q[4],q[9];
u3(pi,6.048194176691069,-6.048194176691069) q[8];
rzz(pi/32) q[4],q[8];
u3(pi,5.223840264389108,-5.223840264389108) q[7];
rzz(0.19640131429629326) q[4],q[7];
u3(pi,0.18849555921538758,-0.18849555921538758) q[6];
rzz(pi/8) q[4],q[6];
u3(pi/2,5.031574793989412,-5.031574793989412) q[5];
rzz(0.7853830994606742) q[4],q[5];
u3(pi/2,pi,-pi) q[3];
rzz(0) q[3],q[18];
u3(pi,5.74345968929286,-5.74345968929286) q[3];
rzz(0) q[3],q[17];
u3(pi,1.5230441184603318,-1.5230441184603318) q[3];
rzz(0) q[3],q[16];
u3(pi,3.5858138548073897,-3.5858138548073897) q[3];
rzz(0) q[3],q[15];
u3(pi,5.648583591154448,-5.648583591154448) q[3];
rzz(0) q[3],q[14];
u3(pi,1.4275397017912022,-1.4275397017912022) q[3];
rzz(0) q[3],q[13];
u3(pi,3.4903094381382602,-3.4903094381382602) q[3];
rzz(0) q[3],q[12];
u3(pi,5.553079174485318,-5.553079174485318) q[3];
u3(pi/2,5.575698641591164,-5.575698641591164) q[11];
rzz(0) q[3],q[11];
u3(pi/2,1.8717609030087987,-1.8717609030087987) q[3];
u3(pi,4.200309377849553,-4.200309377849553) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi,2.614433406317426,-2.614433406317426) q[9];
rzz(pi/128) q[3],q[9];
u3(pi,0.8928406321502192,-0.8928406321502192) q[8];
rzz(0.049100328574073315) q[3],q[8];
u3(pi,2.124344952357418,-2.124344952357418) q[7];
rzz(pi/32) q[3],q[7];
u3(pi,0.668530916683908,-0.668530916683908) q[6];
rzz(0.19640131429629326) q[3],q[6];
u3(pi,1.273601661765302,-1.273601661765302) q[5];
rzz(pi/8) q[3],q[5];
u3(pi/2,1.7561502933566946,-1.7561502933566946) q[4];
rzz(0.7853830994606742) q[3],q[4];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(0) q[2],q[18];
u3(pi,5.1390172627421835,-5.1390172627421835) q[2];
rzz(0) q[2],q[17];
u3(pi,2.851309492398096,-2.851309492398096) q[2];
rzz(0) q[2],q[16];
u3(pi,0.5636017220540089,-0.5636017220540089) q[2];
rzz(0) q[2],q[15];
u3(pi,4.559079258889508,-4.559079258889508) q[2];
rzz(0) q[2],q[14];
u3(pi,2.2713714885454204,-2.2713714885454204) q[2];
rzz(0) q[2],q[13];
u3(pi,6.266849025380919,-6.266849025380919) q[2];
rzz(0) q[2],q[12];
u3(pi,3.979141255036832,-3.979141255036832) q[2];
rzz(0) q[2],q[11];
u3(pi,1.6914334846927446,-1.6914334846927446) q[2];
u3(pi/2,2.8682740927274812,-2.8682740927274812) q[10];
rzz(0) q[2],q[10];
u3(pi/2,5.26028273917075,-5.26028273917075) q[2];
u3(pi,6.209672039085585,-6.209672039085585) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi,1.9113449704440304,-1.9113449704440304) q[8];
rzz(pi/128) q[2],q[8];
u3(pi,0.16273449945595128,-0.16273449945595128) q[7];
rzz(0.049100328574073315) q[2],q[7];
u3(pi,1.3031326327090462,-1.3031326327090462) q[6];
rzz(pi/32) q[2],q[6];
u3(pi,3.1830616766171786,-3.1830616766171786) q[5];
rzz(0.19640131429629326) q[2],q[5];
u3(pi,0.4090353634973911,-0.4090353634973911) q[4];
rzz(pi/8) q[2],q[4];
u3(pi/2,1.8717609030087987,-1.8717609030087987) q[3];
rzz(0.7853830994606742) q[2],q[3];
u3(pi/2,pi,-pi) q[1];
rzz(0) q[1],q[18];
u3(pi,5.74345968929286,-5.74345968929286) q[1];
rzz(0) q[1],q[17];
u3(pi,1.5230441184603318,-1.5230441184603318) q[1];
rzz(0) q[1],q[16];
u3(pi,3.5858138548073897,-3.5858138548073897) q[1];
rzz(0) q[1],q[15];
u3(pi,5.648583591154448,-5.648583591154448) q[1];
rzz(0) q[1],q[14];
u3(pi,1.4275397017912022,-1.4275397017912022) q[1];
rzz(0) q[1],q[13];
u3(pi,3.4903094381382602,-3.4903094381382602) q[1];
rzz(0) q[1],q[12];
u3(pi,5.553079174485318,-5.553079174485318) q[1];
rzz(0) q[1],q[11];
u3(pi,1.3326636036527904,-1.3326636036527904) q[1];
rzz(0) q[1],q[10];
u3(pi,3.3948050214691303,-3.3948050214691303) q[1];
u3(pi/2,3.094468763785946,-3.094468763785946) q[9];
rzz(0) q[1],q[9];
u3(pi/2,5.997300375702915,-5.997300375702915) q[1];
u3(pi,2.731300653030966,-2.731300653030966) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi,4.329114676646735,-4.329114676646735) q[7];
rzz(pi/128) q[1],q[7];
u3(pi,1.7831679901775666,-1.7831679901775666) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi,0.0075398223686155025,-0.0075398223686155025) q[5];
rzz(pi/32) q[1],q[5];
u3(pi,0.9506459369762713,-0.9506459369762713) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi,0.5585751738082653,-0.5585751738082653) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,5.26028273917075,-5.26028273917075) q[2];
rzz(0.7853830994606742) q[1],q[2];
rzz(pi/2) q[0],q[18];
rzz(pi/2) q[0],q[17];
rzz(pi/2) q[0],q[16];
rzz(pi/2) q[0],q[15];
rzz(-pi/2) q[0],q[14];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[10];
rzz(pi/2) q[0],q[9];
u3(pi/2,1.7724865751553613,-1.7724865751553613) q[8];
rzz(pi/2) q[0],q[8];
u3(pi/2,1.6147786239451536,-1.6147786239451536) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,0.3751061628386213,-0.3751061628386213) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,3.730955435403238,-3.730955435403238) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,5.98033577537353,-5.98033577537353) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,5.529203070318036,-5.529203070318036) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,2.1186900855809565,-2.1186900855809565) q[2];
rzz(pi/2) q[0],q[2];
rzz(-pi/2) q[0],q[1];
u3(pi,3*pi/2,-3*pi/2) q[0];
u3(pi/2,5.997300375702915,-5.997300375702915) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,0.54789375878606,-0.54789375878606) q[2];
u3(pi/2,2.904088248978405,-2.904088248978405) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,3.958406743523139,-3.958406743523139) q[3];
u3(pi/2,0.4241150082346221,-0.4241150082346221) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,4.4095394485786334,-4.4095394485786334) q[4];
u3(pi/2,1.0712830948741194,-1.0712830948741194) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,5.301751762198135,-5.301751762198135) q[5];
u3(pi/2,2.06214141781634,-2.06214141781634) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,5.087495143223311,-5.087495143223311) q[6];
u3(pi/2,1.8962653257067992,-1.8962653257067992) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,0.0439822971502571,-0.0439822971502571) q[7];
u3(pi/2,3.1610705280420497,-3.1610705280420497) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,0.2016902483604647,-0.2016902483604647) q[8];
u3(pi/2,3.3313448498666167,-3.3313448498666167) q[8];
rzz(pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[11];
rzz(pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[13];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[15];
rzz(pi/2) q[0],q[16];
rzz(-pi/2) q[0],q[17];
rzz(pi/2) q[0],q[18];
u3(pi/2,1.2849113953182254,-1.2849113953182254) q[1];
u3(pi/2,4.474884575773301,-4.474884575773301) q[2];
rzz(0.7853830994606742) q[1],q[2];
u3(pi/2,1.9949113350295187,-1.9949113350295187) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,2.642079421669016,-2.642079421669016) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi/2,0.4913450910214437,-0.4913450910214437) q[5];
rzz(pi/32) q[1],q[5];
u3(pi/2,3.4670616525016955,-3.4670616525016955) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi/2,4.731866854836946,-4.731866854836946) q[7];
rzz(pi/128) q[1],q[7];
u3(pi/2,4.902141176661513,-4.902141176661513) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi/2,4.4265040489080185,-4.4265040489080185) q[1];
rzz(0) q[1],q[9];
u3(pi,1.7115396776757192,-1.7115396776757192) q[1];
rzz(0) q[1],q[10];
u3(pi,5.707017214511218,-5.707017214511218) q[1];
rzz(0) q[1],q[11];
u3(pi,3.419309444167131,-3.419309444167131) q[1];
rzz(0) q[1],q[12];
u3(pi,1.1316016738230434,-1.1316016738230434) q[1];
rzz(0) q[1],q[13];
u3(pi,5.127079210658542,-5.127079210658542) q[1];
rzz(0) q[1],q[14];
u3(pi,2.839371440314455,-2.839371440314455) q[1];
rzz(0) q[1],q[15];
u3(pi,0.5516636699703676,-0.5516636699703676) q[1];
rzz(0) q[1],q[16];
u3(pi,4.5471412068058665,-4.5471412068058665) q[1];
rzz(0) q[1],q[17];
u3(pi,2.259433436461779,-2.259433436461779) q[1];
rzz(0) q[1],q[18];
u3(pi/2,0.54789375878606,-0.54789375878606) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,6.224123365292098,-6.224123365292098) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi/2,2.1186900855809565,-2.1186900855809565) q[2];
rzz(0) q[2],q[10];
u3(pi,5.6869110215282435,-5.6869110215282435) q[2];
rzz(0) q[2],q[11];
u3(pi,3.3992032511841566,-3.3992032511841566) q[2];
rzz(0) q[2],q[12];
u3(pi,1.1114954808400688,-1.1114954808400688) q[2];
rzz(0) q[2],q[13];
u3(pi,5.106973017675568,-5.106973017675568) q[2];
rzz(0) q[2],q[14];
u3(pi,2.8192652473314803,-2.8192652473314803) q[2];
rzz(0) q[2],q[15];
u3(pi,0.531557476987393,-0.531557476987393) q[2];
rzz(0) q[2],q[16];
u3(pi,4.527035013822892,-4.527035013822892) q[2];
rzz(0) q[2],q[17];
u3(pi,2.2393272434788045,-2.2393272434788045) q[2];
rzz(0) q[2],q[18];
u3(pi/2,1.6022122533307945,-1.6022122533307945) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,2.8525661294595324,-2.8525661294595324) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,pi/100,-pi/100) q[3];
rzz(0) q[3],q[11];
u3(pi,3.5996368624831847,-3.5996368624831847) q[3];
rzz(0) q[3],q[12];
u3(pi,1.3119290921390976,-1.3119290921390976) q[3];
rzz(0) q[3],q[13];
u3(pi,5.307406628974596,-5.307406628974596) q[3];
rzz(0) q[3],q[14];
u3(pi,3.0196988586305094,-3.0196988586305094) q[3];
rzz(0) q[3],q[15];
u3(pi,0.7319910882864218,-0.7319910882864218) q[3];
rzz(0) q[3],q[16];
u3(pi,4.7274686251219205,-4.7274686251219205) q[3];
rzz(0) q[3],q[17];
u3(pi,2.439760854777833,-2.439760854777833) q[3];
rzz(0) q[3],q[18];
u3(pi/2,5.587636693674806,-5.587636693674806) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,5.56187563391537,-5.56187563391537) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,0.8752477132901164,-0.8752477132901164) q[4];
rzz(0) q[4],q[12];
u3(pi,4.443468649237404,-4.443468649237404) q[4];
rzz(0) q[4],q[13];
u3(pi,2.1557608788933162,-2.1557608788933162) q[4];
rzz(0) q[4],q[14];
u3(pi,6.151238415728815,-6.151238415728815) q[4];
rzz(0) q[4],q[15];
u3(pi,3.863530645384728,-3.863530645384728) q[4];
rzz(0) q[4],q[16];
u3(pi,1.5758228750406404,-1.5758228750406404) q[4];
rzz(0) q[4],q[17];
u3(pi,5.571300411876139,-5.571300411876139) q[4];
rzz(0) q[4],q[18];
u3(pi/2,pi/8,-pi/8) q[5];
rzz(0.7853830994606742) q[5],q[6];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,4.517610235862122,-4.517610235862122) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,13*pi/8,-13*pi/8) q[5];
rzz(0) q[5],q[13];
u3(pi,2.3907520093818326,-2.3907520093818326) q[5];
rzz(0) q[5],q[14];
u3(pi,0.10304423903774522,-0.10304423903774522) q[5];
rzz(0) q[5],q[15];
u3(pi,4.098521775873244,-4.098521775873244) q[5];
rzz(0) q[5],q[16];
u3(pi,1.8108140055291568,-1.8108140055291568) q[5];
rzz(0) q[5],q[17];
u3(pi,5.806291542364656,-5.806291542364656) q[5];
rzz(0) q[5],q[18];
u3(pi/2,0.27646015351590175,-0.27646015351590175) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,4.46294652368966,-4.46294652368966) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,1.8472564803107983,-1.8472564803107983) q[6];
rzz(0) q[6],q[14];
u3(pi,5.416105734788803,-5.416105734788803) q[6];
rzz(0) q[6],q[15];
u3(pi,3.128397964444716,-3.128397964444716) q[6];
rzz(0) q[6],q[16];
u3(pi,0.8406901941006286,-0.8406901941006286) q[6];
rzz(0) q[6],q[17];
u3(pi,4.836167730936128,-4.836167730936128) q[6];
rzz(0) q[6],q[18];
u3(pi/2,1.565769778549153,-1.565769778549153) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,4.089725316443193,-4.089725316443193) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,6.278158758933842,-6.278158758933842) q[7];
rzz(0) q[7],q[15];
u3(pi,3.5631943877015435,-3.5631943877015435) q[7];
rzz(0) q[7],q[16];
u3(pi,1.2754866173574562,-1.2754866173574562) q[7];
rzz(0) q[7],q[17];
u3(pi,5.270964154192955,-5.270964154192955) q[7];
rzz(0) q[7],q[18];
u3(pi/2,1.7479821524573609,-1.7479821524573609) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,4.444096967768122,-4.444096967768122) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,3.3187784792522574,-3.3187784792522574) q[8];
rzz(0) q[8],q[16];
u3(pi,5.478309269329881,-5.478309269329881) q[8];
rzz(0) q[8],q[17];
u3(pi,0.37322120724646746,-0.37322120724646746) q[8];
rzz(0) q[8],q[18];
u3(pi/2,3.0699643410879456,-3.0699643410879456) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,5.9558313526755295,-5.9558313526755295) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,1.5054511996002289,-1.5054511996002289) q[9];
rzz(0) q[9],q[17];
u3(pi,5.753512785784347,-5.753512785784347) q[9];
rzz(0) q[9],q[18];
u3(pi/2,5.985362323619274,-5.985362323619274) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,3.812636844396573,-3.812636844396573) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi/2,1.2729733432345842,-1.2729733432345842) q[10];
rzz(0) q[10],q[18];
u3(pi/2,2.412743157956961,-2.412743157956961) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,0.36253979222426214,-0.36253979222426214) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,4.510698732024225,-4.510698732024225) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/128) q[12],q[18];
u3(pi/2,1.3144423662619695,-1.3144423662619695) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
u3(pi/2,0.9418494775462201,-0.9418494775462201) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
u3(pi/2,4.437813782460942,-4.437813782460942) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
u3(pi/2,2.808583832309275,-2.808583832309275) q[16];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/8) q[16],q[18];
u3(pi/2,0.6653893240303181,-0.6653893240303181) q[17];
rzz(0.7853830994606742) q[17],q[18];
u3(pi,pi,-pi) q[0];
u3(pi/2,5.827654372409066,-5.827654372409066) q[1];
u3(pi/2,5.8075481794260915,-5.8075481794260915) q[2];
u3(pi/2,6.007981790725121,-6.007981790725121) q[3];
u3(pi/2,2.85633604064384,-2.85633604064384) q[4];
u3(pi/2,3.0913271711323564,-3.0913271711323564) q[5];
u3(pi/2,2.1212033597038285,-2.1212033597038285) q[6];
u3(pi/2,2.5559997829606558,-2.5559997829606558) q[7];
u3(pi/2,2.5327519973240915,-2.5327519973240915) q[8];
u3(pi/2,3.718389064788879,-3.718389064788879) q[9];
u3(pi/2,4.414565996824377,-4.414565996824377) q[10];
u3(pi/2,3.4978492605068756,-3.4978492605068756) q[18];
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
