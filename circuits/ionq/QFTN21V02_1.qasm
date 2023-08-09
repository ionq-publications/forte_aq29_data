OPENQASM 2.0;
include "qelib1.inc";
qreg q[21];
creg c[21];
u3(pi/2,0.12503538761287378,-0.12503538761287378) q[20];
u3(pi/2,2.4812298778052186,-2.4812298778052186) q[20];
rzz(-pi/2) q[19],q[20];
u3(pi/2,4.052026204600115,-4.052026204600115) q[20];
rzz(pi/8) q[18],q[20];
u3(pi/2,5.1333623959657215,-5.1333623959657215) q[19];
rzz(0.7853830994606742) q[18],q[19];
rzz(0.19640131429629326) q[17],q[20];
rzz(pi/8) q[17],q[19];
u3(pi/2,0.500769868982213,-0.500769868982213) q[18];
rzz(0.7853830994606742) q[17],q[18];
u3(pi,0.37133625165431355,-0.37133625165431355) q[20];
rzz(pi/32) q[16],q[20];
u3(pi,3.937043913478729,-3.937043913478729) q[19];
rzz(0.19640131429629326) q[16],q[19];
u3(pi,3.163583802164921,-3.163583802164921) q[18];
rzz(pi/8) q[16],q[18];
u3(pi/2,2.0514600027941348,-2.0514600027941348) q[17];
rzz(0.7853830994606742) q[16],q[17];
rzz(0.049100328574073315) q[15],q[20];
rzz(pi/32) q[15],q[19];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/8) q[15],q[17];
u3(pi/2,5.144043810987927,-5.144043810987927) q[16];
rzz(0.7853830994606742) q[15],q[16];
u3(pi,4.856273923919102,-4.856273923919102) q[20];
rzz(pi/128) q[14],q[20];
u3(pi,3.338884672235232,-3.338884672235232) q[19];
rzz(0.049100328574073315) q[14],q[19];
u3(pi,1.4319379315062277,-1.4319379315062277) q[18];
rzz(pi/32) q[14],q[18];
u3(pi,0.7389025921243193,-0.7389025921243193) q[17];
rzz(0.19640131429629326) q[14],q[17];
u3(pi,3.7692828657770336,-3.7692828657770336) q[16];
rzz(pi/8) q[14],q[16];
u3(pi/2,3.5927253586452874,-3.5927253586452874) q[15];
rzz(0.7853830994606742) q[14],q[15];
u3(pi,5.407309275358752,-5.407309275358752) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi,0.6716725093374978,-0.6716725093374978) q[19];
rzz(pi/128) q[13],q[19];
u3(pi,1.957212223186441,-1.957212223186441) q[18];
rzz(0.049100328574073315) q[13],q[18];
u3(pi,1.3709910340265858,-1.3709910340265858) q[17];
rzz(pi/32) q[13],q[17];
u3(pi,4.577300496280329,-4.577300496280329) q[16];
rzz(0.19640131429629326) q[13],q[16];
u3(pi,6.194592394348354,-6.194592394348354) q[15];
rzz(pi/8) q[13],q[15];
u3(pi/2,4.867583657472026,-4.867583657472026) q[14];
rzz(0.7853830994606742) q[13],q[14];
u3(pi/2,0.6214070268800611,-0.6214070268800611) q[12];
u3(pi/2,4.0752739902366795,-4.0752739902366795) q[20];
rzz(0) q[12],q[20];
u3(pi/2,3.762999680469854,-3.762999680469854) q[12];
rzz(0.012275082143518329) q[12],q[19];
rzz(pi/128) q[12],q[18];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/32) q[12],q[16];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/8) q[12],q[14];
u3(pi/2,1.8045308202219772,-1.8045308202219772) q[13];
rzz(0.7853830994606742) q[12],q[13];
u3(pi/2,0.9355662922390404,-0.9355662922390404) q[11];
rzz(0) q[11],q[20];
u3(pi,4.504415546717046,-4.504415546717046) q[11];
u3(pi/2,3.690114730906571,-3.690114730906571) q[19];
rzz(0) q[11],q[19];
u3(pi/2,1.7894511754847462,-1.7894511754847462) q[11];
u3(pi,1.1774689265654545,-1.1774689265654545) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi,0.05152211951887261,-0.05152211951887261) q[17];
rzz(pi/128) q[11],q[17];
u3(pi,3.495335986384004,-3.495335986384004) q[16];
rzz(0.049100328574073315) q[11],q[16];
u3(pi,4.293300520395811,-4.293300520395811) q[15];
rzz(pi/32) q[11],q[15];
u3(pi,2.1532476047704443,-2.1532476047704443) q[14];
rzz(0.19640131429629326) q[11],q[14];
u3(pi,0.49197340955216157,-0.49197340955216157) q[13];
rzz(pi/8) q[11],q[13];
u3(pi/2,0.6214070268800611,-0.6214070268800611) q[12];
rzz(0.7853830994606742) q[11],q[12];
u3(pi/2,2.4849997889895263,-2.4849997889895263) q[10];
rzz(0) q[10],q[20];
u3(pi,0.3286105915654924,-0.3286105915654924) q[10];
rzz(0) q[10],q[19];
u3(pi,5.441238476017522,-5.441238476017522) q[10];
u3(pi/2,4.900884539600077,-4.900884539600077) q[18];
rzz(0) q[10],q[18];
u3(pi/2,3.2848492785934877,-3.2848492785934877) q[10];
u3(pi,3.49659262344544,-3.49659262344544) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi,0.13006193585861744,-0.13006193585861744) q[16];
rzz(pi/128) q[10],q[16];
u3(pi,4.635105801106381,-4.635105801106381) q[15];
rzz(0.049100328574073315) q[10],q[15];
u3(pi,6.148725141605944,-6.148725141605944) q[14];
rzz(pi/32) q[10],q[14];
u3(pi,1.124061851454428,-1.124061851454428) q[13];
rzz(0.19640131429629326) q[10],q[13];
u3(pi,3.639220929918417,-3.639220929918417) q[12];
rzz(pi/8) q[10],q[12];
u3(pi/2,1.7894511754847462,-1.7894511754847462) q[11];
rzz(0.7853830994606742) q[10],q[11];
u3(pi/2,6.025574709585223,-6.025574709585223) q[9];
rzz(0) q[9],q[20];
u3(pi,2.7602033054439925,-2.7602033054439925) q[9];
rzz(0) q[9],q[19];
u3(pi,5.655495094992346,-5.655495094992346) q[9];
rzz(0) q[9],q[18];
u3(pi,2.2669732588303946,-2.2669732588303946) q[9];
u3(pi/2,0.7822565707438585,-0.7822565707438585) q[17];
rzz(0) q[9],q[17];
u3(pi/2,5.28478716186875,-5.28478716186875) q[9];
rzz(0.012275082143518329) q[9],q[16];
rzz(pi/128) q[9],q[15];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/32) q[9],q[13];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/8) q[9],q[11];
u3(pi/2,3.2848492785934877,-3.2848492785934877) q[10];
rzz(0.7853830994606742) q[9],q[10];
u3(pi/2,3.656813848778519,-3.656813848778519) q[8];
rzz(0) q[8],q[20];
u3(pi,0.39207076316800615,-0.39207076316800615) q[8];
rzz(0) q[8],q[19];
u3(pi,3.2867342341856416,-3.2867342341856416) q[8];
rzz(0) q[8],q[18];
u3(pi,6.181397705203277,-6.181397705203277) q[8];
rzz(0) q[8],q[17];
u3(pi,2.793504187572044,-2.793504187572044) q[8];
u3(pi/2,3.171751943064255,-3.171751943064255) q[16];
rzz(0) q[8],q[16];
u3(pi/2,5.8113180906104,-5.8113180906104) q[8];
rzz(0.012275082143518329) q[8],q[15];
rzz(pi/128) q[8],q[14];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/32) q[8],q[12];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/8) q[8],q[10];
u3(pi/2,2.143194508278957,-2.143194508278957) q[9];
rzz(0.7853830994606742) q[8],q[9];
u3(pi/2,2.4052033355883453,-2.4052033355883453) q[7];
rzz(0) q[7],q[20];
u3(pi,5.573813685999011,-5.573813685999011) q[7];
rzz(0) q[7],q[19];
u3(pi,2.4856281075202444,-2.4856281075202444) q[7];
rzz(0) q[7],q[18];
u3(pi,5.680627836221064,-5.680627836221064) q[7];
rzz(0) q[7],q[17];
u3(pi,2.5924422577422974,-2.5924422577422974) q[7];
rzz(0) q[7],q[16];
u3(pi,5.787441986443117,-5.787441986443117) q[7];
u3(pi/2,3.1975130028236918,-3.1975130028236918) q[15];
rzz(0) q[7],q[15];
u3(pi/2,2.672238711143478,-2.672238711143478) q[7];
u3(pi,3.768654547246316,-3.768654547246316) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi,2.112406900273777,-2.112406900273777) q[13];
rzz(pi/128) q[7],q[13];
u3(pi,0.7087433026498573,-0.7087433026498573) q[12];
rzz(0.049100328574073315) q[7],q[12];
u3(pi,0.5397256178867265,-0.5397256178867265) q[11];
rzz(pi/32) q[7],q[11];
u3(pi,0.5705132258919065,-0.5705132258919065) q[10];
rzz(0.19640131429629326) q[7],q[10];
u3(pi,5.161636729848031,-5.161636729848031) q[9];
rzz(pi/8) q[7],q[9];
u3(pi/2,2.669725437020606,-2.669725437020606) q[8];
rzz(0.7853830994606742) q[7],q[8];
u3(pi/2,2.061513099285622,-2.061513099285622) q[6];
rzz(0) q[6],q[20];
u3(pi,5.537999529748087,-5.537999529748087) q[6];
rzz(0) q[6],q[19];
u3(pi,3.0655661113729202,-3.0655661113729202) q[6];
rzz(0) q[6],q[18];
u3(pi,0.593132692997753,-0.593132692997753) q[6];
rzz(0) q[6],q[17];
u3(pi,4.4038845818021715,-4.4038845818021715) q[6];
rzz(0) q[6],q[16];
u3(pi,1.930822844896287,-1.930822844896287) q[6];
rzz(0) q[6],q[15];
u3(pi,5.7415747337007055,-5.7415747337007055) q[6];
u3(pi/2,0.9613273519984766,-0.9613273519984766) q[14];
rzz(0) q[6],q[14];
u3(pi/2,2.934875856983585,-2.934875856983585) q[6];
rzz(0.012275082143518329) q[6],q[13];
rzz(pi/128) q[6],q[12];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/32) q[6],q[10];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/8) q[6],q[8];
u3(pi/2,5.813831364733272,-5.813831364733272) q[7];
rzz(0.7853830994606742) q[6],q[7];
u3(pi/2,3.3376280351737964,-3.3376280351737964) q[5];
rzz(0) q[5],q[20];
u3(pi,0.623291982472215,-0.623291982472215) q[5];
rzz(0) q[5],q[19];
u3(pi,4.618769519307714,-4.618769519307714) q[5];
rzz(0) q[5],q[18];
u3(pi,2.3310617489636263,-2.3310617489636263) q[5];
rzz(0) q[5],q[17];
u3(pi,0.043353978619539144,-0.043353978619539144) q[5];
rzz(0) q[5],q[16];
u3(pi,4.038831515455039,-4.038831515455039) q[5];
rzz(0) q[5],q[15];
u3(pi,1.7511237451109507,-1.7511237451109507) q[5];
rzz(0) q[5],q[14];
u3(pi,5.746601281946449,-5.746601281946449) q[5];
u3(pi/2,1.154849459459608,-1.154849459459608) q[13];
rzz(0) q[5],q[13];
u3(pi/2,3.03163691071415,-3.03163691071415) q[5];
rzz(0.012275082143518329) q[5],q[12];
rzz(pi/128) q[5],q[11];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/32) q[5],q[9];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/8) q[5],q[7];
u3(pi/2,6.076468510573378,-6.076468510573378) q[6];
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
u3(pi/2,4.1852297331123225,-4.1852297331123225) q[12];
rzz(0) q[4],q[12];
u3(pi/2,4.828627908567512,-4.828627908567512) q[4];
u3(pi,1.1184069846779663,-1.1184069846779663) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi,5.235149997942032,-5.235149997942032) q[10];
rzz(pi/128) q[4],q[10];
u3(pi,2.231159102579471,-2.231159102579471) q[9];
rzz(0.049100328574073315) q[4],q[9];
u3(pi,5.272220791254391,-5.272220791254391) q[8];
rzz(pi/32) q[4],q[8];
u3(pi,4.5553093477052,-4.5553093477052) q[7];
rzz(0.19640131429629326) q[4],q[7];
u3(pi,2.834973210599429,-2.834973210599429) q[6];
rzz(pi/8) q[4],q[6];
u3(pi/2,3.03163691071415,-3.03163691071415) q[5];
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
u3(pi/2,6.088406562657019,-6.088406562657019) q[11];
rzz(0) q[3],q[11];
u3(pi/2,5.997300375702915,-5.997300375702915) q[3];
u3(pi,4.174548318090117,-4.174548318090117) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi,6.13427381539943,-6.13427381539943) q[9];
rzz(pi/128) q[3],q[9];
u3(pi,2.189061761021368,-2.189061761021368) q[8];
rzz(0.049100328574073315) q[3],q[8];
u3(pi,5.2973535324831085,-5.2973535324831085) q[7];
rzz(pi/32) q[3],q[7];
u3(pi,5.903052596095221,-5.903052596095221) q[6];
rzz(0.19640131429629326) q[3],q[6];
u3(pi,1.7781414419318227,-1.7781414419318227) q[5];
rzz(pi/8) q[3],q[5];
u3(pi/2,1.687035254977719,-1.687035254977719) q[4];
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
u3(pi/2,2.018159120666083,-2.018159120666083) q[10];
rzz(0) q[2],q[10];
u3(pi/2,0.684238879951857,-0.684238879951857) q[2];
rzz(0.012275082143518329) q[2],q[9];
rzz(pi/128) q[2],q[8];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/32) q[2],q[6];
rzz(0.19640131429629326) q[2],q[5];
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
u3(pi/2,3.419937762697849,-3.419937762697849) q[9];
rzz(0) q[1],q[9];
u3(pi/2,3.8390262226867273,-3.8390262226867273) q[1];
u3(pi,2.159530790077624,-2.159530790077624) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi,3.61660146281257,-3.61660146281257) q[7];
rzz(pi/128) q[1],q[7];
u3(pi,4.4648314792818145,-4.4648314792818145) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi,5.825769416816913,-5.825769416816913) q[5];
rzz(pi/32) q[1],q[5];
u3(pi,0.7282211771021141,-0.7282211771021141) q[4];
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
u3(pi/2,1.532468896421101,-1.532468896421101) q[8];
rzz(pi/2) q[0],q[8];
u3(pi/2,6.218468498515636,-6.218468498515636) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,3.0002209841782523,-3.0002209841782523) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,1.7027432182456679,-1.7027432182456679) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,6.051964087875377,-6.051964087875377) q[4];
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
u3(pi/2,4.4811677610804805,-4.4811677610804805) q[4];
u3(pi/2,1.1435397259066846,-1.1435397259066846) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,0.13194689145077132,-0.13194689145077132) q[5];
u3(pi/2,3.1748935357178447,-3.1748935357178447) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,1.429424657383356,-1.429424657383356) q[6];
u3(pi/2,4.522008465577148,-4.522008465577148) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,1.5060795181309468,-1.5060795181309468) q[7];
u3(pi/2,4.623167749022739,-4.623167749022739) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,6.244857876805791,-6.244857876805791) q[8];
u3(pi/2,3.0913271711323564,-3.0913271711323564) q[8];
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
u3(pi/2,2.714336052701581,-2.714336052701581) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi/2,4.745689862512742,-4.745689862512742) q[5];
rzz(pi/32) q[1],q[5];
u3(pi/2,6.0928047923720445,-6.0928047923720445) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi/2,3.052371422227843,-3.052371422227843) q[7];
rzz(pi/128) q[1],q[7];
u3(pi/2,4.662123497927253,-4.662123497927253) q[8];
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
u3(pi/2,3.4073713920834896,-3.4073713920834896) q[9];
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
u3(pi/2,2.5855307539043997,-2.5855307539043997) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,2.003079475928852,-2.003079475928852) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,1.0147344271095031,-1.0147344271095031) q[3];
rzz(0) q[3],q[11];
u3(pi,4.583583681587508,-4.583583681587508) q[3];
rzz(0) q[3],q[12];
u3(pi,2.295875911243421,-2.295875911243421) q[3];
rzz(0) q[3],q[13];
u3(pi,0.008168140899333461,-0.008168140899333461) q[3];
rzz(0) q[3],q[14];
u3(pi,4.003017359204114,-4.003017359204114) q[3];
rzz(0) q[3],q[15];
u3(pi,1.7153095888600272,-1.7153095888600272) q[3];
rzz(0) q[3],q[16];
u3(pi,5.710787125695526,-5.710787125695526) q[3];
rzz(0) q[3],q[17];
u3(pi,3.4230793553514385,-3.4230793553514385) q[3];
rzz(0) q[3],q[18];
u3(pi,1.1353715850073511,-1.1353715850073511) q[3];
rzz(0) q[3],q[19];
u3(pi,5.13084912184285,-5.13084912184285) q[3];
rzz(0) q[3],q[20];
u3(pi/2,2.5176723525868603,-2.5176723525868603) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,6.074583554981224,-6.074583554981224) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,4.088468679381757,-4.088468679381757) q[4];
rzz(0) q[4],q[12];
u3(pi,0.46872562391559713,-0.46872562391559713) q[4];
rzz(0) q[4],q[13];
u3(pi,2.6533891552219395,-2.6533891552219395) q[4];
rzz(0) q[4],q[14];
u3(pi,4.838052686528282,-4.838052686528282) q[4];
rzz(0) q[4],q[15];
u3(pi,0.7395309106550373,-0.7395309106550373) q[4];
rzz(0) q[4],q[16];
u3(pi,2.9235661234306614,-2.9235661234306614) q[4];
rzz(0) q[4],q[17];
u3(pi,5.108229654737003,-5.108229654737003) q[4];
rzz(0) q[4],q[18];
u3(pi,1.0097078788637597,-1.0097078788637597) q[4];
rzz(0) q[4],q[19];
u3(pi,3.1943714101701013,-3.1943714101701013) q[4];
rzz(0) q[4],q[20];
u3(pi/2,1.5060795181309468,-1.5060795181309468) q[5];
rzz(0.7853830994606742) q[5],q[6];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,1.0260441606624264,-1.0260441606624264) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,6.218468498515636,-6.218468498515636) q[5];
rzz(0) q[5],q[13];
u3(pi,2.976973198541688,-2.976973198541688) q[5];
rzz(0) q[5],q[14];
u3(pi,5.918132240832452,-5.918132240832452) q[5];
rzz(0) q[5],q[15];
u3(pi,2.5767342944743485,-2.5767342944743485) q[5];
rzz(0) q[5],q[16];
u3(pi,5.518521655295831,-5.518521655295831) q[5];
rzz(0) q[5],q[17];
u3(pi,2.1764953904070086,-2.1764953904070086) q[5];
rzz(0) q[5],q[18];
u3(pi,5.118282751228491,-5.118282751228491) q[5];
rzz(0) q[5],q[19];
u3(pi,1.776884804870387,-1.776884804870387) q[5];
rzz(0) q[5],q[20];
u3(pi/2,2.902203293386251,-2.902203293386251) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,1.137884859130223,-1.137884859130223) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,4.472999620181147,-4.472999620181147) q[6];
rzz(0) q[6],q[14];
u3(pi,1.913229926036184,-1.913229926036184) q[6];
rzz(0) q[6],q[15];
u3(pi,6.218468498515636,-6.218468498515636) q[6];
rzz(0) q[6],q[16];
u3(pi,4.241150082346221,-4.241150082346221) q[6];
rzz(0) q[6],q[17];
u3(pi,2.263203347646087,-2.263203347646087) q[6];
rzz(0) q[6],q[18];
u3(pi,0.28588493147667116,-0.28588493147667116) q[6];
rzz(0) q[6],q[19];
u3(pi,4.591123503956124,-4.591123503956124) q[6];
rzz(0) q[6],q[20];
u3(pi/2,6.1694596531196355,-6.1694596531196355) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,0.9431061146076559,-0.9431061146076559) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,4.598663326324739,-4.598663326324739) q[7];
rzz(0) q[7],q[15];
u3(pi,0.5560618996853933,-0.5560618996853933) q[7];
rzz(0) q[7],q[16];
u3(pi,1.8956370071760813,-1.8956370071760813) q[7];
rzz(0) q[7],q[17];
u3(pi,3.235212114666769,-3.235212114666769) q[7];
rzz(0) q[7],q[18];
u3(pi,4.574787222157457,-4.574787222157457) q[7];
rzz(0) q[7],q[19];
u3(pi,5.914362329648145,-5.914362329648145) q[7];
rzz(0) q[7],q[20];
u3(pi/2,4.649557127312894,-4.649557127312894) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,3.1792917654328705,-3.1792917654328705) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,6.220353454107791,-6.220353454107791) q[8];
rzz(0) q[8],q[16];
u3(pi,2.096698937005828,-2.096698937005828) q[8];
rzz(0) q[8],q[17];
u3(pi,3.2747961821020004,-3.2747961821020004) q[8];
rzz(0) q[8],q[18];
u3(pi,4.452265108667455,-4.452265108667455) q[8];
rzz(0) q[8],q[19];
u3(pi,5.630362353763627,-5.630362353763627) q[8];
rzz(0) q[8],q[20];
u3(pi/2,3.3954333399998484,-3.3954333399998484) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,3.1535307056734343,-3.1535307056734343) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,1.8246370132049519,-1.8246370132049519) q[9];
rzz(0) q[9],q[17];
u3(pi,5.963999493574864,-5.963999493574864) q[9];
rzz(0) q[9],q[18];
u3(pi,4.818574812076025,-4.818574812076025) q[9];
rzz(0) q[9],q[19];
u3(pi,3.673150130577186,-3.673150130577186) q[9];
rzz(0) q[9],q[20];
u3(pi/2,1.993654697968083,-1.993654697968083) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,3.9049996684121133,-3.9049996684121133) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi/2,3.567592617416569,-3.567592617416569) q[10];
rzz(0) q[10],q[18];
u3(pi,6.170087971650354,-6.170087971650354) q[10];
rzz(0) q[10],q[19];
u3(pi,1.9490440822871076,-1.9490440822871076) q[10];
rzz(0) q[10],q[20];
u3(pi/2,2.9254510790228156,-2.9254510790228156) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,4.882034983678539,-4.882034983678539) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,1.3559113892893546,-1.3559113892893546) q[11];
rzz(0) q[11],q[19];
u3(pi,3.596495269829595,-3.596495269829595) q[11];
rzz(0) q[11],q[20];
u3(pi/2,1.019132656824529,-1.019132656824529) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/128) q[12],q[18];
u3(pi/2,0.5303008399259571,-0.5303008399259571) q[19];
rzz(0.012275082143518329) q[12],q[19];
u3(pi/2,2.5905573021501436,-2.5905573021501436) q[12];
rzz(0) q[12],q[20];
u3(pi/2,1.1316016738230434,-1.1316016738230434) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/128) q[13],q[19];
u3(pi/2,0.9154600992560656,-0.9154600992560656) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi/2,4.078415582890269,-4.078415582890269) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/128) q[14],q[20];
u3(pi/2,3.173008580125691,-3.173008580125691) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/32) q[15],q[19];
rzz(0.049100328574073315) q[15],q[20];
u3(pi/2,3.1472475203662547,-3.1472475203662547) q[16];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/8) q[16],q[18];
rzz(0.19640131429629326) q[16],q[19];
rzz(pi/32) q[16],q[20];
u3(pi/2,0.7577521480458581,-0.7577521480458581) q[17];
rzz(0.7853830994606742) q[17],q[18];
rzz(pi/8) q[17],q[19];
rzz(0.19640131429629326) q[17],q[20];
u3(pi/2,1.7347874633122837,-1.7347874633122837) q[18];
rzz(0.7853830994606742) q[18],q[19];
rzz(pi/8) q[18],q[20];
u3(pi/2,3.6656103082085707,-3.6656103082085707) q[19];
rzz(0.7853830994606742) q[19],q[20];
u3(pi,0,0) q[0];
u3(pi/2,2.236185650825215,-2.236185650825215) q[1];
u3(pi/2,2.9399024052293283,-2.9399024052293283) q[2];
u3(pi/2,2.416513069141269,-2.416513069141269) q[3];
u3(pi/2,5.857813661883529,-5.857813661883529) q[4];
u3(pi/2,4.817946493545307,-4.817946493545307) q[5];
u3(pi/2,2.03135380981116,-2.03135380981116) q[6];
u3(pi/2,1.8717609030087987,-1.8717609030087987) q[7];
u3(pi/2,1.5067078366616649,-1.5067078366616649) q[8];
u3(pi/2,1.5293273037675112,-1.5293273037675112) q[9];
u3(pi/2,4.551539436520892,-4.551539436520892) q[10];
u3(pi/2,5.837079150369836,-5.837079150369836) q[11];
u3(pi/2,5.732149955739937,-5.732149955739937) q[12];
u3(pi/2,4.05076956753868,-4.05076956753868) q[20];
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
