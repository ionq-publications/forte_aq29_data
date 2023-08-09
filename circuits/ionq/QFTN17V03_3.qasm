OPENQASM 2.0;
include "qelib1.inc";
qreg q[17];
creg c[17];
u3(pi/2,0.06220353454107791,-0.06220353454107791) q[16];
u3(pi/2,2.4183980247334227,-2.4183980247334227) q[16];
rzz(pi/2) q[15],q[16];
u3(pi/2,0.8476016979385261,-0.8476016979385261) q[16];
rzz(pi/8) q[14],q[16];
u3(pi/2,1.7253626853515145,-1.7253626853515145) q[15];
rzz(0.7853830994606742) q[14],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/8) q[13],q[15];
u3(pi/2,0.2475575011028757,-0.2475575011028757) q[14];
rzz(0.7853830994606742) q[13],q[14];
u3(pi,4.9743978076940785,-4.9743978076940785) q[16];
rzz(pi/32) q[12],q[16];
u3(pi,0.4668406683234433,-0.4668406683234433) q[15];
rzz(0.19640131429629326) q[12],q[15];
u3(pi,2.8494245368059423,-2.8494245368059423) q[14];
rzz(pi/8) q[12],q[14];
u3(pi/2,0.6176371156957533,-0.6176371156957533) q[13];
rzz(0.7853830994606742) q[12],q[13];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/32) q[11],q[15];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/8) q[11],q[13];
u3(pi/2,4.131194339470578,-4.131194339470578) q[12];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/128) q[10],q[16];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/32) q[10],q[14];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/8) q[10],q[12];
u3(pi/2,2.4699201442522956,-2.4699201442522956) q[11];
rzz(0.7853830994606742) q[10],q[11];
rzz(0.012275082143518329) q[9],q[16];
rzz(pi/128) q[9],q[15];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/32) q[9],q[13];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/8) q[9],q[11];
u3(pi/2,0.8161857714026282,-0.8161857714026282) q[10];
rzz(0.7853830994606742) q[9],q[10];
u3(pi/2,3.264114767079795,-3.264114767079795) q[8];
u3(pi/2,2.8180086102700446,-2.8180086102700446) q[16];
rzz(0) q[8],q[16];
u3(pi/2,0.12252211349000193,-0.12252211349000193) q[8];
u3(pi,5.517893336765113,-5.517893336765113) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi,1.0562034501368884,-1.0562034501368884) q[14];
rzz(pi/128) q[8],q[14];
u3(pi,4.185858051643041,-4.185858051643041) q[13];
rzz(0.049100328574073315) q[8],q[13];
u3(pi,2.8180086102700446,-2.8180086102700446) q[12];
rzz(pi/32) q[8],q[12];
u3(pi,6.0381410801995825,-6.0381410801995825) q[11];
rzz(0.19640131429629326) q[8],q[11];
u3(pi,5.480822543452753,-5.480822543452753) q[10];
rzz(pi/8) q[8],q[10];
u3(pi/2,0.4542742977090841,-0.4542742977090841) q[9];
rzz(0.7853830994606742) q[8],q[9];
u3(pi/2,3.2151059216837945,-3.2151059216837945) q[7];
rzz(0) q[7],q[16];
u3(pi,0.500769868982213,-0.500769868982213) q[7];
u3(pi/2,2.402690061465474,-2.402690061465474) q[15];
rzz(0) q[7],q[15];
u3(pi/2,4.068990804929499,-4.068990804929499) q[7];
rzz(0.012275082143518329) q[7],q[14];
rzz(pi/128) q[7],q[13];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/32) q[7],q[11];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/8) q[7],q[9];
u3(pi/2,3.264114767079795,-3.264114767079795) q[8];
rzz(0.7853830994606742) q[7],q[8];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[6];
rzz(0) q[6],q[16];
u3(pi,4.623167749022739,-4.623167749022739) q[6];
rzz(0) q[6],q[15];
u3(pi,2.6772652593892214,-2.6772652593892214) q[6];
u3(pi/2,6.085893288534147,-6.085893288534147) q[14];
rzz(0) q[6],q[14];
u3(pi/2,0.1338318470429252,-0.1338318470429252) q[6];
rzz(0.012275082143518329) q[6],q[13];
rzz(pi/128) q[6],q[12];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/32) q[6],q[10];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/8) q[6],q[8];
u3(pi/2,0.9273981513397069,-0.9273981513397069) q[7];
rzz(0.7853830994606742) q[6],q[7];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[5];
rzz(0) q[5],q[16];
u3(pi,4.292672201865093,-4.292672201865093) q[5];
rzz(0) q[5],q[15];
u3(pi,6.2021322167169695,-6.2021322167169695) q[5];
rzz(0) q[5],q[14];
u3(pi,1.8284069243892596,-1.8284069243892596) q[5];
u3(pi/2,1.471521998941459,-1.471521998941459) q[13];
rzz(0) q[5],q[13];
u3(pi/2,4.354247417875453,-4.354247417875453) q[5];
rzz(0.012275082143518329) q[5],q[12];
rzz(pi/128) q[5],q[11];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/32) q[5],q[9];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/8) q[5],q[7];
u3(pi/2,3.275424500632718,-3.275424500632718) q[6];
rzz(0.7853830994606742) q[5],q[6];
u3(pi/2,5*pi/8,-5*pi/8) q[4];
rzz(0) q[4],q[16];
u3(pi,4.302096979825863,-4.302096979825863) q[4];
rzz(0) q[4],q[15];
u3(pi,5.838335787431272,-5.838335787431272) q[4];
rzz(0) q[4],q[14];
u3(pi,1.091389287857094,-1.091389287857094) q[4];
rzz(0) q[4],q[13];
u3(pi,2.627628095462503,-2.627628095462503) q[4];
u3(pi/2,1.5048228810695108,-1.5048228810695108) q[12];
rzz(0) q[4],q[12];
u3(pi/2,4.966857985325463,-4.966857985325463) q[4];
u3(pi,5.1333623959657215,-5.1333623959657215) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi,3.432504133312208,-3.432504133312208) q[10];
rzz(pi/128) q[4],q[10];
u3(pi,5.329397777549724,-5.329397777549724) q[9];
rzz(0.049100328574073315) q[4],q[9];
u3(pi,5.154725226010132,-5.154725226010132) q[8];
rzz(pi/32) q[4],q[8];
u3(pi,6.0142649760323,-6.0142649760323) q[7];
rzz(0.19640131429629326) q[4],q[7];
u3(pi,0.7319910882864218,-0.7319910882864218) q[6];
rzz(pi/8) q[4],q[6];
u3(pi/2,4.354247417875453,-4.354247417875453) q[5];
rzz(0.7853830994606742) q[4],q[5];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(0) q[3],q[16];
u3(pi,1.2120264457549421,-1.2120264457549421) q[3];
rzz(0) q[3],q[15];
u3(pi,5.207503982590441,-5.207503982590441) q[3];
rzz(0) q[3],q[14];
u3(pi,2.9197962122463537,-2.9197962122463537) q[3];
rzz(0) q[3],q[13];
u3(pi,0.6320884419022663,-0.6320884419022663) q[3];
rzz(0) q[3],q[12];
u3(pi,4.627565978737765,-4.627565978737765) q[3];
u3(pi/2,3.8013271108436495,-3.8013271108436495) q[11];
rzz(0) q[3],q[11];
u3(pi/2,1.913229926036184,-1.913229926036184) q[3];
rzz(0.012275082143518329) q[3],q[10];
rzz(pi/128) q[3],q[9];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/32) q[3],q[7];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/8) q[3],q[5];
u3(pi/2,4.966857985325463,-4.966857985325463) q[4];
rzz(0.7853830994606742) q[3],q[4];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(0) q[2],q[16];
u3(pi,5.1390172627421835,-5.1390172627421835) q[2];
rzz(0) q[2],q[15];
u3(pi,2.851309492398096,-2.851309492398096) q[2];
rzz(0) q[2],q[14];
u3(pi,0.5636017220540089,-0.5636017220540089) q[2];
rzz(0) q[2],q[13];
u3(pi,4.559079258889508,-4.559079258889508) q[2];
rzz(0) q[2],q[12];
u3(pi,2.2713714885454204,-2.2713714885454204) q[2];
rzz(0) q[2],q[11];
u3(pi,6.266849025380919,-6.266849025380919) q[2];
u3(pi/2,6.144955230421635,-6.144955230421635) q[10];
rzz(0) q[2],q[10];
u3(pi/2,3.552512972679338,-3.552512972679338) q[2];
rzz(0.012275082143518329) q[2],q[9];
rzz(pi/128) q[2],q[8];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/32) q[2],q[6];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/8) q[2],q[4];
u3(pi/2,1.913229926036184,-1.913229926036184) q[3];
rzz(0.7853830994606742) q[2],q[3];
u3(pi/2,pi,-pi) q[1];
rzz(0) q[1],q[16];
u3(pi,5.74345968929286,-5.74345968929286) q[1];
rzz(0) q[1],q[15];
u3(pi,1.5230441184603318,-1.5230441184603318) q[1];
rzz(0) q[1],q[14];
u3(pi,3.5858138548073897,-3.5858138548073897) q[1];
rzz(0) q[1],q[13];
u3(pi,5.648583591154448,-5.648583591154448) q[1];
rzz(0) q[1],q[12];
u3(pi,1.4275397017912022,-1.4275397017912022) q[1];
rzz(0) q[1],q[11];
u3(pi,3.4903094381382602,-3.4903094381382602) q[1];
rzz(0) q[1],q[10];
u3(pi,5.553079174485318,-5.553079174485318) q[1];
u3(pi/2,3.9207076316800618,-3.9207076316800618) q[9];
rzz(0) q[1],q[9];
u3(pi/2,1.8717609030087987,-1.8717609030087987) q[1];
u3(pi,6.159406556628148,-6.159406556628148) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi,0.3022212132753381,-0.3022212132753381) q[7];
rzz(pi/128) q[1],q[7];
u3(pi,3.9320173652329853,-3.9320173652329853) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi,3.397318295592002,-3.397318295592002) q[5];
rzz(pi/32) q[1],q[5];
u3(pi,2.251893614093164,-2.251893614093164) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi,0.6000441968356505,-0.6000441968356505) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,0.41092031908954496,-0.41092031908954496) q[2];
rzz(0.7853830994606742) q[1],q[2];
rzz(pi/2) q[0],q[16];
rzz(pi/2) q[0],q[15];
rzz(-pi/2) q[0],q[14];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[11];
rzz(pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[9];
u3(pi/2,2.1318847747260334,-2.1318847747260334) q[8];
rzz(pi/2) q[0],q[8];
rzz(0.012275082143518329) q[0],q[7];
u3(pi/2,0.2506990937564655,-0.2506990937564655) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,2.4410174918392693,-2.4410174918392693) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,5.820742868571169,-5.820742868571169) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,5.570043774814703,-5.570043774814703) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3.552512972679338,-3.552512972679338) q[2];
rzz(pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi,0,0) q[0];
u3(pi/2,1.8717609030087987,-1.8717609030087987) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,1.9817166458844415,-1.9817166458844415) q[2];
u3(pi/2,4.337911136076786,-4.337911136076786) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,3.999247448019806,-3.999247448019806) q[3];
u3(pi/2,0.4649557127312894,-0.4649557127312894) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,1.108353888186479,-1.108353888186479) q[4];
u3(pi/2,4.053282841661551,-4.053282841661551) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,0.8702211650443727,-0.8702211650443727) q[5];
u3(pi/2,3.9131678093114464,-3.9131678093114464) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,4.963088074141155,-4.963088074141155) q[6];
u3(pi/2,1.7724865751553613,-1.7724865751553613) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi,0.49825659485934115,-0.49825659485934115) q[7];
rzz(0.012275082143518329) q[0],q[7];
u3(pi/2,0.5610884479311371,-0.5610884479311371) q[8];
u3(pi/2,3.6907430494372893,-3.6907430494372893) q[8];
rzz(pi/2) q[0],q[8];
rzz(pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[11];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[14];
rzz(-pi/2) q[0],q[15];
rzz(pi/2) q[0],q[16];
u3(pi/2,3.4425572298036955,-3.4425572298036955) q[1];
u3(pi/2,5.908707462871683,-5.908707462871683) q[2];
rzz(0.7853830994606742) q[1],q[2];
u3(pi/2,5.177344693115979,-5.177344693115979) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,2.4824865148666544,-2.4824865148666544) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi/2,2.34237148251655,-2.34237148251655) q[5];
rzz(pi/32) q[1],q[5];
u3(pi/2,3.343282901950258,-3.343282901950258) q[6];
rzz(0.049100328574073315) q[1],q[6];
rzz(pi/128) q[1],q[7];
u3(pi/2,2.1199467226423923,-2.1199467226423923) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi/2,3.4425572298036955,-3.4425572298036955) q[1];
rzz(0) q[1],q[9];
u3(pi,0.7282211771021141,-0.7282211771021141) q[1];
rzz(0) q[1],q[10];
u3(pi,4.723698713937613,-4.723698713937613) q[1];
rzz(0) q[1],q[11];
u3(pi,2.4359909435935254,-2.4359909435935254) q[1];
rzz(0) q[1],q[12];
u3(pi,0.14828317324943824,-0.14828317324943824) q[1];
rzz(0) q[1],q[13];
u3(pi,4.143132391554219,-4.143132391554219) q[1];
rzz(0) q[1],q[14];
u3(pi,1.855424621210132,-1.855424621210132) q[1];
rzz(0) q[1],q[15];
u3(pi,5.850902158045631,-5.850902158045631) q[1];
rzz(0) q[1],q[16];
u3(pi/2,5.123309299474235,-5.123309299474235) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,3.9087695795964206,-3.9087695795964206) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi/2,0.41092031908954496,-0.41092031908954496) q[2];
rzz(0) q[2],q[10];
u3(pi,3.979141255036832,-3.979141255036832) q[2];
rzz(0) q[2],q[11];
u3(pi,1.6914334846927446,-1.6914334846927446) q[2];
rzz(0) q[2],q[12];
u3(pi,5.6869110215282435,-5.6869110215282435) q[2];
rzz(0) q[2],q[13];
u3(pi,3.3992032511841566,-3.3992032511841566) q[2];
rzz(0) q[2],q[14];
u3(pi,1.1114954808400688,-1.1114954808400688) q[2];
rzz(0) q[2],q[15];
u3(pi,5.106973017675568,-5.106973017675568) q[2];
rzz(0) q[2],q[16];
u3(pi/2,4.784645611417255,-4.784645611417255) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,2.994566117401791,-2.994566117401791) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,3.213849284622358,-3.213849284622358) q[3];
rzz(0) q[3],q[11];
u3(pi,0.4995132319207771,-0.4995132319207771) q[3];
rzz(0) q[3],q[12];
u3(pi,4.494990768756276,-4.494990768756276) q[3];
rzz(0) q[3],q[13];
u3(pi,2.2072829984121887,-2.2072829984121887) q[3];
rzz(0) q[3],q[14];
u3(pi,6.202760535247688,-6.202760535247688) q[3];
rzz(0) q[3],q[15];
u3(pi,3.914424446372882,-3.914424446372882) q[3];
rzz(0) q[3],q[16];
u3(pi/2,2.2864511332826516,-2.2864511332826516) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,3.7843625105142644,-3.7843625105142644) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,3.857247460077548,-3.857247460077548) q[4];
rzz(0) q[4],q[12];
u3(pi,1.3138140477312514,-1.3138140477312514) q[4];
rzz(0) q[4],q[13];
u3(pi,5.65109686527732,-5.65109686527732) q[4];
rzz(0) q[4],q[14];
u3(pi,3.705194375643802,-3.705194375643802) q[4];
rzz(0) q[4],q[15];
u3(pi,1.7599202045410023,-1.7599202045410023) q[4];
rzz(0) q[4],q[16];
u3(pi/2,5.385946445314341,-5.385946445314341) q[5];
rzz(0.7853830994606742) q[5],q[6];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,1.489114917801562,-1.489114917801562) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,3.8151501185194445,-3.8151501185194445) q[5];
rzz(0) q[5],q[13];
u3(pi,0.057805304826052194,-0.057805304826052194) q[5];
rzz(0) q[5],q[14];
u3(pi,1.9672653196779284,-1.9672653196779284) q[5];
rzz(0) q[5],q[15];
u3(pi,3.8767253345298047,-3.8767253345298047) q[5];
rzz(0) q[5],q[16];
u3(pi/2,0.15268140296446395,-0.15268140296446395) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,1.4533007615506384,-1.4533007615506384) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,1.7234777297593604,-1.7234777297593604) q[6];
rzz(0) q[6],q[14];
u3(pi,4.435928826868787,-4.435928826868787) q[6];
rzz(0) q[6],q[15];
u3(pi,0.4360530603182633,-0.4360530603182633) q[6];
rzz(0) q[6],q[16];
u3(pi/2,5.186141152546031,-5.186141152546031) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,2.9267077160842514,-2.9267077160842514) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,3.615344825751134,-3.615344825751134) q[7];
rzz(0) q[7],q[15];
u3(pi,1.533725533482537,-1.533725533482537) q[7];
rzz(0) q[7],q[16];
u3(pi/2,5.2489730056178265,-5.2489730056178265) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,5.526061477664446,-5.526061477664446) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,0.5365840252331366,-0.5365840252331366) q[8];
rzz(0) q[8],q[16];
u3(pi/2,3.8962032089820613,-3.8962032089820613) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,5.941380026469017,-5.941380026469017) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,2.9851413394410216,-2.9851413394410216) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,3.776822688145649,-3.776822688145649) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
u3(pi/2,1.4822034139636644,-1.4822034139636644) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
u3(pi/2,4.588610229833251,-4.588610229833251) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
u3(pi/2,6.062017184366865,-6.062017184366865) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
u3(pi/2,2.3781856387674734,-2.3781856387674734) q[15];
rzz(0.7853830994606742) q[15],q[16];
u3(pi/2,3.1365661053440492,-3.1365661053440492) q[1];
u3(pi/2,2.3920086464432684,-2.3920086464432684) q[2];
u3(pi/2,1.200088393671301,-1.200088393671301) q[3];
u3(pi/2,5.499043780843574,-5.499043780843574) q[4];
u3(pi/2,0.11875220230569418,-0.11875220230569418) q[5];
u3(pi/2,3.1485041574276904,-3.1485041574276904) q[6];
u3(pi/2,5.735291548393526,-5.735291548393526) q[7];
u3(pi/2,3.67817667882293,-3.67817667882293) q[8];
u3(pi/2,5.935096841161837,-5.935096841161837) q[16];
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
