OPENQASM 2.0;
include "qelib1.inc";
qreg q[24];
creg c[24];
u3(pi/2,4.86632702041059,-4.86632702041059) q[23];
u3(pi/2,0.9393362034233481,-0.9393362034233481) q[23];
rzz(-pi/2) q[22],q[23];
u3(pi/2,2.510132530218245,-2.510132530218245) q[23];
rzz(pi/8) q[21],q[23];
u3(pi/2,3.389150154692669,-3.389150154692669) q[22];
rzz(0.7853830994606742) q[21],q[22];
u3(pi,1.1969468010177111,-1.1969468010177111) q[23];
rzz(0.19640131429629326) q[20],q[23];
u3(pi,6.101601251802096,-6.101601251802096) q[22];
rzz(pi/8) q[20],q[22];
u3(pi/2,3.7573448136933925,-3.7573448136933925) q[21];
rzz(0.7853830994606742) q[20],q[21];
rzz(pi/32) q[19],q[23];
rzz(0.19640131429629326) q[19],q[22];
rzz(pi/8) q[19],q[21];
u3(pi/2,0.9896016858807848,-0.9896016858807848) q[20];
rzz(0.7853830994606742) q[19],q[20];
rzz(0.049100328574073315) q[18],q[23];
rzz(pi/32) q[18],q[22];
rzz(0.19640131429629326) q[18],q[21];
rzz(pi/8) q[18],q[20];
u3(pi/2,5.603972975473473,-5.603972975473473) q[19];
rzz(0.7853830994606742) q[18],q[19];
rzz(pi/128) q[17],q[23];
rzz(0.049100328574073315) q[17],q[22];
rzz(pi/32) q[17],q[21];
rzz(0.19640131429629326) q[17],q[20];
rzz(pi/8) q[17],q[19];
u3(pi/2,3.958406743523139,-3.958406743523139) q[18];
rzz(0.7853830994606742) q[17],q[18];
u3(pi,5.522291566480138,-5.522291566480138) q[23];
rzz(0.012275082143518329) q[16],q[23];
u3(pi,4.359273966121196,-4.359273966121196) q[22];
rzz(pi/128) q[16],q[22];
u3(pi,0.18661060362323373,-0.18661060362323373) q[21];
rzz(0.049100328574073315) q[16],q[21];
u3(pi,0.03330088212805181,-0.03330088212805181) q[20];
rzz(pi/32) q[16],q[20];
u3(pi,1.561371548834127,-1.561371548834127) q[19];
rzz(0.19640131429629326) q[16],q[19];
u3(pi,2.442902447431423,-2.442902447431423) q[18];
rzz(pi/8) q[16],q[18];
u3(pi/2,3.566964298885851,-3.566964298885851) q[17];
rzz(0.7853830994606742) q[16],q[17];
u3(pi/2,4.4258757303773,-4.4258757303773) q[15];
u3(pi/2,1.7360441003737197,-1.7360441003737197) q[23];
rzz(0) q[15],q[23];
u3(pi/2,1.2842830767875073,-1.2842830767875073) q[15];
rzz(0.012275082143518329) q[15],q[22];
rzz(pi/128) q[15],q[21];
rzz(0.049100328574073315) q[15],q[20];
rzz(pi/32) q[15],q[19];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/8) q[15],q[17];
u3(pi/2,0.1262920246743097,-0.1262920246743097) q[16];
rzz(0.7853830994606742) q[15],q[16];
u3(pi/2,0.5039114616358028,-0.5039114616358028) q[14];
rzz(0) q[14],q[23];
u3(pi,2.6188316360324517,-2.6188316360324517) q[14];
u3(pi/2,3.0460882369206637,-3.0460882369206637) q[22];
rzz(0) q[14],q[22];
u3(pi/2,4.7337518104291,-4.7337518104291) q[14];
rzz(0.012275082143518329) q[14],q[21];
rzz(pi/128) q[14],q[20];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/32) q[14],q[18];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/8) q[14],q[16];
u3(pi/2,1.2842830767875073,-1.2842830767875073) q[15];
rzz(0.7853830994606742) q[14],q[15];
u3(pi/2,0.5233893360880595,-0.5233893360880595) q[13];
rzz(0) q[13],q[23];
u3(pi,3.048601511043535,-3.048601511043535) q[13];
rzz(0) q[13],q[22];
u3(pi,4.9580615258954115,-4.9580615258954115) q[13];
u3(pi/2,2.899061700732661,-2.899061700732661) q[21];
rzz(0) q[13],q[21];
u3(pi/2,1.2007167122020188,-1.2007167122020188) q[13];
rzz(0.012275082143518329) q[13],q[20];
rzz(pi/128) q[13],q[19];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/32) q[13],q[17];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/8) q[13],q[15];
u3(pi/2,1.5921591568393072,-1.5921591568393072) q[14];
rzz(0.7853830994606742) q[13],q[14];
u3(pi/2,2.015645846543211,-2.015645846543211) q[12];
rzz(0) q[12],q[23];
u3(pi,5.034088068112284,-5.034088068112284) q[12];
rzz(0) q[12],q[22];
u3(pi,1.6455662319503337,-1.6455662319503337) q[12];
rzz(0) q[12],q[21];
u3(pi,4.540229702967969,-4.540229702967969) q[12];
u3(pi/2,5.359557067024187,-5.359557067024187) q[20];
rzz(0) q[12],q[20];
u3(pi/2,1.2754866173574562,-1.2754866173574562) q[12];
u3(pi,5.529831388848754,-5.529831388848754) q[19];
rzz(0.012275082143518329) q[12],q[19];
u3(pi,0.44799111240190453,-0.44799111240190453) q[18];
rzz(pi/128) q[12],q[18];
u3(pi,2.6081502210102463,-2.6081502210102463) q[17];
rzz(0.049100328574073315) q[12],q[17];
u3(pi,2.728159060377376,-2.728159060377376) q[16];
rzz(pi/32) q[12],q[16];
u3(pi,3.9477253285009337,-3.9477253285009337) q[15];
rzz(0.19640131429629326) q[12],q[15];
u3(pi,4.1179996503255,-4.1179996503255) q[14];
rzz(pi/8) q[12],q[14];
u3(pi/2,1.2007167122020188,-1.2007167122020188) q[13];
rzz(0.7853830994606742) q[12],q[13];
u3(pi/2,2.092300707290802,-2.092300707290802) q[11];
rzz(0) q[11],q[23];
u3(pi,4.3328845878310425,-4.3328845878310425) q[11];
rzz(0) q[11],q[22];
u3(pi,5.672459695321731,-5.672459695321731) q[11];
rzz(0) q[11],q[21];
u3(pi,0.7288494956328321,-0.7288494956328321) q[11];
rzz(0) q[11],q[20];
u3(pi,2.06842460312352,-2.06842460312352) q[11];
u3(pi/2,4.115486376202629,-4.115486376202629) q[19];
rzz(0) q[11],q[19];
u3(pi/2,4.30900848366376,-4.30900848366376) q[11];
rzz(0.012275082143518329) q[11],q[18];
rzz(pi/128) q[11],q[17];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/32) q[11],q[15];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/8) q[11],q[13];
u3(pi/2,1.2754866173574562,-1.2754866173574562) q[12];
rzz(0.7853830994606742) q[11],q[12];
u3(pi/2,1.7793980789932589,-1.7793980789932589) q[10];
rzz(0) q[10],q[23];
u3(pi,4.797840300562331,-4.797840300562331) q[10];
rzz(0) q[10],q[22];
u3(pi,1.409318464400381,-1.409318464400381) q[10];
rzz(0) q[10],q[21];
u3(pi,4.303981935418017,-4.303981935418017) q[10];
rzz(0) q[10],q[20];
u3(pi,0.9154600992560656,-0.9154600992560656) q[10];
rzz(0) q[10],q[19];
u3(pi,3.8101235702737015,-3.8101235702737015) q[10];
u3(pi/2,3.111433364115331,-3.111433364115331) q[18];
rzz(0) q[10],q[18];
u3(pi/2,0.5453804846631881,-0.5453804846631881) q[10];
rzz(0.012275082143518329) q[10],q[17];
rzz(pi/128) q[10],q[16];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/32) q[10],q[14];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/8) q[10],q[12];
u3(pi/2,4.30900848366376,-4.30900848366376) q[11];
rzz(0.7853830994606742) q[10],q[11];
u3(pi/2,5.227610175573416,-5.227610175573416) q[9];
rzz(0) q[9],q[23];
u3(pi,1.9628670899629028,-1.9628670899629028) q[9];
rzz(0) q[9],q[22];
u3(pi,4.857530560980538,-4.857530560980538) q[9];
rzz(0) q[9],q[21];
u3(pi,1.4690087248185872,-1.4690087248185872) q[9];
rzz(0) q[9],q[20];
u3(pi,4.36430051436694,-4.36430051436694) q[9];
rzz(0) q[9],q[19];
u3(pi,0.9757786782049896,-0.9757786782049896) q[9];
rzz(0) q[9],q[18];
u3(pi,3.870442149222625,-3.870442149222625) q[9];
u3(pi/2,1.6493361431346414,-1.6493361431346414) q[17];
rzz(0) q[9],q[17];
u3(pi/2,0.6056990636121121,-0.6056990636121121) q[9];
u3(pi,4.369327062612684,-4.369327062612684) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi,0.9091769139488861,-0.9091769139488861) q[15];
rzz(pi/128) q[9],q[15];
u3(pi,1.33894678895997,-1.33894678895997) q[14];
rzz(0.049100328574073315) q[9],q[14];
u3(pi,0.24378758991856794,-0.24378758991856794) q[13];
rzz(pi/32) q[9],q[13];
u3(pi,4.843707553304744,-4.843707553304744) q[12];
rzz(0.19640131429629326) q[9],q[12];
u3(pi,2.919167893715636,-2.919167893715636) q[11];
rzz(pi/8) q[9],q[11];
u3(pi/2,3.686973138252981,-3.686973138252981) q[10];
rzz(0.7853830994606742) q[9],q[10];
u3(pi/2,0.8344070087934491,-0.8344070087934491) q[8];
rzz(0) q[8],q[23];
u3(pi,4.403256263271454,-4.403256263271454) q[8];
rzz(0) q[8],q[22];
u3(pi,2.1155484929273665,-2.1155484929273665) q[8];
rzz(0) q[8],q[21];
u3(pi,6.110397711232148,-6.110397711232148) q[8];
rzz(0) q[8],q[20];
u3(pi,3.822689940888061,-3.822689940888061) q[8];
rzz(0) q[8],q[19];
u3(pi,1.5349821705439728,-1.5349821705439728) q[8];
rzz(0) q[8],q[18];
u3(pi,5.530459707379472,-5.530459707379472) q[8];
rzz(0) q[8],q[17];
u3(pi,3.2427519370353846,-3.2427519370353846) q[8];
u3(pi/2,0.2664070570244145,-0.2664070570244145) q[16];
rzz(0) q[8],q[16];
u3(pi/2,0.5284158843338032,-0.5284158843338032) q[8];
u3(pi,4.659610223804381,-4.659610223804381) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi,5.029689838397259,-5.029689838397259) q[14];
rzz(pi/128) q[8],q[14];
u3(pi,1.1158937105550946,-1.1158937105550946) q[13];
rzz(0.049100328574073315) q[8],q[13];
u3(pi,2.7111944600479916,-2.7111944600479916) q[12];
rzz(pi/32) q[8],q[12];
u3(pi,3.206309462253743,-3.206309462253743) q[11];
rzz(0.19640131429629326) q[8],q[11];
u3(pi,5.927557018793221,-5.927557018793221) q[10];
rzz(pi/8) q[8],q[10];
u3(pi/2,0.6056990636121121,-0.6056990636121121) q[9];
rzz(0.7853830994606742) q[8],q[9];
u3(pi/2,2.061513099285622,-2.061513099285622) q[7];
rzz(0) q[7],q[23];
u3(pi,5.630362353763627,-5.630362353763627) q[7];
rzz(0) q[7],q[22];
u3(pi,3.3426545834195402,-3.3426545834195402) q[7];
rzz(0) q[7],q[21];
u3(pi,1.0549468130754525,-1.0549468130754525) q[7];
rzz(0) q[7],q[20];
u3(pi,5.050424349910951,-5.050424349910951) q[7];
rzz(0) q[7],q[19];
u3(pi,2.762716579566864,-2.762716579566864) q[7];
rzz(0) q[7],q[18];
u3(pi,0.47438049069205873,-0.47438049069205873) q[7];
rzz(0) q[7],q[17];
u3(pi,4.4698580275275575,-4.4698580275275575) q[7];
rzz(0) q[7],q[16];
u3(pi,2.18215025718347,-2.18215025718347) q[7];
u3(pi/2,1.5444069485047422,-1.5444069485047422) q[15];
rzz(0) q[7],q[15];
u3(pi/2,5.750999511661475,-5.750999511661475) q[7];
u3(pi,0.9808052264507333,-0.9808052264507333) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi,1.9873715126609033,-1.9873715126609033) q[13];
rzz(pi/128) q[7],q[13];
u3(pi,0.578053048260522,-0.578053048260522) q[12];
rzz(0.049100328574073315) q[7],q[12];
u3(pi,3.508530675529081,-3.508530675529081) q[11];
rzz(pi/32) q[7],q[11];
u3(pi,1.2692034320502765,-1.2692034320502765) q[10];
rzz(0.19640131429629326) q[7],q[10];
u3(pi,5.635388902009371,-5.635388902009371) q[9];
rzz(pi/8) q[7],q[9];
u3(pi/2,3.670008537923596,-3.670008537923596) q[8];
rzz(0.7853830994606742) q[7],q[8];
u3(pi/2,3.3376280351737964,-3.3376280351737964) q[6];
rzz(0) q[6],q[23];
u3(pi,0.623291982472215,-0.623291982472215) q[6];
rzz(0) q[6],q[22];
u3(pi,4.618769519307714,-4.618769519307714) q[6];
rzz(0) q[6],q[21];
u3(pi,2.3310617489636263,-2.3310617489636263) q[6];
rzz(0) q[6],q[20];
u3(pi,0.043353978619539144,-0.043353978619539144) q[6];
rzz(0) q[6],q[19];
u3(pi,4.038831515455039,-4.038831515455039) q[6];
rzz(0) q[6],q[18];
u3(pi,1.7511237451109507,-1.7511237451109507) q[6];
rzz(0) q[6],q[17];
u3(pi,5.746601281946449,-5.746601281946449) q[6];
rzz(0) q[6],q[16];
u3(pi,3.458893511602362,-3.458893511602362) q[6];
rzz(0) q[6],q[15];
u3(pi,1.171185741258275,-1.171185741258275) q[6];
u3(pi/2,3.6436191596334417,-3.6436191596334417) q[14];
rzz(0) q[6],q[14];
u3(pi/2,4.739406677205562,-4.739406677205562) q[6];
u3(pi,3.212592647560922,-3.212592647560922) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi,3.717760746258161,-3.717760746258161) q[12];
rzz(pi/128) q[6],q[12];
u3(pi,3.943327098785909,-3.943327098785909) q[11];
rzz(0.049100328574073315) q[6],q[11];
u3(pi,4.376238566450582,-4.376238566450582) q[10];
rzz(pi/32) q[6],q[10];
u3(pi,0.6132388859807276,-0.6132388859807276) q[9];
rzz(0.19640131429629326) q[6],q[9];
u3(pi,0.9550441666912971,-0.9550441666912971) q[8];
rzz(pi/8) q[6],q[8];
u3(pi/2,5.750999511661475,-5.750999511661475) q[7];
rzz(0.7853830994606742) q[6],q[7];
u3(pi/2,5*pi/8,-5*pi/8) q[5];
rzz(0) q[5],q[23];
u3(pi,4.302096979825863,-4.302096979825863) q[5];
rzz(0) q[5],q[22];
u3(pi,5.838335787431272,-5.838335787431272) q[5];
rzz(0) q[5],q[21];
u3(pi,1.091389287857094,-1.091389287857094) q[5];
rzz(0) q[5],q[20];
u3(pi,2.627628095462503,-2.627628095462503) q[5];
rzz(0) q[5],q[19];
u3(pi,4.163866903067912,-4.163866903067912) q[5];
rzz(0) q[5],q[18];
u3(pi,5.700105710673321,-5.700105710673321) q[5];
rzz(0) q[5],q[17];
u3(pi,0.9531592110991433,-0.9531592110991433) q[5];
rzz(0) q[5],q[16];
u3(pi,2.489398018704552,-2.489398018704552) q[5];
rzz(0) q[5],q[15];
u3(pi,4.025636826309961,-4.025636826309961) q[5];
rzz(0) q[5],q[14];
u3(pi,5.56187563391537,-5.56187563391537) q[5];
u3(pi/2,2.2537785696853176,-2.2537785696853176) q[13];
rzz(0) q[5],q[13];
u3(pi/2,1.6179202165987434,-1.6179202165987434) q[5];
rzz(0.012275082143518329) q[5],q[12];
rzz(pi/128) q[5],q[11];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/32) q[5],q[9];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/8) q[5],q[7];
u3(pi/2,4.739406677205562,-4.739406677205562) q[6];
rzz(0.7853830994606742) q[5],q[6];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(0) q[4],q[23];
u3(pi,4.353619099344735,-4.353619099344735) q[4];
rzz(0) q[4],q[22];
u3(pi,2.065911329000648,-2.065911329000648) q[4];
rzz(0) q[4],q[21];
u3(pi,6.061388865836147,-6.061388865836147) q[4];
rzz(0) q[4],q[20];
u3(pi,3.7736810954920594,-3.7736810954920594) q[4];
rzz(0) q[4],q[19];
u3(pi,1.485973325147972,-1.485973325147972) q[4];
rzz(0) q[4],q[18];
u3(pi,5.481450861983471,-5.481450861983471) q[4];
rzz(0) q[4],q[17];
u3(pi,3.1937430916393836,-3.1937430916393836) q[4];
rzz(0) q[4],q[16];
u3(pi,0.9060353212952963,-0.9060353212952963) q[4];
rzz(0) q[4],q[15];
u3(pi,4.901512858130795,-4.901512858130795) q[4];
rzz(0) q[4],q[14];
u3(pi,2.613805087786708,-2.613805087786708) q[4];
rzz(0) q[4],q[13];
u3(pi,0.32609731744262055,-0.32609731744262055) q[4];
u3(pi/2,0.14702653618800232,-0.14702653618800232) q[12];
rzz(0) q[4],q[12];
u3(pi/2,3.8943182533899074,-3.8943182533899074) q[4];
u3(pi,3.038548414552048,-3.038548414552048) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi,3.4720882007474394,-3.4720882007474394) q[10];
rzz(pi/128) q[4],q[10];
u3(pi,1.0210176124166828,-1.0210176124166828) q[9];
rzz(0.049100328574073315) q[4],q[9];
u3(pi,0.06911503837897544,-0.06911503837897544) q[8];
rzz(pi/32) q[4],q[8];
u3(pi,2.4856281075202444,-2.4856281075202444) q[7];
rzz(0.19640131429629326) q[4],q[7];
u3(pi,3.485911208423234,-3.485911208423234) q[6];
rzz(pi/8) q[4],q[6];
u3(pi/2,4.759512870188536,-4.759512870188536) q[5];
rzz(0.7853830994606742) q[4],q[5];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(0) q[3],q[23];
u3(pi,5.1390172627421835,-5.1390172627421835) q[3];
rzz(0) q[3],q[22];
u3(pi,2.851309492398096,-2.851309492398096) q[3];
rzz(0) q[3],q[21];
u3(pi,0.5636017220540089,-0.5636017220540089) q[3];
rzz(0) q[3],q[20];
u3(pi,4.559079258889508,-4.559079258889508) q[3];
rzz(0) q[3],q[19];
u3(pi,2.2713714885454204,-2.2713714885454204) q[3];
rzz(0) q[3],q[18];
u3(pi,6.266849025380919,-6.266849025380919) q[3];
rzz(0) q[3],q[17];
u3(pi,3.979141255036832,-3.979141255036832) q[3];
rzz(0) q[3],q[16];
u3(pi,1.6914334846927446,-1.6914334846927446) q[3];
rzz(0) q[3],q[15];
u3(pi,5.6869110215282435,-5.6869110215282435) q[3];
rzz(0) q[3],q[14];
u3(pi,3.3992032511841566,-3.3992032511841566) q[3];
rzz(0) q[3],q[13];
u3(pi,1.1114954808400688,-1.1114954808400688) q[3];
rzz(0) q[3],q[12];
u3(pi,5.106973017675568,-5.106973017675568) q[3];
u3(pi/2,0.32358404331974866,-0.32358404331974866) q[11];
rzz(0) q[3],q[11];
u3(pi/2,2.3920086464432684,-2.3920086464432684) q[3];
rzz(0.012275082143518329) q[3],q[10];
rzz(pi/128) q[3],q[9];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/32) q[3],q[7];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/8) q[3],q[5];
u3(pi/2,0.7527255998001144,-0.7527255998001144) q[4];
rzz(0.7853830994606742) q[3],q[4];
u3(pi/2,pi,-pi) q[2];
rzz(0) q[2],q[23];
u3(pi,5.74345968929286,-5.74345968929286) q[2];
rzz(0) q[2],q[22];
u3(pi,1.5230441184603318,-1.5230441184603318) q[2];
rzz(0) q[2],q[21];
u3(pi,3.5858138548073897,-3.5858138548073897) q[2];
rzz(0) q[2],q[20];
u3(pi,5.648583591154448,-5.648583591154448) q[2];
rzz(0) q[2],q[19];
u3(pi,1.4275397017912022,-1.4275397017912022) q[2];
rzz(0) q[2],q[18];
u3(pi,3.4903094381382602,-3.4903094381382602) q[2];
rzz(0) q[2],q[17];
u3(pi,5.553079174485318,-5.553079174485318) q[2];
rzz(0) q[2],q[16];
u3(pi,1.3326636036527904,-1.3326636036527904) q[2];
rzz(0) q[2],q[15];
u3(pi,3.3948050214691303,-3.3948050214691303) q[2];
rzz(0) q[2],q[14];
u3(pi,5.457574757816189,-5.457574757816189) q[2];
rzz(0) q[2],q[13];
u3(pi,1.2371591869836605,-1.2371591869836605) q[2];
rzz(0) q[2],q[12];
u3(pi,3.2999289233307185,-3.2999289233307185) q[2];
rzz(0) q[2],q[11];
u3(pi,5.362070341147059,-5.362070341147059) q[2];
u3(pi/2,1.986114875599467,-1.986114875599467) q[10];
rzz(0) q[2],q[10];
u3(pi/2,1.6813803882012572,-1.6813803882012572) q[2];
u3(pi,1.0989291102257097,-1.0989291102257097) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi,5.637273857601524,-5.637273857601524) q[8];
rzz(pi/128) q[2],q[8];
u3(pi,0.9481326628533996,-0.9481326628533996) q[7];
rzz(0.049100328574073315) q[2],q[7];
u3(pi,1.2503538761287378,-1.2503538761287378) q[6];
rzz(pi/32) q[2],q[6];
u3(pi,3.800698792312932,-3.800698792312932) q[5];
rzz(0.19640131429629326) q[2],q[5];
u3(pi,3.7711678213691875,-3.7711678213691875) q[4];
rzz(pi/8) q[2],q[4];
u3(pi/2,2.3920086464432684,-2.3920086464432684) q[3];
rzz(0.7853830994606742) q[2],q[3];
u3(pi/2,0,0) q[1];
rzz(0) q[1],q[23];
u3(pi,3.568220935947287,-3.568220935947287) q[1];
rzz(0) q[1],q[22];
u3(pi,1.2805131656031998,-1.2805131656031998) q[1];
rzz(0) q[1],q[21];
u3(pi,5.2759907024386985,-5.2759907024386985) q[1];
rzz(0) q[1],q[20];
u3(pi,2.988282932094611,-2.988282932094611) q[1];
rzz(0) q[1],q[19];
u3(pi,0.7005751617505239,-0.7005751617505239) q[1];
rzz(0) q[1],q[18];
u3(pi,4.696052698586023,-4.696052698586023) q[1];
rzz(0) q[1],q[17];
u3(pi,2.4083449282419354,-2.4083449282419354) q[1];
rzz(0) q[1],q[16];
u3(pi,0.12063715789784804,-0.12063715789784804) q[1];
rzz(0) q[1],q[15];
u3(pi,4.116114694733347,-4.116114694733347) q[1];
rzz(0) q[1],q[14];
u3(pi,1.8284069243892596,-1.8284069243892596) q[1];
rzz(0) q[1],q[13];
u3(pi,5.8238844612247584,-5.8238844612247584) q[1];
rzz(0) q[1],q[12];
u3(pi,3.536176690880671,-3.536176690880671) q[1];
rzz(0) q[1],q[11];
u3(pi,1.2484689205365838,-1.2484689205365838) q[1];
rzz(0) q[1],q[10];
u3(pi,5.243946457372083,-5.243946457372083) q[1];
u3(pi/2,0.14199998794225863,-0.14199998794225863) q[9];
rzz(0) q[1],q[9];
u3(pi/2,2.528982086139784,-2.528982086139784) q[1];
rzz(0.012275082143518329) q[1],q[8];
rzz(pi/128) q[1],q[7];
rzz(0.049100328574073315) q[1],q[6];
rzz(pi/32) q[1],q[5];
rzz(0.19640131429629326) q[1],q[4];
rzz(pi/8) q[1],q[3];
u3(pi/2,1.6813803882012572,-1.6813803882012572) q[2];
rzz(0.7853830994606742) q[1],q[2];
rzz(pi/2) q[0],q[23];
rzz(pi/2) q[0],q[22];
rzz(pi/2) q[0],q[21];
rzz(pi/2) q[0],q[20];
rzz(-pi/2) q[0],q[19];
rzz(pi/2) q[0],q[18];
rzz(pi/2) q[0],q[17];
rzz(pi/2) q[0],q[16];
rzz(-pi/2) q[0],q[15];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[10];
rzz(pi/2) q[0],q[9];
u3(pi/2,3.0938404452552284,-3.0938404452552284) q[8];
rzz(pi/2) q[0],q[8];
u3(pi/2,5.817601275917578,-5.817601275917578) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,3.409884666206361,-3.409884666206361) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,2.841256395906609,-2.841256395906609) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,0.5064247357586746,-0.5064247357586746) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,5.5336013000330615,-5.5336013000330615) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,4.82297304179105,-4.82297304179105) q[2];
rzz(pi/2) q[0],q[2];
rzz(-pi/2) q[0],q[1];
u3(pi,pi/2,-pi/2) q[0];
u3(pi/2,5.6705747397295765,-5.6705747397295765) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,0.11058406140636072,-0.11058406140636072) q[2];
u3(pi/2,2.4667785515987055,-2.4667785515987055) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,0.8212123196483719,-0.8212123196483719) q[3];
u3(pi/2,3.5701058915394412,-3.5701058915394412) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[4];
u3(pi/2,5.022150016028643,-5.022150016028643) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,4.412052722701506,-4.412052722701506) q[5];
u3(pi/2,1.1724423783197107,-1.1724423783197107) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,1.839088339411465,-1.839088339411465) q[6];
u3(pi/2,4.931672147605258,-4.931672147605258) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,1.1052122955328891,-1.1052122955328891) q[7];
u3(pi/2,4.221672207893964,-4.221672207893964) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,4.664636772050124,-4.664636772050124) q[8];
u3(pi/2,1.5104777478459726,-1.5104777478459726) q[8];
rzz(-pi/2) q[0],q[8];
rzz(pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[12];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[15];
rzz(pi/2) q[0],q[16];
rzz(pi/2) q[0],q[17];
rzz(pi/2) q[0],q[18];
rzz(pi/2) q[0],q[19];
rzz(pi/2) q[0],q[20];
rzz(-pi/2) q[0],q[21];
rzz(pi/2) q[0],q[22];
rzz(pi/2) q[0],q[23];
u3(pi/2,4.09977841293468,-4.09977841293468) q[1];
u3(pi/2,0.8959822248038091,-0.8959822248038091) q[2];
rzz(0.7853830994606742) q[1],q[2];
u3(pi/2,1.9993095647445442,-1.9993095647445442) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,0.3097610356439536,-0.3097610356439536) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi/2,5.8848313587044006,-5.8848313587044006) q[5];
rzz(pi/32) q[1],q[5];
u3(pi/2,3.360875820810361,-3.360875820810361) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi/2,2.6508758810990676,-2.6508758810990676) q[7];
rzz(pi/128) q[1],q[7];
u3(pi/2,3.081274074640869,-3.081274074640869) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi/2,4.09977841293468,-4.09977841293468) q[1];
rzz(0) q[1],q[9];
u3(pi,0.4190884599888784,-0.4190884599888784) q[1];
rzz(0) q[1],q[10];
u3(pi,2.4812298778052186,-2.4812298778052186) q[1];
rzz(0) q[1],q[11];
u3(pi,4.543999614152276,-4.543999614152276) q[1];
rzz(0) q[1],q[12];
u3(pi,0.32358404331974866,-0.32358404331974866) q[1];
rzz(0) q[1],q[13];
u3(pi,2.386353779666807,-2.386353779666807) q[1];
rzz(0) q[1],q[14];
u3(pi,4.448495197483147,-4.448495197483147) q[1];
rzz(0) q[1],q[15];
u3(pi,0.22807962665061898,-0.22807962665061898) q[1];
rzz(0) q[1],q[16];
u3(pi,2.290849362997677,-2.290849362997677) q[1];
rzz(0) q[1],q[17];
u3(pi,4.353619099344735,-4.353619099344735) q[1];
rzz(0) q[1],q[18];
u3(pi,0.13257520998148928,-0.13257520998148928) q[1];
rzz(0) q[1],q[19];
u3(pi,2.1953449463285475,-2.1953449463285475) q[1];
rzz(0) q[1],q[20];
u3(pi,4.258114682675606,-4.258114682675606) q[1];
rzz(0) q[1],q[21];
u3(pi,0.03769911184307752,-0.03769911184307752) q[1];
rzz(0) q[1],q[22];
u3(pi,2.0998405296594176,-2.0998405296594176) q[1];
rzz(0) q[1],q[23];
u3(pi/2,0.11058406140636072,-0.11058406140636072) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,3.2710262709176923,-3.2710262709176923) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi/2,4.82297304179105,-4.82297304179105) q[2];
rzz(0) q[2],q[10];
u3(pi,2.1080086705587515,-2.1080086705587515) q[2];
rzz(0) q[2],q[11];
u3(pi,6.10348620739425,-6.10348620739425) q[2];
rzz(0) q[2],q[12];
u3(pi,3.8157784370501626,-3.8157784370501626) q[2];
rzz(0) q[2],q[13];
u3(pi,1.5280706667060753,-1.5280706667060753) q[2];
rzz(0) q[2],q[14];
u3(pi,5.523548203541575,-5.523548203541575) q[2];
rzz(0) q[2],q[15];
u3(pi,3.235840433197487,-3.235840433197487) q[2];
rzz(0) q[2],q[16];
u3(pi,0.9481326628533996,-0.9481326628533996) q[2];
rzz(0) q[2],q[17];
u3(pi,4.943610199688899,-4.943610199688899) q[2];
rzz(0) q[2],q[18];
u3(pi,2.655902429344811,-2.655902429344811) q[2];
rzz(0) q[2],q[19];
u3(pi,0.36819465900072373,-0.36819465900072373) q[2];
rzz(0) q[2],q[20];
u3(pi,4.363672195836223,-4.363672195836223) q[2];
rzz(0) q[2],q[21];
u3(pi,2.0759644254921357,-2.0759644254921357) q[2];
rzz(0) q[2],q[22];
u3(pi,6.071441962327635,-6.071441962327635) q[2];
rzz(0) q[2],q[23];
u3(pi/2,4.748203136635613,-4.748203136635613) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,5.118282751228491,-5.118282751228491) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,0.03581415625092364,-0.03581415625092364) q[3];
rzz(0) q[3],q[11];
u3(pi,3.6046634107289286,-3.6046634107289286) q[3];
rzz(0) q[3],q[12];
u3(pi,1.3169556403848413,-1.3169556403848413) q[3];
rzz(0) q[3],q[13];
u3(pi,5.31243317722034,-5.31243317722034) q[3];
rzz(0) q[3],q[14];
u3(pi,3.024725406876253,-3.024725406876253) q[3];
rzz(0) q[3],q[15];
u3(pi,0.7370176365321655,-0.7370176365321655) q[3];
rzz(0) q[3],q[16];
u3(pi,4.732495173367664,-4.732495173367664) q[3];
rzz(0) q[3],q[17];
u3(pi,2.444787403023577,-2.444787403023577) q[3];
rzz(0) q[3],q[18];
u3(pi,pi/20,-pi/20) q[3];
rzz(0) q[3],q[19];
u3(pi,4.151928850984271,-4.151928850984271) q[3];
rzz(0) q[3],q[20];
u3(pi,1.8642210806401833,-1.8642210806401833) q[3];
rzz(0) q[3],q[21];
u3(pi,5.859698617475682,-5.859698617475682) q[3];
rzz(0) q[3],q[22];
u3(pi,3.5719908471315946,-3.5719908471315946) q[3];
rzz(0) q[3],q[23];
u3(pi/2,3.255318307649744,-3.255318307649744) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,0.3066194429903638,-0.3066194429903638) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,1.684521980854847,-1.684521980854847) q[4];
rzz(0) q[4],q[12];
u3(pi,4.852504012734794,-4.852504012734794) q[4];
rzz(0) q[4],q[13];
u3(pi,1.7643184342560279,-1.7643184342560279) q[4];
rzz(0) q[4],q[14];
u3(pi,4.959318162956848,-4.959318162956848) q[4];
rzz(0) q[4],q[15];
u3(pi,1.8711325844780808,-1.8711325844780808) q[4];
rzz(0) q[4],q[16];
u3(pi,5.0661323131789,-5.0661323131789) q[4];
rzz(0) q[4],q[17];
u3(pi,1.9779467347001338,-1.9779467347001338) q[4];
rzz(0) q[4],q[18];
u3(pi,5.172946463400954,-5.172946463400954) q[4];
rzz(0) q[4],q[19];
u3(pi,2.084760884922187,-2.084760884922187) q[4];
rzz(0) q[4],q[20];
u3(pi,5.279760613623006,-5.279760613623006) q[4];
rzz(0) q[4],q[21];
u3(pi,2.19157503514424,-2.19157503514424) q[4];
rzz(0) q[4],q[22];
u3(pi,5.386574763845059,-5.386574763845059) q[4];
rzz(0) q[4],q[23];
u3(pi/2,5.7868136679124,-5.7868136679124) q[5];
rzz(0.7853830994606742) q[5],q[6];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,0.13069025438933538,-0.13069025438933538) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,1.0744246875277093,-1.0744246875277093) q[5];
rzz(0) q[5],q[13];
u3(pi,3.737238620710418,-3.737238620710418) q[5];
rzz(0) q[5],q[14];
u3(pi,5.92190215201676,-5.92190215201676) q[5];
rzz(0) q[5],q[15];
u3(pi,1.8233803761435161,-1.8233803761435161) q[5];
rzz(0) q[5],q[16];
u3(pi,4.008043907449858,-4.008043907449858) q[5];
rzz(0) q[5],q[17];
u3(pi,6.192707438756201,-6.192707438756201) q[5];
rzz(0) q[5],q[18];
u3(pi,2.094185662882956,-2.094185662882956) q[5];
rzz(0) q[5],q[19];
u3(pi,4.2788491941892985,-4.2788491941892985) q[5];
rzz(0) q[5],q[20];
u3(pi,0.18032741831605412,-0.18032741831605412) q[5];
rzz(0) q[5],q[21];
u3(pi,2.3649909496223964,-2.3649909496223964) q[5];
rzz(0) q[5],q[22];
u3(pi,4.54902616239802,-4.54902616239802) q[5];
rzz(0) q[5],q[23];
u3(pi/2,3.31186697541436,-3.31186697541436) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,5.37714998588429,-5.37714998588429) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,1.7410706486194634,-1.7410706486194634) q[6];
rzz(0) q[6],q[14];
u3(pi,5.881061447520093,-5.881061447520093) q[6];
rzz(0) q[6],q[15];
u3(pi,4.735008447490537,-4.735008447490537) q[6];
rzz(0) q[6],q[16];
u3(pi,3.589583765991698,-3.589583765991698) q[6];
rzz(0) q[6],q[17];
u3(pi,2.444159084492859,-2.444159084492859) q[6];
rzz(0) q[6],q[18];
u3(pi,1.2987344029940204,-1.2987344029940204) q[6];
rzz(0) q[6],q[19];
u3(pi,0.1533097214951819,-0.1533097214951819) q[6];
rzz(0) q[6],q[20];
u3(pi,5.290442028645211,-5.290442028645211) q[6];
rzz(0) q[6],q[21];
u3(pi,4.145017347146373,-4.145017347146373) q[6];
rzz(0) q[6],q[22];
u3(pi,2.9995926656475342,-2.9995926656475342) q[6];
rzz(0) q[6],q[23];
u3(pi/2,5.76796411199086,-5.76796411199086) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,3.626026240773339,-3.626026240773339) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,1.0555751316061706,-1.0555751316061706) q[7];
rzz(0) q[7],q[15];
u3(pi,4.624424386084176,-4.624424386084176) q[7];
rzz(0) q[7],q[16];
u3(pi,2.3367166157400883,-2.3367166157400883) q[7];
rzz(0) q[7],q[17];
u3(pi,0.04900884539600077,-0.04900884539600077) q[7];
rzz(0) q[7],q[18];
u3(pi,4.0444863822315,-4.0444863822315) q[7];
rzz(0) q[7],q[19];
u3(pi,1.7567786118874125,-1.7567786118874125) q[7];
rzz(0) q[7],q[20];
u3(pi,5.752256148722911,-5.752256148722911) q[7];
rzz(0) q[7],q[21];
u3(pi,3.464548378378824,-3.464548378378824) q[7];
rzz(0) q[7],q[22];
u3(pi,1.1768406080347364,-1.1768406080347364) q[7];
rzz(0) q[7],q[23];
u3(pi/2,3.06870770402651,-3.06870770402651) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,1.5261857111139214,-1.5261857111139214) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,1.4979113772316133,-1.4979113772316133) q[8];
rzz(0) q[8],q[16];
u3(pi,5.066760631709618,-5.066760631709618) q[8];
rzz(0) q[8],q[17];
u3(pi,2.779052861365531,-2.779052861365531) q[8];
rzz(0) q[8],q[18];
u3(pi,0.4913450910214437,-0.4913450910214437) q[8];
rzz(0) q[8],q[19];
u3(pi,4.4868226278569425,-4.4868226278569425) q[8];
rzz(0) q[8],q[20];
u3(pi,7*pi/10,-7*pi/10) q[8];
rzz(0) q[8],q[21];
u3(pi,6.194592394348354,-6.194592394348354) q[8];
rzz(0) q[8],q[22];
u3(pi,3.9068846240042667,-3.9068846240042667) q[8];
rzz(0) q[8],q[23];
u3(pi/2,0.11749556524425828,-0.11749556524425828) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,0.24818581963359365,-0.24818581963359365) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,1.6882918920391548,-1.6882918920391548) q[9];
rzz(0) q[9],q[17];
u3(pi,5.256512827986442,-5.256512827986442) q[9];
rzz(0) q[9],q[18];
u3(pi,2.968805057642354,-2.968805057642354) q[9];
rzz(0) q[9],q[19];
u3(pi,0.6810972872982671,-0.6810972872982671) q[9];
rzz(0) q[9],q[20];
u3(pi,4.676574824133766,-4.676574824133766) q[9];
rzz(0) q[9],q[21];
u3(pi,2.388867053789679,-2.388867053789679) q[9];
rzz(0) q[9],q[22];
u3(pi,0.10115928344559133,-0.10115928344559133) q[9];
rzz(0) q[9],q[23];
u3(pi/2,5.108857973267722,-5.108857973267722) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,4.772079240802896,-4.772079240802896) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi/2,3.5380616464728254,-3.5380616464728254) q[10];
rzz(0) q[10],q[18];
u3(pi,1.3948671381938682,-1.3948671381938682) q[10];
rzz(0) q[10],q[19];
u3(pi,0.24944245669502957,-0.24944245669502957) q[10];
rzz(0) q[10],q[20];
u3(pi,5.386574763845059,-5.386574763845059) q[10];
rzz(0) q[10],q[21];
u3(pi,4.241150082346221,-4.241150082346221) q[10];
rzz(0) q[10],q[22];
u3(pi,3.0957254008473822,-3.0957254008473822) q[10];
rzz(0) q[10],q[23];
u3(pi/2,0.2990796206217483,-0.2990796206217483) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,6.234804780314303,-6.234804780314303) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,1.8711325844780808,-1.8711325844780808) q[11];
rzz(0) q[11],q[19];
u3(pi,5.998557012764351,-5.998557012764351) q[11];
rzz(0) q[11],q[20];
u3(pi,4.827371271506076,-4.827371271506076) q[11];
rzz(0) q[11],q[21];
u3(pi,3.656813848778519,-3.656813848778519) q[11];
rzz(0) q[11],q[22];
u3(pi,2.486256426050962,-2.486256426050962) q[11];
rzz(0) q[11],q[23];
u3(pi/2,0.12377875055143785,-0.12377875055143785) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/128) q[12],q[18];
u3(pi/2,4.097265138811808,-4.097265138811808) q[19];
rzz(0.012275082143518329) q[12],q[19];
u3(pi/2,4.836167730936128,-4.836167730936128) q[12];
rzz(0) q[12],q[20];
u3(pi,1.155477777990326,-1.155477777990326) q[12];
rzz(0) q[12],q[21];
u3(pi,3.218247514337384,-3.218247514337384) q[12];
rzz(0) q[12],q[22];
u3(pi,5.280388932153724,-5.280388932153724) q[12];
rzz(0) q[12],q[23];
u3(pi/2,2.229274146987317,-2.229274146987317) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/128) q[13],q[19];
u3(pi/2,5.341335829633366,-5.341335829633366) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi/2,3.8000704737822137,-3.8000704737822137) q[13];
rzz(0) q[13],q[21];
u3(pi,0.9927432785343746,-0.9927432785343746) q[13];
rzz(0) q[13],q[22];
u3(pi,4.803495167338793,-4.803495167338793) q[13];
rzz(0) q[13],q[23];
u3(pi/2,0.47815040187636654,-0.47815040187636654) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/128) q[14],q[20];
u3(pi/2,6.022433116931634,-6.022433116931634) q[21];
rzz(0.012275082143518329) q[14],q[21];
u3(pi/2,5.190539382261056,-5.190539382261056) q[14];
rzz(0) q[14],q[22];
u3(pi,1.5092211107845366,-1.5092211107845366) q[14];
rzz(0) q[14],q[23];
u3(pi/2,4.661495179396535,-4.661495179396535) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/32) q[15],q[19];
rzz(0.049100328574073315) q[15],q[20];
rzz(pi/128) q[15],q[21];
u3(pi/2,3.0278669995298424,-3.0278669995298424) q[22];
rzz(0.012275082143518329) q[15],q[22];
u3(pi/2,6.232291506191432,-6.232291506191432) q[15];
rzz(0) q[15],q[23];
u3(pi/2,3.3834952879162072,-3.3834952879162072) q[16];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/8) q[16],q[18];
rzz(0.19640131429629326) q[16],q[19];
rzz(pi/32) q[16],q[20];
rzz(0.049100328574073315) q[16],q[21];
rzz(pi/128) q[16],q[22];
u3(pi/2,1.7178228629828987,-1.7178228629828987) q[23];
rzz(0.012275082143518329) q[16],q[23];
u3(pi/2,4.7657960554957155,-4.7657960554957155) q[17];
rzz(0.7853830994606742) q[17],q[18];
rzz(pi/8) q[17],q[19];
rzz(0.19640131429629326) q[17],q[20];
rzz(pi/32) q[17],q[21];
rzz(0.049100328574073315) q[17],q[22];
rzz(pi/128) q[17],q[23];
u3(pi/2,6.228521595007123,-6.228521595007123) q[18];
rzz(0.7853830994606742) q[18],q[19];
rzz(pi/8) q[18],q[20];
rzz(0.19640131429629326) q[18],q[21];
rzz(pi/32) q[18],q[22];
rzz(0.049100328574073315) q[18],q[23];
u3(pi/2,4.090981953504628,-4.090981953504628) q[19];
rzz(0.7853830994606742) q[19],q[20];
rzz(pi/8) q[19],q[21];
rzz(0.19640131429629326) q[19],q[22];
rzz(pi/32) q[19],q[23];
u3(pi/2,5.335052644326186,-5.335052644326186) q[20];
rzz(0.7853830994606742) q[20],q[21];
rzz(pi/8) q[20],q[22];
rzz(0.19640131429629326) q[20],q[23];
u3(pi/2,2.874557278034661,-2.874557278034661) q[21];
rzz(0.7853830994606742) q[21],q[22];
rzz(pi/8) q[21],q[23];
u3(pi/2,6.163176467812456,-6.163176467812456) q[22];
rzz(0.7853830994606742) q[22],q[23];
u3(pi/2,4.702335883893202,-4.702335883893202) q[1];
u3(pi/2,3.356477591095335,-3.356477591095335) q[2];
u3(pi/2,0.8576547944300136,-0.8576547944300136) q[3];
u3(pi/2,2.271999807076138,-2.271999807076138) q[4];
u3(pi/2,0.9292831069318608,-0.9292831069318608) q[5];
u3(pi/2,0.8563981573685776,-0.8563981573685776) q[6];
u3(pi/2,4.745061543982024,-4.745061543982024) q[7];
u3(pi/2,1.1919202527719674,-1.1919202527719674) q[8];
u3(pi/2,3.670008537923596,-3.670008537923596) q[9];
u3(pi/2,0.9519025740377073,-0.9519025740377073) q[10];
u3(pi/2,0.32986722862692824,-0.32986722862692824) q[11];
u3(pi/2,1.5996989792079226,-1.5996989792079226) q[12];
u3(pi/2,1.9967962906216727,-1.9967962906216727) q[13];
u3(pi/2,4.111716465018321,-4.111716465018321) q[14];
u3(pi/2,3.0906988526016383,-3.0906988526016383) q[15];
u3(pi/2,1.7115396776757192,-1.7115396776757192) q[23];
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
