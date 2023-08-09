OPENQASM 2.0;
include "qelib1.inc";
qreg q[19];
creg c[19];
u3(pi/2,4.178946547805143,-4.178946547805143) q[18];
u3(pi/2,0.2519557308179014,-0.2519557308179014) q[18];
rzz(-pi/2) q[17],q[18];
u3(pi/2,4.964344711202592,-4.964344711202592) q[18];
rzz(pi/8) q[16],q[18];
u3(pi/2,1.0147344271095031,-1.0147344271095031) q[17];
rzz(0.7853830994606742) q[16],q[17];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/8) q[15],q[17];
u3(pi/2,4.150043895392116,-4.150043895392116) q[16];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/32) q[14],q[18];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/8) q[14],q[16];
u3(pi/2,0.9160884177867837,-0.9160884177867837) q[15];
rzz(0.7853830994606742) q[14],q[15];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/32) q[13],q[17];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/8) q[13],q[15];
u3(pi/2,4.03506160427073,-4.03506160427073) q[14];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/128) q[12],q[18];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/32) q[12],q[16];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/8) q[12],q[14];
u3(pi/2,0.5240176546187775,-0.5240176546187775) q[13];
rzz(0.7853830994606742) q[12],q[13];
rzz(0.012275082143518329) q[11],q[18];
rzz(pi/128) q[11],q[17];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/32) q[11],q[15];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/8) q[11],q[13];
u3(pi/2,3.5726191656623127,-3.5726191656623127) q[12];
rzz(0.7853830994606742) q[11],q[12];
u3(pi/2,1.7241060482900783,-1.7241060482900783) q[10];
u3(pi/2,1.822752057612798,-1.822752057612798) q[18];
rzz(0) q[10],q[18];
u3(pi/2,4.865698701879872,-4.865698701879872) q[10];
rzz(0.012275082143518329) q[10],q[17];
rzz(pi/128) q[10],q[16];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/32) q[10],q[14];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/8) q[10],q[12];
u3(pi/2,5.237034953534185,-5.237034953534185) q[11];
rzz(0.7853830994606742) q[10],q[11];
u3(pi/2,2.098583892597982,-2.098583892597982) q[9];
rzz(0) q[9],q[18];
u3(pi,4.437185463930224,-4.437185463930224) q[9];
u3(pi/2,4.156327080699296,-4.156327080699296) q[17];
rzz(0) q[9],q[17];
u3(pi/2,0.49323004661359754,-0.49323004661359754) q[9];
rzz(0.012275082143518329) q[9],q[16];
rzz(pi/128) q[9],q[15];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/32) q[9],q[13];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/8) q[9],q[11];
u3(pi/2,1.7241060482900783,-1.7241060482900783) q[10];
rzz(0.7853830994606742) q[9],q[10];
u3(pi/2,0.6138672045114455,-0.6138672045114455) q[8];
rzz(0) q[8],q[18];
u3(pi,3.7893890587600083,-3.7893890587600083) q[8];
rzz(0) q[8],q[17];
u3(pi,0.7169114435491908,-0.7169114435491908) q[8];
u3(pi/2,1.0084512418023237,-1.0084512418023237) q[16];
rzz(0) q[8],q[16];
u3(pi/2,3.892433297797754,-3.892433297797754) q[8];
u3(pi,5.867238439844297,-5.867238439844297) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi,1.227734409022891,-1.227734409022891) q[14];
rzz(pi/128) q[8],q[14];
u3(pi,5.445008387201829,-5.445008387201829) q[13];
rzz(0.049100328574073315) q[8],q[13];
u3(pi,5.732149955739937,-5.732149955739937) q[12];
rzz(pi/32) q[8],q[12];
u3(pi,3.9835394847518577,-3.9835394847518577) q[11];
rzz(0.19640131429629326) q[8],q[11];
u3(pi,3.9049996684121133,-3.9049996684121133) q[10];
rzz(pi/8) q[8],q[10];
u3(pi/2,0.49323004661359754,-0.49323004661359754) q[9];
rzz(0.7853830994606742) q[8],q[9];
u3(pi/2,2.110521944681623,-2.110521944681623) q[7];
rzz(0) q[7],q[18];
u3(pi,5.679371199159628,-5.679371199159628) q[7];
rzz(0) q[7],q[17];
u3(pi,3.3916634288155403,-3.3916634288155403) q[7];
rzz(0) q[7],q[16];
u3(pi,1.1039556584714532,-1.1039556584714532) q[7];
u3(pi/2,4.5358314732529434,-4.5358314732529434) q[15];
rzz(0) q[7],q[15];
u3(pi/2,4.67217659441874,-4.67217659441874) q[7];
u3(pi,3.901229757227805,-3.901229757227805) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi,5.861583573067835,-5.861583573067835) q[13];
rzz(pi/128) q[7],q[13];
u3(pi,2.2066546798814706,-2.2066546798814706) q[12];
rzz(0.049100328574073315) q[7],q[12];
u3(pi,4.407026174455762,-4.407026174455762) q[11];
rzz(pi/32) q[7],q[11];
u3(pi,5.10320310649126,-5.10320310649126) q[10];
rzz(0.19640131429629326) q[7],q[10];
u3(pi,5.311804858689622,-5.311804858689622) q[9];
rzz(pi/8) q[7],q[9];
u3(pi/2,0.7508406442079605,-0.7508406442079605) q[8];
rzz(0.7853830994606742) q[7],q[8];
u3(pi/2,2.4542121809843462,-2.4542121809843462) q[6];
rzz(0) q[6],q[18];
u3(pi,6.023061435462352,-6.023061435462352) q[6];
rzz(0) q[6],q[17];
u3(pi,3.7353536651182644,-3.7353536651182644) q[6];
rzz(0) q[6],q[16];
u3(pi,1.4476458947741766,-1.4476458947741766) q[6];
rzz(0) q[6],q[15];
u3(pi,5.443123431609675,-5.443123431609675) q[6];
u3(pi/2,6.239831328560047,-6.239831328560047) q[14];
rzz(0) q[6],q[14];
u3(pi/2,2.728159060377376,-2.728159060377376) q[6];
u3(pi,0.09927432785343747,-0.09927432785343747) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi,0.08984954989266808,-0.08984954989266808) q[12];
rzz(pi/128) q[6],q[12];
u3(pi,4.579185451872482,-4.579185451872482) q[11];
rzz(0.049100328574073315) q[6],q[11];
u3(pi,1.0065662862101699,-1.0065662862101699) q[10];
rzz(pi/32) q[6],q[10];
u3(pi,5.7868136679124,-5.7868136679124) q[9];
rzz(0.19640131429629326) q[6],q[9];
u3(pi,3.276681137694154,-3.276681137694154) q[8];
rzz(pi/8) q[6],q[8];
u3(pi/2,4.67217659441874,-4.67217659441874) q[7];
rzz(0.7853830994606742) q[6],q[7];
u3(pi/2,9*pi/8,-9*pi/8) q[5];
rzz(0) q[5],q[18];
u3(pi,0.26954864967800424,-0.26954864967800424) q[5];
rzz(0) q[5],q[17];
u3(pi,3.1642121206956397,-3.1642121206956397) q[5];
rzz(0) q[5],q[16];
u3(pi,6.058875591713275,-6.058875591713275) q[5];
rzz(0) q[5],q[15];
u3(pi,2.670353755551324,-2.670353755551324) q[5];
rzz(0) q[5],q[14];
u3(pi,5.56501722656896,-5.56501722656896) q[5];
u3(pi/2,5.123937618004953,-5.123937618004953) q[13];
rzz(0) q[5],q[13];
u3(pi/2,2.300274140958446,-2.300274140958446) q[5];
rzz(0.012275082143518329) q[5],q[12];
rzz(pi/128) q[5],q[11];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/32) q[5],q[9];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/8) q[5],q[7];
u3(pi/2,2.728159060377376,-2.728159060377376) q[6];
rzz(0.7853830994606742) q[5],q[6];
u3(pi/2,9*pi/8,-9*pi/8) q[4];
rzz(0) q[4],q[18];
u3(pi,0.26954864967800424,-0.26954864967800424) q[4];
rzz(0) q[4],q[17];
u3(pi,3.1642121206956397,-3.1642121206956397) q[4];
rzz(0) q[4],q[16];
u3(pi,6.058875591713275,-6.058875591713275) q[4];
rzz(0) q[4],q[15];
u3(pi,2.670353755551324,-2.670353755551324) q[4];
rzz(0) q[4],q[14];
u3(pi,5.56501722656896,-5.56501722656896) q[4];
rzz(0) q[4],q[13];
u3(pi,2.1771237089377267,-2.1771237089377267) q[4];
u3(pi/2,3.6586988043706734,-3.6586988043706734) q[12];
rzz(0) q[4],q[12];
u3(pi/2,5.194937611976082,-5.194937611976082) q[4];
u3(pi,2.189061761021368,-2.189061761021368) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi,5.942008344999735,-5.942008344999735) q[10];
rzz(pi/128) q[4],q[10];
u3(pi,4.462318205158942,-4.462318205158942) q[9];
rzz(0.049100328574073315) q[4],q[9];
u3(pi,1.4068051902775094,-1.4068051902775094) q[8];
rzz(pi/32) q[4],q[8];
u3(pi,1.5576016376498194,-1.5576016376498194) q[7];
rzz(0.19640131429629326) q[4],q[7];
u3(pi,1.4149733311768429,-1.4149733311768429) q[6];
rzz(pi/8) q[4],q[6];
u3(pi/2,2.300274140958446,-2.300274140958446) q[5];
rzz(0.7853830994606742) q[4],q[5];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(0) q[3],q[18];
u3(pi,5.924415426139632,-5.924415426139632) q[3];
rzz(0) q[3],q[17];
u3(pi,3.6367076557955444,-3.6367076557955444) q[3];
rzz(0) q[3],q[16];
u3(pi,1.348999885451457,-1.348999885451457) q[3];
rzz(0) q[3],q[15];
u3(pi,5.344477422286956,-5.344477422286956) q[3];
rzz(0) q[3],q[14];
u3(pi,3.0567696519428686,-3.0567696519428686) q[3];
rzz(0) q[3],q[13];
u3(pi,0.7690618815987813,-0.7690618815987813) q[3];
rzz(0) q[3],q[12];
u3(pi,4.76453941843428,-4.76453941843428) q[3];
u3(pi/2,4.44472528629884,-4.44472528629884) q[11];
rzz(0) q[3],q[11];
u3(pi/2,2.050203365732699,-2.050203365732699) q[3];
u3(pi,0.4693539424463151,-0.4693539424463151) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi,0.6572211831309848,-0.6572211831309848) q[9];
rzz(pi/128) q[3],q[9];
u3(pi,1.8868405477460297,-1.8868405477460297) q[8];
rzz(0.049100328574073315) q[3],q[8];
u3(pi,5.307406628974596,-5.307406628974596) q[7];
rzz(pi/32) q[3],q[7];
u3(pi,1.9848582385380313,-1.9848582385380313) q[6];
rzz(0.19640131429629326) q[3],q[6];
u3(pi,0.9104335510103221,-0.9104335510103221) q[5];
rzz(pi/8) q[3],q[5];
u3(pi/2,5.194937611976082,-5.194937611976082) q[4];
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
u3(pi/2,5.79561012734245,-5.79561012734245) q[10];
rzz(0) q[2],q[10];
u3(pi/2,5.26028273917075,-5.26028273917075) q[2];
u3(pi,3.159813890980614,-3.159813890980614) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi,2.3071856447963444,-2.3071856447963444) q[8];
rzz(pi/128) q[2],q[8];
u3(pi,2.6640705702441445,-2.6640705702441445) q[7];
rzz(0.049100328574073315) q[2],q[7];
u3(pi,2.614433406317426,-2.614433406317426) q[6];
rzz(pi/32) q[2],q[6];
u3(pi,1.2534954687823274,-1.2534954687823274) q[5];
rzz(0.19640131429629326) q[2],q[5];
u3(pi,3.8050970220279576,-3.8050970220279576) q[4];
rzz(pi/8) q[2],q[4];
u3(pi/2,5.191796019322492,-5.191796019322492) q[3];
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
u3(pi/2,6.200875579655533,-6.200875579655533) q[9];
rzz(0) q[1],q[9];
u3(pi/2,5.997300375702915,-5.997300375702915) q[1];
rzz(0.012275082143518329) q[1],q[8];
rzz(pi/128) q[1],q[7];
rzz(0.049100328574073315) q[1],q[6];
rzz(pi/32) q[1],q[5];
rzz(0.19640131429629326) q[1],q[4];
rzz(pi/8) q[1],q[3];
u3(pi/2,2.1186900855809565,-2.1186900855809565) q[2];
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
u3(pi/2,0.9946282341265285,-0.9946282341265285) q[8];
rzz(pi/2) q[0],q[8];
u3(pi/2,5.721468540717731,-5.721468540717731) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,1.3609379375350983,-1.3609379375350983) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,6.128618948622969,-6.128618948622969) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,2.415256432079833,-2.415256432079833) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,2.050203365732699,-2.050203365732699) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,5.26028273917075,-5.26028273917075) q[2];
rzz(pi/2) q[0],q[2];
rzz(-pi/2) q[0],q[1];
u3(pi,3*pi/2,-3*pi/2) q[0];
u3(pi/2,5.997300375702915,-5.997300375702915) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,0.54789375878606,-0.54789375878606) q[2];
u3(pi/2,2.904088248978405,-2.904088248978405) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,3.6209996925275956,-3.6209996925275956) q[3];
u3(pi/2,0.08670795723907829,-0.08670795723907829) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,3.986052758874729,-3.986052758874729) q[4];
u3(pi/2,0.6477964051702153,-0.6477964051702153) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,4.557822621828072,-4.557822621828072) q[5];
u3(pi/2,1.318212277446277,-1.318212277446277) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,2.931734264329995,-2.931734264329995) q[6];
u3(pi/2,6.024318072523787,-6.024318072523787) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,1.0090795603330416,-1.0090795603330416) q[7];
u3(pi/2,4.126167791224834,-4.126167791224834) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,2.565424560921425,-2.565424560921425) q[8];
u3(pi/2,5.694450843896859,-5.694450843896859) q[8];
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
u3(pi/2,4.4265040489080185,-4.4265040489080185) q[1];
u3(pi/2,4.474884575773301,-4.474884575773301) q[2];
rzz(0.7853830994606742) q[1],q[2];
u3(pi/2,1.6575042840339747,-1.6575042840339747) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,2.218592731965112,-2.218592731965112) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi/2,6.030601257830967,-6.030601257830967) q[5];
rzz(pi/32) q[1],q[5];
u3(pi/2,1.3119290921390976,-1.3119290921390976) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi/2,5.69696411801973,-5.69696411801973) q[7];
rzz(pi/128) q[1],q[7];
u3(pi/2,0.9820618635121693,-0.9820618635121693) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi/2,4.4265040489080185,-4.4265040489080185) q[1];
rzz(0) q[1],q[9];
u3(pi,0.7451857774314989,-0.7451857774314989) q[1];
rzz(0) q[1],q[10];
u3(pi,2.8079555137785572,-2.8079555137785572) q[1];
rzz(0) q[1],q[11];
u3(pi,4.8707252501256155,-4.8707252501256155) q[1];
rzz(0) q[1],q[12];
u3(pi,0.6496813607623693,-0.6496813607623693) q[1];
rzz(0) q[1],q[13];
u3(pi,2.7124510971094273,-2.7124510971094273) q[1];
rzz(0) q[1],q[14];
u3(pi,4.775220833456485,-4.775220833456485) q[1];
rzz(0) q[1],q[15];
u3(pi,0.5548052626239575,-0.5548052626239575) q[1];
rzz(0) q[1],q[16];
u3(pi,2.6169466804402974,-2.6169466804402974) q[1];
rzz(0) q[1],q[17];
u3(pi,4.679716416787356,-4.679716416787356) q[1];
rzz(0) q[1],q[18];
u3(pi/2,0.54789375878606,-0.54789375878606) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,3.0473448739820994,-3.0473448739820994) q[9];
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
u3(pi/2,1.2648052023352507,-1.2648052023352507) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,5.7805304826052195,-5.7805304826052195) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,5.97719418271994,-5.97719418271994) q[3];
rzz(0) q[3],q[11];
u3(pi,3.262229811487641,-3.262229811487641) q[3];
rzz(0) q[3],q[12];
u3(pi,0.9745220411435538,-0.9745220411435538) q[3];
rzz(0) q[3],q[13];
u3(pi,4.969999577979053,-4.969999577979053) q[3];
rzz(0) q[3],q[14];
u3(pi,2.6822918076349653,-2.6822918076349653) q[3];
rzz(0) q[3],q[15];
u3(pi,0.394584037290878,-0.394584037290878) q[3];
rzz(0) q[3],q[16];
u3(pi,4.390061574126377,-4.390061574126377) q[3];
rzz(0) q[3],q[17];
u3(pi,2.1023538037822895,-2.1023538037822895) q[3];
rzz(0) q[3],q[18];
u3(pi/2,5.164150003970902,-5.164150003970902) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,4.430902278623044,-4.430902278623044) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,0.4517610235862123,-0.4517610235862123) q[4];
rzz(0) q[4],q[12];
u3(pi,4.019981959533499,-4.019981959533499) q[4];
rzz(0) q[4],q[13];
u3(pi,1.732274189189412,-1.732274189189412) q[4];
rzz(0) q[4],q[14];
u3(pi,5.727751726024911,-5.727751726024911) q[4];
rzz(0) q[4],q[15];
u3(pi,3.4400439556808236,-3.4400439556808236) q[4];
rzz(0) q[4],q[16];
u3(pi,1.1523361853367362,-1.1523361853367362) q[4];
rzz(0) q[4],q[17];
u3(pi,5.147813722172235,-5.147813722172235) q[4];
rzz(0) q[4],q[18];
u3(pi/2,5.932583567038965,-5.932583567038965) q[5];
rzz(0.7853830994606742) q[5],q[6];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,0.4995132319207771,-0.4995132319207771) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,4.361787240244069,-4.361787240244069) q[5];
rzz(0) q[5],q[13];
u3(pi,0.41720350439672454,-0.41720350439672454) q[5];
rzz(0) q[5],q[14];
u3(pi,1.9534423120021334,-1.9534423120021334) q[5];
rzz(0) q[5],q[15];
u3(pi,3.489681119607542,-3.489681119607542) q[5];
rzz(0) q[5],q[16];
u3(pi,5.025919927212952,-5.025919927212952) q[5];
rzz(0) q[5],q[17];
u3(pi,0.27897342763877364,-0.27897342763877364) q[5];
rzz(0) q[5],q[18];
u3(pi/2,1.2629202467430969,-1.2629202467430969) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,5.10634469914485,-5.10634469914485) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,2.8337165735379934,-2.8337165735379934) q[6];
rzz(0) q[6],q[14];
u3(pi,5.359557067024187,-5.359557067024187) q[6];
rzz(0) q[6],q[15];
u3(pi,0.9858317746964772,-0.9858317746964772) q[6];
rzz(0) q[6],q[16];
u3(pi,2.895291789548353,-2.895291789548353) q[6];
rzz(0) q[6],q[17];
u3(pi,4.8047518044002295,-4.8047518044002295) q[6];
rzz(0) q[6],q[18];
u3(pi/2,2.530867041731937,-2.530867041731937) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,3.0800174375794334,-3.0800174375794334) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,0.9600707149370408,-0.9600707149370408) q[7];
rzz(0) q[7],q[15];
u3(pi,4.528291650884328,-4.528291650884328) q[7];
rzz(0) q[7],q[16];
u3(pi,2.2405838805402403,-2.2405838805402403) q[7];
rzz(0) q[7],q[17];
u3(pi,6.23606141737574,-6.23606141737574) q[7];
rzz(0) q[7],q[18];
u3(pi/2,0.9701238114285282,-0.9701238114285282) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,1.3760175822723293,-1.3760175822723293) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,2.5409201382234246,-2.5409201382234246) q[8];
rzz(0) q[8],q[16];
u3(pi,6.264335751258048,-6.264335751258048) q[8];
rzz(0) q[8],q[17];
u3(pi,4.286389016557914,-4.286389016557914) q[8];
rzz(0) q[8],q[18];
u3(pi/2,6.176371156957533,-6.176371156957533) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,0.9902300044115028,-0.9902300044115028) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,4.611858015469816,-4.611858015469816) q[9];
rzz(0) q[9],q[17];
u3(pi,2.4560971365765005,-2.4560971365765005) q[9];
rzz(0) q[9],q[18];
u3(pi/2,2.6295130510546567,-2.6295130510546567) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,0.9958848711879644,-0.9958848711879644) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi/2,4.200309377849553,-4.200309377849553) q[10];
rzz(0) q[10],q[18];
u3(pi/2,1.2817698026646356,-1.2817698026646356) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,1.8045308202219772,-1.8045308202219772) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,3.634194381672673,-3.634194381672673) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/128) q[12],q[18];
u3(pi/2,1.9584688602478768,-1.9584688602478768) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
u3(pi/2,6.215326905862047,-6.215326905862047) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
u3(pi/2,4.511327050554943,-4.511327050554943) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
u3(pi/2,4.125539472694116,-4.125539472694116) q[16];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/8) q[16],q[18];
u3(pi/2,0.9902300044115028,-0.9902300044115028) q[17];
rzz(0.7853830994606742) q[17],q[18];
u3(pi,pi,-pi) q[0];
u3(pi/2,0.9983981453108364,-0.9983981453108364) q[1];
u3(pi/2,5.8075481794260915,-5.8075481794260915) q[2];
u3(pi/2,5.6705747397295765,-5.6705747397295765) q[3];
u3(pi/2,2.432849350939936,-2.432849350939936) q[4];
u3(pi/2,2.6182033175017336,-2.6182033175017336) q[5];
u3(pi/2,1.046778672176119,-1.046778672176119) q[6];
u3(pi/2,3.521725364674158,-3.521725364674158) q[7];
u3(pi/2,1.7266193224129502,-1.7266193224129502) q[8];
u3(pi/2,0.29970793915246624,-0.29970793915246624) q[9];
u3(pi/2,1.0587167242597604,-1.0587167242597604) q[10];
u3(pi/2,4.939840288504591,-4.939840288504591) q[18];
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
