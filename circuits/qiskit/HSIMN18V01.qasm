OPENQASM 2.0;
include "qelib1.inc";
qreg q[18];
creg c[18];
u3(3.115983980606456,-1.5616042974456645,pi/2) q[0];
u3(0.0025855305740396513,-0.05786776311630515,pi/2) q[1];
cx q[1],q[0];
u3(0,0,-1.560328) q[0];
u3(1.560328,0.0,0.0) q[1];
cx q[0],q[1];
u3(-1.560328,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.025608672983337577,1.5799883561441286,3.1415889803846895) q[0];
u3(3.1058196320513964,1.6225444161149305,-pi/2) q[2];
u3(0.038562764822781315,0.057291430161645174,pi/2) q[3];
cx q[3],q[2];
u3(0,0,-1.560328) q[2];
u3(1.560328,0.0,0.0) q[3];
cx q[2],q[3];
u3(-1.560328,0.0,0.0) q[3];
cx q[3],q[2];
cx q[2],q[1];
u3(0,0,-1.560328) q[1];
u3(1.560328,0.0,0.0) q[2];
cx q[1],q[2];
u3(-1.560328,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.0025855305740396517,-0.057867763116304705,-3.673205103638111e-06) q[1];
cx q[1],q[0];
u3(0,0,-1.560328) q[0];
u3(1.560328,0.0,0.0) q[1];
cx q[0],q[1];
u3(-1.560328,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.025608672983337577,1.5799883561441286,3.1415889803846895) q[0];
u3(0.03577302153839667,-1.519048237474863,pi/2) q[2];
u3(3.092184378593017,1.5947552656450465,-pi/2) q[4];
u3(0.04161048111705812,-3.1356855381947124,-pi/2) q[5];
cx q[5],q[4];
u3(0,0,-1.560328) q[4];
u3(1.560328,0.0,0.0) q[5];
cx q[4],q[5];
u3(-1.560328,0.0,0.0) q[5];
cx q[5],q[4];
cx q[4],q[3];
u3(0,0,-1.560328) q[3];
u3(1.560328,0.0,0.0) q[4];
cx q[3],q[4];
u3(-1.560328,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.03856276482278132,0.057291430161645174,-3.6732051031940216e-06) q[3];
cx q[3],q[2];
u3(0,0,-1.560328) q[2];
u3(1.560328,0.0,0.0) q[3];
cx q[2],q[3];
u3(-1.560328,0.0,0.0) q[3];
cx q[3],q[2];
cx q[2],q[1];
u3(0,0,-1.560328) q[1];
u3(1.560328,0.0,0.0) q[2];
cx q[1],q[2];
u3(-1.560328,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.0025855305740396517,-0.057867763116304705,-3.673205103638111e-06) q[1];
cx q[1],q[0];
u3(0,0,-1.560328) q[0];
u3(1.560328,0.0,0.0) q[1];
cx q[0],q[1];
u3(-1.560328,0.0,0.0) q[1];
cx q[1],q[0];
u3(0,0,-1.5708) q[0];
u3(0.03577302153839667,-1.519048237474863,pi/2) q[2];
u3(0.049408274996776136,-1.5468373879447461,pi/2) q[4];
u3(3.110860019672049,-1.6124164145071633,pi/2) q[6];
u3(0.0005077145349519084,-3.1258205187600208,-pi/2) q[7];
cx q[7],q[6];
u3(0,0,-1.560328) q[6];
u3(1.560328,0.0,0.0) q[7];
cx q[6],q[7];
u3(-1.560328,0.0,0.0) q[7];
cx q[7],q[6];
cx q[6],q[5];
u3(0,0,-1.560328) q[5];
u3(1.560328,0.0,0.0) q[6];
cx q[5],q[6];
u3(-1.560328,0.0,0.0) q[6];
cx q[6],q[5];
u3(0.041610481117058125,-3.1356855381947124,3.1415889803846895) q[5];
cx q[5],q[4];
u3(0,0,-1.560328) q[4];
u3(1.560328,0.0,0.0) q[5];
cx q[4],q[5];
u3(-1.560328,0.0,0.0) q[5];
cx q[5],q[4];
cx q[4],q[3];
u3(0,0,-1.560328) q[3];
u3(1.560328,0.0,0.0) q[4];
cx q[3],q[4];
u3(-1.560328,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.03856276482278132,0.057291430161645174,-3.6732051031940216e-06) q[3];
cx q[3],q[2];
u3(0,0,-1.560328) q[2];
u3(1.560328,0.0,0.0) q[3];
cx q[2],q[3];
u3(-1.560328,0.0,0.0) q[3];
cx q[3],q[2];
cx q[2],q[1];
u3(0,0,-1.560328) q[1];
u3(1.560328,0.0,0.0) q[2];
cx q[1],q[2];
u3(-1.560328,0.0,0.0) q[2];
cx q[2],q[1];
u3(0,0,-1.5708) q[1];
u3(0.049408274996776136,-1.5468373879447461,pi/2) q[4];
u3(0.030732633917744504,1.5291762390826298,-pi/2) q[6];
u3(3.109709753974632,1.5540133625796138,-pi/2) q[8];
u3(0.04067306778111446,-3.0837780246151008,-pi/2) q[9];
cx q[9],q[8];
u3(0,0,-1.560328) q[8];
u3(1.560328,0.0,0.0) q[9];
cx q[8],q[9];
u3(-1.560328,0.0,0.0) q[9];
cx q[9],q[8];
cx q[8],q[7];
u3(0,0,-1.560328) q[7];
u3(1.560328,0.0,0.0) q[8];
cx q[7],q[8];
u3(-1.560328,0.0,0.0) q[8];
cx q[8],q[7];
u3(0.0005077145349519085,-3.1258205187600216,3.1415889803846895) q[7];
cx q[7],q[6];
u3(0,0,-1.560328) q[6];
u3(1.560328,0.0,0.0) q[7];
cx q[6],q[7];
u3(-1.560328,0.0,0.0) q[7];
cx q[7],q[6];
cx q[6],q[5];
u3(0,0,-1.560328) q[5];
u3(1.560328,0.0,0.0) q[6];
cx q[5],q[6];
u3(-1.560328,0.0,0.0) q[6];
cx q[6],q[5];
u3(0.041610481117058125,-3.1356855381947124,3.1415889803846895) q[5];
cx q[5],q[4];
u3(0,0,-1.560328) q[4];
u3(1.560328,0.0,0.0) q[5];
cx q[4],q[5];
u3(-1.560328,0.0,0.0) q[5];
cx q[5],q[4];
cx q[4],q[3];
u3(0,0,-1.560328) q[3];
u3(1.560328,0.0,0.0) q[4];
cx q[3],q[4];
u3(-1.560328,0.0,0.0) q[4];
cx q[4],q[3];
u3(0,0,-1.5708) q[3];
u3(0.030732633917744504,1.5291762390826298,-pi/2) q[6];
u3(0.031882899615161,-1.5875792910101798,pi/2) q[8];
u3(3.098892019819915,-1.551044166471155,pi/2) q[10];
u3(0.005176834476686487,0.02290916352230976,pi/2) q[11];
cx q[11],q[10];
u3(0,0,-1.560328) q[10];
u3(1.560328,0.0,0.0) q[11];
cx q[10],q[11];
u3(-1.560328,0.0,0.0) q[11];
cx q[11],q[10];
cx q[10],q[9];
u3(1.560328,0.0,0.0) q[10];
u3(0,0,-1.560328) q[9];
cx q[9],q[10];
u3(-1.560328,0.0,0.0) q[10];
cx q[10],q[9];
u3(0.04270063376987876,1.590548487118638,-pi/2) q[10];
u3(0.04067306778111445,-3.083778024615101,3.1415889803846895) q[9];
cx q[9],q[8];
u3(0,0,-1.560328) q[8];
u3(1.560328,0.0,0.0) q[9];
cx q[8],q[9];
u3(-1.560328,0.0,0.0) q[9];
cx q[9],q[8];
cx q[8],q[7];
u3(0,0,-1.560328) q[7];
u3(1.560328,0.0,0.0) q[8];
cx q[7],q[8];
u3(-1.560328,0.0,0.0) q[8];
cx q[8],q[7];
u3(0.0005077145349519085,-3.1258205187600216,3.1415889803846895) q[7];
cx q[7],q[6];
u3(0,0,-1.560328) q[6];
u3(1.560328,0.0,0.0) q[7];
cx q[6],q[7];
u3(-1.560328,0.0,0.0) q[7];
cx q[7],q[6];
cx q[6],q[5];
u3(0,0,-1.560328) q[5];
u3(1.560328,0.0,0.0) q[6];
cx q[5],q[6];
u3(-1.560328,0.0,0.0) q[6];
cx q[6],q[5];
u3(0,0,-1.5708) q[5];
u3(0.031882899615161,-1.5875792910101798,pi/2) q[8];
u3(3.112628067537924,-1.5132531844183195,pi/2) q[12];
u3(0.03479198259494653,0.022668751260862763,pi/2) q[13];
cx q[13],q[12];
u3(0,0,-1.560328) q[12];
u3(1.560328,0.0,0.0) q[13];
cx q[12],q[13];
u3(-1.560328,0.0,0.0) q[13];
cx q[13],q[12];
cx q[12],q[11];
u3(0,0,-1.560328) q[11];
u3(1.560328,0.0,0.0) q[12];
cx q[11],q[12];
u3(-1.560328,0.0,0.0) q[12];
cx q[12],q[11];
u3(0.005176834476686487,0.022909163522309317,-3.6732051031940216e-06) q[11];
cx q[11],q[10];
u3(0,0,-1.560328) q[10];
u3(1.560328,0.0,0.0) q[11];
cx q[10],q[11];
u3(-1.560328,0.0,0.0) q[11];
cx q[11],q[10];
cx q[10],q[9];
u3(1.560328,0.0,0.0) q[10];
u3(0.02896458605186947,1.6283394691714737,-pi/2) q[12];
u3(0,0,-1.560328) q[9];
cx q[9],q[10];
u3(-1.560328,0.0,0.0) q[10];
cx q[10],q[9];
u3(0.04270063376987876,1.590548487118638,-pi/2) q[10];
u3(0.04067306778111445,-3.083778024615101,3.1415889803846895) q[9];
cx q[9],q[8];
u3(0,0,-1.560328) q[8];
u3(1.560328,0.0,0.0) q[9];
cx q[8],q[9];
u3(-1.560328,0.0,0.0) q[9];
cx q[9],q[8];
cx q[8],q[7];
u3(0,0,-1.560328) q[7];
u3(1.560328,0.0,0.0) q[8];
cx q[7],q[8];
u3(-1.560328,0.0,0.0) q[8];
cx q[8],q[7];
u3(0,0,-1.5708) q[7];
u3(3.1108850458092534,-1.6065130586118777,pi/2) q[14];
u3(0.02753073209146705,3.109040906934305,-pi/2) q[15];
cx q[15],q[14];
u3(0,0,-1.560328) q[14];
u3(1.560328,0.0,0.0) q[15];
cx q[14],q[15];
u3(-1.560328,0.0,0.0) q[15];
cx q[15],q[14];
cx q[14],q[13];
u3(0,0,-1.560328) q[13];
u3(1.560328,0.0,0.0) q[14];
cx q[13],q[14];
u3(-1.560328,0.0,0.0) q[14];
cx q[14],q[13];
u3(0.03479198259494653,0.022668751260862763,-3.6732051031940216e-06) q[13];
cx q[13],q[12];
u3(0,0,-1.560328) q[12];
u3(1.560328,0.0,0.0) q[13];
cx q[12],q[13];
u3(-1.560328,0.0,0.0) q[13];
cx q[13],q[12];
cx q[12],q[11];
u3(0,0,-1.560328) q[11];
u3(1.560328,0.0,0.0) q[12];
cx q[11],q[12];
u3(-1.560328,0.0,0.0) q[12];
cx q[12],q[11];
u3(0.005176834476686487,0.022909163522309317,-3.6732051031940216e-06) q[11];
cx q[11],q[10];
u3(0,0,-1.560328) q[10];
u3(1.560328,0.0,0.0) q[11];
cx q[10],q[11];
u3(-1.560328,0.0,0.0) q[11];
cx q[11],q[10];
cx q[10],q[9];
u3(1.560328,0.0,0.0) q[10];
u3(0.02896458605186947,1.6283394691714737,-pi/2) q[12];
u3(0.030707607780540165,1.5350795949779155,-pi/2) q[14];
u3(0,0,-1.560328) q[9];
cx q[9],q[10];
u3(-1.560328,0.0,0.0) q[10];
cx q[10],q[9];
u3(0,0,-1.5708) q[9];
u3(3.079938585374757,1.5341285699522773,-pi/2) q[16];
u3(0.05509524348948113,-0.05590290607771742,pi/2) q[17];
cx q[17],q[16];
u3(0,0,-1.560328) q[16];
u3(1.560328,0.0,0.0) q[17];
cx q[16],q[17];
u3(-1.560328,0.0,0.0) q[17];
cx q[17],q[16];
cx q[16],q[15];
u3(0,0,-1.560328) q[15];
u3(1.560328,0.0,0.0) q[16];
cx q[15],q[16];
u3(-1.560328,0.0,0.0) q[16];
cx q[16],q[15];
u3(0.027530732091467055,3.109040906934304,3.1415889803846895) q[15];
cx q[15],q[14];
u3(0,0,-1.560328) q[14];
u3(1.560328,0.0,0.0) q[15];
cx q[14],q[15];
u3(-1.560328,0.0,0.0) q[15];
cx q[15],q[14];
cx q[14],q[13];
u3(0,0,-1.560328) q[13];
u3(1.560328,0.0,0.0) q[14];
cx q[13],q[14];
u3(-1.560328,0.0,0.0) q[14];
cx q[14],q[13];
u3(0.03479198259494653,0.022668751260862763,-3.6732051031940216e-06) q[13];
cx q[13],q[12];
u3(0,0,-1.560328) q[12];
u3(1.560328,0.0,0.0) q[13];
cx q[12],q[13];
u3(-1.560328,0.0,0.0) q[13];
cx q[13],q[12];
cx q[12],q[11];
u3(0,0,-1.560328) q[11];
u3(1.560328,0.0,0.0) q[12];
cx q[11],q[12];
u3(-1.560328,0.0,0.0) q[12];
cx q[12],q[11];
u3(0,0,-1.5708) q[11];
u3(0.030707607780540165,1.5350795949779155,-pi/2) q[14];
u3(0.06165406821503616,-1.6074640836375158,pi/2) q[16];
u3(0.05509524348948113,-0.05590290607771742,pi/2) q[17];
cx q[17],q[16];
u3(0,0,-1.560328) q[16];
u3(1.560328,0.0,0.0) q[17];
cx q[16],q[17];
u3(-1.560328,0.0,0.0) q[17];
cx q[17],q[16];
cx q[16],q[15];
u3(0,0,-1.560328) q[15];
u3(1.560328,0.0,0.0) q[16];
cx q[15],q[16];
u3(-1.560328,0.0,0.0) q[16];
cx q[16],q[15];
u3(0.027530732091467055,3.109040906934304,3.1415889803846895) q[15];
cx q[15],q[14];
u3(0,0,-1.560328) q[14];
u3(1.560328,0.0,0.0) q[15];
cx q[14],q[15];
u3(-1.560328,0.0,0.0) q[15];
cx q[15],q[14];
cx q[14],q[13];
u3(0,0,-1.560328) q[13];
u3(1.560328,0.0,0.0) q[14];
cx q[13],q[14];
u3(-1.560328,0.0,0.0) q[14];
cx q[14],q[13];
u3(0,0,-1.5708) q[13];
u3(0.06165406821503616,-1.6074640836375158,pi/2) q[16];
u3(0.05509524348948113,-0.05590290607771742,pi/2) q[17];
cx q[17],q[16];
u3(0,0,-1.560328) q[16];
u3(1.560328,0.0,0.0) q[17];
cx q[16],q[17];
u3(-1.560328,0.0,0.0) q[17];
cx q[17],q[16];
cx q[16],q[15];
u3(0,0,-1.560328) q[15];
u3(1.560328,0.0,0.0) q[16];
cx q[15],q[16];
u3(-1.560328,0.0,0.0) q[16];
cx q[16],q[15];
u3(0,0,-1.5708) q[15];
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
