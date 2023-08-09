OPENQASM 2.0;
include "qelib1.inc";
qreg q[17];
creg c[17];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[15];
u3(pi/2,6.099716296209943,-6.099716296209943) q[16];
rzz(pi/2) q[15],q[16];
u3(pi/2,1.3873273158252526,-1.3873273158252526) q[16];
u3(pi/2,4.008672225980576,-4.008672225980576) q[16];
rzz(pi/2) q[15],q[16];
u3(pi/2,3.666238626739289,-3.666238626739289) q[14];
rzz(-pi/2) q[14],q[16];
u3(pi/2,0.8670795723907829,-0.8670795723907829) q[16];
u3(pi/2,2.9681767391116365,-2.9681767391116365) q[16];
rzz(pi/2) q[14],q[16];
u3(pi/2,5.789326942035271,-5.789326942035271) q[13];
rzz(pi/2) q[13],q[16];
u3(pi/2,2.9681767391116365,-2.9681767391116365) q[16];
u3(pi/2,4.028778418963551,-4.028778418963551) q[16];
rzz(pi/2) q[13],q[16];
u3(pi/2,3.666238626739289,-3.666238626739289) q[12];
rzz(-pi/2) q[12],q[16];
u3(pi/2,4.028778418963551,-4.028778418963551) q[16];
u3(pi/2,5.049167712849515,-5.049167712849515) q[16];
rzz(pi/2) q[12],q[16];
u3(pi/2,5.867238439844297,-5.867238439844297) q[11];
rzz(pi/2) q[11],q[16];
u3(pi/2,1.9075750592597223,-1.9075750592597223) q[16];
u3(pi/2,3.007760806546868,-3.007760806546868) q[16];
rzz(pi/2) q[11],q[16];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[10];
rzz(pi/2) q[10],q[16];
u3(pi/2,6.149353460136661,-6.149353460136661) q[16];
u3(pi/2,0.8080176305032948,-0.8080176305032948) q[16];
rzz(-pi/2) q[10],q[16];
u3(pi/2,6.185167616387585,-6.185167616387585) q[9];
rzz(pi/2) q[9],q[16];
u3(pi/2,0.8080176305032948,-0.8080176305032948) q[16];
u3(pi/2,2.065911329000648,-2.065911329000648) q[16];
rzz(pi/2) q[9],q[16];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[8];
rzz(pi/2) q[8],q[16];
u3(pi/2,5.207503982590441,-5.207503982590441) q[16];
u3(pi/2,5.833309239185528,-5.833309239185528) q[16];
rzz(-pi/2) q[8],q[16];
u3(pi/2,0.36819465900072373,-0.36819465900072373) q[7];
rzz(pi/2) q[7],q[16];
u3(pi/2,5.833309239185528,-5.833309239185528) q[16];
u3(pi/2,1.440106072405561,-1.440106072405561) q[16];
rzz(pi/2) q[7],q[16];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
rzz(pi/2) q[6],q[16];
u3(pi/2,1.440106072405561,-1.440106072405561) q[16];
u3(pi/2,2.077849381084289,-2.077849381084289) q[16];
rzz(-pi/2) q[6],q[16];
u3(pi/2,0,0) q[5];
rzz(pi/2) q[5],q[16];
u3(pi/2,2.077849381084289,-2.077849381084289) q[16];
u3(pi/2,3.943327098785909,-3.943327098785909) q[16];
rzz(pi/2) q[5],q[16];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[4],q[16];
u3(pi/2,3.943327098785909,-3.943327098785909) q[16];
u3(pi/2,4.532689880599354,-4.532689880599354) q[16];
rzz(-pi/2) q[4],q[16];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[3],q[16];
u3(pi/2,4.532689880599354,-4.532689880599354) q[16];
u3(pi/2,0.21299998191338798,-0.21299998191338798) q[16];
rzz(pi/2) q[3],q[16];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[2],q[16];
u3(pi/2,0.21299998191338798,-0.21299998191338798) q[16];
u3(pi/2,0.9983981453108364,-0.9983981453108364) q[16];
rzz(pi/2) q[2],q[16];
u3(pi/2,0,0) q[1];
u3(pi/2,5.710787125695526,-5.710787125695526) q[16];
rzz(pi/2) q[1],q[16];
u3(pi/2,pi,-pi) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,pi/2,-pi/2) q[1];
u3(pi/2,5*pi/4,-5*pi/4) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,15*pi/8,-15*pi/8) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,pi,-pi) q[3];
u3(pi/2,6.086521607064865,-6.086521607064865) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,pi,-pi) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,pi/2,-pi/2) q[5];
u3(pi/2,4.663380134988689,-4.663380134988689) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,5.301751762198135,-5.301751762198135) q[6];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,3.509787312590517,-3.509787312590517) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,5.080583639385413,-5.080583639385413) q[7];
u3(pi/2,1.926424615181261,-1.926424615181261) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[8];
rzz(-pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[8];
u3(pi/2,3.0435749627977917,-3.0435749627977917) q[9];
rzz(pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[9];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[10];
rzz(-pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[10];
u3(pi/2,2.7256457862545047,-2.7256457862545047) q[11];
rzz(pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[11];
u3(pi/2,3.666238626739289,-3.666238626739289) q[12];
rzz(-pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[12];
u3(pi/2,2.6477342884454775,-2.6477342884454775) q[13];
rzz(pi/2) q[0],q[13];
rzz(-pi/2) q[0],q[13];
u3(pi/2,3.666238626739289,-3.666238626739289) q[14];
rzz(-pi/2) q[0],q[14];
rzz(-pi/2) q[0],q[14];
u3(pi/2,3.669380219392878,-3.669380219392878) q[15];
rzz(-pi/2) q[0],q[15];
rzz(pi/2) q[0],q[15];
u3(pi/2,3*pi/4,-3*pi/4) q[1];
u3(pi/2,pi,-pi) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
u3(pi/2,13*pi/8,-13*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,6.087149925595583,-6.087149925595583) q[3];
u3(pi/2,2.552229871776348,-2.552229871776348) q[3];
rzz(-pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,5.399769452990137,-5.399769452990137) q[4];
u3(pi/2,2.061513099285622,-2.061513099285622) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,4.663380134988689,-4.663380134988689) q[5];
u3(pi/2,1.4237697906068942,-1.4237697906068942) q[5];
rzz(pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[6];
u3(pi/2,5.276619020969417,-5.276619020969417) q[6];
u3(pi/2,2.0860175219836226,-2.0860175219836226) q[6];
rzz(pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u3(pi/2,1.926424615181261,-1.926424615181261) q[7];
u3(pi/2,5.043512846073054,-5.043512846073054) q[7];
rzz(pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[8];
u3(pi/2,5.246459731494954,-5.246459731494954) q[8];
u3(pi/2,2.092300707290802,-2.092300707290802) q[8];
rzz(pi/2) q[1],q[8];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[10];
rzz(pi/2) q[1],q[10];
rzz(pi/2) q[1],q[11];
rzz(pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[12];
rzz(pi/2) q[1],q[12];
rzz(pi/2) q[1],q[13];
rzz(pi/2) q[1],q[13];
rzz(pi/2) q[1],q[14];
rzz(pi/2) q[1],q[14];
rzz(pi/2) q[1],q[15];
rzz(pi/2) q[1],q[15];
u3(pi/2,pi/8,-pi/8) q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[3];
rzz(-pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,2.061513099285622,-2.061513099285622) q[4];
u3(pi/2,4.810406671176691,-4.810406671176691) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,1.4237697906068942,-1.4237697906068942) q[5];
u3(pi/2,4.368698744081967,-4.368698744081967) q[5];
rzz(-pi/2) q[2],q[5];
rzz(pi/2) q[2],q[6];
u3(pi/2,2.0860175219836226,-2.0860175219836226) q[6];
u3(pi/2,5.129592484781415,-5.129592484781415) q[6];
rzz(pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,5.043512846073054,-5.043512846073054) q[7];
u3(pi/2,1.85291134708726,-1.85291134708726) q[7];
rzz(pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[8];
u3(pi/2,5.233893360880595,-5.233893360880595) q[8];
u3(pi/2,2.0677962845928017,-2.0677962845928017) q[8];
rzz(pi/2) q[2],q[8];
rzz(pi/2) q[2],q[9];
u3(pi/2,4.604946511631919,-4.604946511631919) q[9];
u3(pi/2,1.4514158059584845,-1.4514158059584845) q[9];
rzz(pi/2) q[2],q[9];
rzz(-pi/2) q[2],q[10];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[12];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[13];
rzz(-pi/2) q[2],q[14];
rzz(pi/2) q[2],q[14];
rzz(pi/2) q[2],q[15];
rzz(pi/2) q[2],q[15];
u3(pi/2,3.3376280351737964,-3.3376280351737964) q[3];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,4.810406671176691,-4.810406671176691) q[4];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[4];
rzz(-pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,1.2271060904921731,-1.2271060904921731) q[5];
u3(pi/2,3.9759996623832423,-3.9759996623832423) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,5.129592484781415,-5.129592484781415) q[6];
u3(pi/2,1.7919644496076181,-1.7919644496076181) q[6];
rzz(-pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,1.85291134708726,-1.85291134708726) q[7];
u3(pi/2,4.896486309885051,-4.896486309885051) q[7];
rzz(pi/2) q[3],q[7];
rzz(pi/2) q[3],q[8];
u3(pi/2,2.0677962845928017,-2.0677962845928017) q[8];
u3(pi/2,5.160380092786594,-5.160380092786594) q[8];
rzz(pi/2) q[3],q[8];
rzz(-pi/2) q[3],q[9];
u3(pi/2,4.593008459548278,-4.593008459548278) q[9];
u3(pi/2,1.426911383260484,-1.426911383260484) q[9];
rzz(pi/2) q[3],q[9];
rzz(pi/2) q[3],q[10];
u3(pi/2,5.23263672381916,-5.23263672381916) q[10];
u3(pi/2,2.078477699615007,-2.078477699615007) q[10];
rzz(pi/2) q[3],q[10];
rzz(-pi/2) q[3],q[11];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[12];
rzz(pi/2) q[3],q[12];
rzz(pi/2) q[3],q[13];
rzz(pi/2) q[3],q[13];
rzz(pi/2) q[3],q[14];
rzz(pi/2) q[3],q[14];
rzz(pi/2) q[3],q[15];
rzz(-pi/2) q[3],q[15];
u3(pi/2,5.595804834574139,-5.595804834574139) q[4];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,3.9759996623832423,-3.9759996623832423) q[5];
u3(pi/2,0.04900884539600077,-0.04900884539600077) q[5];
rzz(pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,4.933557103197411,-4.933557103197411) q[6];
u3(pi/2,1.399265367908894,-1.399265367908894) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,4.896486309885051,-4.896486309885051) q[7];
u3(pi/2,1.5582299561805373,-1.5582299561805373) q[7];
rzz(pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,5.160380092786594,-5.160380092786594) q[8];
u3(pi/2,1.9207697484047996,-1.9207697484047996) q[8];
rzz(-pi/2) q[4],q[8];
rzz(pi/2) q[4],q[9];
u3(pi/2,1.426911383260484,-1.426911383260484) q[9];
u3(pi/2,4.518866872923558,-4.518866872923558) q[9];
rzz(pi/2) q[4],q[9];
rzz(pi/2) q[4],q[10];
u3(pi/2,2.078477699615007,-2.078477699615007) q[10];
u3(pi/2,5.1955659305068,-5.1955659305068) q[10];
rzz(-pi/2) q[4],q[10];
rzz(pi/2) q[4],q[11];
u3(pi/2,1.1460530000295566,-1.1460530000295566) q[11];
u3(pi/2,4.275707601535708,-4.275707601535708) q[11];
rzz(pi/2) q[4],q[11];
rzz(pi/2) q[4],q[12];
rzz(pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[13];
rzz(pi/2) q[4],q[13];
rzz(pi/2) q[4],q[14];
rzz(pi/2) q[4],q[14];
rzz(-pi/2) q[4],q[15];
rzz(pi/2) q[4],q[15];
u3(pi/2,4.761397825780691,-4.761397825780691) q[5];
u3(pi/2,13*pi/8,-13*pi/8) q[5];
rzz(pi/2) q[5],q[6];
u3(pi/2,1.399265367908894,-1.399265367908894) q[6];
u3(pi/2,3.7554598581012386,-3.7554598581012386) q[6];
rzz(pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,1.5582299561805373,-1.5582299561805373) q[7];
u3(pi/2,4.307123528071607,-4.307123528071607) q[7];
rzz(-pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,5.0623624019945925,-5.0623624019945925) q[8];
u3(pi/2,1.7241060482900783,-1.7241060482900783) q[8];
rzz(pi/2) q[5],q[8];
rzz(pi/2) q[5],q[9];
u3(pi/2,4.518866872923558,-4.518866872923558) q[9];
u3(pi/2,1.2792565285417639,-1.2792565285417639) q[9];
rzz(-pi/2) q[5],q[9];
rzz(pi/2) q[5],q[10];
u3(pi/2,2.0539732769170067,-2.0539732769170067) q[10];
u3(pi/2,5.1465570851108,-5.1465570851108) q[10];
rzz(pi/2) q[5],q[10];
rzz(pi/2) q[5],q[11];
u3(pi/2,4.275707601535708,-4.275707601535708) q[11];
u3(pi/2,1.109610525247915,-1.109610525247915) q[11];
rzz(-pi/2) q[5],q[11];
rzz(pi/2) q[5],q[12];
u3(pi/2,5.228866812634852,-5.228866812634852) q[12];
u3(pi/2,2.0747077884306995,-2.0747077884306995) q[12];
rzz(pi/2) q[5],q[12];
rzz(pi/2) q[5],q[13];
rzz(-pi/2) q[5],q[13];
rzz(pi/2) q[5],q[14];
rzz(pi/2) q[5],q[14];
rzz(pi/2) q[5],q[15];
rzz(pi/2) q[5],q[15];
u3(pi/2,5.326256184896136,-5.326256184896136) q[6];
u3(pi/2,1.2271060904921731,-1.2271060904921731) q[6];
rzz(-pi/2) q[6],q[7];
u3(pi/2,4.307123528071607,-4.307123528071607) q[7];
u3(pi/2,0.38013271108436497,-0.38013271108436497) q[7];
rzz(pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,1.7241060482900783,-1.7241060482900783) q[8];
u3(pi/2,4.472999620181147,-4.472999620181147) q[8];
rzz(pi/2) q[6],q[8];
rzz(-pi/2) q[6],q[9];
u3(pi/2,1.2792565285417639,-1.2792565285417639) q[9];
u3(pi/2,4.224813800547554,-4.224813800547554) q[9];
rzz(pi/2) q[6],q[9];
rzz(pi/2) q[6],q[10];
u3(pi/2,5.1465570851108,-5.1465570851108) q[10];
u3(pi/2,1.9069467407290044,-1.9069467407290044) q[10];
rzz(pi/2) q[6],q[10];
rzz(pi/2) q[6],q[11];
u3(pi/2,4.251203178837708,-4.251203178837708) q[11];
u3(pi/2,1.0606016798519142,-1.0606016798519142) q[11];
rzz(pi/2) q[6],q[11];
rzz(pi/2) q[6],q[12];
u3(pi/2,2.0747077884306995,-2.0747077884306995) q[12];
u3(pi/2,5.191796019322492,-5.191796019322492) q[12];
rzz(pi/2) q[6],q[12];
rzz(pi/2) q[6],q[13];
u3(pi/2,4.209734155810323,-4.209734155810323) q[13];
u3(pi/2,1.0562034501368884,-1.0562034501368884) q[13];
rzz(-pi/2) q[6],q[13];
rzz(pi/2) q[6],q[14];
rzz(pi/2) q[6],q[14];
rzz(pi/2) q[6],q[15];
rzz(-pi/2) q[6],q[15];
u3(pi/2,5.092521691469055,-5.092521691469055) q[7];
u3(pi/2,1.6443095948888977,-1.6443095948888977) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,4.472999620181147,-4.472999620181147) q[8];
u3(pi/2,0.5460088031939061,-0.5460088031939061) q[8];
rzz(pi/2) q[7],q[8];
rzz(pi/2) q[7],q[9];
u3(pi/2,4.224813800547554,-4.224813800547554) q[9];
u3(pi/2,0.6905220652590365,-0.6905220652590365) q[9];
rzz(-pi/2) q[7],q[9];
rzz(pi/2) q[7],q[10];
u3(pi/2,1.9069467407290044,-1.9069467407290044) q[10];
u3(pi/2,4.851875694204076,-4.851875694204076) q[10];
rzz(pi/2) q[7],q[10];
rzz(pi/2) q[7],q[11];
u3(pi/2,1.0606016798519142,-1.0606016798519142) q[11];
u3(pi/2,4.103548324118988,-4.103548324118988) q[11];
rzz(-pi/2) q[7],q[11];
rzz(pi/2) q[7],q[12];
u3(pi/2,5.191796019322492,-5.191796019322492) q[12];
u3(pi/2,2.001194520336698,-2.001194520336698) q[12];
rzz(pi/2) q[7],q[12];
rzz(pi/2) q[7],q[13];
u3(pi/2,4.197796103726682,-4.197796103726682) q[13];
u3(pi/2,1.0316990274388882,-1.0316990274388882) q[13];
rzz(pi/2) q[7],q[13];
rzz(-pi/2) q[7],q[14];
u3(pi/2,5.228866812634852,-5.228866812634852) q[14];
u3(pi/2,2.0747077884306995,-2.0747077884306995) q[14];
rzz(pi/2) q[7],q[14];
rzz(pi/2) q[7],q[15];
rzz(pi/2) q[7],q[15];
u3(pi/2,2.1168051299888027,-2.1168051299888027) q[8];
u3(pi/2,3.055513014881433,-3.055513014881433) q[8];
rzz(-pi/2) q[8],q[9];
u3(pi/2,0.6905220652590365,-0.6905220652590365) q[9];
u3(pi/2,3.0467165554513813,-3.0467165554513813) q[9];
rzz(pi/2) q[8],q[9];
rzz(pi/2) q[8],q[10];
u3(pi/2,4.851875694204076,-4.851875694204076) q[10];
u3(pi/2,1.3175839589155591,-1.3175839589155591) q[10];
rzz(pi/2) q[8],q[10];
rzz(pi/2) q[8],q[11];
u3(pi/2,0.9619556705291947,-0.9619556705291947) q[11];
u3(pi/2,3.907512942534985,-3.907512942534985) q[11];
rzz(-pi/2) q[8],q[11];
rzz(pi/2) q[8],q[12];
u3(pi/2,2.001194520336698,-2.001194520336698) q[12];
u3(pi/2,5.044769483134489,-5.044769483134489) q[12];
rzz(pi/2) q[8],q[12];
rzz(pi/2) q[8],q[13];
u3(pi/2,1.0316990274388882,-1.0316990274388882) q[13];
u3(pi/2,4.123654517101962,-4.123654517101962) q[13];
rzz(-pi/2) q[8],q[13];
rzz(pi/2) q[8],q[14];
u3(pi/2,2.0747077884306995,-2.0747077884306995) q[14];
u3(pi/2,5.191796019322492,-5.191796019322492) q[14];
rzz(pi/2) q[8],q[14];
rzz(pi/2) q[8],q[15];
u3(pi/2,2.0891591146372126,-2.0891591146372126) q[15];
u3(pi/2,5.218813716143364,-5.218813716143364) q[15];
rzz(pi/2) q[8],q[15];
u3(pi/2,4.617512882246278,-4.617512882246278) q[9];
u3(pi/2,1.1598760077053516,-1.1598760077053516) q[9];
rzz(-pi/2) q[9],q[10];
u3(pi/2,4.4591766125053525,-4.4591766125053525) q[10];
u3(pi/2,0.532185795518111,-0.532185795518111) q[10];
rzz(pi/2) q[9],q[10];
rzz(pi/2) q[9],q[11];
u3(pi/2,0.7659202889451915,-0.7659202889451915) q[11];
u3(pi/2,3.5148138608362607,-3.5148138608362607) q[11];
rzz(-pi/2) q[9],q[11];
rzz(pi/2) q[9],q[12];
u3(pi/2,5.044769483134489,-5.044769483134489) q[12];
u3(pi/2,1.7065131294299756,-1.7065131294299756) q[12];
rzz(pi/2) q[9],q[12];
rzz(pi/2) q[9],q[13];
u3(pi/2,0.9820618635121693,-0.9820618635121693) q[13];
u3(pi/2,4.025636826309961,-4.025636826309961) q[13];
rzz(pi/2) q[9],q[13];
rzz(-pi/2) q[9],q[14];
u3(pi/2,2.050203365732699,-2.050203365732699) q[14];
u3(pi/2,5.142787173926491,-5.142787173926491) q[14];
rzz(pi/2) q[9],q[14];
rzz(pi/2) q[9],q[15];
u3(pi/2,5.218813716143364,-5.218813716143364) q[15];
u3(pi/2,2.052716639855571,-2.052716639855571) q[15];
rzz(pi/2) q[9],q[15];
u3(pi/2,2.1029821223130076,-2.1029821223130076) q[10];
u3(pi/2,2.730672334500248,-2.730672334500248) q[10];
rzz(-pi/2) q[10],q[11];
u3(pi/2,3.5148138608362607,-3.5148138608362607) q[11];
u3(pi/2,5.871008351028605,-5.871008351028605) q[11];
rzz(pi/2) q[10],q[11];
rzz(pi/2) q[10],q[12];
u3(pi/2,1.7065131294299756,-1.7065131294299756) q[12];
u3(pi/2,4.455406701321044,-4.455406701321044) q[12];
rzz(pi/2) q[10],q[12];
rzz(pi/2) q[10],q[13];
u3(pi/2,4.025636826309961,-4.025636826309961) q[13];
u3(pi/2,0.6880087911361646,-0.6880087911361646) q[13];
rzz(pi/2) q[10],q[13];
rzz(pi/2) q[10],q[14];
u3(pi/2,5.142787173926491,-5.142787173926491) q[14];
u3(pi/2,1.9031768295446967,-1.9031768295446967) q[14];
rzz(pi/2) q[10],q[14];
rzz(pi/2) q[10],q[15];
u3(pi/2,2.052716639855571,-2.052716639855571) q[15];
u3(pi/2,5.145300448049363,-5.145300448049363) q[15];
rzz(-pi/2) q[10],q[15];
u3(pi/2,1.1586193706439158,-1.1586193706439158) q[11];
u3(pi/2,3.8283448076645215,-3.8283448076645215) q[11];
rzz(pi/2) q[11],q[12];
u3(pi/2,4.455406701321044,-4.455406701321044) q[12];
u3(pi/2,0.5284158843338032,-0.5284158843338032) q[12];
rzz(pi/2) q[11],q[12];
rzz(pi/2) q[11],q[13];
u3(pi/2,0.6880087911361646,-0.6880087911361646) q[13];
u3(pi/2,3.436902363027234,-3.436902363027234) q[13];
rzz(-pi/2) q[11],q[13];
rzz(pi/2) q[11],q[14];
u3(pi/2,1.9031768295446967,-1.9031768295446967) q[14];
u3(pi/2,4.8481057830197685,-4.8481057830197685) q[14];
rzz(pi/2) q[11],q[14];
rzz(pi/2) q[11],q[15];
u3(pi/2,2.00370779445957,-2.00370779445957) q[15];
u3(pi/2,5.046654438726644,-5.046654438726644) q[15];
rzz(pi/2) q[11],q[15];
u3(pi/2,2.0992122111287,-2.0992122111287) q[12];
u3(pi/2,2.6489909255069133,-2.6489909255069133) q[12];
rzz(-pi/2) q[12],q[13];
u3(pi/2,3.436902363027234,-3.436902363027234) q[13];
u3(pi/2,5.793096853219579,-5.793096853219579) q[13];
rzz(pi/2) q[12],q[13];
rzz(pi/2) q[12],q[14];
u3(pi/2,4.8481057830197685,-4.8481057830197685) q[14];
u3(pi/2,1.3138140477312514,-1.3138140477312514) q[14];
rzz(-pi/2) q[12],q[14];
rzz(pi/2) q[12],q[15];
u3(pi/2,5.046654438726644,-5.046654438726644) q[15];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[15];
rzz(pi/2) q[12],q[15];
u3(pi/2,1.0807078728348887,-1.0807078728348887) q[13];
u3(pi/2,3.7114775609509816,-3.7114775609509816) q[13];
rzz(pi/2) q[13],q[14];
u3(pi/2,4.455406701321044,-4.455406701321044) q[14];
u3(pi/2,0.5284158843338032,-0.5284158843338032) q[14];
rzz(pi/2) q[13],q[14];
rzz(-pi/2) q[13],q[15];
u3(pi/2,4.850619057142641,-4.850619057142641) q[15];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[15];
rzz(pi/2) q[13],q[15];
u3(pi/2,5.240804864718492,-5.240804864718492) q[14];
u3(pi/2,1.5695396897334606,-1.5695396897334606) q[14];
rzz(pi/2) q[14],q[15];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[15];
u3(pi/2,3.6725218120464684,-3.6725218120464684) q[15];
rzz(pi/2) q[14],q[15];
u3(pi,3*pi/2,-3*pi/2) q[1];
u3(pi,3*pi/4,-3*pi/4) q[2];
u3(pi,5*pi/4,-5*pi/4) q[3];
u3(pi,0.5893627818134451,-0.5893627818134451) q[4];
u3(pi,1.3747609452108935,-1.3747609452108935) q[5];
u3(pi,3.2886191897777954,-3.2886191897777954) q[6];
u3(pi,1.1045839770021713,-1.1045839770021713) q[7];
u3(pi,5.927557018793221,-5.927557018793221) q[8];
u3(pi,0.23938936020354223,-0.23938936020354223) q[9];
u3(pi,6.197105668471226,-6.197105668471226) q[10];
u3(pi,4.002389040673396,-4.002389040673396) q[11];
u3(pi,5.875406580743632,-5.875406580743632) q[12];
u3(pi,2.1947166277978294,-2.1947166277978294) q[13];
u3(pi,3.164840439226358,-3.164840439226358) q[14];
u3(pi/2,5.2433181388413646,-5.2433181388413646) q[15];
u3(pi/2,0.010681415022205296,-0.010681415022205296) q[15];
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
