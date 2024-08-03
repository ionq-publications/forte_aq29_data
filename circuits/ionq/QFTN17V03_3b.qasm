OPENQASM 2.0;
include "qelib1.inc";
qreg q[17];
creg c[17];
u(pi/2,-1.5085927922538187,1.5085927922538187) q[16];
u(pi/2,3.9891943515283197,-3.9891943515283197) q[16];
rzz(-pi/2) q[15],q[16];
rzz(-pi/2) q[14],q[16];
u(pi/2,0.8476016979385261,-0.8476016979385261) q[16];
u(pi/2,4.381893433227043,-4.381893433227043) q[16];
rzz(-pi/2) q[14],q[16];
rzz(-pi/2) q[14],q[15];
u(pi/2,-1.4162299682382788,1.4162299682382788) q[15];
u(pi/2,4.081557175543859,-4.081557175543859) q[15];
rzz(-pi/2) q[14],q[15];
rzz(-pi/2) q[13],q[16];
u(pi/2,4.381893433227043,-4.381893433227043) q[16];
u(pi/2,1.4363361612212535,-1.4363361612212535) q[16];
rzz(-pi/2) q[13],q[16];
rzz(pi/2) q[13],q[15];
u(pi/2,4.081557175543859,-4.081557175543859) q[15];
u(pi/2,1.3326636036527901,-1.3326636036527901) q[15];
rzz(pi/2) q[13],q[15];
rzz(pi/2) q[13],q[14];
u(pi/2,3.3891501546926683,-3.3891501546926683) q[14];
u(pi/2,2.60375199129522,-2.60375199129522) q[14];
rzz(-pi/2) q[13],q[14];
rzz(pi/2) q[12],q[16];
u(pi/2,1.4363361612212535,-1.4363361612212535) q[16];
u(pi/2,4.675946505603048,-4.675946505603048) q[16];
rzz(-pi/2) q[12],q[16];
rzz(pi/2) q[12],q[15];
u(pi/2,4.474256257242583,-4.474256257242583) q[15];
u(pi/2,1.5286989852367934,-1.5286989852367934) q[15];
rzz(pi/2) q[12],q[15];
rzz(-pi/2) q[12],q[14];
u(pi/2,2.60375199129522,-2.60375199129522) q[14];
u(pi/2,-0.14514158059584847,0.14514158059584847) q[14];
rzz(-pi/2) q[12],q[14];
rzz(pi/2) q[12],q[13];
u(pi/2,3.7592297692855468,-3.7592297692855468) q[13];
u(pi/2,2.9738316058880985,-2.9738316058880985) q[13];
rzz(pi/2) q[12],q[13];
rzz(-pi/2) q[11],q[16];
u(pi/2,4.675946505603048,-4.675946505603048) q[16];
u(pi/2,1.5833626974092558,-1.5833626974092558) q[16];
rzz(-pi/2) q[11],q[16];
rzz(pi/2) q[11],q[15];
u(pi/2,1.5286989852367934,-1.5286989852367934) q[15];
u(pi/2,-1.5148759775609983,1.5148759775609983) q[15];
rzz(pi/2) q[11],q[15];
rzz(pi/2) q[11],q[14];
u(pi/2,2.9964510729939446,-2.9964510729939446) q[14];
u(pi/2,0.05089380098815455,-0.05089380098815455) q[14];
rzz(pi/2) q[11],q[14];
rzz(-pi/2) q[11],q[13];
u(pi/2,2.9738316058880985,-2.9738316058880985) q[13];
u(pi/2,0.22493803399702927,-0.22493803399702927) q[13];
rzz(pi/2) q[11],q[13];
rzz(pi/2) q[11],q[12];
u(pi/2,4.131194339470578,-4.131194339470578) q[12];
u(pi/2,3.3457961760731294,-3.3457961760731294) q[12];
rzz(pi/2) q[11],q[12];
rzz(-pi/2) q[10],q[16];
u(pi/2,1.5833626974092558,-1.5833626974092558) q[16];
u(pi/2,-1.533097214951819,1.533097214951819) q[16];
rzz(-pi/2) q[10],q[16];
rzz(-pi/2) q[10],q[15];
u(pi/2,1.6267166760287952,-1.6267166760287952) q[15];
u(pi/2,-1.4652388136342795,1.4652388136342795) q[15];
rzz(-pi/2) q[10],q[15];
rzz(pi/2) q[10],q[14];
u(pi/2,0.05089380098815455,-0.05089380098815455) q[14];
u(pi/2,3.2905041453699493,-3.2905041453699493) q[14];
rzz(-pi/2) q[10],q[14];
rzz(pi/2) q[10],q[13];
u(pi/2,0.22493803399702927,-0.22493803399702927) q[13];
u(pi/2,3.562566069170825,-3.562566069170825) q[13];
rzz(-pi/2) q[10],q[13];
rzz(-pi/2) q[10],q[12];
u(pi/2,3.3457961760731294,-3.3457961760731294) q[12];
u(pi/2,0.5969026041820604,-0.5969026041820604) q[12];
rzz(pi/2) q[10],q[12];
rzz(pi/2) q[10],q[11];
u(pi/2,-0.6716725093374978,0.6716725093374978) q[11];
u(pi/2,-1.457070672734946,1.457070672734946) q[11];
rzz(-pi/2) q[10],q[11];
rzz(-pi/2) q[9],q[16];
u(pi/2,-1.533097214951819,1.533097214951819) q[16];
u(pi/2,1.6204334907216156,-1.6204334907216156) q[16];
rzz(-pi/2) q[9],q[16];
rzz(pi/2) q[9],q[15];
u(pi/2,1.6763538399555138,-1.6763538399555138) q[15];
u(pi/2,-1.440734390936279,1.440734390936279) q[15];
rzz(-pi/2) q[9],q[15];
rzz(-pi/2) q[9],q[14];
u(pi/2,3.2905041453699493,-3.2905041453699493) q[14];
u(pi/2,0.19854865570687497,-0.19854865570687497) q[14];
rzz(pi/2) q[9],q[14];
rzz(-pi/2) q[9],q[13];
u(pi/2,3.562566069170825,-3.562566069170825) q[13];
u(pi/2,0.5189911063730337,-0.5189911063730337) q[13];
rzz(pi/2) q[9],q[13];
rzz(-pi/2) q[9],q[12];
u(pi/2,3.738495257771854,-3.738495257771854) q[12];
u(pi/2,0.7929379857660637,-0.7929379857660637) q[12];
rzz(pi/2) q[9],q[12];
rzz(-pi/2) q[9],q[11];
u(pi/2,1.6845219808548473,-1.6845219808548473) q[11];
u(pi/2,-1.064371591036222,1.064371591036222) q[11];
rzz(pi/2) q[9],q[11];
rzz(-pi/2) q[9],q[10];
u(pi/2,0.8161857714026284,-0.8161857714026284) q[10];
u(pi/2,0.0307876080051801,-0.0307876080051801) q[10];
rzz(-pi/2) q[9],q[10];
rzz(-pi/2) q[8],q[16];
rzz(-pi/2) q[8],q[16];
rzz(-pi/2) q[8],q[15];
u(pi/2,1.7008582626535143,-1.7008582626535143) q[15];
u(pi/2,-1.428796338852638,1.428796338852638) q[15];
rzz(pi/2) q[8],q[15];
rzz(-pi/2) q[8],q[14];
u(pi/2,0.19854865570687497,-0.19854865570687497) q[14];
u(pi/2,3.364645731994668,-3.364645731994668) q[14];
rzz(-pi/2) q[8],q[14];
rzz(-pi/2) q[8],q[13];
u(pi/2,0.5189911063730337,-0.5189911063730337) q[13];
u(pi/2,3.710220923889546,-3.710220923889546) q[13];
rzz(-pi/2) q[8],q[13];
rzz(pi/2) q[8],q[12];
u(pi/2,3.934530639355857,-3.934530639355857) q[12];
u(pi/2,0.8909556765580651,-0.8909556765580651) q[12];
rzz(-pi/2) q[8],q[12];
rzz(-pi/2) q[8],q[11];
u(pi/2,-1.064371591036222,1.064371591036222) q[11];
u(pi/2,2.273256444137574,-2.273256444137574) q[11];
rzz(pi/2) q[8],q[11];
rzz(pi/2) q[8],q[10];
u(pi/2,3.1723802615949728,-3.1723802615949728) q[10];
u(pi/2,0.42348668970390424,-0.42348668970390424) q[10];
rzz(pi/2) q[8],q[10];
rzz(pi/2) q[8],q[9];
u(pi/2,3.595866951298877,-3.595866951298877) q[9];
u(pi/2,2.8104687879014287,-2.8104687879014287) q[9];
rzz(pi/2) q[8],q[9];
rzz(-pi/2) q[7],q[16];
rzz(-pi/2) q[7],q[16];
rzz(pi/2) q[7],q[15];
rzz(-pi/2) q[7],q[15];
rzz(-pi/2) q[7],q[14];
u(pi/2,3.364645731994668,-3.364645731994668) q[14];
u(pi/2,0.2349911304885164,-0.2349911304885164) q[14];
rzz(-pi/2) q[7],q[14];
rzz(pi/2) q[7],q[13];
u(pi/2,0.5686282702997527,-0.5686282702997527) q[13];
u(pi/2,3.7347253465875463,-3.7347253465875463) q[13];
rzz(pi/2) q[7],q[13];
rzz(-pi/2) q[7],q[12];
u(pi/2,0.8909556765580651,-0.8909556765580651) q[12];
u(pi/2,4.082185494074578,-4.082185494074578) q[12];
rzz(-pi/2) q[7],q[12];
rzz(-pi/2) q[7],q[11];
u(pi/2,-0.8683362094522188,0.8683362094522188) q[11];
u(pi/2,2.3712741349295756,-2.3712741349295756) q[11];
rzz(-pi/2) q[7],q[11];
rzz(pi/2) q[7],q[10];
u(pi/2,0.42348668970390424,-0.42348668970390424) q[10];
u(pi/2,3.7611147248777,-3.7611147248777) q[10];
rzz(-pi/2) q[7],q[10];
rzz(pi/2) q[7],q[9];
u(pi/2,-0.33112386568836416,0.33112386568836416) q[9];
u(pi/2,3.2031678696001533,-3.2031678696001533) q[9];
rzz(-pi/2) q[7],q[9];
rzz(pi/2) q[7],q[8];
u(pi/2,3.2641147670797945,-3.2641147670797945) q[8];
u(pi/2,2.4787166036823463,-2.4787166036823463) q[8];
rzz(pi/2) q[7],q[8];
rzz(-pi/2) q[6],q[16];
rzz(pi/2) q[6],q[16];
rzz(pi/2) q[6],q[15];
rzz(-pi/2) q[6],q[15];
rzz(pi/2) q[6],q[14];
rzz(-pi/2) q[6],q[14];
rzz(-pi/2) q[6],q[13];
u(pi/2,0.5931326929977527,-0.5931326929977527) q[13];
u(pi/2,3.7466633986711875,-3.7466633986711875) q[13];
rzz(-pi/2) q[6],q[13];
rzz(-pi/2) q[6],q[12];
u(pi/2,4.082185494074578,-4.082185494074578) q[12];
u(pi/2,0.9650972631827845,-0.9650972631827845) q[12];
rzz(pi/2) q[6],q[12];
rzz(-pi/2) q[6],q[11];
u(pi/2,2.3712741349295756,-2.3712741349295756) q[11];
u(pi/2,-0.7206813547334985,0.7206813547334985) q[11];
rzz(-pi/2) q[6],q[11];
rzz(-pi/2) q[6],q[10];
u(pi/2,3.7611147248777,-3.7611147248777) q[10];
u(pi/2,0.718168080610627,-0.718168080610627) q[10];
rzz(-pi/2) q[6],q[10];
rzz(pi/2) q[6],q[9];
u(pi/2,0.06157521601035976,-0.06157521601035976) q[9];
u(pi/2,3.399203251184156,-3.399203251184156) q[9];
rzz(pi/2) q[6],q[9];
rzz(-pi/2) q[6],q[8];
u(pi/2,2.4787166036823463,-2.4787166036823463) q[8];
u(pi/2,-0.2701769682087223,0.2701769682087223) q[8];
rzz(-pi/2) q[6],q[8];
rzz(pi/2) q[6],q[7];
u(pi/2,3.2151059216837945,-3.2151059216837945) q[7];
u(pi/2,2.429707758286346,-2.429707758286346) q[7];
rzz(pi/2) q[6],q[7];
rzz(pi/2) q[5],q[16];
rzz(-pi/2) q[5],q[16];
rzz(-pi/2) q[5],q[15];
rzz(pi/2) q[5],q[15];
rzz(pi/2) q[5],q[14];
rzz(-pi/2) q[5],q[14];
rzz(-pi/2) q[5],q[13];
rzz(pi/2) q[5],q[13];
rzz(-pi/2) q[5],q[12];
u(pi/2,4.106689916772577,-4.106689916772577) q[12];
u(pi/2,0.9770353152664257,-0.9770353152664257) q[12];
rzz(pi/2) q[5],q[12];
rzz(-pi/2) q[5],q[11];
u(pi/2,-0.7206813547334985,0.7206813547334985) q[11];
u(pi/2,2.445415721554295,-2.445415721554295) q[11];
rzz(-pi/2) q[5],q[11];
rzz(-pi/2) q[5],q[10];
u(pi/2,0.718168080610627,-0.718168080610627) q[10];
u(pi/2,3.90876957959642,-3.90876957959642) q[10];
rzz(-pi/2) q[5],q[10];
rzz(pi/2) q[5],q[9];
u(pi/2,3.399203251184156,-3.399203251184156) q[9];
u(pi/2,0.3556282883863644,-0.3556282883863644) q[9];
rzz(-pi/2) q[5],q[9];
rzz(-pi/2) q[5],q[8];
u(pi/2,-0.2701769682087223,0.2701769682087223) q[8];
u(pi/2,3.0680793854957917,-3.0680793854957917) q[8];
rzz(-pi/2) q[5],q[8];
rzz(-pi/2) q[5],q[7];
u(pi/2,2.429707758286346,-2.429707758286346) q[7];
u(pi/2,-0.319185813604723,0.319185813604723) q[7];
rzz(pi/2) q[5],q[7];
rzz(pi/2) q[5],q[6];
u(pi/2,4.025008507779242,-4.025008507779242) q[6];
u(pi/2,3.239610344381794,-3.239610344381794) q[6];
rzz(-pi/2) q[5],q[6];
rzz(-pi/2) q[4],q[16];
rzz(-pi/2) q[4],q[16];
rzz(pi/2) q[4],q[15];
rzz(-pi/2) q[4],q[15];
rzz(-pi/2) q[4],q[14];
rzz(-pi/2) q[4],q[14];
rzz(-pi/2) q[4],q[13];
rzz(pi/2) q[4],q[13];
rzz(-pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[11];
u(pi/2,-0.6961769320354982,0.6961769320354982) q[11];
u(pi/2,2.4573537736379363,-2.4573537736379363) q[11];
rzz(-pi/2) q[4],q[11];
rzz(pi/2) q[4],q[10];
u(pi/2,3.90876957959642,-3.90876957959642) q[10];
u(pi/2,0.7916813487046279,-0.7916813487046279) q[10];
rzz(-pi/2) q[4],q[10];
rzz(-pi/2) q[4],q[9];
u(pi/2,3.497220941976158,-3.497220941976158) q[9];
u(pi/2,0.4052654523130834,-0.4052654523130834) q[9];
rzz(-pi/2) q[4],q[9];
rzz(-pi/2) q[4],q[8];
u(pi/2,-0.07351326809400116,0.07351326809400116) q[8];
u(pi/2,3.1660970762877936,-3.1660970762877936) q[8];
rzz(pi/2) q[4],q[8];
rzz(-pi/2) q[4],q[7];
u(pi/2,-0.319185813604723,0.319185813604723) q[7];
u(pi/2,3.0190705400997917,-3.0190705400997917) q[7];
rzz(pi/2) q[4],q[7];
rzz(-pi/2) q[4],q[6];
u(pi/2,3.239610344381794,-3.239610344381794) q[6];
u(pi/2,0.49071677249072554,-0.49071677249072554) q[6];
rzz(-pi/2) q[4],q[6];
rzz(-pi/2) q[4],q[5];
u(pi/2,-1.3747609452108935,1.3747609452108935) q[5];
u(pi/2,4.123026198571244,-4.123026198571244) q[5];
rzz(-pi/2) q[4],q[5];
rzz(-pi/2) q[3],q[16];
rzz(pi/2) q[3],q[16];
rzz(-pi/2) q[3],q[15];
rzz(pi/2) q[3],q[15];
rzz(-pi/2) q[3],q[14];
rzz(-pi/2) q[3],q[14];
rzz(pi/2) q[3],q[13];
rzz(-pi/2) q[3],q[13];
rzz(pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[11];
rzz(-pi/2) q[3],q[11];
rzz(pi/2) q[3],q[10];
u(pi/2,3.9332740022944206,-3.9332740022944206) q[10];
u(pi/2,0.8036194007882691,-0.8036194007882691) q[10];
rzz(pi/2) q[3],q[10];
rzz(-pi/2) q[3],q[9];
u(pi/2,0.4052654523130834,-0.4052654523130834) q[9];
u(pi/2,3.5713625286008766,-3.5713625286008766) q[9];
rzz(-pi/2) q[3],q[9];
rzz(pi/2) q[3],q[8];
u(pi/2,3.1660970762877936,-3.1660970762877936) q[8];
u(pi/2,0.07351326809400116,-0.07351326809400116) q[8];
rzz(-pi/2) q[3],q[8];
rzz(pi/2) q[3],q[7];
u(pi/2,3.0190705400997917,-3.0190705400997917) q[7];
u(pi/2,-0.02450442269800024,0.02450442269800024) q[7];
rzz(pi/2) q[3],q[7];
rzz(-pi/2) q[3],q[6];
u(pi/2,0.49071677249072554,-0.49071677249072554) q[6];
u(pi/2,3.8289731261952404,-3.8289731261952404) q[6];
rzz(pi/2) q[3],q[6];
rzz(-pi/2) q[3],q[5];
u(pi/2,0.9814335449814515,-0.9814335449814515) q[5];
u(pi/2,4.516353598800687,-4.516353598800687) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[4];
u(pi/2,-3*pi/8,3*pi/8) q[4];
u(pi/2,11*pi/8,-11*pi/8) q[4];
rzz(pi/2) q[3],q[4];
rzz(-pi/2) q[2],q[16];
rzz(-pi/2) q[2],q[16];
rzz(pi/2) q[2],q[15];
rzz(-pi/2) q[2],q[15];
rzz(-pi/2) q[2],q[14];
rzz(-pi/2) q[2],q[14];
rzz(-pi/2) q[2],q[13];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[11];
rzz(-pi/2) q[2],q[11];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[10];
rzz(-pi/2) q[2],q[9];
u(pi/2,3.5713625286008766,-3.5713625286008766) q[9];
u(pi/2,0.4417079270947246,-0.4417079270947246) q[9];
rzz(-pi/2) q[2],q[9];
rzz(-pi/2) q[2],q[8];
u(pi/2,0.07351326809400116,-0.07351326809400116) q[8];
u(pi/2,3.239610344381794,-3.239610344381794) q[8];
rzz(pi/2) q[2],q[8];
rzz(pi/2) q[2],q[7];
u(pi/2,-0.02450442269800024,0.02450442269800024) q[7];
u(pi/2,3.1660970762877936,-3.1660970762877936) q[7];
rzz(-pi/2) q[2],q[7];
rzz(pi/2) q[2],q[6];
u(pi/2,3.8289731261952404,-3.8289731261952404) q[6];
u(pi/2,pi/4,-pi/4) q[6];
rzz(-pi/2) q[2],q[6];
rzz(pi/2) q[2],q[5];
u(pi/2,4.516353598800687,-4.516353598800687) q[5];
u(pi/2,pi/2,-pi/2) q[5];
rzz(-pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[4];
u(pi/2,11*pi/8,-11*pi/8) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[2],q[4];
rzz(pi/2) q[2],q[3];
u(pi/2,pi/4,-pi/4) q[3];
u(pi/2,0,0) q[3];
rzz(pi/2) q[2],q[3];
rzz(-pi/2) q[1],q[16];
rzz(pi/2) q[1],q[16];
rzz(-pi/2) q[1],q[15];
rzz(-pi/2) q[1],q[15];
rzz(pi/2) q[1],q[14];
rzz(-pi/2) q[1],q[14];
rzz(-pi/2) q[1],q[13];
rzz(-pi/2) q[1],q[13];
rzz(-pi/2) q[1],q[12];
rzz(pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[10];
rzz(-pi/2) q[1],q[10];
rzz(pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[8];
u(pi/2,3.239610344381794,-3.239610344381794) q[8];
u(pi/2,0.11058406140636068,-0.11058406140636068) q[8];
rzz(pi/2) q[1],q[8];
rzz(-pi/2) q[1],q[7];
u(pi/2,0.02450442269800046,-0.02450442269800046) q[7];
u(pi/2,3.190601498985794,-3.190601498985794) q[7];
rzz(pi/2) q[1],q[7];
rzz(pi/2) q[1],q[6];
u(pi/2,pi/4,-pi/4) q[6];
u(pi/2,3.9759996623832423,-3.9759996623832423) q[6];
rzz(-pi/2) q[1],q[6];
rzz(-pi/2) q[1],q[5];
u(pi/2,-pi/2,pi/2) q[5];
u(pi/2,1.668814017586898,-1.668814017586898) q[5];
rzz(-pi/2) q[1],q[5];
rzz(pi/2) q[1],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-1.3747609452108935,1.3747609452108935) q[4];
rzz(-pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[3];
u(pi/2,pi,-pi) q[3];
u(pi/2,pi/8,-pi/8) q[3];
rzz(-pi/2) q[1],q[3];
rzz(pi/2) q[1],q[2];
u(pi/2,pi/2,-pi/2) q[2];
u(pi/2,pi/4,-pi/4) q[2];
rzz(pi/2) q[1],q[2];
rzz(-pi/2) q[0],q[16];
rzz(-pi/2) q[0],q[15];
rzz(-pi/2) q[0],q[14];
rzz(-pi/2) q[0],q[13];
rzz(pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[8];
rzz(pi/2) q[0],q[7];
u(pi/2,3.190601498985794,-3.190601498985794) q[7];
u(pi/2,0.06157521601035976,-0.06157521601035976) q[7];
rzz(-pi/2) q[0],q[7];
rzz(-pi/2) q[0],q[6];
rzz(-pi/2) q[0],q[5];
rzz(pi/2) q[0],q[4];
rzz(pi/2) q[0],q[3];
rzz(-pi/2) q[0],q[2];
rzz(-pi/2) q[0],q[1];
u(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[0],q[1];
u(pi/2,pi/4,-pi/4) q[2];
u(pi/2,-pi/2,pi/2) q[2];
rzz(pi/2) q[0],q[2];
u(pi/2,9*pi/8,-9*pi/8) q[3];
u(pi/2,pi/4,-pi/4) q[3];
rzz(-pi/2) q[0],q[3];
u(pi/2,1.7668317083788998,-1.7668317083788998) q[4];
u(pi/2,-3*pi/8,3*pi/8) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,1.668814017586898,-1.668814017586898) q[5];
u(pi/2,-1.3747609452108935,1.3747609452108935) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,3.9759996623832423,-3.9759996623832423) q[6];
u(pi/2,0.8834158541894497,-0.8834158541894497) q[6];
rzz(pi/2) q[0],q[6];
rzz(pi/2) q[0],q[7];
u(pi/2,3.2031678696001533,-3.2031678696001533) q[7];
u(pi/2,0.07351326809400116,-0.07351326809400116) q[7];
rzz(-pi/2) q[0],q[7];
u(pi/2,3.2521767149961534,-3.2521767149961534) q[8];
u(pi/2,0.12252211349000208,-0.12252211349000208) q[8];
rzz(-pi/2) q[0],q[8];
rzz(pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[13];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[15];
rzz(pi/2) q[0],q[16];
u(pi/2,pi,-pi) q[1];
rzz(-pi/2) q[1],q[2];
u(pi/2,pi/2,-pi/2) q[2];
u(pi/2,-pi/4,pi/4) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u(pi/2,5*pi/4,-5*pi/4) q[3];
u(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u(pi/2,-3*pi/8,3*pi/8) q[4];
u(pi/2,2.1601591086083416,-2.1601591086083416) q[4];
rzz(-pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u(pi/2,-1.3747609452108935,1.3747609452108935) q[5];
u(pi/2,1.8654777177016193,-1.8654777177016193) q[5];
rzz(pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[6];
u(pi/2,4.025008507779242,-4.025008507779242) q[6];
u(pi/2,0.9324246995854506,-0.9324246995854506) q[6];
rzz(pi/2) q[1],q[6];
rzz(-pi/2) q[1],q[7];
u(pi/2,0.07351326809400116,-0.07351326809400116) q[7];
u(pi/2,3.239610344381794,-3.239610344381794) q[7];
rzz(pi/2) q[1],q[7];
rzz(pi/2) q[1],q[8];
u(pi/2,3.2641147670797945,-3.2641147670797945) q[8];
u(pi/2,0.13508848410436114,-0.13508848410436114) q[8];
rzz(pi/2) q[1],q[8];
rzz(pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[9];
rzz(pi/2) q[1],q[10];
rzz(-pi/2) q[1],q[10];
rzz(pi/2) q[1],q[11];
rzz(pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[12];
rzz(pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[13];
rzz(pi/2) q[1],q[13];
rzz(pi/2) q[1],q[14];
rzz(-pi/2) q[1],q[14];
rzz(pi/2) q[1],q[15];
rzz(pi/2) q[1],q[15];
rzz(pi/2) q[1],q[16];
rzz(pi/2) q[1],q[16];
u(pi/2,5*pi/4,-5*pi/4) q[2];
u(pi/2,pi,-pi) q[2];
rzz(-pi/2) q[2],q[3];
u(pi/2,11*pi/8,-11*pi/8) q[3];
u(pi/2,5*pi/8,-5*pi/8) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u(pi/2,-0.9814335449814514,0.9814335449814514) q[4];
u(pi/2,2.5522298717763476,-2.5522298717763476) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u(pi/2,1.8654777177016193,-1.8654777177016193) q[5];
u(pi/2,-1.0800795543041708,1.0800795543041708) q[5];
rzz(-pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[6];
u(pi/2,4.074017353175243,-4.074017353175243) q[6];
u(pi/2,1.03107070890817,-1.03107070890817) q[6];
rzz(pi/2) q[2],q[6];
rzz(-pi/2) q[2],q[7];
u(pi/2,0.09801769079200162,-0.09801769079200162) q[7];
u(pi/2,3.288619189777795,-3.288619189777795) q[7];
rzz(pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[8];
u(pi/2,3.276681137694154,-3.276681137694154) q[8];
u(pi/2,0.15959290680236138,-0.15959290680236138) q[8];
rzz(-pi/2) q[2],q[8];
rzz(pi/2) q[2],q[9];
u(pi/2,0.45992916448554544,-0.45992916448554544) q[9];
u(pi/2,3.6140881886896983,-3.6140881886896983) q[9];
rzz(pi/2) q[2],q[9];
rzz(pi/2) q[2],q[10];
rzz(-pi/2) q[2],q[10];
rzz(-pi/2) q[2],q[11];
rzz(pi/2) q[2],q[11];
rzz(-pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[13];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[14];
rzz(-pi/2) q[2],q[14];
rzz(pi/2) q[2],q[15];
rzz(-pi/2) q[2],q[15];
rzz(-pi/2) q[2],q[16];
rzz(-pi/2) q[2],q[16];
u(pi/2,pi/8,-pi/8) q[3];
u(pi/2,0,0) q[3];
rzz(pi/2) q[3],q[4];
u(pi/2,2.5522298717763476,-2.5522298717763476) q[4];
u(pi/2,0.19603538158400302,-0.19603538158400302) q[4];
rzz(pi/2) q[3],q[4];
rzz(-pi/2) q[3],q[5];
u(pi/2,-1.0800795543041708,1.0800795543041708) q[5];
u(pi/2,2.454212180984346,-2.454212180984346) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u(pi/2,1.03107070890817,-1.03107070890817) q[6];
u(pi/2,4.368698744081967,-4.368698744081967) q[6];
rzz(pi/2) q[3],q[6];
rzz(-pi/2) q[3],q[7];
u(pi/2,0.1470265361880021,-0.1470265361880021) q[7];
u(pi/2,3.387265199100515,-3.387265199100515) q[7];
rzz(-pi/2) q[3],q[7];
rzz(pi/2) q[3],q[8];
u(pi/2,3.3011855603921543,-3.3011855603921543) q[8];
u(pi/2,0.2086017521983623,-0.2086017521983623) q[8];
rzz(-pi/2) q[3],q[8];
rzz(-pi/2) q[3],q[9];
u(pi/2,0.4724955350999047,-0.4724955350999047) q[9];
u(pi/2,3.638592611387698,-3.638592611387698) q[9];
rzz(pi/2) q[3],q[9];
rzz(-pi/2) q[3],q[10];
u(pi/2,0.8130441787490383,-0.8130441787490383) q[10];
u(pi/2,3.9665748844224726,-3.9665748844224726) q[10];
rzz(pi/2) q[3],q[10];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[12];
rzz(pi/2) q[3],q[13];
rzz(pi/2) q[3],q[13];
rzz(pi/2) q[3],q[14];
rzz(pi/2) q[3],q[14];
rzz(-pi/2) q[3],q[15];
rzz(pi/2) q[3],q[15];
rzz(pi/2) q[3],q[16];
rzz(pi/2) q[3],q[16];
u(pi/2,-1.3747609452108935,1.3747609452108935) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[4],q[5];
u(pi/2,2.454212180984346,-2.454212180984346) q[5];
u(pi/2,0.09801769079200162,-0.09801769079200162) q[5];
rzz(-pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u(pi/2,4.368698744081967,-4.368698744081967) q[6];
u(pi/2,1.6198051721908975,-1.6198051721908975) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u(pi/2,0.24567254551072204,-0.24567254551072204) q[7];
u(pi/2,3.5833005806845186,-3.5833005806845186) q[7];
rzz(pi/2) q[4],q[7];
rzz(-pi/2) q[4],q[8];
u(pi/2,0.2086017521983623,-0.2086017521983623) q[8];
u(pi/2,3.448212096580157,-3.448212096580157) q[8];
rzz(pi/2) q[4],q[8];
rzz(pi/2) q[4],q[9];
u(pi/2,3.638592611387698,-3.638592611387698) q[9];
u(pi/2,0.5460088031939061,-0.5460088031939061) q[9];
rzz(pi/2) q[4],q[9];
rzz(pi/2) q[4],q[10];
u(pi/2,3.9665748844224726,-3.9665748844224726) q[10];
u(pi/2,0.8501149720613981,-0.8501149720613981) q[10];
rzz(-pi/2) q[4],q[10];
rzz(-pi/2) q[4],q[11];
u(pi/2,-0.6641326869688822,0.6641326869688822) q[11];
u(pi/2,2.489398018704552,-2.489398018704552) q[11];
rzz(pi/2) q[4],q[11];
rzz(pi/2) q[4],q[12];
rzz(pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[13];
rzz(pi/2) q[4],q[13];
rzz(pi/2) q[4],q[14];
rzz(-pi/2) q[4],q[14];
rzz(pi/2) q[4],q[15];
rzz(-pi/2) q[4],q[15];
rzz(pi/2) q[4],q[16];
rzz(pi/2) q[4],q[16];
u(pi/2,-1.472778636002895,1.472778636002895) q[5];
u(pi/2,1.7668317083788998,-1.7668317083788998) q[5];
rzz(-pi/2) q[5],q[6];
u(pi/2,-1.5217874813988959,1.5217874813988959) q[6];
u(pi/2,2.4052033355883458,-2.4052033355883458) q[6];
rzz(pi/2) q[5],q[6];
rzz(-pi/2) q[5],q[7];
u(pi/2,0.4417079270947246,-0.4417079270947246) q[7];
u(pi/2,3.9759996623832423,-3.9759996623832423) q[7];
rzz(pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u(pi/2,3.448212096580157,-3.448212096580157) q[8];
u(pi/2,0.5032831431050848,-0.5032831431050848) q[8];
rzz(pi/2) q[5],q[8];
rzz(pi/2) q[5],q[9];
u(pi/2,0.5460088031939061,-0.5460088031939061) q[9];
u(pi/2,3.7856191475757006,-3.7856191475757006) q[9];
rzz(-pi/2) q[5],q[9];
rzz(pi/2) q[5],q[10];
u(pi/2,3.991707625651191,-3.991707625651191) q[10];
u(pi/2,0.899123817457399,-0.899123817457399) q[10];
rzz(pi/2) q[5],q[10];
rzz(-pi/2) q[5],q[11];
u(pi/2,-0.6521946348852411,0.6521946348852411) q[11];
u(pi/2,2.5139024414025526,-2.5139024414025526) q[11];
rzz(pi/2) q[5],q[11];
rzz(-pi/2) q[5],q[12];
u(pi/2,0.9933715970650927,-0.9933715970650927) q[12];
u(pi/2,4.1469023027385274,-4.1469023027385274) q[12];
rzz(pi/2) q[5],q[12];
rzz(-pi/2) q[5],q[13];
rzz(-pi/2) q[5],q[13];
rzz(pi/2) q[5],q[14];
rzz(-pi/2) q[5],q[14];
rzz(pi/2) q[5],q[15];
rzz(pi/2) q[5],q[15];
rzz(pi/2) q[5],q[16];
rzz(pi/2) q[5],q[16];
u(pi/2,3.9759996623832423,-3.9759996623832423) q[6];
u(pi/2,0.8834158541894497,-0.8834158541894497) q[6];
rzz(-pi/2) q[6],q[7];
u(pi/2,0.8344070087934488,-0.8344070087934488) q[7];
u(pi/2,-1.5217874813988959,1.5217874813988959) q[7];
rzz(pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u(pi/2,0.5032831431050848,-0.5032831431050848) q[8];
u(pi/2,4.037574878393602,-4.037574878393602) q[8];
rzz(pi/2) q[6],q[8];
rzz(-pi/2) q[6],q[9];
u(pi/2,3.7856191475757006,-3.7856191475757006) q[9];
u(pi/2,0.8406901941006284,-0.8406901941006284) q[9];
rzz(-pi/2) q[6],q[9];
rzz(pi/2) q[6],q[10];
u(pi/2,0.899123817457399,-0.899123817457399) q[10];
u(pi/2,4.138734161839193,-4.138734161839193) q[10];
rzz(pi/2) q[6],q[10];
rzz(pi/2) q[6],q[11];
u(pi/2,2.5139024414025526,-2.5139024414025526) q[11];
u(pi/2,-0.578053048260522,0.578053048260522) q[11];
rzz(pi/2) q[6],q[11];
rzz(-pi/2) q[6],q[12];
u(pi/2,1.0053096491487334,-1.0053096491487334) q[12];
u(pi/2,4.171406725436528,-4.171406725436528) q[12];
rzz(pi/2) q[6],q[12];
rzz(pi/2) q[6],q[13];
u(pi/2,pi/5,-pi/5) q[13];
u(pi/2,3.782477554922111,-3.782477554922111) q[13];
rzz(pi/2) q[6],q[13];
rzz(pi/2) q[6],q[14];
rzz(-pi/2) q[6],q[14];
rzz(pi/2) q[6],q[15];
rzz(pi/2) q[6],q[15];
rzz(pi/2) q[6],q[16];
rzz(pi/2) q[6],q[16];
u(pi/2,0.0490088453960007,-0.0490088453960007) q[7];
u(pi/2,3.2151059216837945,-3.2151059216837945) q[7];
rzz(-pi/2) q[7],q[8];
u(pi/2,0.895982224803809,-0.895982224803809) q[8];
u(pi/2,-1.4602122653885359,1.4602122653885359) q[8];
rzz(-pi/2) q[7],q[8];
rzz(pi/2) q[7],q[9];
u(pi/2,3.9822828476904215,-3.9822828476904215) q[9];
u(pi/2,1.2333892757993525,-1.2333892757993525) q[9];
rzz(pi/2) q[7],q[9];
rzz(pi/2) q[7],q[10];
u(pi/2,4.138734161839193,-4.138734161839193) q[10];
u(pi/2,1.1931768898334036,-1.1931768898334036) q[10];
rzz(-pi/2) q[7],q[10];
rzz(pi/2) q[7],q[11];
u(pi/2,-0.578053048260522,0.578053048260522) q[11];
u(pi/2,2.6615572961212726,-2.6615572961212726) q[11];
rzz(-pi/2) q[7],q[11];
rzz(pi/2) q[7],q[12];
u(pi/2,4.171406725436528,-4.171406725436528) q[12];
u(pi/2,1.079451235773453,-1.079451235773453) q[12];
rzz(pi/2) q[7],q[12];
rzz(-pi/2) q[7],q[13];
u(pi/2,0.6408849013323175,-0.6408849013323175) q[13];
u(pi/2,3.8069819776201115,-3.8069819776201115) q[13];
rzz(-pi/2) q[7],q[13];
rzz(pi/2) q[7],q[14];
u(pi/2,3.3891501546926683,-3.3891501546926683) q[14];
u(pi/2,0.2601238717172347,-0.2601238717172347) q[14];
rzz(pi/2) q[7],q[14];
rzz(pi/2) q[7],q[15];
rzz(-pi/2) q[7],q[15];
rzz(pi/2) q[7],q[16];
rzz(-pi/2) q[7],q[16];
u(pi/2,3.2521767149961534,-3.2521767149961534) q[8];
u(pi/2,0.12252211349000208,-0.12252211349000208) q[8];
rzz(pi/2) q[8],q[9];
u(pi/2,1.2333892757993525,-1.2333892757993525) q[9];
u(pi/2,-1.122805214392992,1.122805214392992) q[9];
rzz(pi/2) q[8],q[9];
rzz(-pi/2) q[8],q[10];
u(pi/2,1.1931768898334036,-1.1931768898334036) q[10];
u(pi/2,-1.5557166820576656,1.5557166820576656) q[10];
rzz(pi/2) q[8],q[10];
rzz(pi/2) q[8],q[11];
u(pi/2,-0.48003535746852033,0.48003535746852033) q[11];
u(pi/2,2.8575926777052754,-2.8575926777052754) q[11];
rzz(pi/2) q[8],q[11];
rzz(pi/2) q[8],q[12];
u(pi/2,1.079451235773453,-1.079451235773453) q[12];
u(pi/2,4.319061580155248,-4.319061580155248) q[12];
rzz(-pi/2) q[8],q[12];
rzz(pi/2) q[8],q[13];
u(pi/2,0.6653893240303184,-0.6653893240303184) q[13];
u(pi/2,3.8559908230161124,-3.8559908230161124) q[13];
rzz(pi/2) q[8],q[13];
rzz(pi/2) q[8],q[14];
u(pi/2,0.2601238717172347,-0.2601238717172347) q[14];
u(pi/2,3.426220948005028,-3.426220948005028) q[14];
rzz(pi/2) q[8],q[14];
rzz(-pi/2) q[8],q[15];
u(pi/2,1.7285042780051043,-1.7285042780051043) q[15];
u(pi/2,-1.4011503235010476,1.4011503235010476) q[15];
rzz(pi/2) q[8],q[15];
rzz(-pi/2) q[8],q[16];
rzz(pi/2) q[8],q[16];
u(pi,3.589583765991698,-3.589583765991698) q[9];
rzz(pi/2) q[9],q[10];
u(pi/2,-1.5557166820576656,1.5557166820576656) q[10];
u(pi/2,2.3712741349295756,-2.3712741349295756) q[10];
rzz(-pi/2) q[9],q[10];
rzz(pi/2) q[9],q[11];
u(pi/2,2.8575926777052754,-2.8575926777052754) q[11];
u(pi/2,0.10869910581420683,-0.10869910581420683) q[11];
rzz(pi/2) q[9],q[11];
rzz(pi/2) q[9],q[12];
u(pi/2,1.1774689265654543,-1.1774689265654543) q[12];
u(pi/2,4.515096961739251,-4.515096961739251) q[12];
rzz(pi/2) q[9],q[12];
rzz(-pi/2) q[9],q[13];
u(pi/2,0.7143981694263193,-0.7143981694263193) q[13];
u(pi/2,3.9540085138081134,-3.9540085138081134) q[13];
rzz(pi/2) q[9],q[13];
rzz(pi/2) q[9],q[14];
u(pi/2,3.426220948005028,-3.426220948005028) q[14];
u(pi/2,0.3336371398112359,-0.3336371398112359) q[14];
rzz(pi/2) q[9],q[14];
rzz(pi/2) q[9],q[15];
u(pi/2,-1.4011503235010476,1.4011503235010476) q[15];
u(pi/2,1.764946752786746,-1.764946752786746) q[15];
rzz(-pi/2) q[9],q[15];
rzz(pi/2) q[9],q[16];
u(pi/2,1.6298582686823844,-1.6298582686823844) q[16];
u(pi/2,-1.4997963328237673,1.4997963328237673) q[16];
rzz(pi/2) q[9],q[16];
u(pi,3.942070461724472,-3.942070461724472) q[10];
rzz(pi/2) q[10],q[11];
u(pi/2,0.10869910581420683,-0.10869910581420683) q[11];
u(pi/2,4.035689922801448,-4.035689922801448) q[11];
rzz(pi/2) q[10],q[11];
rzz(-pi/2) q[10],q[12];
u(pi/2,1.3735043081494576,-1.3735043081494576) q[12];
u(pi/2,-1.3753892637416114,1.3753892637416114) q[12];
rzz(pi/2) q[10],q[12];
rzz(pi/2) q[10],q[13];
u(pi/2,3.9540085138081134,-3.9540085138081134) q[13];
u(pi/2,1.0090795603330416,-1.0090795603330416) q[13];
rzz(-pi/2) q[10],q[13];
rzz(pi/2) q[10],q[14];
u(pi/2,0.3336371398112359,-0.3336371398112359) q[14];
u(pi/2,3.573247484193031,-3.573247484193031) q[14];
rzz(-pi/2) q[10],q[14];
rzz(pi/2) q[10],q[15];
u(pi/2,-1.3766459008030474,1.3766459008030474) q[15];
u(pi/2,1.8139555981827464,-1.8139555981827464) q[15];
rzz(pi/2) q[10],q[15];
rzz(-pi/2) q[10],q[16];
u(pi/2,1.6417963207660256,-1.6417963207660256) q[16];
u(pi/2,-1.4752919101257669,1.4752919101257669) q[16];
rzz(pi/2) q[10],q[16];
u(pi,2.4648935960065517,-2.4648935960065517) q[11];
rzz(-pi/2) q[11],q[12];
u(pi/2,1.7662033898481817,-1.7662033898481817) q[12];
u(pi/2,-0.5899911003441632,0.5899911003441632) q[12];
rzz(pi/2) q[11],q[12];
rzz(pi/2) q[11],q[13];
u(pi/2,4.150672213922834,-4.150672213922834) q[13];
u(pi/2,1.4017786420317657,-1.4017786420317657) q[13];
rzz(-pi/2) q[11],q[13];
rzz(pi/2) q[11],q[14];
u(pi/2,0.4316548306032373,-0.4316548306032373) q[14];
u(pi/2,3.7692828657770336,-3.7692828657770336) q[14];
rzz(-pi/2) q[11],q[14];
rzz(pi/2) q[11],q[15];
u(pi/2,1.8139555981827464,-1.8139555981827464) q[15];
u(pi/2,-1.2289910460843272,1.2289910460843272) q[15];
rzz(pi/2) q[11],q[15];
rzz(pi/2) q[11],q[16];
u(pi/2,-1.4752919101257669,1.4752919101257669) q[16];
u(pi/2,1.7153095888600274,-1.7153095888600274) q[16];
rzz(pi/2) q[11],q[16];
u(pi,4.122397880040527,-4.122397880040527) q[12];
rzz(-pi/2) q[12],q[13];
u(pi/2,1.4017786420317657,-1.4017786420317657) q[13];
u(pi/2,-0.9544158481605791,0.9544158481605791) q[13];
rzz(pi/2) q[12],q[13];
rzz(-pi/2) q[12],q[14];
u(pi/2,3.7692828657770336,-3.7692828657770336) q[14];
u(pi/2,1.0203892938859647,-1.0203892938859647) q[14];
rzz(-pi/2) q[12],q[14];
rzz(pi/2) q[12],q[15];
u(pi/2,-1.2289910460843272,1.2289910460843272) q[15];
u(pi/2,2.108636989089469,-2.108636989089469) q[15];
rzz(-pi/2) q[12],q[15];
rzz(pi/2) q[12],q[16];
u(pi/2,1.7153095888600274,-1.7153095888600274) q[16];
u(pi/2,-1.3276370554070467,1.3276370554070467) q[16];
rzz(pi/2) q[12],q[16];
rzz(pi/2) q[13],q[14];
u(pi/2,4.161981947475758,-4.161981947475758) q[14];
u(pi/2,1.805787457283413,-1.805787457283413) q[14];
rzz(pi/2) q[13],q[14];
rzz(-pi/2) q[13],q[15];
u(pi/2,2.108636989089469,-2.108636989089469) q[15];
u(pi/2,-0.6402565828015998,0.6402565828015998) q[15];
rzz(pi/2) q[13],q[15];
rzz(pi/2) q[13],q[16];
u(pi/2,-1.3276370554070467,1.3276370554070467) q[16];
u(pi/2,2.0099909797667492,-2.0099909797667492) q[16];
rzz(pi/2) q[13],q[16];
rzz(pi/2) q[14],q[15];
u(pi/2,-0.6402565828015998,0.6402565828015998) q[15];
u(pi/2,3.2867342341856416,-3.2867342341856416) q[15];
rzz(-pi/2) q[14],q[15];
rzz(pi/2) q[14],q[16];
u(pi/2,2.0099909797667492,-2.0099909797667492) q[16];
u(pi/2,-0.7389025921243194,0.7389025921243194) q[16];
rzz(pi/2) q[14],q[16];
rzz(-pi/2) q[15],q[16];
u(pi/2,2.4026900614654734,-2.4026900614654734) q[16];
u(pi/2,0.04649557127312898,-0.04649557127312898) q[16];
rzz(pi/2) q[15],q[16];
u(pi,pi/2,-pi/2) q[0];
u(pi,pi/2,-pi/2) q[1];
u(pi,pi,-pi) q[2];
u(pi,3*pi/4,-3*pi/4) q[3];
u(pi,11*pi/8,-11*pi/8) q[4];
u(pi,5*pi/4,-5*pi/4) q[5];
u(pi,-3*pi/8,3*pi/8) q[6];
u(pi,2.6508758810990676,-2.6508758810990676) q[7];
u(pi,3.239610344381794,-3.239610344381794) q[8];
u(pi,-1.141026451783813,1.141026451783813) q[9];
u(pi,1.5770795121020762,-1.5770795121020762) q[10];
u(pi,2.8531944479902496,-2.8531944479902496) q[11];
u(pi,4.31843326162453,-4.31843326162453) q[12];
u(pi,-0.8576547944300135,0.8576547944300135) q[13];
u(pi,2.6401944660768626,-2.6401944660768626) q[14];
u(pi,1.3477432483900214,-1.3477432483900214) q[15];
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
