OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u(pi/2,-0.03958406743523146,0.03958406743523146) q[15];
u(pi/2,-0.8249822308326796,0.8249822308326796) q[15];
rzz(-pi/2) q[14],q[15];
rzz(pi/2) q[13],q[15];
u(pi/2,2.3166104227571136,-2.3166104227571136) q[15];
u(pi/2,-0.4322831491339556,0.4322831491339556) q[15];
rzz(-pi/2) q[13],q[15];
rzz(pi/2) q[13],q[14];
u(pi/2,0.9229999216246809,-0.9229999216246809) q[14];
u(pi/2,0.13760175822723286,-0.13760175822723286) q[14];
rzz(pi/2) q[13],q[14];
rzz(-pi/2) q[12],q[15];
u(pi/2,2.709309504455838,-2.709309504455838) q[15];
u(pi/2,-0.23561944901923448,0.23561944901923448) q[15];
rzz(pi/2) q[12],q[15];
rzz(-pi/2) q[12],q[14];
u(pi/2,3.279194411817026,-3.279194411817026) q[14];
u(pi/2,0.5303008399259568,-0.5303008399259568) q[14];
rzz(-pi/2) q[12],q[14];
rzz(-pi/2) q[12],q[13];
u(pi/2,2.4925396113581426,-2.4925396113581426) q[13];
u(pi/2,1.707141447960694,-1.707141447960694) q[13];
rzz(-pi/2) q[12],q[13];
rzz(pi/2) q[11],q[15];
u(pi/2,-0.23561944901923448,0.23561944901923448) q[15];
u(pi/2,3.00399089536256,-3.00399089536256) q[15];
rzz(-pi/2) q[11],q[15];
rzz(pi/2) q[11],q[14];
u(pi/2,3.671893493515751,-3.671893493515751) q[14];
u(pi/2,0.72633622150996,-0.72633622150996) q[14];
rzz(pi/2) q[11],q[14];
rzz(-pi/2) q[11],q[13];
u(pi/2,-1.4344512056290994,1.4344512056290994) q[13];
u(pi/2,2.099840529659418,-2.099840529659418) q[13];
rzz(pi/2) q[11],q[13];
rzz(pi/2) q[11],q[12];
u(pi/2,-0.2607521902479528,0.2607521902479528) q[12];
u(pi/2,-1.046150353645401,1.046150353645401) q[12];
rzz(-pi/2) q[11],q[12];
rzz(-pi/2) q[10],q[15];
u(pi/2,-0.13760175822723286,0.13760175822723286) q[15];
u(pi/2,3.052999740758561,-3.052999740758561) q[15];
rzz(-pi/2) q[10],q[15];
rzz(pi/2) q[10],q[14];
u(pi/2,3.8679288750997536,-3.8679288750997536) q[14];
u(pi/2,0.8243539123019614,-0.8243539123019614) q[14];
rzz(-pi/2) q[10],q[14];
rzz(-pi/2) q[10],q[13];
u(pi/2,2.099840529659418,-2.099840529659418) q[13];
u(pi/2,-0.8450884238156543,0.8450884238156543) q[13];
rzz(-pi/2) q[10],q[13];
rzz(-pi/2) q[10],q[12];
u(pi/2,-1.046150353645401,1.046150353645401) q[12];
u(pi/2,2.488141381643116,-2.488141381643116) q[12];
rzz(pi/2) q[10],q[12];
rzz(pi/2) q[10],q[11];
u(pi/2,3.6876014567836988,-3.6876014567836988) q[11];
u(pi/2,2.9022032933862505,-2.9022032933862505) q[11];
rzz(-pi/2) q[10],q[11];
rzz(-pi/2) q[9],q[15];
u(pi/2,3.052999740758561,-3.052999740758561) q[15];
u(pi/2,-0.0640884901332317,0.0640884901332317) q[15];
rzz(-pi/2) q[9],q[15];
rzz(pi/2) q[9],q[14];
u(pi/2,3.9659465658917545,-3.9659465658917545) q[14];
u(pi/2,0.8733627576979623,-0.8733627576979623) q[14];
rzz(pi/2) q[9],q[14];
rzz(pi/2) q[9],q[13];
u(pi/2,2.296504229774139,-2.296504229774139) q[13];
u(pi/2,-0.7470707330236528,0.7470707330236528) q[13];
rzz(pi/2) q[9],q[13];
rzz(pi/2) q[9],q[12];
u(pi/2,2.488141381643116,-2.488141381643116) q[12];
u(pi/2,-0.45741589036267394,0.45741589036267394) q[12];
rzz(pi/2) q[9],q[12];
rzz(-pi/2) q[9],q[11];
u(pi/2,-0.23938936020354218,0.23938936020354218) q[11];
u(pi/2,3.294902375084975,-3.294902375084975) q[11];
rzz(-pi/2) q[9],q[11];
rzz(pi/2) q[9],q[10];
u(pi/2,2.1966015833899837,-2.1966015833899837) q[10];
u(pi/2,1.4112034199925354,-1.4112034199925354) q[10];
rzz(pi/2) q[9],q[10];
rzz(pi/2) q[8],q[15];
u(pi/2,3.0775041634565614,-3.0775041634565614) q[15];
u(pi/2,-0.051522119518872644,0.051522119518872644) q[15];
rzz(-pi/2) q[8],q[15];
rzz(-pi/2) q[8],q[14];
u(pi/2,4.0149554112877555,-4.0149554112877555) q[14];
u(pi/2,0.8978671803959632,-0.8978671803959632) q[14];
rzz(-pi/2) q[8],q[14];
rzz(pi/2) q[8],q[13];
u(pi/2,-0.7470707330236528,0.7470707330236528) q[13];
u(pi/2,2.443530765962141,-2.443530765962141) q[13];
rzz(pi/2) q[8],q[13];
rzz(-pi/2) q[8],q[12];
u(pi/2,2.6841767632271196,-2.6841767632271196) q[12];
u(pi/2,-0.35876988103995444,0.35876988103995444) q[12];
rzz(-pi/2) q[8],q[12];
rzz(pi/2) q[8],q[11];
u(pi/2,0.15330972149518174,-0.15330972149518174) q[11];
u(pi/2,3.491566075199696,-3.491566075199696) q[11];
rzz(-pi/2) q[8],q[11];
rzz(pi/2) q[8],q[10];
u(pi/2,4.5527960735823285,-4.5527960735823285) q[10];
u(pi/2,1.8039025016912595,-1.8039025016912595) q[10];
rzz(-pi/2) q[8],q[10];
rzz(pi/2) q[8],q[9];
u(pi/2,-0.9569291222834511,0.9569291222834511) q[9];
u(pi/2,4.540858021498687,-4.540858021498687) q[9];
rzz(-pi/2) q[8],q[9];
rzz(pi/2) q[7],q[15];
rzz(pi/2) q[7],q[15];
rzz(pi/2) q[7],q[14];
u(pi/2,0.8978671803959632,-0.8978671803959632) q[14];
u(pi/2,4.052026204600115,-4.052026204600115) q[14];
rzz(pi/2) q[7],q[14];
rzz(-pi/2) q[7],q[13];
u(pi/2,2.443530765962141,-2.443530765962141) q[13];
u(pi/2,-0.6735574649296516,0.6735574649296516) q[13];
rzz(-pi/2) q[7],q[13];
rzz(pi/2) q[7],q[12];
u(pi/2,-0.35876988103995444,0.35876988103995444) q[12];
u(pi/2,2.8318316179458396,-2.8318316179458396) q[12];
rzz(pi/2) q[7],q[12];
rzz(-pi/2) q[7],q[11];
u(pi/2,0.3499734216099031,-0.3499734216099031) q[11];
u(pi/2,3.589583765991698,-3.589583765991698) q[11];
rzz(pi/2) q[7],q[11];
rzz(-pi/2) q[7],q[10];
u(pi/2,-1.337690151898534,1.337690151898534) q[10];
u(pi/2,2.0005662018059804,-2.0005662018059804) q[10];
rzz(pi/2) q[7],q[10];
rzz(-pi/2) q[7],q[9];
u(pi/2,4.540858021498687,-4.540858021498687) q[9];
u(pi/2,1.791964449607618,-1.791964449607618) q[9];
rzz(-pi/2) q[7],q[9];
rzz(pi/2) q[7],q[8];
u(pi/2,2.5032210263803467,-2.5032210263803467) q[8];
u(pi/2,1.7178228629828989,-1.7178228629828989) q[8];
rzz(-pi/2) q[7],q[8];
rzz(pi/2) q[6],q[15];
rzz(-pi/2) q[6],q[15];
rzz(pi/2) q[6],q[14];
rzz(-pi/2) q[6],q[14];
rzz(-pi/2) q[6],q[13];
u(pi/2,2.468035188660142,-2.468035188660142) q[13];
u(pi/2,-0.6609910943152923,0.6609910943152923) q[13];
rzz(pi/2) q[6],q[13];
rzz(-pi/2) q[6],q[12];
u(pi/2,2.8318316179458396,-2.8318316179458396) q[12];
u(pi/2,-0.2852566129459533,0.2852566129459533) q[12];
rzz(-pi/2) q[6],q[12];
rzz(-pi/2) q[6],q[11];
u(pi/2,3.589583765991698,-3.589583765991698) q[11];
u(pi/2,0.4969999577979052,-0.4969999577979052) q[11];
rzz(-pi/2) q[6],q[11];
rzz(pi/2) q[6],q[10];
u(pi/2,-1.141026451783813,1.141026451783813) q[10];
u(pi/2,2.0985838925979814,-2.0985838925979814) q[10];
rzz(-pi/2) q[6],q[10];
rzz(-pi/2) q[6],q[9];
u(pi/2,-1.3496282039821752,1.3496282039821752) q[9];
u(pi/2,1.9879998311916212,-1.9879998311916212) q[9];
rzz(-pi/2) q[6],q[9];
rzz(-pi/2) q[6],q[8];
u(pi/2,1.7178228629828989,-1.7178228629828989) q[8];
u(pi/2,-1.03107070890817,1.03107070890817) q[8];
rzz(pi/2) q[6],q[8];
rzz(pi/2) q[6],q[7];
u(pi/2,2.454212180984346,-2.454212180984346) q[7];
u(pi/2,1.668814017586898,-1.668814017586898) q[7];
rzz(-pi/2) q[6],q[7];
rzz(-pi/2) q[5],q[15];
rzz(-pi/2) q[5],q[15];
rzz(pi/2) q[5],q[14];
rzz(-pi/2) q[5],q[14];
rzz(pi/2) q[5],q[13];
rzz(-pi/2) q[5],q[13];
rzz(-pi/2) q[5],q[12];
u(pi/2,2.85633604064384,-2.85633604064384) q[12];
u(pi/2,-0.2733185608623121,0.2733185608623121) q[12];
rzz(pi/2) q[5],q[12];
rzz(-pi/2) q[5],q[11];
u(pi/2,3.638592611387698,-3.638592611387698) q[11];
u(pi/2,0.5215043804959056,-0.5215043804959056) q[11];
rzz(-pi/2) q[5],q[11];
rzz(-pi/2) q[5],q[10];
u(pi/2,-1.0430087609918113,1.0430087609918113) q[10];
u(pi/2,2.1475927379939823,-2.1475927379939823) q[10];
rzz(-pi/2) q[5],q[10];
rzz(pi/2) q[5],q[9];
u(pi/2,1.9879998311916212,-1.9879998311916212) q[9];
u(pi/2,-1.0555751316061706,1.0555751316061706) q[9];
rzz(-pi/2) q[5],q[9];
rzz(-pi/2) q[5],q[8];
u(pi/2,-1.03107070890817,1.03107070890817) q[8];
u(pi/2,2.307185644796344,-2.307185644796344) q[8];
rzz(-pi/2) q[5],q[8];
rzz(-pi/2) q[5],q[7];
u(pi/2,1.668814017586898,-1.668814017586898) q[7];
u(pi/2,-1.0800795543041708,1.0800795543041708) q[7];
rzz(pi/2) q[5],q[7];
rzz(pi/2) q[5],q[6];
u(pi/2,3.7309554354032386,-3.7309554354032386) q[6];
u(pi/2,2.9455572720057903,-2.9455572720057903) q[6];
rzz(pi/2) q[5],q[6];
rzz(-pi/2) q[4],q[15];
rzz(-pi/2) q[4],q[15];
rzz(pi/2) q[4],q[14];
rzz(-pi/2) q[4],q[14];
rzz(-pi/2) q[4],q[13];
rzz(pi/2) q[4],q[13];
rzz(-pi/2) q[4],q[12];
rzz(pi/2) q[4],q[12];
rzz(pi/2) q[4],q[11];
u(pi/2,0.5215043804959056,-0.5215043804959056) q[11];
u(pi/2,3.6756634047000576,-3.6756634047000576) q[11];
rzz(-pi/2) q[4],q[11];
rzz(-pi/2) q[4],q[10];
u(pi/2,-0.9939999155958105,0.9939999155958105) q[10];
u(pi/2,2.172097160691983,-2.172097160691983) q[10];
rzz(pi/2) q[4],q[10];
rzz(pi/2) q[4],q[9];
u(pi/2,-1.0555751316061706,1.0555751316061706) q[9];
u(pi/2,2.1350263673796235,-2.1350263673796235) q[9];
rzz(-pi/2) q[4],q[9];
rzz(-pi/2) q[4],q[8];
u(pi/2,-0.8344070087934491,0.8344070087934491) q[8];
u(pi/2,2.4052033355883458,-2.4052033355883458) q[8];
rzz(-pi/2) q[4],q[8];
rzz(-pi/2) q[4],q[7];
u(pi/2,-1.0800795543041708,1.0800795543041708) q[7];
u(pi/2,2.2581767994003434,-2.2581767994003434) q[7];
rzz(pi/2) q[4],q[7];
rzz(-pi/2) q[4],q[6];
u(pi/2,-0.19603538158400302,0.19603538158400302) q[6];
u(pi/2,3.337628035173796,-3.337628035173796) q[6];
rzz(-pi/2) q[4],q[6];
rzz(pi/2) q[4],q[5];
u(pi/2,9*pi/8,-9*pi/8) q[5];
u(pi/2,7*pi/8,-7*pi/8) q[5];
rzz(-pi/2) q[4],q[5];
rzz(pi/2) q[3],q[15];
rzz(-pi/2) q[3],q[15];
rzz(pi/2) q[3],q[14];
rzz(-pi/2) q[3],q[14];
rzz(-pi/2) q[3],q[13];
rzz(pi/2) q[3],q[13];
rzz(-pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[12];
rzz(pi/2) q[3],q[11];
rzz(-pi/2) q[3],q[11];
rzz(pi/2) q[3],q[10];
u(pi/2,2.172097160691983,-2.172097160691983) q[10];
u(pi/2,-0.9569291222834511,0.9569291222834511) q[10];
rzz(-pi/2) q[3],q[10];
rzz(-pi/2) q[3],q[9];
u(pi/2,2.1350263673796235,-2.1350263673796235) q[9];
u(pi/2,-0.9814335449814514,0.9814335449814514) q[9];
rzz(-pi/2) q[3],q[9];
rzz(-pi/2) q[3],q[8];
u(pi/2,2.4052033355883458,-2.4052033355883458) q[8];
u(pi/2,-0.6873804726054468,0.6873804726054468) q[8];
rzz(pi/2) q[3],q[8];
rzz(-pi/2) q[3],q[7];
u(pi/2,-0.8834158541894498,0.8834158541894498) q[7];
u(pi/2,3*pi/4,-3*pi/4) q[7];
rzz(-pi/2) q[3],q[7];
rzz(-pi/2) q[3],q[6];
u(pi/2,3.337628035173796,-3.337628035173796) q[6];
u(pi/2,pi/8,-pi/8) q[6];
rzz(-pi/2) q[3],q[6];
rzz(pi/2) q[3],q[5];
u(pi/2,7*pi/8,-7*pi/8) q[5];
u(pi/2,0,0) q[5];
rzz(-pi/2) q[3],q[5];
rzz(pi/2) q[3],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(-pi/2) q[3],q[4];
rzz(-pi/2) q[2],q[15];
rzz(pi/2) q[2],q[15];
rzz(pi/2) q[2],q[14];
rzz(-pi/2) q[2],q[14];
rzz(-pi/2) q[2],q[13];
rzz(-pi/2) q[2],q[13];
rzz(pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[11];
rzz(-pi/2) q[2],q[11];
rzz(-pi/2) q[2],q[10];
rzz(pi/2) q[2],q[10];
rzz(-pi/2) q[2],q[9];
u(pi/2,-0.9814335449814514,0.9814335449814514) q[9];
u(pi/2,2.172097160691983,-2.172097160691983) q[9];
rzz(-pi/2) q[2],q[9];
rzz(pi/2) q[2],q[8];
u(pi/2,-0.6873804726054468,0.6873804726054468) q[8];
u(pi/2,2.4787166036823463,-2.4787166036823463) q[8];
rzz(pi/2) q[2],q[8];
rzz(pi/2) q[2],q[7];
u(pi/2,-pi/4,pi/4) q[7];
u(pi/2,2.4052033355883458,-2.4052033355883458) q[7];
rzz(pi/2) q[2],q[7];
rzz(pi/2) q[2],q[6];
u(pi/2,9*pi/8,-9*pi/8) q[6];
u(pi/2,0.49071677249072554,-0.49071677249072554) q[6];
rzz(-pi/2) q[2],q[6];
rzz(-pi/2) q[2],q[5];
u(pi/2,0,0) q[5];
u(pi/2,3.337628035173796,-3.337628035173796) q[5];
rzz(pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-3*pi/8,3*pi/8) q[4];
rzz(-pi/2) q[2],q[4];
rzz(pi/2) q[2],q[3];
u(pi/2,3*pi/4,-3*pi/4) q[3];
u(pi/2,pi/2,-pi/2) q[3];
rzz(-pi/2) q[2],q[3];
rzz(pi/2) q[1],q[15];
rzz(pi/2) q[1],q[15];
rzz(-pi/2) q[1],q[14];
rzz(pi/2) q[1],q[14];
rzz(-pi/2) q[1],q[13];
rzz(pi/2) q[1],q[13];
rzz(-pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[11];
rzz(pi/2) q[1],q[10];
rzz(-pi/2) q[1],q[10];
rzz(-pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[8];
u(pi/2,2.4787166036823463,-2.4787166036823463) q[8];
u(pi/2,-0.6503096792930873,0.6503096792930873) q[8];
rzz(pi/2) q[1],q[8];
rzz(-pi/2) q[1],q[7];
u(pi/2,2.4052033355883458,-2.4052033355883458) q[7];
u(pi/2,-0.7118848953034472,0.7118848953034472) q[7];
rzz(-pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[6];
u(pi/2,3.6323094260805187,-3.6323094260805187) q[6];
u(pi/2,0.5397256178867265,-0.5397256178867265) q[6];
rzz(-pi/2) q[1],q[6];
rzz(pi/2) q[1],q[5];
u(pi/2,0.19603538158400302,-0.19603538158400302) q[5];
u(pi/2,3.436274044496516,-3.436274044496516) q[5];
rzz(pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[4];
u(pi/2,5*pi/8,-5*pi/8) q[4];
u(pi/2,-0.9814335449814514,0.9814335449814514) q[4];
rzz(-pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[3];
u(pi/2,pi/2,-pi/2) q[3];
u(pi/2,-3*pi/8,3*pi/8) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[2];
u(pi/2,0,0) q[2];
u(pi/2,-pi/4,pi/4) q[2];
rzz(-pi/2) q[1],q[2];
rzz(-pi/2) q[0],q[15];
rzz(-pi/2) q[0],q[14];
rzz(pi/2) q[0],q[13];
rzz(-pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[9];
rzz(pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[7];
rzz(-pi/2) q[0],q[6];
rzz(-pi/2) q[0],q[5];
rzz(-pi/2) q[0],q[4];
rzz(pi/2) q[0],q[3];
rzz(pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[0],q[1];
u(pi/2,3*pi/4,-3*pi/4) q[2];
u(pi/2,0,0) q[2];
rzz(pi/2) q[0],q[2];
u(pi/2,5*pi/8,-5*pi/8) q[3];
u(pi/2,-pi/4,pi/4) q[3];
rzz(-pi/2) q[0],q[3];
u(pi/2,2.1601591086083416,-2.1601591086083416) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3.436274044496516,-3.436274044496516) q[5];
u(pi/2,pi/8,-pi/8) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,3.6813182714765196,-3.6813182714765196) q[6];
u(pi/2,0.589362781813445,-0.589362781813445) q[6];
rzz(pi/2) q[0],q[6];
u(pi/2,2.429707758286346,-2.429707758286346) q[7];
u(pi/2,-0.6873804726054468,0.6873804726054468) q[7];
rzz(-pi/2) q[0],q[7];
u(pi/2,2.4912829742967055,-2.4912829742967055) q[8];
u(pi/2,-0.638371627209446,0.638371627209446) q[8];
rzz(-pi/2) q[0],q[8];
rzz(pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[13];
rzz(-pi/2) q[0],q[14];
rzz(-pi/2) q[0],q[15];
u(pi/2,0,0) q[1];
rzz(-pi/2) q[1],q[2];
u(pi/2,pi,-pi) q[2];
u(pi/2,pi/4,-pi/4) q[2];
rzz(-pi/2) q[1],q[2];
rzz(-pi/2) q[1],q[3];
u(pi/2,-pi/4,pi/4) q[3];
u(pi/2,7*pi/8,-7*pi/8) q[3];
rzz(-pi/2) q[1],q[3];
rzz(-pi/2) q[1],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,-0.5893627818134451,0.5893627818134451) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u(pi/2,pi/8,-pi/8) q[5];
u(pi/2,3.6323094260805187,-3.6323094260805187) q[5];
rzz(-pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u(pi/2,0.589362781813445,-0.589362781813445) q[6];
u(pi/2,3.7799642807992395,-3.7799642807992395) q[6];
rzz(-pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u(pi/2,2.454212180984346,-2.454212180984346) q[7];
u(pi/2,-0.6628760499074464,0.6628760499074464) q[7];
rzz(-pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[8];
u(pi/2,-0.638371627209446,0.638371627209446) q[8];
u(pi/2,2.515787396994706,-2.515787396994706) q[8];
rzz(pi/2) q[1],q[8];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[10];
rzz(-pi/2) q[1],q[10];
rzz(pi/2) q[1],q[11];
rzz(pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[12];
rzz(pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[13];
rzz(-pi/2) q[1],q[13];
rzz(-pi/2) q[1],q[14];
rzz(pi/2) q[1],q[14];
rzz(pi/2) q[1],q[15];
rzz(-pi/2) q[1],q[15];
u(pi/2,-pi/4,pi/4) q[2];
u(pi/2,pi,-pi) q[2];
rzz(pi/2) q[2],q[3];
u(pi/2,-pi/8,pi/8) q[3];
u(pi/2,9*pi/8,-9*pi/8) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u(pi/2,-0.5893627818134451,0.5893627818134451) q[4];
u(pi/2,2.9455572720057903,-2.9455572720057903) q[4];
rzz(pi/2) q[2],q[4];
rzz(-pi/2) q[2],q[5];
u(pi/2,3.6323094260805187,-3.6323094260805187) q[5];
u(pi/2,0.6873804726054469,-0.6873804726054469) q[5];
rzz(pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[6];
u(pi/2,3.7799642807992395,-3.7799642807992395) q[6];
u(pi/2,0.7363893180014478,-0.7363893180014478) q[6];
rzz(pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u(pi/2,2.4787166036823463,-2.4787166036823463) q[7];
u(pi/2,-0.6138672045114456,0.6138672045114456) q[7];
rzz(-pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[8];
u(pi/2,-0.6258052565950868,0.6258052565950868) q[8];
u(pi/2,2.5402918196927065,-2.5402918196927065) q[8];
rzz(-pi/2) q[2],q[8];
rzz(-pi/2) q[2],q[9];
u(pi/2,2.190318398082804,-2.190318398082804) q[9];
u(pi/2,-0.9387078848926302,0.9387078848926302) q[9];
rzz(pi/2) q[2],q[9];
rzz(-pi/2) q[2],q[10];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[11];
rzz(-pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[13];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[14];
rzz(pi/2) q[2],q[14];
rzz(-pi/2) q[2],q[15];
rzz(pi/2) q[2],q[15];
u(pi/2,-3*pi/8,3*pi/8) q[3];
u(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(pi/2) q[3],q[4];
u(pi/2,2.9455572720057903,-2.9455572720057903) q[4];
u(pi/2,0.589362781813445,-0.589362781813445) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u(pi/2,0.6873804726054469,-0.6873804726054469) q[5];
u(pi/2,4.221672207893964,-4.221672207893964) q[5];
rzz(-pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u(pi/2,0.7363893180014478,-0.7363893180014478) q[6];
u(pi/2,4.074017353175243,-4.074017353175243) q[6];
rzz(pi/2) q[3],q[6];
rzz(-pi/2) q[3],q[7];
u(pi/2,-0.6138672045114456,0.6138672045114456) q[7];
u(pi/2,2.626371458401067,-2.626371458401067) q[7];
rzz(pi/2) q[3],q[7];
rzz(-pi/2) q[3],q[8];
u(pi/2,2.5402918196927065,-2.5402918196927065) q[8];
u(pi/2,-0.5522919885010857,0.5522919885010857) q[8];
rzz(-pi/2) q[3],q[8];
rzz(pi/2) q[3],q[9];
u(pi/2,-0.9387078848926302,0.9387078848926302) q[9];
u(pi/2,2.2273891913951633,-2.2273891913951633) q[9];
rzz(pi/2) q[3],q[9];
rzz(pi/2) q[3],q[10];
u(pi/2,2.199743176043573,-2.199743176043573) q[10];
u(pi/2,-0.9292831069318608,0.9292831069318608) q[10];
rzz(-pi/2) q[3],q[10];
rzz(-pi/2) q[3],q[11];
rzz(-pi/2) q[3],q[11];
rzz(-pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[13];
rzz(pi/2) q[3],q[13];
rzz(pi/2) q[3],q[14];
rzz(-pi/2) q[3],q[14];
rzz(pi/2) q[3],q[15];
rzz(-pi/2) q[3],q[15];
u(pi/2,2.1601591086083416,-2.1601591086083416) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[4],q[5];
u(pi/2,4.221672207893964,-4.221672207893964) q[5];
u(pi/2,1.8654777177016193,-1.8654777177016193) q[5];
rzz(pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u(pi/2,4.074017353175243,-4.074017353175243) q[6];
u(pi/2,1.3251237812841747,-1.3251237812841747) q[6];
rzz(pi/2) q[4],q[6];
rzz(-pi/2) q[4],q[7];
u(pi/2,-0.515221195188726,0.515221195188726) q[7];
u(pi/2,2.822406839985071,-2.822406839985071) q[7];
rzz(pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u(pi/2,2.5893006650887074,-2.5893006650887074) q[8];
u(pi/2,-0.4542742977090841,0.4542742977090841) q[8];
rzz(pi/2) q[4],q[8];
rzz(-pi/2) q[4],q[9];
u(pi/2,-0.9142034621946298,0.9142034621946298) q[9];
u(pi/2,2.276398036791164,-2.276398036791164) q[9];
rzz(-pi/2) q[4],q[9];
rzz(pi/2) q[4],q[10];
u(pi/2,2.212309546657932,-2.212309546657932) q[10];
u(pi/2,-0.9047786842338604,0.9047786842338604) q[10];
rzz(pi/2) q[4],q[10];
rzz(pi/2) q[4],q[11];
u(pi/2,3.6800616344150834,-3.6800616344150834) q[11];
u(pi/2,0.5504070329089319,-0.5504070329089319) q[11];
rzz(pi/2) q[4],q[11];
rzz(-pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[12];
rzz(pi/2) q[4],q[13];
rzz(pi/2) q[4],q[13];
rzz(pi/2) q[4],q[14];
rzz(-pi/2) q[4],q[14];
rzz(-pi/2) q[4],q[15];
rzz(pi/2) q[4],q[15];
u(pi/2,0.2946813909067225,-0.2946813909067225) q[5];
u(pi/2,0.19603538158400302,-0.19603538158400302) q[5];
rzz(pi/2) q[5],q[6];
u(pi/2,1.3251237812841747,-1.3251237812841747) q[6];
u(pi/2,-1.03107070890817,1.03107070890817) q[6];
rzz(pi/2) q[5],q[6];
rzz(-pi/2) q[5],q[7];
u(pi/2,-0.319185813604723,0.319185813604723) q[7];
u(pi/2,3.2151059216837945,-3.2151059216837945) q[7];
rzz(pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u(pi/2,-0.4542742977090841,0.4542742977090841) q[8];
u(pi/2,2.88398205599543,-2.88398205599543) q[8];
rzz(-pi/2) q[5],q[8];
rzz(pi/2) q[5],q[9];
u(pi/2,-0.8651946167986291,0.8651946167986291) q[9];
u(pi/2,2.3744157275831657,-2.3744157275831657) q[9];
rzz(-pi/2) q[5],q[9];
rzz(pi/2) q[5],q[10];
u(pi/2,-0.9047786842338604,0.9047786842338604) q[10];
u(pi/2,2.2858228147519335,-2.2858228147519335) q[10];
rzz(pi/2) q[5],q[10];
rzz(pi/2) q[5],q[11];
u(pi/2,0.5504070329089319,-0.5504070329089319) q[11];
u(pi/2,3.717132427727443,-3.717132427727443) q[11];
rzz(pi/2) q[5],q[11];
rzz(-pi/2) q[5],q[12];
u(pi/2,-0.2550973234714913,0.2550973234714913) q[12];
u(pi/2,2.8984333822019437,-2.8984333822019437) q[12];
rzz(pi/2) q[5],q[12];
rzz(pi/2) q[5],q[13];
rzz(pi/2) q[5],q[13];
rzz(pi/2) q[5],q[14];
rzz(-pi/2) q[5],q[14];
rzz(pi/2) q[5],q[15];
rzz(pi/2) q[5],q[15];
u(pi/2,0.5397256178867265,-0.5397256178867265) q[6];
u(pi/2,3.7309554354032386,-3.7309554354032386) q[6];
rzz(pi/2) q[6],q[7];
u(pi/2,3.2151059216837945,-3.2151059216837945) q[7];
u(pi/2,0.8589114314914492,-0.8589114314914492) q[7];
rzz(pi/2) q[6],q[7];
rzz(-pi/2) q[6],q[8];
u(pi/2,2.88398205599543,-2.88398205599543) q[8];
u(pi/2,0.13508848410436114,-0.13508848410436114) q[8];
rzz(pi/2) q[6],q[8];
rzz(-pi/2) q[6],q[9];
u(pi/2,2.3744157275831657,-2.3744157275831657) q[9];
u(pi/2,-0.5705132258919063,0.5705132258919063) q[9];
rzz(pi/2) q[6],q[9];
rzz(-pi/2) q[6],q[10];
u(pi/2,-0.8557698388378596,0.8557698388378596) q[10];
u(pi/2,2.383840505543935,-2.383840505543935) q[10];
rzz(-pi/2) q[6],q[10];
rzz(pi/2) q[6],q[11];
u(pi/2,3.717132427727443,-3.717132427727443) q[11];
u(pi/2,0.6245486195336509,-0.6245486195336509) q[11];
rzz(-pi/2) q[6],q[11];
rzz(pi/2) q[6],q[12];
u(pi/2,2.8984333822019437,-2.8984333822019437) q[12];
u(pi/2,-0.21865484868984963,0.21865484868984963) q[12];
rzz(pi/2) q[6],q[12];
rzz(-pi/2) q[6],q[13];
u(pi/2,2.4900263372352702,-2.4900263372352702) q[13];
u(pi/2,-0.6389999457401639,0.6389999457401639) q[13];
rzz(pi/2) q[6],q[13];
rzz(pi/2) q[6],q[14];
rzz(pi/2) q[6],q[14];
rzz(pi/2) q[6],q[15];
rzz(-pi/2) q[6],q[15];
u(pi/2,-0.7118848953034472,0.7118848953034472) q[7];
u(pi/2,-0.7363893180014475,0.7363893180014475) q[7];
rzz(pi/2) q[7],q[8];
u(pi/2,0.13508848410436114,-0.13508848410436114) q[8];
u(pi/2,4.062079301091602,-4.062079301091602) q[8];
rzz(pi/2) q[7],q[8];
rzz(pi/2) q[7],q[9];
u(pi/2,-0.5705132258919063,0.5705132258919063) q[9];
u(pi/2,2.9637785093966107,-2.9637785093966107) q[9];
rzz(pi/2) q[7],q[9];
rzz(-pi/2) q[7],q[10];
u(pi/2,2.383840505543935,-2.383840505543935) q[10];
u(pi/2,-0.561716766461855,0.561716766461855) q[10];
rzz(pi/2) q[7],q[10];
rzz(pi/2) q[7],q[11];
u(pi/2,3.766141273123444,-3.766141273123444) q[11];
u(pi/2,0.7225663103256523,-0.7225663103256523) q[11];
rzz(pi/2) q[7],q[11];
rzz(pi/2) q[7],q[12];
u(pi/2,-0.21865484868984963,0.21865484868984963) q[12];
u(pi/2,2.971946650295944,-2.971946650295944) q[12];
rzz(-pi/2) q[7],q[12];
rzz(pi/2) q[7],q[13];
u(pi/2,-0.6389999457401639,0.6389999457401639) q[13];
u(pi/2,2.52709713054763,-2.52709713054763) q[13];
rzz(pi/2) q[7],q[13];
rzz(pi/2) q[7],q[14];
u(pi/2,4.068990804929499,-4.068990804929499) q[14];
u(pi/2,0.9393362034233483,-0.9393362034233483) q[14];
rzz(pi/2) q[7],q[14];
rzz(-pi/2) q[7],q[15];
rzz(pi/2) q[7],q[15];
u(pi/2,-0.6503096792930873,0.6503096792930873) q[8];
u(pi/2,2.5032210263803467,-2.5032210263803467) q[8];
rzz(pi/2) q[8],q[9];
u(pi/2,2.9637785093966107,-2.9637785093966107) q[9];
u(pi/2,0.6075840192042659,-0.6075840192042659) q[9];
rzz(pi/2) q[8],q[9];
rzz(-pi/2) q[8],q[10];
u(pi/2,2.5798758871279377,-2.5798758871279377) q[10];
u(pi/2,-0.16901768476313084,0.16901768476313084) q[10];
rzz(-pi/2) q[8],q[10];
rzz(pi/2) q[8],q[11];
u(pi/2,0.7225663103256523,-0.7225663103256523) q[11];
u(pi/2,4.060194345499449,-4.060194345499449) q[11];
rzz(pi/2) q[8],q[11];
rzz(pi/2) q[8],q[12];
u(pi/2,-0.1696460032938487,0.1696460032938487) q[12];
u(pi/2,3.069964341087946,-3.069964341087946) q[12];
rzz(pi/2) q[8],q[12];
rzz(-pi/2) q[8],q[13];
u(pi/2,-0.6144955230421635,0.6144955230421635) q[13];
u(pi/2,2.576105975943631,-2.576105975943631) q[13];
rzz(pi/2) q[8],q[13];
rzz(pi/2) q[8],q[14];
u(pi/2,0.9393362034233483,-0.9393362034233483) q[14];
u(pi/2,4.105433279711142,-4.105433279711142) q[14];
rzz(pi/2) q[8],q[14];
rzz(pi/2) q[8],q[15];
u(pi/2,3.103265223215998,-3.103265223215998) q[15];
u(pi/2,-0.02638937829015431,0.02638937829015431) q[15];
rzz(-pi/2) q[8],q[15];
rzz(pi/2) q[9],q[10];
u(pi/2,2.9725749688266623,-2.9725749688266623) q[10];
u(pi/2,0.6163804786343174,-0.6163804786343174) q[10];
rzz(pi/2) q[9],q[10];
rzz(pi/2) q[9],q[11];
u(pi/2,4.060194345499449,-4.060194345499449) q[11];
u(pi/2,1.3113007736083797,-1.3113007736083797) q[11];
rzz(pi/2) q[9],q[11];
rzz(-pi/2) q[9],q[12];
u(pi/2,-0.07162831250184731,0.07162831250184731) q[12];
u(pi/2,3.266628041202667,-3.266628041202667) q[12];
rzz(-pi/2) q[9],q[12];
rzz(pi/2) q[9],q[13];
u(pi/2,2.576105975943631,-2.576105975943631) q[13];
u(pi/2,-0.4674689868541613,0.4674689868541613) q[13];
rzz(-pi/2) q[9],q[13];
rzz(-pi/2) q[9],q[14];
u(pi/2,0.9638406261213484,-0.9638406261213484) q[14];
u(pi/2,4.154442125107143,-4.154442125107143) q[14];
rzz(-pi/2) q[9],q[14];
rzz(-pi/2) q[9],q[15];
u(pi/2,-0.02638937829015431,0.02638937829015431) q[15];
u(pi/2,3.1397076979976397,-3.1397076979976397) q[15];
rzz(pi/2) q[9],q[15];
rzz(pi/2) q[10],q[11];
u(pi/2,1.3113007736083797,-1.3113007736083797) q[11];
u(pi/2,-1.0448937165839651,1.0448937165839651) q[11];
rzz(-pi/2) q[10],q[11];
rzz(-pi/2) q[10],q[12];
u(pi/2,3.266628041202667,-3.266628041202667) q[12];
u(pi/2,0.517734469311598,-0.517734469311598) q[12];
rzz(pi/2) q[10],q[12];
rzz(pi/2) q[10],q[13];
u(pi/2,2.674123666735632,-2.674123666735632) q[13];
u(pi/2,-0.27143360527015803,0.27143360527015803) q[13];
rzz(pi/2) q[10],q[13];
rzz(pi/2) q[10],q[14];
u(pi/2,1.0128494715173493,-1.0128494715173493) q[14];
u(pi/2,4.252459815899144,-4.252459815899144) q[14];
rzz(-pi/2) q[10],q[14];
rzz(pi/2) q[10],q[15];
u(pi/2,3.1397076979976397,-3.1397076979976397) q[15];
u(pi/2,0.04775220833456473,-0.04775220833456473) q[15];
rzz(pi/2) q[10],q[15];
u(pi,0.5259026102109314,-0.5259026102109314) q[11];
rzz(-pi/2) q[11],q[12];
u(pi/2,3.6593271229013915,-3.6593271229013915) q[12];
u(pi/2,1.3031326327090462,-1.3031326327090462) q[12];
rzz(pi/2) q[11],q[12];
rzz(-pi/2) q[11],q[13];
u(pi/2,2.8701590483196346,-2.8701590483196346) q[13];
u(pi/2,0.12126547642856589,-0.12126547642856589) q[13];
rzz(-pi/2) q[11],q[13];
rzz(pi/2) q[11],q[14];
u(pi/2,1.1108671623093511,-1.1108671623093511) q[14];
u(pi/2,4.449123516013865,-4.449123516013865) q[14];
rzz(pi/2) q[11],q[14];
rzz(pi/2) q[11],q[15];
u(pi/2,0.04775220833456473,-0.04775220833456473) q[15];
u(pi/2,3.2873625527163597,-3.2873625527163597) q[15];
rzz(-pi/2) q[11],q[15];
u(pi,-0.26766369408585033,0.26766369408585033) q[12];
rzz(pi/2) q[12],q[13];
u(pi/2,3.2628581300183592,-3.2628581300183592) q[13];
u(pi/2,0.9066636398260144,-0.9066636398260144) q[13];
rzz(pi/2) q[12],q[13];
rzz(pi/2) q[12],q[14];
u(pi/2,4.449123516013865,-4.449123516013865) q[14];
u(pi/2,1.7002299441227957,-1.7002299441227957) q[14];
rzz(pi/2) q[12],q[14];
rzz(-pi/2) q[12],q[15];
u(pi/2,3.2873625527163597,-3.2873625527163597) q[15];
u(pi/2,0.3418052807105696,-0.3418052807105696) q[15];
rzz(pi/2) q[12],q[15];
rzz(-pi/2) q[13],q[14];
u(pi/2,-1.441362709466997,1.441362709466997) q[14];
u(pi/2,2.4856281075202444,-2.4856281075202444) q[14];
rzz(pi/2) q[13],q[14];
rzz(-pi/2) q[13],q[15];
u(pi/2,3.4833979343003625,-3.4833979343003625) q[15];
u(pi/2,0.7345043624092935,-0.7345043624092935) q[15];
rzz(-pi/2) q[13],q[15];
u(pi,0.9148317807253479,-0.9148317807253479) q[14];
rzz(pi/2) q[14],q[15];
u(pi/2,3.876097015999087,-3.876097015999087) q[15];
u(pi/2,1.5199025258067418,-1.5199025258067418) q[15];
rzz(pi/2) q[14],q[15];
u(pi,pi,-pi) q[0];
u(pi,0,0) q[1];
u(pi,0,0) q[2];
u(pi,-pi/4,pi/4) q[3];
u(pi,pi/4,-pi/4) q[4];
u(pi,0.9814335449814515,-0.9814335449814515) q[5];
u(pi,-0.5893627818134451,0.5893627818134451) q[6];
u(pi,1.03107070890817,-1.03107070890817) q[7];
u(pi,-1.3251237812841747,1.3251237812841747) q[8];
u(pi,2.626371458401067,-2.626371458401067) q[9];
u(pi,3.9759996623832423,-3.9759996623832423) q[10];
u(pi,-0.15016812884159214,0.15016812884159214) q[11];
u(pi,3.321291753375129,-3.321291753375129) q[12];
u(pi,1.1316016738230434,-1.1316016738230434) q[13];
u(pi,2.597468805988041,-2.597468805988041) q[14];
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
