OPENQASM 2.0;
include "qelib1.inc";
qreg q[18];
creg c[18];
u3(pi/2,3.666238626739289,-3.666238626739289) q[16];
u3(pi/2,1.670070654648334,-1.670070654648334) q[17];
rzz(-pi/2) q[16],q[17];
u3(pi/2,0.09927432785343747,-0.09927432785343747) q[17];
u3(pi/2,1.9546989490635691,-1.9546989490635691) q[17];
rzz(pi/2) q[16],q[17];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[15];
rzz(pi/2) q[15],q[17];
u3(pi/2,1.9546989490635691,-1.9546989490635691) q[17];
u3(pi/2,2.523327219363322,-2.523327219363322) q[17];
rzz(pi/2) q[15],q[17];
u3(pi/2,2.7991590543485056,-2.7991590543485056) q[14];
rzz(pi/2) q[14],q[17];
u3(pi/2,5.6649198729531145,-5.6649198729531145) q[17];
u3(pi/2,1.3860706787638166,-1.3860706787638166) q[17];
rzz(-pi/2) q[14],q[17];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[13];
rzz(pi/2) q[13],q[17];
u3(pi/2,4.52766333235361,-4.52766333235361) q[17];
u3(pi/2,5.394114586213675,-5.394114586213675) q[17];
rzz(pi/2) q[13],q[17];
u3(pi/2,3.340141309296668,-3.340141309296668) q[12];
rzz(pi/2) q[12],q[17];
u3(pi/2,2.2525219326238815,-2.2525219326238815) q[17];
u3(pi/2,3.660583759962827,-3.660583759962827) q[17];
rzz(-pi/2) q[12],q[17];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[11];
rzz(pi/2) q[11],q[17];
u3(pi/2,3.660583759962827,-3.660583759962827) q[17];
u3(pi/2,3.986052758874729,-3.986052758874729) q[17];
rzz(pi/2) q[11],q[17];
u3(pi/2,2.365619268153114,-2.365619268153114) q[10];
rzz(pi/2) q[10],q[17];
u3(pi/2,0.8444601052849363,-0.8444601052849363) q[17];
u3(pi/2,3.3357430795816425,-3.3357430795816425) q[17];
rzz(pi/2) q[10],q[17];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
rzz(-pi/2) q[9],q[17];
u3(pi/2,0.19415042599184923,-0.19415042599184923) q[17];
u3(pi/2,2.0344954024647497,-2.0344954024647497) q[17];
rzz(pi/2) q[9],q[17];
u3(pi/2,4.749459773697049,-4.749459773697049) q[8];
rzz(pi/2) q[8],q[17];
u3(pi/2,2.0344954024647497,-2.0344954024647497) q[17];
u3(pi/2,2.5748493388821942,-2.5748493388821942) q[17];
rzz(pi/2) q[8],q[17];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
rzz(-pi/2) q[7],q[17];
u3(pi/2,2.5748493388821942,-2.5748493388821942) q[17];
u3(pi/2,4.636362438167817,-4.636362438167817) q[17];
rzz(pi/2) q[7],q[17];
u3(pi/2,1.7178228629828987,-1.7178228629828987) q[6];
rzz(pi/2) q[6],q[17];
u3(pi/2,4.636362438167817,-4.636362438167817) q[17];
u3(pi/2,5.617795983149268,-5.617795983149268) q[17];
rzz(pi/2) q[6],q[17];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(-pi/2) q[5],q[17];
u3(pi/2,5.617795983149268,-5.617795983149268) q[17];
u3(pi/2,0.5127079210658543,-0.5127079210658543) q[17];
rzz(pi/2) q[5],q[17];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[4];
rzz(pi/2) q[4],q[17];
u3(pi/2,3.654300574655647,-3.654300574655647) q[17];
u3(pi/2,4.439698738053096,-4.439698738053096) q[17];
rzz(pi/2) q[4],q[17];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
u3(pi/2,2.868902411258199,-2.868902411258199) q[17];
rzz(-pi/2) q[3],q[17];
u3(pi/2,pi,-pi) q[1];
u3(pi/2,5*pi/4,-5*pi/4) q[1];
rzz(-pi/2) q[0],q[1];
rzz(-pi/2) q[0],q[2];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
u3(pi/2,5*pi/8,-5*pi/8) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,15*pi/8,-15*pi/8) q[3];
u3(pi/2,2.552229871776348,-2.552229871776348) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,5.301751762198135,-5.301751762198135) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,3.730955435403238,-3.730955435403238) q[4];
u3(pi/2,0.4907167724907257,-0.4907167724907257) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,5.399769452990137,-5.399769452990137) q[5];
u3(pi/2,2.2091679540043425,-2.2091679540043425) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,4.8594155165726916,-4.8594155165726916) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,3.2886191897777954,-3.2886191897777954) q[6];
u3(pi/2,0.12252211349000193,-0.12252211349000193) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,4.466716434873968,-4.466716434873968) q[7];
u3(pi/2,1.3131857292005336,-1.3131857292005336) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,1.6078671201072563,-1.6078671201072563) q[8];
rzz(-pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[8];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
rzz(-pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[9];
u3(pi/2,5.507211921742907,-5.507211921742907) q[10];
rzz(pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[10];
u3(pi/2,3.669380219392878,-3.669380219392878) q[11];
rzz(-pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[11];
u3(pi/2,3.340141309296668,-3.340141309296668) q[12];
rzz(pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[12];
u3(pi/2,3.669380219392878,-3.669380219392878) q[13];
rzz(-pi/2) q[0],q[13];
rzz(-pi/2) q[0],q[13];
u3(pi/2,2.7991590543485056,-2.7991590543485056) q[14];
rzz(pi/2) q[0],q[14];
rzz(-pi/2) q[0],q[14];
u3(pi/2,3.669380219392878,-3.669380219392878) q[15];
rzz(-pi/2) q[0],q[15];
rzz(-pi/2) q[0],q[15];
u3(pi/2,3.666238626739289,-3.666238626739289) q[16];
rzz(pi/2) q[0],q[16];
rzz(-pi/2) q[0],q[16];
u3(pi/2,7*pi/4,-7*pi/4) q[1];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,5*pi/8,-5*pi/8) q[2];
u3(pi/2,11*pi/8,-11*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,2.552229871776348,-2.552229871776348) q[3];
u3(pi/2,5.301751762198135,-5.301751762198135) q[3];
rzz(-pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,0.4907167724907257,-0.4907167724907257) q[4];
u3(pi/2,3.436274044496516,-3.436274044496516) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
u3(pi/2,2.110521944681623,-2.110521944681623) q[5];
rzz(-pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,0.12252211349000193,-0.12252211349000193) q[6];
u3(pi/2,3.2151059216837945,-3.2151059216837945) q[6];
rzz(pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u3(pi/2,4.454778382790327,-4.454778382790327) q[7];
u3(pi/2,1.288681306502533,-1.288681306502533) q[7];
rzz(pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[8];
u3(pi/2,0.030787608005179972,-0.030787608005179972) q[8];
u3(pi/2,3.159813890980614,-3.159813890980614) q[8];
rzz(pi/2) q[1],q[8];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[10];
rzz(pi/2) q[1],q[10];
rzz(pi/2) q[1],q[11];
rzz(pi/2) q[1],q[11];
rzz(pi/2) q[1],q[12];
rzz(pi/2) q[1],q[12];
rzz(pi/2) q[1],q[13];
rzz(pi/2) q[1],q[13];
rzz(-pi/2) q[1],q[14];
rzz(pi/2) q[1],q[14];
rzz(pi/2) q[1],q[15];
rzz(pi/2) q[1],q[15];
rzz(pi/2) q[1],q[16];
rzz(pi/2) q[1],q[16];
u3(pi/2,15*pi/8,-15*pi/8) q[2];
u3(pi/2,pi/4,-pi/4) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[3];
u3(pi/2,4.516353598800687,-4.516353598800687) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,3.436274044496516,-3.436274044496516) q[4];
u3(pi/2,6.185167616387585,-6.185167616387585) q[4];
rzz(-pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,5.252114598271416,-5.252114598271416) q[5];
u3(pi/2,1.91448656309762,-1.91448656309762) q[5];
rzz(pi/2) q[2],q[5];
rzz(pi/2) q[2],q[6];
u3(pi/2,3.2151059216837945,-3.2151059216837945) q[6];
u3(pi/2,6.258680884481586,-6.258680884481586) q[6];
rzz(-pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,1.288681306502533,-1.288681306502533) q[7];
u3(pi/2,4.381265114696325,-4.381265114696325) q[7];
rzz(pi/2) q[2],q[7];
rzz(pi/2) q[2],q[8];
u3(pi/2,3.159813890980614,-3.159813890980614) q[8];
u3(pi/2,6.276902121872407,-6.276902121872407) q[8];
rzz(pi/2) q[2],q[8];
rzz(-pi/2) q[2],q[9];
u3(pi/2,2.0954422999443922,-2.0954422999443922) q[9];
u3(pi/2,5.224468582919826,-5.224468582919826) q[9];
rzz(pi/2) q[2],q[9];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[10];
rzz(-pi/2) q[2],q[11];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[12];
rzz(pi/2) q[2],q[12];
rzz(pi/2) q[2],q[13];
rzz(-pi/2) q[2],q[13];
rzz(pi/2) q[2],q[14];
rzz(pi/2) q[2],q[14];
rzz(-pi/2) q[2],q[15];
rzz(pi/2) q[2],q[15];
rzz(pi/2) q[2],q[16];
rzz(pi/2) q[2],q[16];
u3(pi/2,6.087149925595583,-6.087149925595583) q[3];
u3(pi/2,7*pi/8,-7*pi/8) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.0435749627977917,-3.0435749627977917) q[4];
u3(pi/2,5.399769452990137,-5.399769452990137) q[4];
rzz(-pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,1.91448656309762,-1.91448656309762) q[5];
u3(pi/2,4.663380134988689,-4.663380134988689) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,3.1170882308917927,-3.1170882308917927) q[6];
u3(pi/2,6.062017184366865,-6.062017184366865) q[6];
rzz(-pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,4.381265114696325,-4.381265114696325) q[7];
u3(pi/2,1.141026451783813,-1.141026451783813) q[7];
rzz(pi/2) q[3],q[7];
rzz(pi/2) q[3],q[8];
u3(pi/2,6.276902121872407,-6.276902121872407) q[8];
u3(pi/2,3.086300622886613,-3.086300622886613) q[8];
rzz(pi/2) q[3],q[8];
rzz(-pi/2) q[3],q[9];
u3(pi/2,2.082875929330033,-2.082875929330033) q[9];
u3(pi/2,5.199964160221826,-5.199964160221826) q[9];
rzz(pi/2) q[3],q[9];
rzz(pi/2) q[3],q[10];
u3(pi/2,0.7841415263360123,-0.7841415263360123) q[10];
u3(pi/2,3.9131678093114464,-3.9131678093114464) q[10];
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
rzz(pi/2) q[3],q[16];
rzz(pi/2) q[3],q[16];
u3(pi/2,3.82897312619524,-3.82897312619524) q[4];
u3(pi/2,4.516353598800687,-4.516353598800687) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,4.663380134988689,-4.663380134988689) q[5];
u3(pi/2,0.7363893180014475,-0.7363893180014475) q[5];
rzz(pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,2.9204245307770718,-2.9204245307770718) q[6];
u3(pi/2,5.66931810266814,-5.66931810266814) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,1.141026451783813,-1.141026451783813) q[7];
u3(pi/2,4.0865837237896026,-4.0865837237896026) q[7];
rzz(-pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,3.086300622886613,-3.086300622886613) q[8];
u3(pi/2,6.129875585684404,-6.129875585684404) q[8];
rzz(pi/2) q[4],q[8];
rzz(pi/2) q[4],q[9];
u3(pi/2,5.199964160221826,-5.199964160221826) q[9];
u3(pi/2,2.0093626612360316,-2.0093626612360316) q[9];
rzz(-pi/2) q[4],q[9];
rzz(pi/2) q[4],q[10];
u3(pi/2,3.9131678093114464,-3.9131678093114464) q[10];
u3(pi/2,0.7470707330236528,-0.7470707330236528) q[10];
rzz(pi/2) q[4],q[10];
rzz(pi/2) q[4],q[11];
u3(pi/2,2.0897874331679303,-2.0897874331679303) q[11];
u3(pi/2,5.218813716143364,-5.218813716143364) q[11];
rzz(pi/2) q[4],q[11];
rzz(-pi/2) q[4],q[12];
rzz(pi/2) q[4],q[12];
rzz(pi/2) q[4],q[13];
rzz(pi/2) q[4],q[13];
rzz(-pi/2) q[4],q[14];
rzz(pi/2) q[4],q[14];
rzz(pi/2) q[4],q[15];
rzz(pi/2) q[4],q[15];
rzz(pi/2) q[4],q[16];
rzz(-pi/2) q[4],q[16];
u3(pi/2,2.3071856447963444,-2.3071856447963444) q[5];
u3(pi/2,5.007070371291412,-5.007070371291412) q[5];
rzz(pi/2) q[5],q[6];
u3(pi/2,5.66931810266814,-5.66931810266814) q[6];
u3(pi/2,1.7423272856808991,-1.7423272856808991) q[6];
rzz(pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,0.9449910701998098,-0.9449910701998098) q[7];
u3(pi/2,3.6938846420908784,-3.6938846420908784) q[7];
rzz(-pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,6.129875585684404,-6.129875585684404) q[8];
u3(pi/2,2.79161923197989,-2.79161923197989) q[8];
rzz(pi/2) q[5],q[8];
rzz(pi/2) q[5],q[9];
u3(pi/2,5.150955314825825,-5.150955314825825) q[9];
u3(pi/2,1.9113449704440304,-1.9113449704440304) q[9];
rzz(-pi/2) q[5],q[9];
rzz(pi/2) q[5],q[10];
u3(pi/2,0.7470707330236528,-0.7470707330236528) q[10];
u3(pi/2,3.839654541217445,-3.839654541217445) q[10];
rzz(pi/2) q[5],q[10];
rzz(pi/2) q[5],q[11];
u3(pi/2,5.218813716143364,-5.218813716143364) q[11];
u3(pi/2,2.052716639855571,-2.052716639855571) q[11];
rzz(-pi/2) q[5],q[11];
rzz(pi/2) q[5],q[12];
u3(pi/2,4.902769495192231,-4.902769495192231) q[12];
u3(pi/2,1.7492387895187966,-1.7492387895187966) q[12];
rzz(pi/2) q[5],q[12];
rzz(pi/2) q[5],q[13];
rzz(pi/2) q[5],q[13];
rzz(-pi/2) q[5],q[14];
rzz(pi/2) q[5],q[14];
rzz(pi/2) q[5],q[15];
rzz(pi/2) q[5],q[15];
rzz(-pi/2) q[5],q[16];
rzz(pi/2) q[5],q[16];
u3(pi/2,3.313123612475796,-3.313123612475796) q[6];
u3(pi/2,3.8779819715912405,-3.8779819715912405) q[6];
rzz(pi/2) q[6],q[7];
u3(pi/2,0.5522919885010856,-0.5522919885010856) q[7];
u3(pi/2,2.90848647869343,-2.90848647869343) q[7];
rzz(pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,2.79161923197989,-2.79161923197989) q[8];
u3(pi/2,5.54051280387096,-5.54051280387096) q[8];
rzz(pi/2) q[6],q[8];
rzz(pi/2) q[6],q[9];
u3(pi/2,5.052937624033824,-5.052937624033824) q[9];
u3(pi/2,1.7146812703293088,-1.7146812703293088) q[9];
rzz(pi/2) q[6],q[9];
rzz(pi/2) q[6],q[10];
u3(pi/2,3.839654541217445,-3.839654541217445) q[10];
u3(pi/2,0.6000441968356505,-0.6000441968356505) q[10];
rzz(-pi/2) q[6],q[10];
rzz(pi/2) q[6],q[11];
u3(pi/2,5.194309293445364,-5.194309293445364) q[11];
u3(pi/2,2.00370779445957,-2.00370779445957) q[11];
rzz(pi/2) q[6],q[11];
rzz(pi/2) q[6],q[12];
u3(pi/2,1.7492387895187966,-1.7492387895187966) q[12];
u3(pi/2,4.865698701879872,-4.865698701879872) q[12];
rzz(-pi/2) q[6],q[12];
rzz(pi/2) q[6],q[13];
u3(pi/2,5.230751768227005,-5.230751768227005) q[13];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[13];
rzz(pi/2) q[6],q[13];
rzz(pi/2) q[6],q[14];
rzz(-pi/2) q[6],q[14];
rzz(pi/2) q[6],q[15];
rzz(pi/2) q[6],q[15];
rzz(pi/2) q[6],q[16];
rzz(-pi/2) q[6],q[16];
u3(pi/2,1.337690151898534,-1.337690151898534) q[7];
u3(pi/2,1.8409732950036186,-1.8409732950036186) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,5.54051280387096,-5.54051280387096) q[8];
u3(pi/2,1.6135219868837176,-1.6135219868837176) q[8];
rzz(pi/2) q[7],q[8];
rzz(pi/2) q[7],q[9];
u3(pi/2,1.7146812703293088,-1.7146812703293088) q[9];
u3(pi/2,4.463574842220378,-4.463574842220378) q[9];
rzz(pi/2) q[7],q[9];
rzz(-pi/2) q[7],q[10];
u3(pi/2,0.6000441968356505,-0.6000441968356505) q[10];
u3(pi/2,3.5449731503107227,-3.5449731503107227) q[10];
rzz(pi/2) q[7],q[10];
rzz(pi/2) q[7],q[11];
u3(pi/2,2.00370779445957,-2.00370779445957) q[11];
u3(pi/2,5.047282757257362,-5.047282757257362) q[11];
rzz(pi/2) q[7],q[11];
rzz(-pi/2) q[7],q[12];
u3(pi/2,4.865698701879872,-4.865698701879872) q[12];
u3(pi/2,1.6750972028940778,-1.6750972028940778) q[12];
rzz(pi/2) q[7],q[12];
rzz(pi/2) q[7],q[13];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[13];
u3(pi/2,5.194309293445364,-5.194309293445364) q[13];
rzz(pi/2) q[7],q[13];
rzz(pi/2) q[7],q[14];
u3(pi/2,1.216424675469968,-1.216424675469968) q[14];
u3(pi/2,4.34607927697612,-4.34607927697612) q[14];
rzz(-pi/2) q[7],q[14];
rzz(pi/2) q[7],q[15];
rzz(pi/2) q[7],q[15];
rzz(pi/2) q[7],q[16];
rzz(-pi/2) q[7],q[16];
u3(pi/2,3.1843183136786144,-3.1843183136786144) q[8];
u3(pi/2,5.289185391583776,-5.289185391583776) q[8];
rzz(pi/2) q[8],q[9];
u3(pi/2,4.463574842220378,-4.463574842220378) q[9];
u3(pi/2,0.5365840252331366,-0.5365840252331366) q[9];
rzz(pi/2) q[8],q[9];
rzz(pi/2) q[8],q[10];
u3(pi/2,3.5449731503107227,-3.5449731503107227) q[10];
u3(pi/2,0.010681415022205296,-0.010681415022205296) q[10];
rzz(pi/2) q[8],q[10];
rzz(-pi/2) q[8],q[11];
u3(pi/2,1.9056901036675686,-1.9056901036675686) q[11];
u3(pi/2,4.850619057142641,-4.850619057142641) q[11];
rzz(pi/2) q[8],q[11];
rzz(pi/2) q[8],q[12];
u3(pi/2,1.6750972028940778,-1.6750972028940778) q[12];
u3(pi/2,4.718672165691869,-4.718672165691869) q[12];
rzz(-pi/2) q[8],q[12];
rzz(pi/2) q[8],q[13];
u3(pi/2,5.194309293445364,-5.194309293445364) q[13];
u3(pi/2,2.00370779445957,-2.00370779445957) q[13];
rzz(pi/2) q[8],q[13];
rzz(pi/2) q[8],q[14];
u3(pi/2,1.2044866233863267,-1.2044866233863267) q[14];
u3(pi/2,4.321574854278119,-4.321574854278119) q[14];
rzz(pi/2) q[8],q[14];
rzz(-pi/2) q[8],q[15];
u3(pi/2,5.230751768227005,-5.230751768227005) q[15];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[15];
rzz(pi/2) q[8],q[15];
rzz(pi/2) q[8],q[16];
rzz(pi/2) q[8],q[16];
u3(pi/2,5.2489730056178265,-5.2489730056178265) q[9];
u3(pi/2,1.8409732950036186,-1.8409732950036186) q[9];
rzz(-pi/2) q[9],q[10];
u3(pi/2,3.1522740686119985,-3.1522740686119985) q[10];
u3(pi/2,5.508468558804344,-5.508468558804344) q[10];
rzz(pi/2) q[9],q[10];
rzz(pi/2) q[9],q[11];
u3(pi/2,4.850619057142641,-4.850619057142641) q[11];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[11];
rzz(pi/2) q[9],q[11];
rzz(pi/2) q[9],q[12];
u3(pi/2,1.5770795121020762,-1.5770795121020762) q[12];
u3(pi/2,4.522636784107866,-4.522636784107866) q[12];
rzz(pi/2) q[9],q[12];
rzz(pi/2) q[9],q[13];
u3(pi/2,2.00370779445957,-2.00370779445957) q[13];
u3(pi/2,5.047282757257362,-5.047282757257362) q[13];
rzz(pi/2) q[9],q[13];
rzz(pi/2) q[9],q[14];
u3(pi/2,4.321574854278119,-4.321574854278119) q[14];
u3(pi/2,1.1309733552923256,-1.1309733552923256) q[14];
rzz(-pi/2) q[9],q[14];
rzz(pi/2) q[9],q[15];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[15];
u3(pi/2,5.194309293445364,-5.194309293445364) q[15];
rzz(pi/2) q[9],q[15];
rzz(pi/2) q[9],q[16];
u3(pi/2,2.0872741590450583,-2.0872741590450583) q[16];
u3(pi/2,5.216300442020493,-5.216300442020493) q[16];
rzz(-pi/2) q[9],q[16];
u3(pi/2,3.937672232009447,-3.937672232009447) q[10];
u3(pi/2,6.160663193689585,-6.160663193689585) q[10];
rzz(pi/2) q[10],q[11];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[11];
u3(pi/2,3.6725218120464684,-3.6725218120464684) q[11];
rzz(pi/2) q[10],q[11];
rzz(pi/2) q[10],q[12];
u3(pi/2,4.522636784107866,-4.522636784107866) q[12];
u3(pi/2,0.9883450488193489,-0.9883450488193489) q[12];
rzz(pi/2) q[10],q[12];
rzz(-pi/2) q[10],q[13];
u3(pi/2,1.9056901036675686,-1.9056901036675686) q[13];
u3(pi/2,4.850619057142641,-4.850619057142641) q[13];
rzz(pi/2) q[10],q[13];
rzz(pi/2) q[10],q[14];
u3(pi/2,4.272566008882119,-4.272566008882119) q[14];
u3(pi/2,1.032955664500324,-1.032955664500324) q[14];
rzz(-pi/2) q[10],q[14];
rzz(pi/2) q[10],q[15];
u3(pi/2,5.194309293445364,-5.194309293445364) q[15];
u3(pi/2,2.00370779445957,-2.00370779445957) q[15];
rzz(pi/2) q[10],q[15];
rzz(pi/2) q[10],q[16];
u3(pi/2,2.0747077884306995,-2.0747077884306995) q[16];
u3(pi/2,5.191796019322492,-5.191796019322492) q[16];
rzz(pi/2) q[10],q[16];
u3(pi/2,5.2433181388413646,-5.2433181388413646) q[11];
u3(pi/2,0.20546015954477248,-0.20546015954477248) q[11];
rzz(-pi/2) q[11],q[12];
u3(pi/2,4.129937702409142,-4.129937702409142) q[12];
u3(pi/2,0.20294688542190065,-0.20294688542190065) q[12];
rzz(pi/2) q[11],q[12];
rzz(pi/2) q[11],q[13];
u3(pi/2,4.850619057142641,-4.850619057142641) q[13];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[13];
rzz(pi/2) q[11],q[13];
rzz(-pi/2) q[11],q[14];
u3(pi/2,1.032955664500324,-1.032955664500324) q[14];
u3(pi/2,3.977884617975396,-3.977884617975396) q[14];
rzz(pi/2) q[11],q[14];
rzz(pi/2) q[11],q[15];
u3(pi/2,2.00370779445957,-2.00370779445957) q[15];
u3(pi/2,5.046654438726644,-5.046654438726644) q[15];
rzz(pi/2) q[11],q[15];
rzz(pi/2) q[11],q[16];
u3(pi/2,5.191796019322492,-5.191796019322492) q[16];
u3(pi/2,2.001194520336698,-2.001194520336698) q[16];
rzz(-pi/2) q[11],q[16];
u3(pi/2,1.773743212216797,-1.773743212216797) q[12];
u3(pi/2,4.751973047819921,-4.751973047819921) q[12];
rzz(pi/2) q[12],q[13];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[13];
u3(pi/2,3.6725218120464684,-3.6725218120464684) q[13];
rzz(pi/2) q[12],q[13];
rzz(pi/2) q[12],q[14];
u3(pi/2,3.977884617975396,-3.977884617975396) q[14];
u3(pi/2,0.44359288268687874,-0.44359288268687874) q[14];
rzz(-pi/2) q[12],q[14];
rzz(pi/2) q[12],q[15];
u3(pi/2,5.046654438726644,-5.046654438726644) q[15];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[15];
rzz(pi/2) q[12],q[15];
rzz(pi/2) q[12],q[16];
u3(pi/2,5.142787173926491,-5.142787173926491) q[16];
u3(pi/2,1.9031768295446967,-1.9031768295446967) q[16];
rzz(pi/2) q[12],q[16];
u3(pi/2,5.2433181388413646,-5.2433181388413646) q[13];
u3(pi/2,5.947034893245479,-5.947034893245479) q[13];
rzz(-pi/2) q[13],q[14];
u3(pi/2,0.44359288268687874,-0.44359288268687874) q[14];
u3(pi/2,2.7997873728792237,-2.7997873728792237) q[14];
rzz(pi/2) q[13],q[14];
rzz(pi/2) q[13],q[15];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[15];
u3(pi/2,4.457919975443916,-4.457919975443916) q[15];
rzz(pi/2) q[13],q[15];
rzz(-pi/2) q[13],q[16];
u3(pi/2,5.044769483134489,-5.044769483134489) q[16];
u3(pi/2,1.7065131294299756,-1.7065131294299756) q[16];
rzz(pi/2) q[13],q[16];
u3(pi/2,1.228991046084327,-1.228991046084327) q[14];
u3(pi/2,1.6625308322797185,-1.6625308322797185) q[14];
rzz(pi/2) q[14],q[15];
u3(pi/2,4.457919975443916,-4.457919975443916) q[15];
u3(pi/2,0.5309291584566751,-0.5309291584566751) q[15];
rzz(pi/2) q[14],q[15];
rzz(-pi/2) q[14],q[16];
u3(pi/2,4.8481057830197685,-4.8481057830197685) q[16];
u3(pi/2,1.3138140477312514,-1.3138140477312514) q[16];
rzz(pi/2) q[14],q[16];
u3(pi/2,2.101725485251572,-2.101725485251572) q[15];
u3(pi/2,4.241150082346221,-4.241150082346221) q[15];
rzz(pi/2) q[15],q[16];
u3(pi/2,1.3138140477312514,-1.3138140477312514) q[16];
u3(pi/2,3.670008537923596,-3.670008537923596) q[16];
rzz(pi/2) q[15],q[16];
u3(pi,pi,-pi) q[1];
u3(pi,3*pi/2,-3*pi/2) q[2];
u3(pi,5*pi/4,-5*pi/4) q[3];
u3(pi,15*pi/8,-15*pi/8) q[4];
u3(pi,4.123026198571244,-4.123026198571244) q[5];
u3(pi,1.0800795543041708,-1.0800795543041708) q[6];
u3(pi,3.583300580684518,-3.583300580684518) q[7];
u3(pi,5.3752650302921365,-5.3752650302921365) q[8];
u3(pi,4.2398934452847845,-4.2398934452847845) q[9];
u3(pi,5.789326942035271,-5.789326942035271) q[10];
u3(pi,0.8168140899333463,-0.8168140899333463) q[11];
u3(pi,5.8408490615541435,-5.8408490615541435) q[12];
u3(pi,1.790707812546182,-1.790707812546182) q[13];
u3(pi,0.3675663404700058,-0.3675663404700058) q[14];
u3(pi,2.819893565862198,-2.819893565862198) q[15];
u3(pi/2,5.240804864718492,-5.240804864718492) q[16];
u3(pi/2,5.525433159133728,-5.525433159133728) q[16];
u3(pi,0.3097610356439536,-0.3097610356439536) q[17];
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
