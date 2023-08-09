OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[11];
u3(pi/2,4.368698744081967,-4.368698744081967) q[12];
rzz(pi/2) q[11],q[12];
u3(pi/2,5.9394950708768635,-5.9394950708768635) q[12];
u3(pi/2,0.5705132258919065,-0.5705132258919065) q[12];
rzz(pi/2) q[11],q[12];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[10];
rzz(-pi/2) q[10],q[12];
u3(pi/2,0.5705132258919065,-0.5705132258919065) q[12];
u3(pi/2,1.8836989550924401,-1.8836989550924401) q[12];
rzz(pi/2) q[10],q[12];
u3(pi/2,4.69982260977033,-4.69982260977033) q[9];
rzz(pi/2) q[9],q[12];
u3(pi/2,5.025291608682233,-5.025291608682233) q[12];
u3(pi/2,5.54051280387096,-5.54051280387096) q[12];
rzz(pi/2) q[9],q[12];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[8];
rzz(-pi/2) q[8],q[12];
u3(pi/2,5.54051280387096,-5.54051280387096) q[12];
u3(pi/2,1.368477759903714,-1.368477759903714) q[12];
rzz(pi/2) q[8],q[12];
u3(pi/2,0.7118848953034471,-0.7118848953034471) q[7];
rzz(pi/2) q[7],q[12];
u3(pi/2,1.368477759903714,-1.368477759903714) q[12];
u3(pi/2,2.4479289956771666,-2.4479289956771666) q[12];
rzz(pi/2) q[7],q[12];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
rzz(pi/2) q[6],q[12];
u3(pi/2,5.58952164926696,-5.58952164926696) q[12];
u3(pi/2,0.28839820559954304,-0.28839820559954304) q[12];
rzz(-pi/2) q[6],q[12];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[5];
rzz(pi/2) q[5],q[12];
u3(pi/2,0.28839820559954304,-0.28839820559954304) q[12];
u3(pi/2,1.4664954506957153,-1.4664954506957153) q[12];
rzz(pi/2) q[5],q[12];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[4],q[12];
u3(pi/2,4.608088104285509,-4.608088104285509) q[12];
u3(pi/2,5.393486267682957,-5.393486267682957) q[12];
rzz(-pi/2) q[4],q[12];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,0.6810972872982671,-0.6810972872982671) q[12];
rzz(pi/2) q[3],q[12];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
u3(pi/2,7*pi/4,-7*pi/4) q[1];
rzz(-pi/2) q[0],q[1];
rzz(-pi/2) q[0],q[2];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,15*pi/8,-15*pi/8) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,4.516353598800687,-4.516353598800687) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,6.087149925595583,-6.087149925595583) q[5];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,5.301751762198135,-5.301751762198135) q[6];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,3.85347754889324,-3.85347754889324) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,5.424273875688137,-5.424273875688137) q[7];
u3(pi/2,2.2701148514839846,-2.2701148514839846) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[8];
rzz(-pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[8];
u3(pi/2,1.5582299561805373,-1.5582299561805373) q[9];
rzz(pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[9];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[10];
rzz(-pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[10];
u3(pi/2,3.669380219392878,-3.669380219392878) q[11];
rzz(-pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[11];
u3(pi/2,pi/4,-pi/4) q[1];
u3(pi/2,pi,-pi) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,15*pi/8,-15*pi/8) q[2];
u3(pi/2,5*pi/8,-5*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[3];
u3(pi/2,4.9084243619686925,-4.9084243619686925) q[3];
rzz(-pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[4];
u3(pi/2,5.203105752875415,-5.203105752875415) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[5];
u3(pi/2,5.9394950708768635,-5.9394950708768635) q[5];
rzz(-pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[6];
u3(pi/2,5.227610175573416,-5.227610175573416) q[6];
rzz(pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u3(pi/2,2.2701148514839846,-2.2701148514839846) q[7];
u3(pi/2,5.387203082375778,-5.387203082375778) q[7];
rzz(pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[8];
u3(pi/2,5.246459731494954,-5.246459731494954) q[8];
u3(pi/2,2.092300707290802,-2.092300707290802) q[8];
rzz(pi/2) q[1],q[8];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[10];
rzz(pi/2) q[1],q[10];
rzz(pi/2) q[1],q[11];
rzz(pi/2) q[1],q[11];
u3(pi/2,9*pi/8,-9*pi/8) q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(-pi/2) q[2],q[3];
u3(pi/2,4.9084243619686925,-4.9084243619686925) q[3];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,5.203105752875415,-5.203105752875415) q[4];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[4];
rzz(pi/2) q[2],q[4];
rzz(-pi/2) q[2],q[5];
u3(pi/2,5.9394950708768635,-5.9394950708768635) q[5];
u3(pi/2,2.6018670357030667,-2.6018670357030667) q[5];
rzz(pi/2) q[2],q[5];
rzz(pi/2) q[2],q[6];
u3(pi/2,5.227610175573416,-5.227610175573416) q[6];
u3(pi/2,1.9879998311916212,-1.9879998311916212) q[6];
rzz(pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,5.387203082375778,-5.387203082375778) q[7];
u3(pi/2,2.1966015833899837,-2.1966015833899837) q[7];
rzz(-pi/2) q[2],q[7];
rzz(pi/2) q[2],q[8];
u3(pi/2,2.092300707290802,-2.092300707290802) q[8];
u3(pi/2,5.209388938182594,-5.209388938182594) q[8];
rzz(pi/2) q[2],q[8];
rzz(pi/2) q[2],q[9];
u3(pi/2,3.1202298235453823,-3.1202298235453823) q[9];
u3(pi/2,6.249256106520817,-6.249256106520817) q[9];
rzz(-pi/2) q[2],q[9];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[11];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
u3(pi/2,15*pi/8,-15*pi/8) q[3];
rzz(-pi/2) q[3],q[4];
u3(pi/2,4.810406671176691,-4.810406671176691) q[4];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,2.6018670357030667,-2.6018670357030667) q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
rzz(pi/2) q[3],q[5];
rzz(-pi/2) q[3],q[6];
u3(pi/2,5.129592484781415,-5.129592484781415) q[6];
u3(pi/2,1.7919644496076181,-1.7919644496076181) q[6];
rzz(pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,5.338194236979777,-5.338194236979777) q[7];
u3(pi/2,2.098583892597982,-2.098583892597982) q[7];
rzz(pi/2) q[3],q[7];
rzz(pi/2) q[3],q[8];
u3(pi/2,5.209388938182594,-5.209388938182594) q[8];
u3(pi/2,2.018787439196801,-2.018787439196801) q[8];
rzz(pi/2) q[3],q[8];
rzz(pi/2) q[3],q[9];
u3(pi/2,3.1076634529310234,-3.1076634529310234) q[9];
u3(pi/2,6.2247516838228165,-6.2247516838228165) q[9];
rzz(pi/2) q[3],q[9];
rzz(-pi/2) q[3],q[10];
u3(pi/2,2.091044070229366,-2.091044070229366) q[10];
u3(pi/2,5.2200703532048,-5.2200703532048) q[10];
rzz(pi/2) q[3],q[10];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[11];
u3(pi/2,2.4542121809843462,-2.4542121809843462) q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
u3(pi/2,1.4237697906068942,-1.4237697906068942) q[5];
rzz(pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,1.7919644496076181,-1.7919644496076181) q[6];
u3(pi/2,4.540858021498687,-4.540858021498687) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,2.098583892597982,-2.098583892597982) q[7];
u3(pi/2,5.043512846073054,-5.043512846073054) q[7];
rzz(-pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,2.018787439196801,-2.018787439196801) q[8];
u3(pi/2,5.0623624019945925,-5.0623624019945925) q[8];
rzz(pi/2) q[4],q[8];
rzz(pi/2) q[4],q[9];
u3(pi/2,6.2247516838228165,-6.2247516838228165) q[9];
u3(pi/2,3.034150184837022,-3.034150184837022) q[9];
rzz(-pi/2) q[4],q[9];
rzz(pi/2) q[4],q[10];
u3(pi/2,5.2200703532048,-5.2200703532048) q[10];
u3(pi/2,2.0539732769170067,-2.0539732769170067) q[10];
rzz(pi/2) q[4],q[10];
rzz(pi/2) q[4],q[11];
u3(pi/2,2.0897874331679303,-2.0897874331679303) q[11];
u3(pi/2,5.218813716143364,-5.218813716143364) q[11];
rzz(pi/2) q[4],q[11];
u3(pi/2,2.994566117401791,-2.994566117401791) q[5];
u3(pi/2,3.3376280351737964,-3.3376280351737964) q[5];
rzz(-pi/2) q[5],q[6];
u3(pi/2,1.399265367908894,-1.399265367908894) q[6];
u3(pi/2,3.7554598581012386,-3.7554598581012386) q[6];
rzz(pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,1.901920192483261,-1.901920192483261) q[7];
u3(pi/2,4.6508137643743295,-4.6508137643743295) q[7];
rzz(pi/2) q[5],q[7];
rzz(-pi/2) q[5],q[8];
u3(pi/2,1.9207697484047996,-1.9207697484047996) q[8];
u3(pi/2,4.865698701879872,-4.865698701879872) q[8];
rzz(pi/2) q[5],q[8];
rzz(pi/2) q[5],q[9];
u3(pi/2,6.175742838426816,-6.175742838426816) q[9];
u3(pi/2,2.9361324940450206,-2.9361324940450206) q[9];
rzz(pi/2) q[5],q[9];
rzz(pi/2) q[5],q[10];
u3(pi/2,2.0539732769170067,-2.0539732769170067) q[10];
u3(pi/2,5.1465570851108,-5.1465570851108) q[10];
rzz(-pi/2) q[5],q[10];
rzz(pi/2) q[5],q[11];
u3(pi/2,5.218813716143364,-5.218813716143364) q[11];
u3(pi/2,2.052716639855571,-2.052716639855571) q[11];
rzz(pi/2) q[5],q[11];
u3(pi/2,5.326256184896136,-5.326256184896136) q[6];
u3(pi/2,pi/2,-pi/2) q[6];
rzz(-pi/2) q[6],q[7];
u3(pi/2,1.5092211107845366,-1.5092211107845366) q[7];
u3(pi/2,3.865415600976881,-3.865415600976881) q[7];
rzz(pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,4.865698701879872,-4.865698701879872) q[8];
u3(pi/2,1.3314069665913544,-1.3314069665913544) q[8];
rzz(pi/2) q[6],q[8];
rzz(pi/2) q[6],q[9];
u3(pi/2,2.9361324940450206,-2.9361324940450206) q[9];
u3(pi/2,5.881061447520093,-5.881061447520093) q[9];
rzz(-pi/2) q[6],q[9];
rzz(pi/2) q[6],q[10];
u3(pi/2,2.004964431521006,-2.004964431521006) q[10];
u3(pi/2,5.048539394318797,-5.048539394318797) q[10];
rzz(pi/2) q[6],q[10];
rzz(pi/2) q[6],q[11];
u3(pi/2,2.052716639855571,-2.052716639855571) q[11];
u3(pi/2,5.145300448049363,-5.145300448049363) q[11];
rzz(-pi/2) q[6],q[11];
u3(pi/2,5.436211927771778,-5.436211927771778) q[7];
u3(pi/2,5.914990648178863,-5.914990648178863) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,1.3314069665913544,-1.3314069665913544) q[8];
u3(pi/2,3.687601456783699,-3.687601456783699) q[8];
rzz(pi/2) q[7],q[8];
rzz(pi/2) q[7],q[9];
u3(pi/2,2.7394687939302997,-2.7394687939302997) q[9];
u3(pi/2,5.488362365821369,-5.488362365821369) q[9];
rzz(pi/2) q[7],q[9];
rzz(-pi/2) q[7],q[10];
u3(pi/2,1.9069467407290044,-1.9069467407290044) q[10];
u3(pi/2,4.851875694204076,-4.851875694204076) q[10];
rzz(pi/2) q[7],q[10];
rzz(pi/2) q[7],q[11];
u3(pi/2,2.00370779445957,-2.00370779445957) q[11];
u3(pi/2,5.047282757257362,-5.047282757257362) q[11];
rzz(pi/2) q[7],q[11];
u3(pi/2,2.1168051299888027,-2.1168051299888027) q[8];
u3(pi/2,2.6628139331827088,-2.6628139331827088) q[8];
rzz(-pi/2) q[8],q[9];
u3(pi/2,2.3467697122315756,-2.3467697122315756) q[9];
u3(pi/2,4.702964202423921,-4.702964202423921) q[9];
rzz(pi/2) q[8],q[9];
rzz(pi/2) q[8],q[10];
u3(pi/2,4.851875694204076,-4.851875694204076) q[10];
u3(pi/2,1.3175839589155591,-1.3175839589155591) q[10];
rzz(pi/2) q[8],q[10];
rzz(pi/2) q[8],q[11];
u3(pi/2,5.047282757257362,-5.047282757257362) q[11];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[11];
rzz(pi/2) q[8],q[11];
u3(pi/2,6.273760529218817,-6.273760529218817) q[9];
u3(pi/2,2.0740794698999814,-2.0740794698999814) q[9];
rzz(pi/2) q[9],q[10];
u3(pi/2,1.3175839589155591,-1.3175839589155591) q[10];
u3(pi/2,3.673778449107904,-3.673778449107904) q[10];
rzz(pi/2) q[9],q[10];
rzz(pi/2) q[9],q[11];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[11];
u3(pi/2,4.457919975443916,-4.457919975443916) q[11];
rzz(-pi/2) q[9],q[11];
u3(pi/2,5.244574775902801,-5.244574775902801) q[10];
u3(pi/2,5.5009287364357276,-5.5009287364357276) q[10];
rzz(pi/2) q[10],q[11];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[11];
u3(pi/2,3.6725218120464684,-3.6725218120464684) q[11];
rzz(pi/2) q[10],q[11];
u3(pi,pi/2,-pi/2) q[1];
u3(pi,3*pi/4,-3*pi/4) q[2];
u3(pi,pi/4,-pi/4) q[3];
u3(pi,2.94555727200579,-2.94555727200579) q[4];
u3(pi,1.6688140175868982,-1.6688140175868982) q[5];
u3(pi,1.5217874813988959,-1.5217874813988959) q[6];
u3(pi,15*pi/8,-15*pi/8) q[7];
u3(pi,4.221672207893964,-4.221672207893964) q[8];
u3(pi,3.638592611387698,-3.638592611387698) q[9];
u3(pi,pi,-pi) q[10];
u3(pi/2,2.101725485251572,-2.101725485251572) q[11];
u3(pi/2,2.759574986913274,-2.759574986913274) q[11];
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
