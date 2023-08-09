OPENQASM 2.0;
include "qelib1.inc";
qreg q[15];
creg c[15];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[13];
u3(pi/2,4.577300496280329,-4.577300496280329) q[14];
rzz(-pi/2) q[13],q[14];
u3(pi/2,3.006504169485432,-3.006504169485432) q[14];
u3(pi/2,5.413592460665932,-5.413592460665932) q[14];
rzz(pi/2) q[13],q[14];
u3(pi/2,3.666238626739289,-3.666238626739289) q[12];
rzz(pi/2) q[12],q[14];
u3(pi/2,5.413592460665932,-5.413592460665932) q[14];
u3(pi/2,0.8023627637268332,-0.8023627637268332) q[14];
rzz(pi/2) q[12],q[14];
u3(pi/2,4.072760716113808,-4.072760716113808) q[11];
rzz(pi/2) q[11],q[14];
u3(pi/2,0.8023627637268332,-0.8023627637268332) q[14];
u3(pi/2,1.0046813306180158,-1.0046813306180158) q[14];
rzz(-pi/2) q[11],q[14];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[10];
rzz(pi/2) q[10],q[14];
u3(pi/2,1.0046813306180158,-1.0046813306180158) q[14];
u3(pi/2,3.7416368504254436,-3.7416368504254436) q[14];
rzz(pi/2) q[10],q[14];
u3(pi/2,5.289185391583776,-5.289185391583776) q[9];
rzz(pi/2) q[9],q[14];
u3(pi/2,3.7416368504254436,-3.7416368504254436) q[14];
u3(pi/2,6.073326917919788,-6.073326917919788) q[14];
rzz(-pi/2) q[9],q[14];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[8];
rzz(pi/2) q[8],q[14];
u3(pi/2,2.931734264329995,-2.931734264329995) q[14];
u3(pi/2,4.452893427198172,-4.452893427198172) q[14];
rzz(pi/2) q[8],q[14];
u3(pi/2,3.068079385495792,-3.068079385495792) q[7];
rzz(pi/2) q[7],q[14];
u3(pi/2,1.3113007736083797,-1.3113007736083797) q[14];
u3(pi/2,1.409946782931099,-1.409946782931099) q[14];
rzz(pi/2) q[7],q[14];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
rzz(-pi/2) q[6],q[14];
u3(pi/2,1.409946782931099,-1.409946782931099) q[14];
u3(pi/2,4.354875736406171,-4.354875736406171) q[14];
rzz(pi/2) q[6],q[14];
u3(pi/2,4.516353598800687,-4.516353598800687) q[5];
rzz(pi/2) q[5],q[14];
u3(pi/2,4.354875736406171,-4.354875736406171) q[14];
u3(pi/2,0.820584001117654,-0.820584001117654) q[14];
rzz(pi/2) q[5],q[14];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[4],q[14];
u3(pi/2,3.9621766547074473,-3.9621766547074473) q[14];
u3(pi/2,0.035185837720205684,-0.035185837720205684) q[14];
rzz(pi/2) q[4],q[14];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,1.6059821645151022,-1.6059821645151022) q[14];
rzz(pi/2) q[3],q[14];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
u3(pi/2,7*pi/4,-7*pi/4) q[1];
rzz(-pi/2) q[0],q[1];
rzz(pi/2) q[0],q[2];
u3(pi/2,0,0) q[2];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
u3(pi/2,5.301751762198135,-5.301751762198135) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
u3(pi/2,5.399769452990137,-5.399769452990137) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,6.087149925595583,-6.087149925595583) q[5];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,5.301751762198135,-5.301751762198135) q[6];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,6.209672039085585,-6.209672039085585) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,4.638875712290688,-4.638875712290688) q[7];
u3(pi/2,1.4847166880865363,-1.4847166880865363) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[8];
rzz(-pi/2) q[0],q[8];
rzz(pi/2) q[0],q[8];
u3(pi/2,5.289185391583776,-5.289185391583776) q[9];
rzz(-pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[9];
u3(pi/2,3.669380219392878,-3.669380219392878) q[10];
rzz(-pi/2) q[0],q[10];
rzz(pi/2) q[0],q[10];
u3(pi/2,4.072760716113808,-4.072760716113808) q[11];
rzz(-pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[11];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[12];
rzz(-pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[12];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[13];
rzz(pi/2) q[0],q[13];
rzz(-pi/2) q[0],q[13];
u3(pi/2,pi/4,-pi/4) q[1];
u3(pi/2,pi,-pi) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
u3(pi/2,13*pi/8,-13*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(-pi/2) q[1],q[3];
u3(pi/2,2.159530790077624,-2.159530790077624) q[3];
u3(pi/2,4.9084243619686925,-4.9084243619686925) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,5.399769452990137,-5.399769452990137) q[4];
u3(pi/2,2.061513099285622,-2.061513099285622) q[4];
rzz(pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[5];
u3(pi/2,6.037512761668864,-6.037512761668864) q[5];
u3(pi/2,2.7979024172870695,-2.7979024172870695) q[5];
rzz(pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[6];
u3(pi/2,5.227610175573416,-5.227610175573416) q[6];
rzz(pi/2) q[1],q[6];
rzz(-pi/2) q[1],q[7];
u3(pi/2,4.626309341676329,-4.626309341676329) q[7];
u3(pi/2,1.4602122653885359,-1.4602122653885359) q[7];
rzz(pi/2) q[1],q[7];
rzz(pi/2) q[1],q[8];
u3(pi/2,2.1048670779051615,-2.1048670779051615) q[8];
u3(pi/2,5.233893360880595,-5.233893360880595) q[8];
rzz(pi/2) q[1],q[8];
rzz(pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[9];
rzz(pi/2) q[1],q[10];
rzz(pi/2) q[1],q[10];
rzz(pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[11];
rzz(pi/2) q[1],q[12];
rzz(pi/2) q[1],q[12];
rzz(pi/2) q[1],q[13];
rzz(pi/2) q[1],q[13];
u3(pi/2,pi/8,-pi/8) q[2];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(-pi/2) q[2],q[3];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[3];
u3(pi/2,4.123026198571244,-4.123026198571244) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,2.061513099285622,-2.061513099285622) q[4];
u3(pi/2,4.810406671176691,-4.810406671176691) q[4];
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
u3(pi/2,1.4602122653885359,-1.4602122653885359) q[7];
u3(pi/2,4.5527960735823285,-4.5527960735823285) q[7];
rzz(pi/2) q[2],q[7];
rzz(pi/2) q[2],q[8];
u3(pi/2,5.233893360880595,-5.233893360880595) q[8];
u3(pi/2,2.0677962845928017,-2.0677962845928017) q[8];
rzz(pi/2) q[2],q[8];
rzz(-pi/2) q[2],q[9];
u3(pi/2,3.70896428682811,-3.70896428682811) q[9];
u3(pi/2,0.5554335811546754,-0.5554335811546754) q[9];
rzz(pi/2) q[2],q[9];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[12];
rzz(pi/2) q[2],q[12];
rzz(pi/2) q[2],q[13];
rzz(-pi/2) q[2],q[13];
u3(pi/2,2.552229871776348,-2.552229871776348) q[3];
u3(pi/2,7*pi/8,-7*pi/8) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,4.810406671176691,-4.810406671176691) q[4];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,2.6018670357030667,-2.6018670357030667) q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
rzz(-pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,1.9879998311916212,-1.9879998311916212) q[6];
u3(pi/2,4.933557103197411,-4.933557103197411) q[6];
rzz(pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,4.5527960735823285,-4.5527960735823285) q[7];
u3(pi/2,1.3131857292005336,-1.3131857292005336) q[7];
rzz(pi/2) q[3],q[7];
rzz(-pi/2) q[3],q[8];
u3(pi/2,5.209388938182594,-5.209388938182594) q[8];
u3(pi/2,2.018787439196801,-2.018787439196801) q[8];
rzz(pi/2) q[3],q[8];
rzz(pi/2) q[3],q[9];
u3(pi/2,0.5554335811546754,-0.5554335811546754) q[9];
u3(pi/2,3.6725218120464684,-3.6725218120464684) q[9];
rzz(pi/2) q[3],q[9];
rzz(-pi/2) q[3],q[10];
u3(pi/2,5.23263672381916,-5.23263672381916) q[10];
u3(pi/2,2.078477699615007,-2.078477699615007) q[10];
rzz(pi/2) q[3],q[10];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[12];
rzz(pi/2) q[3],q[13];
rzz(pi/2) q[3],q[13];
u3(pi/2,5.595804834574139,-5.595804834574139) q[4];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[4];
rzz(-pi/2) q[4],q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
u3(pi/2,1.4237697906068942,-1.4237697906068942) q[5];
rzz(pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,4.933557103197411,-4.933557103197411) q[6];
u3(pi/2,1.399265367908894,-1.399265367908894) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,1.3131857292005336,-1.3131857292005336) q[7];
u3(pi/2,4.258114682675606,-4.258114682675606) q[7];
rzz(-pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,2.018787439196801,-2.018787439196801) q[8];
u3(pi/2,5.0623624019945925,-5.0623624019945925) q[8];
rzz(pi/2) q[4],q[8];
rzz(pi/2) q[4],q[9];
u3(pi/2,3.6725218120464684,-3.6725218120464684) q[9];
u3(pi/2,0.4819203130606743,-0.4819203130606743) q[9];
rzz(-pi/2) q[4],q[9];
rzz(pi/2) q[4],q[10];
u3(pi/2,2.078477699615007,-2.078477699615007) q[10];
u3(pi/2,5.1955659305068,-5.1955659305068) q[10];
rzz(pi/2) q[4],q[10];
rzz(pi/2) q[4],q[11];
u3(pi/2,5.634760583478653,-5.634760583478653) q[11];
u3(pi/2,2.4806015592745005,-2.4806015592745005) q[11];
rzz(pi/2) q[4],q[11];
rzz(-pi/2) q[4],q[12];
rzz(pi/2) q[4],q[12];
rzz(pi/2) q[4],q[13];
rzz(pi/2) q[4],q[13];
u3(pi/2,6.136158770991584,-6.136158770991584) q[5];
u3(pi/2,1.0800795543041708,-1.0800795543041708) q[5];
rzz(-pi/2) q[5],q[6];
u3(pi/2,4.540858021498687,-4.540858021498687) q[6];
u3(pi/2,0.6138672045114455,-0.6138672045114455) q[6];
rzz(pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,1.1165220290858124,-1.1165220290858124) q[7];
u3(pi/2,3.865415600976881,-3.865415600976881) q[7];
rzz(pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,5.0623624019945925,-5.0623624019945925) q[8];
u3(pi/2,1.7241060482900783,-1.7241060482900783) q[8];
rzz(pi/2) q[5],q[8];
rzz(pi/2) q[5],q[9];
u3(pi/2,3.6235129666504675,-3.6235129666504675) q[9];
u3(pi/2,0.38327430373795474,-0.38327430373795474) q[9];
rzz(pi/2) q[5],q[9];
rzz(pi/2) q[5],q[10];
u3(pi/2,5.1955659305068,-5.1955659305068) q[10];
u3(pi/2,2.004964431521006,-2.004964431521006) q[10];
rzz(-pi/2) q[5],q[10];
rzz(pi/2) q[5],q[11];
u3(pi/2,2.4806015592745005,-2.4806015592745005) q[11];
u3(pi/2,5.597689790166293,-5.597689790166293) q[11];
rzz(pi/2) q[5],q[11];
rzz(pi/2) q[5],q[12];
u3(pi/2,2.0872741590450583,-2.0872741590450583) q[12];
u3(pi/2,5.216300442020493,-5.216300442020493) q[12];
rzz(pi/2) q[5],q[12];
rzz(pi/2) q[5],q[13];
rzz(pi/2) q[5],q[13];
u3(pi/2,5.326256184896136,-5.326256184896136) q[6];
u3(pi/2,0.4417079270947249,-0.4417079270947249) q[6];
rzz(pi/2) q[6],q[7];
u3(pi/2,3.865415600976881,-3.865415600976881) q[7];
u3(pi/2,6.221610091169226,-6.221610091169226) q[7];
rzz(-pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,1.7241060482900783,-1.7241060482900783) q[8];
u3(pi/2,4.472999620181147,-4.472999620181147) q[8];
rzz(pi/2) q[6],q[8];
rzz(pi/2) q[6],q[9];
u3(pi/2,0.38327430373795474,-0.38327430373795474) q[9];
u3(pi/2,3.328831575743745,-3.328831575743745) q[9];
rzz(-pi/2) q[6],q[9];
rzz(pi/2) q[6],q[10];
u3(pi/2,5.1465570851108,-5.1465570851108) q[10];
u3(pi/2,1.9069467407290044,-1.9069467407290044) q[10];
rzz(pi/2) q[6],q[10];
rzz(pi/2) q[6],q[11];
u3(pi/2,5.597689790166293,-5.597689790166293) q[11];
u3(pi/2,2.4070882911804996,-2.4070882911804996) q[11];
rzz(pi/2) q[6],q[11];
rzz(-pi/2) q[6],q[12];
u3(pi/2,2.0747077884306995,-2.0747077884306995) q[12];
u3(pi/2,5.191796019322492,-5.191796019322492) q[12];
rzz(pi/2) q[6],q[12];
rzz(pi/2) q[6],q[13];
u3(pi/2,5.230751768227005,-5.230751768227005) q[13];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[13];
rzz(pi/2) q[6],q[13];
u3(pi/2,4.6508137643743295,-4.6508137643743295) q[7];
u3(pi/2,0.024504422698000385,-0.024504422698000385) q[7];
rzz(-pi/2) q[7],q[8];
u3(pi/2,1.3314069665913544,-1.3314069665913544) q[8];
u3(pi/2,3.687601456783699,-3.687601456783699) q[8];
rzz(pi/2) q[7],q[8];
rzz(pi/2) q[7],q[9];
u3(pi/2,0.18723892215395166,-0.18723892215395166) q[9];
u3(pi/2,2.9361324940450206,-2.9361324940450206) q[9];
rzz(pi/2) q[7],q[9];
rzz(pi/2) q[7],q[10];
u3(pi/2,1.9069467407290044,-1.9069467407290044) q[10];
u3(pi/2,4.851875694204076,-4.851875694204076) q[10];
rzz(-pi/2) q[7],q[10];
rzz(pi/2) q[7],q[11];
u3(pi/2,2.4070882911804996,-2.4070882911804996) q[11];
u3(pi/2,5.450663253978291,-5.450663253978291) q[11];
rzz(pi/2) q[7],q[11];
rzz(pi/2) q[7],q[12];
u3(pi/2,5.191796019322492,-5.191796019322492) q[12];
u3(pi/2,2.001194520336698,-2.001194520336698) q[12];
rzz(-pi/2) q[7],q[12];
rzz(pi/2) q[7],q[13];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[13];
u3(pi/2,5.194309293445364,-5.194309293445364) q[13];
rzz(pi/2) q[7],q[13];
u3(pi/2,5.258397783578595,-5.258397783578595) q[8];
u3(pi/2,5.301123443667417,-5.301123443667417) q[8];
rzz(pi/2) q[8],q[9];
u3(pi/2,2.9361324940450206,-2.9361324940450206) q[9];
u3(pi/2,5.2923269842373655,-5.2923269842373655) q[9];
rzz(-pi/2) q[8],q[9];
rzz(pi/2) q[8],q[10];
u3(pi/2,1.7102830406142833,-1.7102830406142833) q[10];
u3(pi/2,4.4591766125053525,-4.4591766125053525) q[10];
rzz(pi/2) q[8],q[10];
rzz(pi/2) q[8],q[11];
u3(pi/2,5.450663253978291,-5.450663253978291) q[11];
u3(pi/2,2.113035218804495,-2.113035218804495) q[11];
rzz(-pi/2) q[8],q[11];
rzz(pi/2) q[8],q[12];
u3(pi/2,5.142787173926491,-5.142787173926491) q[12];
u3(pi/2,1.9031768295446967,-1.9031768295446967) q[12];
rzz(pi/2) q[8],q[12];
rzz(pi/2) q[8],q[13];
u3(pi/2,5.194309293445364,-5.194309293445364) q[13];
u3(pi/2,2.00370779445957,-2.00370779445957) q[13];
rzz(pi/2) q[8],q[13];
u3(pi/2,0.5799380038526758,-0.5799380038526758) q[9];
u3(pi/2,1.3439733372057134,-1.3439733372057134) q[9];
rzz(-pi/2) q[9],q[10];
u3(pi/2,1.3175839589155591,-1.3175839589155591) q[10];
u3(pi/2,3.673778449107904,-3.673778449107904) q[10];
rzz(pi/2) q[9],q[10];
rzz(pi/2) q[9],q[11];
u3(pi/2,5.2546278723942885,-5.2546278723942885) q[11];
u3(pi/2,1.7203361371057706,-1.7203361371057706) q[11];
rzz(pi/2) q[9],q[11];
rzz(-pi/2) q[9],q[12];
u3(pi/2,5.044769483134489,-5.044769483134489) q[12];
u3(pi/2,1.7065131294299756,-1.7065131294299756) q[12];
rzz(pi/2) q[9],q[12];
rzz(pi/2) q[9],q[13];
u3(pi/2,2.00370779445957,-2.00370779445957) q[13];
u3(pi/2,5.047282757257362,-5.047282757257362) q[13];
rzz(pi/2) q[9],q[13];
u3(pi/2,2.1029821223130076,-2.1029821223130076) q[10];
u3(pi/2,3.2703979523869746,-3.2703979523869746) q[10];
rzz(pi/2) q[10],q[11];
u3(pi/2,1.7203361371057706,-1.7203361371057706) q[11];
u3(pi/2,4.076530627298116,-4.076530627298116) q[11];
rzz(pi/2) q[10],q[11];
rzz(pi/2) q[10],q[12];
u3(pi/2,1.7065131294299756,-1.7065131294299756) q[12];
u3(pi/2,4.455406701321044,-4.455406701321044) q[12];
rzz(pi/2) q[10],q[12];
rzz(pi/2) q[10],q[13];
u3(pi/2,5.047282757257362,-5.047282757257362) q[13];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[13];
rzz(-pi/2) q[10],q[13];
u3(pi/2,5.647326954093012,-5.647326954093012) q[11];
u3(pi/2,1.136628222068787,-1.136628222068787) q[11];
rzz(pi/2) q[11],q[12];
u3(pi/2,4.455406701321044,-4.455406701321044) q[12];
u3(pi/2,0.5284158843338032,-0.5284158843338032) q[12];
rzz(pi/2) q[11],q[12];
rzz(pi/2) q[11],q[13];
u3(pi/2,4.850619057142641,-4.850619057142641) q[13];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[13];
rzz(-pi/2) q[11],q[13];
u3(pi/2,5.240804864718492,-5.240804864718492) q[12];
u3(pi/2,1.9986812462138264,-1.9986812462138264) q[12];
rzz(pi/2) q[12],q[13];
u3(pi/2,4.457919975443916,-4.457919975443916) q[13];
u3(pi/2,0.5309291584566751,-0.5309291584566751) q[13];
rzz(pi/2) q[12],q[13];
u3(pi,3*pi/2,-3*pi/2) q[1];
u3(pi,5*pi/4,-5*pi/4) q[2];
u3(pi,pi/4,-pi/4) q[3];
u3(pi,3*pi/2,-3*pi/2) q[4];
u3(pi,0.9814335449814514,-0.9814335449814514) q[5];
u3(pi,15*pi/8,-15*pi/8) q[6];
u3(pi,9*pi/8,-9*pi/8) q[7];
u3(pi,2.3436281195779856,-2.3436281195779856) q[8];
u3(pi,2.2211060060879837,-2.2211060060879837) q[9];
u3(pi,1.3527697966357648,-1.3527697966357648) q[10];
u3(pi,3.3307165313358986,-3.3307165313358986) q[11];
u3(pi,3.878610290121958,-3.878610290121958) q[12];
u3(pi/2,2.101725485251572,-2.101725485251572) q[13];
u3(pi/2,2.937389131106457,-2.937389131106457) q[13];
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
