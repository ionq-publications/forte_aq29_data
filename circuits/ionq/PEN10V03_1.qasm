OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[8];
u3(pi/2,5.693822525366141,-5.693822525366141) q[9];
rzz(-pi/2) q[8],q[9];
u3(pi/2,4.123026198571244,-4.123026198571244) q[9];
u3(pi/2,5.326256184896136,-5.326256184896136) q[9];
rzz(pi/2) q[8],q[9];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
rzz(pi/2) q[7],q[9];
u3(pi/2,2.184663531306342,-2.184663531306342) q[9];
u3(pi/2,2.9204245307770718,-2.9204245307770718) q[9];
rzz(pi/2) q[7],q[9];
u3(pi/2,0.7363893180014475,-0.7363893180014475) q[6];
rzz(pi/2) q[6],q[9];
u3(pi/2,6.062017184366865,-6.062017184366865) q[9];
u3(pi/2,1.4482742133048947,-1.4482742133048947) q[9];
rzz(-pi/2) q[6],q[9];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(pi/2) q[5],q[9];
u3(pi/2,4.589866866894688,-4.589866866894688) q[9];
u3(pi/2,4.785902248478691,-4.785902248478691) q[9];
rzz(pi/2) q[5],q[9];
u3(pi/2,4.516353598800687,-4.516353598800687) q[4];
rzz(pi/2) q[4],q[9];
u3(pi/2,1.6443095948888977,-1.6443095948888977) q[9];
u3(pi/2,4.393203166779967,-4.393203166779967) q[9];
rzz(-pi/2) q[4],q[9];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(pi/2) q[3],q[9];
u3(pi/2,1.2516105131901736,-1.2516105131901736) q[9];
u3(pi/2,3.6078050033825186,-3.6078050033825186) q[9];
rzz(pi/2) q[3],q[9];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,5.178601330177416,-5.178601330177416) q[9];
rzz(pi/2) q[2],q[9];
u3(pi/2,pi,-pi) q[1];
u3(pi/2,7*pi/4,-7*pi/4) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
u3(pi/2,5*pi/8,-5*pi/8) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,7*pi/8,-7*pi/8) q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,4.516353598800687,-4.516353598800687) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,6.087149925595583,-6.087149925595583) q[4];
u3(pi/2,2.8469112626830704,-2.8469112626830704) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,3.82897312619524,-3.82897312619524) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,0.7363893180014475,-0.7363893180014475) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,5.448778298386137,-5.448778298386137) q[6];
u3(pi/2,2.282681222098344,-2.282681222098344) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,1.3251237812841747,-1.3251237812841747) q[7];
u3(pi/2,4.454778382790327,-4.454778382790327) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[8];
rzz(-pi/2) q[0],q[8];
rzz(pi/2) q[0],q[8];
u3(pi/2,pi/4,-pi/4) q[1];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,5*pi/8,-5*pi/8) q[2];
u3(pi/2,11*pi/8,-11*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[3];
rzz(-pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,2.8469112626830704,-2.8469112626830704) q[4];
u3(pi/2,5.792468534688861,-5.792468534688861) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
u3(pi/2,2.110521944681623,-2.110521944681623) q[5];
rzz(-pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,5.424273875688137,-5.424273875688137) q[6];
u3(pi/2,2.233672376702343,-2.233672376702343) q[6];
rzz(pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u3(pi/2,4.454778382790327,-4.454778382790327) q[7];
u3(pi/2,1.288681306502533,-1.288681306502533) q[7];
rzz(pi/2) q[1],q[7];
rzz(pi/2) q[1],q[8];
u3(pi/2,5.246459731494954,-5.246459731494954) q[8];
u3(pi/2,2.092300707290802,-2.092300707290802) q[8];
rzz(pi/2) q[1],q[8];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,5.301751762198135,-5.301751762198135) q[3];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[3];
rzz(pi/2) q[2],q[3];
rzz(-pi/2) q[2],q[4];
u3(pi/2,2.6508758810990676,-2.6508758810990676) q[4];
u3(pi/2,5.399769452990137,-5.399769452990137) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,5.252114598271416,-5.252114598271416) q[5];
u3(pi/2,1.91448656309762,-1.91448656309762) q[5];
rzz(pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[6];
u3(pi/2,5.3752650302921365,-5.3752650302921365) q[6];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[6];
rzz(pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,1.288681306502533,-1.288681306502533) q[7];
u3(pi/2,4.381265114696325,-4.381265114696325) q[7];
rzz(pi/2) q[2],q[7];
rzz(pi/2) q[2],q[8];
u3(pi/2,2.092300707290802,-2.092300707290802) q[8];
u3(pi/2,5.209388938182594,-5.209388938182594) q[8];
rzz(-pi/2) q[2],q[8];
u3(pi/2,6.087149925595583,-6.087149925595583) q[3];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,5.399769452990137,-5.399769452990137) q[4];
u3(pi/2,1.472778636002895,-1.472778636002895) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,1.91448656309762,-1.91448656309762) q[5];
u3(pi/2,4.663380134988689,-4.663380134988689) q[5];
rzz(-pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[6];
u3(pi/2,5.080583639385413,-5.080583639385413) q[6];
rzz(pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,4.381265114696325,-4.381265114696325) q[7];
u3(pi/2,1.141026451783813,-1.141026451783813) q[7];
rzz(-pi/2) q[3],q[7];
rzz(pi/2) q[3],q[8];
u3(pi/2,2.0677962845928017,-2.0677962845928017) q[8];
u3(pi/2,5.160380092786594,-5.160380092786594) q[8];
rzz(pi/2) q[3],q[8];
u3(pi/2,6.185167616387585,-6.185167616387585) q[4];
u3(pi/2,3*pi/8,-3*pi/8) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[5];
u3(pi/2,3.8779819715912405,-3.8779819715912405) q[5];
rzz(-pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,5.080583639385413,-5.080583639385413) q[6];
u3(pi/2,1.5462919040968963,-1.5462919040968963) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,4.282619105373606,-4.282619105373606) q[7];
u3(pi/2,0.9449910701998098,-0.9449910701998098) q[7];
rzz(pi/2) q[4],q[7];
rzz(-pi/2) q[4],q[8];
u3(pi/2,2.018787439196801,-2.018787439196801) q[8];
u3(pi/2,5.0623624019945925,-5.0623624019945925) q[8];
rzz(pi/2) q[4],q[8];
u3(pi/2,2.3071856447963444,-2.3071856447963444) q[5];
u3(pi/2,4.025008507779242,-4.025008507779242) q[5];
rzz(pi/2) q[5],q[6];
u3(pi/2,1.5462919040968963,-1.5462919040968963) q[6];
u3(pi/2,3.902486394289241,-3.902486394289241) q[6];
rzz(pi/2) q[5],q[6];
rzz(-pi/2) q[5],q[7];
u3(pi/2,4.0865837237896026,-4.0865837237896026) q[7];
u3(pi/2,0.5522919885010856,-0.5522919885010856) q[7];
rzz(pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,5.0623624019945925,-5.0623624019945925) q[8];
u3(pi/2,1.7241060482900783,-1.7241060482900783) q[8];
rzz(pi/2) q[5],q[8];
u3(pi/2,2.3316900674943444,-2.3316900674943444) q[6];
u3(pi/2,5.399769452990137,-5.399769452990137) q[6];
rzz(pi/2) q[6],q[7];
u3(pi/2,0.5522919885010856,-0.5522919885010856) q[7];
u3(pi/2,2.90848647869343,-2.90848647869343) q[7];
rzz(-pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,1.7241060482900783,-1.7241060482900783) q[8];
u3(pi/2,4.472999620181147,-4.472999620181147) q[8];
rzz(pi/2) q[6],q[8];
u3(pi/2,1.337690151898534,-1.337690151898534) q[7];
u3(pi/2,2.159530790077624,-2.159530790077624) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,4.472999620181147,-4.472999620181147) q[8];
u3(pi/2,0.5460088031939061,-0.5460088031939061) q[8];
rzz(-pi/2) q[7],q[8];
u3(pi,3*pi/2,-3*pi/2) q[1];
u3(pi,pi/4,-pi/4) q[2];
u3(pi,7*pi/8,-7*pi/8) q[3];
u3(pi,4.516353598800687,-4.516353598800687) q[4];
u3(pi,0.19603538158400308,-0.19603538158400308) q[5];
u3(pi,5.841477380084861,-5.841477380084861) q[6];
u3(pi,1.595300749492897,-1.595300749492897) q[7];
u3(pi/2,2.1168051299888027,-2.1168051299888027) q[8];
u3(pi/2,2.491282974296706,-2.491282974296706) q[8];
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
