OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
u3(pi/2,4.221672207893964,-4.221672207893964) q[7];
rzz(-pi/2) q[6],q[7];
u3(pi/2,2.6508758810990676,-2.6508758810990676) q[7];
u3(pi/2,4.368698744081967,-4.368698744081967) q[7];
rzz(pi/2) q[6],q[7];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(pi/2) q[5],q[7];
u3(pi/2,4.368698744081967,-4.368698744081967) q[7];
u3(pi/2,4.663380134988689,-4.663380134988689) q[7];
rzz(pi/2) q[5],q[7];
u3(pi/2,4.9084243619686925,-4.9084243619686925) q[4];
rzz(pi/2) q[4],q[7];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[7];
u3(pi/2,4.074017353175243,-4.074017353175243) q[7];
rzz(-pi/2) q[4],q[7];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(pi/2) q[3],q[7];
u3(pi/2,0.9324246995854506,-0.9324246995854506) q[7];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[7];
rzz(pi/2) q[3],q[7];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(pi/2) q[2],q[7];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[7];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[7];
rzz(-pi/2) q[2],q[7];
u3(pi/2,pi/2,-pi/2) q[1];
u3(pi/2,5.252114598271416,-5.252114598271416) q[7];
rzz(pi/2) q[1],q[7];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,pi,-pi) q[1];
u3(pi/2,7*pi/4,-7*pi/4) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,9*pi/8,-9*pi/8) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,7*pi/8,-7*pi/8) q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,4.9084243619686925,-4.9084243619686925) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,3.338256353704514,-3.338256353704514) q[4];
u3(pi/2,0.09801769079200154,-0.09801769079200154) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,3.82897312619524,-3.82897312619524) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,5.301751762198135,-5.301751762198135) q[6];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,pi/4,-pi/4) q[1];
u3(pi/2,pi,-pi) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,pi/8,-pi/8) q[2];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[3];
rzz(-pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,3.2396103443817945,-3.2396103443817945) q[4];
u3(pi/2,6.185167616387585,-6.185167616387585) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
u3(pi/2,2.110521944681623,-2.110521944681623) q[5];
rzz(pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,5.276619020969417,-5.276619020969417) q[6];
u3(pi/2,2.0860175219836226,-2.0860175219836226) q[6];
rzz(pi/2) q[1],q[6];
u3(pi/2,11*pi/8,-11*pi/8) q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,5.301751762198135,-5.301751762198135) q[3];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[3];
rzz(pi/2) q[2],q[3];
rzz(-pi/2) q[2],q[4];
u3(pi/2,3.0435749627977917,-3.0435749627977917) q[4];
u3(pi/2,5.792468534688861,-5.792468534688861) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,2.110521944681623,-2.110521944681623) q[5];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[5];
rzz(pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[6];
u3(pi/2,5.227610175573416,-5.227610175573416) q[6];
u3(pi/2,1.9879998311916212,-1.9879998311916212) q[6];
rzz(pi/2) q[2],q[6];
u3(pi/2,6.087149925595583,-6.087149925595583) q[3];
u3(pi/2,pi/8,-pi/8) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,5.792468534688861,-5.792468534688861) q[4];
u3(pi/2,1.865477717701619,-1.865477717701619) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[5];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[5];
rzz(-pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,1.9879998311916212,-1.9879998311916212) q[6];
u3(pi/2,4.933557103197411,-4.933557103197411) q[6];
rzz(pi/2) q[3],q[6];
u3(pi/2,0.2946813909067226,-0.2946813909067226) q[4];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,4.663380134988689,-4.663380134988689) q[5];
u3(pi/2,0.7363893180014475,-0.7363893180014475) q[5];
rzz(-pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,4.933557103197411,-4.933557103197411) q[6];
u3(pi/2,1.399265367908894,-1.399265367908894) q[6];
rzz(pi/2) q[4],q[6];
u3(pi/2,5.448778298386137,-5.448778298386137) q[5];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[5];
rzz(pi/2) q[5],q[6];
u3(pi/2,1.399265367908894,-1.399265367908894) q[6];
u3(pi/2,3.7554598581012386,-3.7554598581012386) q[6];
rzz(-pi/2) q[5],q[6];
u3(pi,pi/2,-pi/2) q[1];
u3(pi,0,0) q[2];
u3(pi,pi/8,-pi/8) q[3];
u3(pi,2.1601591086083416,-2.1601591086083416) q[4];
u3(pi,2.1601591086083416,-2.1601591086083416) q[5];
u3(pi/2,2.184663531306342,-2.184663531306342) q[6];
u3(pi/2,2.3071856447963444,-2.3071856447963444) q[6];
u3(pi,3*pi/4,-3*pi/4) q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
