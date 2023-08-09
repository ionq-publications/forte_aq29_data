OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(pi/2,1.7674600269096177,-1.7674600269096177) q[6];
u3(pi/2,2.5528581903070657,-2.5528581903070657) q[6];
rzz(-pi/2) q[5],q[6];
rzz(-pi/2) q[4],q[6];
u3(pi/2,2.552229871776348,-2.552229871776348) q[6];
u3(pi/2,5.301751762198135,-5.301751762198135) q[6];
rzz(-pi/2) q[4],q[6];
rzz(-pi/2) q[4],q[5];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[5];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[5];
rzz(-pi/2) q[4],q[5];
rzz(-pi/2) q[3],q[6];
u3(pi/2,5.301123443667417,-5.301123443667417) q[6];
u3(pi/2,5*pi/8,-5*pi/8) q[6];
rzz(-pi/2) q[3],q[6];
rzz(-pi/2) q[3],q[5];
u3(pi/2,4.9084243619686925,-4.9084243619686925) q[5];
u3(pi/2,1.3741326266801754,-1.3741326266801754) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[4];
u3(pi/2,5*pi/8,-5*pi/8) q[4];
u3(pi/2,7*pi/8,-7*pi/8) q[4];
rzz(-pi/2) q[3],q[4];
rzz(-pi/2) q[2],q[6];
u3(pi/2,13*pi/8,-13*pi/8) q[6];
u3(pi/2,1.865477717701619,-1.865477717701619) q[6];
rzz(pi/2) q[2],q[6];
rzz(-pi/2) q[2],q[5];
u3(pi/2,1.3741326266801754,-1.3741326266801754) q[5];
u3(pi/2,11*pi/8,-11*pi/8) q[5];
rzz(-pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[4];
u3(pi/2,7*pi/8,-7*pi/8) q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
rzz(-pi/2) q[2],q[4];
rzz(-pi/2) q[2],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
u3(pi/2,pi,-pi) q[3];
rzz(-pi/2) q[2],q[3];
rzz(-pi/2) q[1],q[6];
u3(pi/2,5.007070371291412,-5.007070371291412) q[6];
u3(pi/2,1.8164688723056186,-1.8164688723056186) q[6];
rzz(-pi/2) q[1],q[6];
rzz(pi/2) q[1],q[5];
u3(pi/2,3*pi/8,-3*pi/8) q[5];
u3(pi/2,4.221672207893964,-4.221672207893964) q[5];
rzz(-pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[4];
rzz(-pi/2) q[1],q[4];
rzz(pi/2) q[1],q[3];
u3(pi/2,pi,-pi) q[3];
u3(pi/2,15*pi/8,-15*pi/8) q[3];
rzz(-pi/2) q[1],q[3];
rzz(pi/2) q[1],q[2];
u3(pi/2,0,0) q[2];
u3(pi/2,pi/4,-pi/4) q[2];
rzz(-pi/2) q[1],q[2];
rzz(pi/2) q[0],q[6];
rzz(-pi/2) q[0],q[5];
rzz(-pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[3];
rzz(-pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi,-pi) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
u3(pi/2,0,0) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,15*pi/8,-15*pi/8) q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[4];
u3(pi/2,13*pi/8,-13*pi/8) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,4.221672207893964,-4.221672207893964) q[5];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,4.9580615258954115,-4.9580615258954115) q[6];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,0,0) q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(-pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
u3(pi/2,13*pi/8,-13*pi/8) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,5*pi/8,-5*pi/8) q[4];
u3(pi/2,4.9084243619686925,-4.9084243619686925) q[4];
rzz(pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[5];
u3(pi/2,4.123026198571244,-4.123026198571244) q[5];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[5];
rzz(pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[6];
u3(pi/2,4.8594155165726916,-4.8594155165726916) q[6];
rzz(pi/2) q[1],q[6];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,pi,-pi) q[2];
rzz(-pi/2) q[2],q[3];
u3(pi/2,5*pi/8,-5*pi/8) q[3];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,4.9084243619686925,-4.9084243619686925) q[4];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[4];
rzz(pi/2) q[2],q[4];
rzz(-pi/2) q[2],q[5];
u3(pi/2,4.025008507779242,-4.025008507779242) q[5];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(pi/2) q[2],q[5];
rzz(pi/2) q[2],q[6];
u3(pi/2,4.8594155165726916,-4.8594155165726916) q[6];
u3(pi/2,1.6198051721908973,-1.6198051721908973) q[6];
rzz(pi/2) q[2],q[6];
u3(pi/2,7*pi/8,-7*pi/8) q[3];
u3(pi/2,pi,-pi) q[3];
rzz(-pi/2) q[3],q[4];
u3(pi/2,4.516353598800687,-4.516353598800687) q[4];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
u3(pi/2,3.436274044496516,-3.436274044496516) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,1.6198051721908973,-1.6198051721908973) q[6];
u3(pi/2,4.565362444196688,-4.565362444196688) q[6];
rzz(-pi/2) q[3],q[6];
u3(pi/2,5.301751762198135,-5.301751762198135) q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,3.436274044496516,-3.436274044496516) q[5];
u3(pi/2,5.792468534688861,-5.792468534688861) q[5];
rzz(pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,1.4237697906068942,-1.4237697906068942) q[6];
u3(pi/2,4.172663362497963,-4.172663362497963) q[6];
rzz(-pi/2) q[4],q[6];
u3(pi/2,1.0800795543041708,-1.0800795543041708) q[5];
u3(pi/2,4.123026198571244,-4.123026198571244) q[5];
rzz(pi/2) q[5],q[6];
u3(pi/2,1.03107070890817,-1.03107070890817) q[6];
u3(pi/2,3.387265199100515,-3.387265199100515) q[6];
rzz(pi/2) q[5],q[6];
u3(pi,0,0) q[0];
u3(pi,0,0) q[1];
u3(pi,pi/2,-pi/2) q[2];
u3(pi,3*pi/4,-3*pi/4) q[3];
u3(pi,3*pi/8,-3*pi/8) q[4];
u3(pi,11*pi/8,-11*pi/8) q[5];
u3(pi/2,1.8164688723056186,-1.8164688723056186) q[6];
u3(pi/2,1.865477717701619,-1.865477717701619) q[6];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
