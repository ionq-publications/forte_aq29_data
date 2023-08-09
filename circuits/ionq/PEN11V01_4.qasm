OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
u3(pi/2,1.1290883997001717,-1.1290883997001717) q[10];
rzz(-pi/2) q[9],q[10];
u3(pi/2,5.841477380084861,-5.841477380084861) q[10];
u3(pi/2,2.473061736905885,-2.473061736905885) q[10];
rzz(pi/2) q[9],q[10];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[8];
rzz(pi/2) q[8],q[10];
u3(pi/2,2.473061736905885,-2.473061736905885) q[10];
u3(pi/2,5.160380092786594,-5.160380092786594) q[10];
rzz(pi/2) q[8],q[10];
u3(pi/2,1.0555751316061706,-1.0555751316061706) q[7];
rzz(pi/2) q[7],q[10];
u3(pi/2,5.160380092786594,-5.160380092786594) q[10];
u3(pi/2,1.110867162309351,-1.110867162309351) q[10];
rzz(-pi/2) q[7],q[10];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
rzz(pi/2) q[6],q[10];
u3(pi/2,4.252459815899144,-4.252459815899144) q[10];
u3(pi/2,5.577583597183319,-5.577583597183319) q[10];
rzz(pi/2) q[6],q[10];
u3(pi/2,7*pi/8,-7*pi/8) q[5];
rzz(pi/2) q[5],q[10];
u3(pi/2,2.4359909435935254,-2.4359909435935254) q[10];
u3(pi/2,2.9267077160842514,-2.9267077160842514) q[10];
rzz(-pi/2) q[5],q[10];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[4],q[10];
u3(pi/2,2.9267077160842514,-2.9267077160842514) q[10];
u3(pi/2,5.086866824692593,-5.086866824692593) q[10];
rzz(pi/2) q[4],q[10];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[3],q[10];
u3(pi/2,5.086866824692593,-5.086866824692593) q[10];
u3(pi/2,6.264964069788765,-6.264964069788765) q[10];
rzz(pi/2) q[3],q[10];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(-pi/2) q[2],q[10];
u3(pi/2,6.264964069788765,-6.264964069788765) q[10];
u3(pi/2,0.7671769260066275,-0.7671769260066275) q[10];
rzz(pi/2) q[2],q[10];
u3(pi/2,0,0) q[1];
u3(pi/2,5.479565906391317,-5.479565906391317) q[10];
rzz(pi/2) q[1],q[10];
u3(pi/2,pi,-pi) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,pi/2,-pi/2) q[1];
u3(pi/2,5*pi/4,-5*pi/4) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,0,0) q[2];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,pi,-pi) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,4.515725280269969,-4.515725280269969) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
u3(pi/2,5.399769452990137,-5.399769452990137) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,7*pi/8,-7*pi/8) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,3*pi/8,-3*pi/8) q[5];
u3(pi/2,4.270681053289964,-4.270681053289964) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,2.159530790077624,-2.159530790077624) q[6];
u3(pi/2,5.276619020969417,-5.276619020969417) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,1.0555751316061706,-1.0555751316061706) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,5.76796411199086,-5.76796411199086) q[7];
u3(pi/2,2.613805087786708,-2.613805087786708) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[8];
rzz(-pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[8];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
rzz(-pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[9];
u3(pi/2,3*pi/4,-3*pi/4) q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
rzz(-pi/2) q[1],q[2];
u3(pi/2,15*pi/8,-15*pi/8) q[2];
u3(pi/2,5*pi/8,-5*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[3];
u3(pi/2,4.123026198571244,-4.123026198571244) q[3];
rzz(pi/2) q[1],q[3];
rzz(-pi/2) q[1],q[4];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[4];
u3(pi/2,5.203105752875415,-5.203105752875415) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,1.1290883997001717,-1.1290883997001717) q[5];
u3(pi/2,4.172663362497963,-4.172663362497963) q[5];
rzz(pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,5.276619020969417,-5.276619020969417) q[6];
u3(pi/2,2.0860175219836226,-2.0860175219836226) q[6];
rzz(-pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u3(pi/2,2.613805087786708,-2.613805087786708) q[7];
u3(pi/2,5.7308933186785005,-5.7308933186785005) q[7];
rzz(pi/2) q[1],q[7];
rzz(pi/2) q[1],q[8];
u3(pi/2,5.246459731494954,-5.246459731494954) q[8];
u3(pi/2,2.092300707290802,-2.092300707290802) q[8];
rzz(-pi/2) q[1],q[8];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[9];
u3(pi/2,9*pi/8,-9*pi/8) q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,4.123026198571244,-4.123026198571244) q[3];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[3];
rzz(-pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,5.203105752875415,-5.203105752875415) q[4];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,4.172663362497963,-4.172663362497963) q[5];
u3(pi/2,0.8344070087934491,-0.8344070087934491) q[5];
rzz(-pi/2) q[2],q[5];
rzz(pi/2) q[2],q[6];
u3(pi/2,5.227610175573416,-5.227610175573416) q[6];
u3(pi/2,1.9879998311916212,-1.9879998311916212) q[6];
rzz(pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,5.7308933186785005,-5.7308933186785005) q[7];
u3(pi/2,2.5402918196927065,-2.5402918196927065) q[7];
rzz(pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[8];
u3(pi/2,2.092300707290802,-2.092300707290802) q[8];
u3(pi/2,5.209388938182594,-5.209388938182594) q[8];
rzz(pi/2) q[2],q[8];
rzz(pi/2) q[2],q[9];
u3(pi/2,5.237034953534185,-5.237034953534185) q[9];
u3(pi/2,2.082875929330033,-2.082875929330033) q[9];
rzz(pi/2) q[2],q[9];
u3(pi/2,4.9084243619686925,-4.9084243619686925) q[3];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(-pi/2) q[3],q[4];
u3(pi/2,4.810406671176691,-4.810406671176691) q[4];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,3.9759996623832423,-3.9759996623832423) q[5];
u3(pi/2,0.4417079270947249,-0.4417079270947249) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,1.9879998311916212,-1.9879998311916212) q[6];
u3(pi/2,4.933557103197411,-4.933557103197411) q[6];
rzz(-pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,2.5402918196927065,-2.5402918196927065) q[7];
u3(pi/2,5.583866782490499,-5.583866782490499) q[7];
rzz(pi/2) q[3],q[7];
rzz(pi/2) q[3],q[8];
u3(pi/2,5.209388938182594,-5.209388938182594) q[8];
u3(pi/2,2.018787439196801,-2.018787439196801) q[8];
rzz(-pi/2) q[3],q[8];
rzz(pi/2) q[3],q[9];
u3(pi/2,2.082875929330033,-2.082875929330033) q[9];
u3(pi/2,5.199964160221826,-5.199964160221826) q[9];
rzz(pi/2) q[3],q[9];
u3(pi/2,5.595804834574139,-5.595804834574139) q[4];
u3(pi/2,5*pi/8,-5*pi/8) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,0.4417079270947249,-0.4417079270947249) q[5];
u3(pi/2,2.7979024172870695,-2.7979024172870695) q[5];
rzz(pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,1.7919644496076181,-1.7919644496076181) q[6];
u3(pi/2,4.540858021498687,-4.540858021498687) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,5.583866782490499,-5.583866782490499) q[7];
u3(pi/2,2.245610428785984,-2.245610428785984) q[7];
rzz(-pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,5.160380092786594,-5.160380092786594) q[8];
u3(pi/2,1.9207697484047996,-1.9207697484047996) q[8];
rzz(pi/2) q[4],q[8];
rzz(pi/2) q[4],q[9];
u3(pi/2,5.199964160221826,-5.199964160221826) q[9];
u3(pi/2,2.0093626612360316,-2.0093626612360316) q[9];
rzz(pi/2) q[4],q[9];
u3(pi/2,4.368698744081967,-4.368698744081967) q[5];
u3(pi/2,5.399769452990137,-5.399769452990137) q[5];
rzz(pi/2) q[5],q[6];
u3(pi/2,4.540858021498687,-4.540858021498687) q[6];
u3(pi/2,0.6138672045114455,-0.6138672045114455) q[6];
rzz(pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,5.387203082375778,-5.387203082375778) q[7];
u3(pi/2,1.85291134708726,-1.85291134708726) q[7];
rzz(pi/2) q[5],q[7];
rzz(-pi/2) q[5],q[8];
u3(pi/2,5.0623624019945925,-5.0623624019945925) q[8];
u3(pi/2,1.7241060482900783,-1.7241060482900783) q[8];
rzz(pi/2) q[5],q[8];
rzz(pi/2) q[5],q[9];
u3(pi/2,2.0093626612360316,-2.0093626612360316) q[9];
u3(pi/2,5.052937624033824,-5.052937624033824) q[9];
rzz(pi/2) q[5],q[9];
u3(pi/2,2.184663531306342,-2.184663531306342) q[6];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[6];
rzz(-pi/2) q[6],q[7];
u3(pi/2,4.994504000677053,-4.994504000677053) q[7];
u3(pi/2,1.0675131836898117,-1.0675131836898117) q[7];
rzz(pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,1.7241060482900783,-1.7241060482900783) q[8];
u3(pi/2,4.472999620181147,-4.472999620181147) q[8];
rzz(pi/2) q[6],q[8];
rzz(pi/2) q[6],q[9];
u3(pi/2,5.052937624033824,-5.052937624033824) q[9];
u3(pi/2,1.7146812703293088,-1.7146812703293088) q[9];
rzz(-pi/2) q[6],q[9];
u3(pi/2,5.779902164074501,-5.779902164074501) q[7];
u3(pi/2,1.9879998311916212,-1.9879998311916212) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,4.472999620181147,-4.472999620181147) q[8];
u3(pi/2,0.5460088031939061,-0.5460088031939061) q[8];
rzz(pi/2) q[7],q[8];
rzz(pi/2) q[7],q[9];
u3(pi/2,4.856273923919102,-4.856273923919102) q[9];
u3(pi/2,1.321982188630585,-1.321982188630585) q[9];
rzz(-pi/2) q[7],q[9];
u3(pi/2,5.258397783578595,-5.258397783578595) q[8];
u3(pi/2,1.0065662862101699,-1.0065662862101699) q[8];
rzz(pi/2) q[8],q[9];
u3(pi/2,4.463574842220378,-4.463574842220378) q[9];
u3(pi/2,0.5365840252331366,-0.5365840252331366) q[9];
rzz(pi/2) q[8],q[9];
u3(pi,0,0) q[1];
u3(pi,3*pi/4,-3*pi/4) q[2];
u3(pi,3*pi/8,-3*pi/8) q[3];
u3(pi,11*pi/8,-11*pi/8) q[4];
u3(pi,5.792468534688861,-5.792468534688861) q[5];
u3(pi,2.8959201080790713,-2.8959201080790713) q[6];
u3(pi,0.9079202768874501,-0.9079202768874501) q[7];
u3(pi,1.2516105131901736,-1.2516105131901736) q[8];
u3(pi/2,2.107380352028033,-2.107380352028033) q[9];
u3(pi/2,3.4482120965801566,-3.4482120965801566) q[9];
u3(pi,4.381265114696325,-4.381265114696325) q[10];
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
