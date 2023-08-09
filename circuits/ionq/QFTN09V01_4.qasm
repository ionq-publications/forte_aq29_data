OPENQASM 2.0;
include "qelib1.inc";
qreg q[9];
creg c[9];
u3(pi/2,2.061513099285622,-2.061513099285622) q[8];
u3(pi/2,2.8469112626830704,-2.8469112626830704) q[8];
rzz(pi/2) q[7],q[8];
rzz(-pi/2) q[6],q[8];
u3(pi/2,2.8469112626830704,-2.8469112626830704) q[8];
u3(pi/2,5.595804834574139,-5.595804834574139) q[8];
rzz(-pi/2) q[6],q[8];
rzz(pi/2) q[6],q[7];
u3(pi/2,3.583300580684518,-3.583300580684518) q[7];
u3(pi/2,4.368698744081967,-4.368698744081967) q[7];
rzz(-pi/2) q[6],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,2.4542121809843462,-2.4542121809843462) q[8];
u3(pi/2,5.399769452990137,-5.399769452990137) q[8];
rzz(-pi/2) q[5],q[8];
rzz(-pi/2) q[5],q[7];
u3(pi/2,1.2271060904921731,-1.2271060904921731) q[7];
u3(pi/2,3.9759996623832423,-3.9759996623832423) q[7];
rzz(-pi/2) q[5],q[7];
rzz(-pi/2) q[5],q[6];
u3(pi/2,2.061513099285622,-2.061513099285622) q[6];
u3(pi/2,2.8469112626830704,-2.8469112626830704) q[6];
rzz(-pi/2) q[5],q[6];
rzz(-pi/2) q[4],q[8];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[8];
u3(pi/2,5.301123443667417,-5.301123443667417) q[8];
rzz(-pi/2) q[4],q[8];
rzz(-pi/2) q[4],q[7];
u3(pi/2,0.8344070087934491,-0.8344070087934491) q[7];
u3(pi/2,3.779964280799239,-3.779964280799239) q[7];
rzz(pi/2) q[4],q[7];
rzz(-pi/2) q[4],q[6];
u3(pi/2,2.8469112626830704,-2.8469112626830704) q[6];
u3(pi/2,5.595804834574139,-5.595804834574139) q[6];
rzz(-pi/2) q[4],q[6];
rzz(pi/2) q[4],q[5];
u3(pi/2,3.3376280351737964,-3.3376280351737964) q[5];
u3(pi/2,4.123026198571244,-4.123026198571244) q[5];
rzz(pi/2) q[4],q[5];
rzz(-pi/2) q[3],q[8];
u3(pi/2,5.301123443667417,-5.301123443667417) q[8];
u3(pi/2,2.110521944681623,-2.110521944681623) q[8];
rzz(-pi/2) q[3],q[8];
rzz(-pi/2) q[3],q[7];
u3(pi/2,0.638371627209446,-0.638371627209446) q[7];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[7];
rzz(pi/2) q[3],q[7];
rzz(-pi/2) q[3],q[6];
u3(pi/2,5.595804834574139,-5.595804834574139) q[6];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[6];
rzz(-pi/2) q[3],q[6];
rzz(-pi/2) q[3],q[5];
u3(pi/2,4.123026198571244,-4.123026198571244) q[5];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[4];
u3(pi/2,9*pi/8,-9*pi/8) q[4];
u3(pi/2,11*pi/8,-11*pi/8) q[4];
rzz(-pi/2) q[3],q[4];
rzz(-pi/2) q[2],q[8];
u3(pi/2,5.252114598271416,-5.252114598271416) q[8];
u3(pi/2,2.0860175219836226,-2.0860175219836226) q[8];
rzz(-pi/2) q[2],q[8];
rzz(-pi/2) q[2],q[7];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[7];
u3(pi/2,0.4907167724907257,-0.4907167724907257) q[7];
rzz(-pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[6];
u3(pi/2,5.399769452990137,-5.399769452990137) q[6];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[6];
rzz(-pi/2) q[2],q[6];
rzz(pi/2) q[2],q[5];
u3(pi/2,3.730955435403238,-3.730955435403238) q[5];
u3(pi/2,pi/8,-pi/8) q[5];
rzz(-pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[4];
u3(pi/2,11*pi/8,-11*pi/8) q[4];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[2],q[4];
rzz(-pi/2) q[2],q[3];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(-pi/2) q[2],q[3];
rzz(-pi/2) q[1],q[8];
u3(pi/2,5.227610175573416,-5.227610175573416) q[8];
u3(pi/2,2.0740794698999814,-2.0740794698999814) q[8];
rzz(-pi/2) q[1],q[8];
rzz(-pi/2) q[1],q[7];
u3(pi/2,3.6323094260805187,-3.6323094260805187) q[7];
u3(pi/2,0.4662123497927253,-0.4662123497927253) q[7];
rzz(pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[6];
u3(pi/2,5.301751762198135,-5.301751762198135) q[6];
u3(pi/2,2.110521944681623,-2.110521944681623) q[6];
rzz(-pi/2) q[1],q[6];
rzz(-pi/2) q[1],q[5];
u3(pi/2,9*pi/8,-9*pi/8) q[5];
u3(pi/2,0.2946813909067226,-0.2946813909067226) q[5];
rzz(pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[4];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[4];
rzz(-pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
u3(pi/2,13*pi/8,-13*pi/8) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[2];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(-pi/2) q[1],q[2];
rzz(-pi/2) q[0],q[8];
rzz(pi/2) q[0],q[7];
rzz(-pi/2) q[0],q[6];
rzz(-pi/2) q[0],q[5];
rzz(-pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[3];
rzz(pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi,-pi) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
u3(pi/2,0,0) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,5*pi/8,-5*pi/8) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[4];
u3(pi/2,9*pi/8,-9*pi/8) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,3.436274044496516,-3.436274044496516) q[5];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,2.110521944681623,-2.110521944681623) q[6];
u3(pi/2,5.203105752875415,-5.203105752875415) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,0.4662123497927253,-0.4662123497927253) q[7];
u3(pi/2,3.583300580684518,-3.583300580684518) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,2.0740794698999814,-2.0740794698999814) q[8];
u3(pi/2,5.203105752875415,-5.203105752875415) q[8];
rzz(-pi/2) q[0],q[8];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(-pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,pi/8,-pi/8) q[4];
u3(pi/2,3.3376280351737964,-3.3376280351737964) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,3.3376280351737964,-3.3376280351737964) q[5];
u3(pi/2,0.09801769079200154,-0.09801769079200154) q[5];
rzz(pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,2.061513099285622,-2.061513099285622) q[6];
u3(pi/2,5.154096907479415,-5.154096907479415) q[6];
rzz(pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u3(pi/2,0.4417079270947249,-0.4417079270947249) q[7];
u3(pi/2,3.5587961579865177,-3.5587961579865177) q[7];
rzz(-pi/2) q[1],q[7];
rzz(pi/2) q[1],q[8];
u3(pi/2,5.203105752875415,-5.203105752875415) q[8];
u3(pi/2,2.049575047201981,-2.049575047201981) q[8];
rzz(pi/2) q[1],q[8];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
u3(pi/2,pi/8,-pi/8) q[3];
rzz(pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,3.3376280351737964,-3.3376280351737964) q[4];
u3(pi/2,6.087149925595583,-6.087149925595583) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,0.09801769079200154,-0.09801769079200154) q[5];
u3(pi/2,3.0435749627977917,-3.0435749627977917) q[5];
rzz(pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[6];
u3(pi/2,2.012504253889621,-2.012504253889621) q[6];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[6];
rzz(pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,0.41720350439672454,-0.41720350439672454) q[7];
u3(pi/2,3.509787312590517,-3.509787312590517) q[7];
rzz(pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[8];
u3(pi/2,5.191167700791774,-5.191167700791774) q[8];
u3(pi/2,2.0250706245039805,-2.0250706245039805) q[8];
rzz(pi/2) q[2],q[8];
u3(pi/2,13*pi/8,-13*pi/8) q[3];
u3(pi/2,7*pi/4,-7*pi/4) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,6.087149925595583,-6.087149925595583) q[4];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,3.0435749627977917,-3.0435749627977917) q[5];
u3(pi/2,5.792468534688861,-5.792468534688861) q[5];
rzz(-pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[6];
u3(pi/2,1.7178228629828987,-1.7178228629828987) q[6];
rzz(pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,3.509787312590517,-3.509787312590517) q[7];
u3(pi/2,0.27017696820872217,-0.27017696820872217) q[7];
rzz(-pi/2) q[3],q[7];
rzz(pi/2) q[3],q[8];
u3(pi/2,2.0250706245039805,-2.0250706245039805) q[8];
u3(pi/2,5.117654432697773,-5.117654432697773) q[8];
rzz(pi/2) q[3],q[8];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[4];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,2.6508758810990676,-2.6508758810990676) q[5];
u3(pi/2,5.007070371291412,-5.007070371291412) q[5];
rzz(pi/2) q[4],q[5];
rzz(-pi/2) q[4],q[6];
u3(pi/2,4.8594155165726916,-4.8594155165726916) q[6];
u3(pi/2,1.3251237812841747,-1.3251237812841747) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,3.4117696217985154,-3.4117696217985154) q[7];
u3(pi/2,0.07351326809400116,-0.07351326809400116) q[7];
rzz(-pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,5.117654432697773,-5.117654432697773) q[8];
u3(pi/2,1.8774157697852605,-1.8774157697852605) q[8];
rzz(pi/2) q[4],q[8];
u3(pi/2,0.2946813909067226,-0.2946813909067226) q[5];
u3(pi/2,3.3376280351737964,-3.3376280351737964) q[5];
rzz(pi/2) q[5],q[6];
u3(pi/2,1.3251237812841747,-1.3251237812841747) q[6];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[6];
rzz(pi/2) q[5],q[6];
rzz(-pi/2) q[5],q[7];
u3(pi/2,0.07351326809400116,-0.07351326809400116) q[7];
u3(pi/2,2.82240683998507,-2.82240683998507) q[7];
rzz(pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,1.8774157697852605,-1.8774157697852605) q[8];
u3(pi/2,4.82297304179105,-4.82297304179105) q[8];
rzz(pi/2) q[5],q[8];
u3(pi/2,5.252114598271416,-5.252114598271416) q[6];
u3(pi/2,2.061513099285622,-2.061513099285622) q[6];
rzz(-pi/2) q[6],q[7];
u3(pi/2,5.963999493574864,-5.963999493574864) q[7];
u3(pi/2,2.0370086765876216,-2.0370086765876216) q[7];
rzz(pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,4.82297304179105,-4.82297304179105) q[8];
u3(pi/2,1.288681306502533,-1.288681306502533) q[8];
rzz(pi/2) q[6],q[8];
u3(pi/2,3.6078050033825186,-3.6078050033825186) q[7];
u3(pi/2,0.4417079270947249,-0.4417079270947249) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,1.288681306502533,-1.288681306502533) q[8];
u3(pi/2,3.6448757966948775,-3.6448757966948775) q[8];
rzz(-pi/2) q[7],q[8];
u3(pi,pi,-pi) q[1];
u3(pi,pi/2,-pi/2) q[2];
u3(pi,pi/2,-pi/2) q[3];
u3(pi,7*pi/8,-7*pi/8) q[4];
u3(pi,11*pi/8,-11*pi/8) q[5];
u3(pi,2.552229871776348,-2.552229871776348) q[6];
u3(pi,5.399769452990137,-5.399769452990137) q[7];
u3(pi/2,2.0740794698999814,-2.0740794698999814) q[8];
u3(pi/2,5.203105752875415,-5.203105752875415) q[8];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
measure q[8] -> c[8];
