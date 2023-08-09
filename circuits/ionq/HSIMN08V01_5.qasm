OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(pi/2,5.085610187631157,-5.085610187631157) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.540574920595697,-3.540574920595697) q[0];
u3(pi/2,3.286105915654924,-3.286105915654924) q[1];
u3(pi/2,0.1338318470429252,-0.1338318470429252) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,1.7046281738378217,-1.7046281738378217) q[1];
u3(pi/2,1.7153095888600272,-1.7153095888600272) q[1];
u3(pi/2,0.409663682028109,-0.409663682028109) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,1.532468896421101,-1.532468896421101) q[3];
u3(pi/2,4.743176588389869,-4.743176588389869) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,3.134681149751896,-3.134681149751896) q[2];
u3(pi/2,3.0454599183899456,-3.0454599183899456) q[3];
u3(pi/2,6.176999475488251,-6.176999475488251) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,1.4646104951035617,-1.4646104951035617) q[3];
u3(pi/2,1.4746635915950488,-1.4746635915950488) q[3];
u3(pi/2,0.0037699111843077513,-0.0037699111843077513) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,4.716158891568997,-4.716158891568997) q[2];
u3(pi/2,1.7153095888600272,-1.7153095888600272) q[1];
rzz(-pi/2) q[2],q[1];
u3(pi/2,3.286105915654924,-3.286105915654924) q[1];
u3(pi/2,4.716158891568997,-4.716158891568997) q[2];
u3(pi/2,1.563884822956999,-1.563884822956999) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,3.134681149751896,-3.134681149751896) q[2];
u3(pi/2,3.1453625647741013,-3.1453625647741013) q[2];
u3(pi/2,0.15456635855661782,-0.15456635855661782) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,pi/2,-pi/2) q[5];
u3(pi/2,4.670919957357304,-4.670919957357304) q[5];
u3(pi/2,4.715530573038279,-4.715530573038279) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,3.0938404452552284,-3.0938404452552284) q[4];
u3(pi/2,3.0938404452552284,-3.0938404452552284) q[5];
u3(pi/2,6.225380002353534,-6.225380002353534) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,1.5129910219688443,-1.5129910219688443) q[5];
u3(pi/2,1.5230441184603318,-1.5230441184603318) q[5];
u3(pi/2,6.2461145138672265,-6.2461145138672265) q[4];
rzz(-pi/2) q[5],q[4];
u3(pi/2,1.533725533482537,-1.533725533482537) q[4];
u3(pi/2,4.616256245184842,-4.616256245184842) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3.0454599183899456,-3.0454599183899456) q[3];
u3(pi/2,4.67531818707233,-4.67531818707233) q[4];
u3(pi/2,1.5230441184603318,-1.5230441184603318) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.0938404452552284,-3.0938404452552284) q[4];
u3(pi/2,3.1045218602774334,-3.1045218602774334) q[4];
u3(pi/2,6.1977339870019446,-6.1977339870019446) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,4.649557127312894,-4.649557127312894) q[6];
rzz(-pi/2) q[7],q[6];
u3(pi/2,6.249884425051534,-6.249884425051534) q[6];
u3(pi/2,1.156106096521044,-1.156106096521044) q[7];
u3(pi/2,4.287017335088632,-4.287017335088632) q[7];
rzz(pi/2) q[6],q[7];
u3(pi/2,5.857813661883529,-5.857813661883529) q[7];
u3(pi/2,5.868495076905734,-5.868495076905734) q[7];
u3(pi/2,3.1189731864839465,-3.1189731864839465) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,1.5481768596890502,-1.5481768596890502) q[6];
u3(pi/2,1.5230441184603318,-1.5230441184603318) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,6.235433098845021,-6.235433098845021) q[5];
u3(pi/2,4.689769513278843,-4.689769513278843) q[6];
u3(pi/2,1.5374954446668447,-1.5374954446668447) q[6];
rzz(-pi/2) q[5],q[6];
u3(pi/2,6.249884425051534,-6.249884425051534) q[6];
u3(pi/2,6.260565840073739,-6.260565840073739) q[6];
u3(pi/2,6.2461145138672265,-6.2461145138672265) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,4.866955338941307,-4.866955338941307) q[1];
u3(pi/2,5.121424343882081,-5.121424343882081) q[0];
u3(pi/2,0.41846014145816046,-0.41846014145816046) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,2.015017528012493,-2.015017528012493) q[0];
u3(pi/2,1.7806547160546946,-1.7806547160546946) q[1];
u3(pi/2,4.911565954622282,-4.911565954622282) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,0.19917697423759287,-0.19917697423759287) q[1];
u3(pi/2,3.3307165313358986,-3.3307165313358986) q[1];
u3(pi/2,2.004336112990288,-2.004336112990288) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,3.0561413334121506,-3.0561413334121506) q[3];
u3(pi/2,4.588610229833251,-4.588610229833251) q[3];
u3(pi/2,0.05529203070318036,-0.05529203070318036) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,4.729981899244793,-4.729981899244793) q[2];
u3(pi/2,6.102229570332814,-6.102229570332814) q[3];
u3(pi/2,2.9499555017208157,-2.9499555017208157) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,4.520751828515713,-4.520751828515713) q[3];
u3(pi/2,4.531433243537918,-4.531433243537918) q[3];
u3(pi/2,1.5990706606772047,-1.5990706606772047) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,0.028274333882308135,-0.028274333882308135) q[2];
u3(pi/2,3.3307165313358986,-3.3307165313358986) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,1.7599202045410023,-1.7599202045410023) q[1];
u3(pi/2,3.169866987472101,-3.169866987472101) q[2];
u3(pi/2,0.017592918860102842,-0.017592918860102842) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,1.5883892456549995,-1.5883892456549995) q[2];
u3(pi/2,4.719928802753305,-4.719928802753305) q[2];
u3(pi/2,1.7492387895187966,-1.7492387895187966) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3.1045218602774334,-3.1045218602774334) q[5];
u3(pi/2,4.7167872100997155,-4.7167872100997155) q[5];
u3(pi/2,6.269990618034509,-6.269990618034509) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,4.648928808782176,-4.648928808782176) q[4];
u3(pi/2,6.28192867011815,-6.28192867011815) q[5];
u3(pi/2,3.129654601506152,-3.129654601506152) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,4.7004509283010485,-4.7004509283010485) q[5];
u3(pi/2,4.7111323433232535,-4.7111323433232535) q[5];
u3(pi/2,1.51738925168387,-1.51738925168387) q[4];
rzz(-pi/2) q[5],q[4];
u3(pi/2,3.0881855784787664,-3.0881855784787664) q[4];
u3(pi/2,1.3898405899481245,-1.3898405899481245) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,6.102229570332814,-6.102229570332814) q[3];
u3(pi/2,6.22977823206856,-6.22977823206856) q[4];
u3(pi/2,3.0781324819872795,-3.0781324819872795) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,4.648928808782176,-4.648928808782176) q[4];
u3(pi/2,4.658981905273664,-4.658981905273664) q[4];
u3(pi/2,2.9706900132345084,-2.9706900132345084) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3.0775041634565614,-3.0775041634565614) q[6];
rzz(-pi/2) q[7],q[6];
u3(pi/2,4.677831461195202,-4.677831461195202) q[6];
u3(pi/2,1.171185741258275,-1.171185741258275) q[7];
u3(pi/2,1.1818671562804801,-1.1818671562804801) q[7];
rzz(pi/2) q[6],q[7];
u3(pi/2,2.7526634830753767,-2.7526634830753767) q[7];
u3(pi/2,2.763344898097582,-2.763344898097582) q[7];
u3(pi/2,4.6671500461729964,-4.6671500461729964) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,3.0963537193781003,-3.0963537193781003) q[6];
u3(pi/2,4.7111323433232535,-4.7111323433232535) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,3.1403360165283574,-3.1403360165283574) q[5];
u3(pi/2,6.237946372967893,-6.237946372967893) q[6];
u3(pi/2,6.248627787990099,-6.248627787990099) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,1.5362388076054088,-1.5362388076054088) q[6];
u3(pi/2,1.5462919040968963,-1.5462919040968963) q[6];
u3(pi/2,3.129654601506152,-3.129654601506152) q[5];
rzz(-pi/2) q[6],q[5];
u3(pi/2,0.17844246272390027,-0.17844246272390027) q[1];
u3(pi/2,0.43353978619539146,-0.43353978619539146) q[0];
u3(pi/2,2.0131325724203397,-2.0131325724203397) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.6096899589746725,-3.6096899589746725) q[0];
u3(pi/2,3.3753271470168738,-3.3753271470168738) q[1];
u3(pi/2,0.2230530784048753,-0.2230530784048753) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,1.7938494051997718,-1.7938494051997718) q[1];
u3(pi/2,4.9253889622980775,-4.9253889622980775) q[1];
u3(pi/2,3.599008543952467,-3.599008543952467) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,6.112282666824301,-6.112282666824301) q[3];
u3(pi/2,1.3615662560658164,-1.3615662560658164) q[3];
u3(pi/2,1.5261857111139214,-1.5261857111139214) q[2];
rzz(-pi/2) q[3],q[2];
u3(pi/2,3.0630528372500483,-3.0630528372500483) q[2];
u3(pi/2,6.016778250155172,-6.016778250155172) q[3];
u3(pi/2,2.8645041815431735,-2.8645041815431735) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,4.43530050833807,-4.43530050833807) q[3];
u3(pi/2,1.2830264397260716,-1.2830264397260716) q[3];
u3(pi/2,3.052999740758561,-3.052999740758561) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,1.4822034139636644,-1.4822034139636644) q[2];
u3(pi/2,1.7837963087082844,-1.7837963087082844) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,0.21299998191338798,-0.21299998191338798) q[1];
u3(pi/2,4.6237960675534575,-4.6237960675534575) q[2];
u3(pi/2,1.471521998941459,-1.471521998941459) q[2];
rzz(-pi/2) q[1],q[2];
u3(pi/2,6.183910979326148,-6.183910979326148) q[2];
u3(pi/2,3.03163691071415,-3.03163691071415) q[2];
u3(pi/2,3.343911220480976,-3.343911220480976) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3.129654601506152,-3.129654601506152) q[5];
u3(pi/2,4.658981905273664,-4.658981905273664) q[5];
u3(pi/2,1.5412653558511524,-1.5412653558511524) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,6.203388853778405,-6.203388853778405) q[4];
u3(pi/2,6.235433098845021,-6.235433098845021) q[5];
u3(pi/2,6.2461145138672265,-6.2461145138672265) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,1.533725533482537,-1.533725533482537) q[5];
u3(pi/2,1.5444069485047422,-1.5444069485047422) q[5];
u3(pi/2,6.192707438756201,-6.192707438756201) q[4];
rzz(-pi/2) q[5],q[4];
u3(pi/2,1.4803184583715105,-1.4803184583715105) q[4];
u3(pi/2,4.424619093315865,-4.424619093315865) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,2.853822766520968,-2.853822766520968) q[3];
u3(pi/2,4.621911111961304,-4.621911111961304) q[4];
u3(pi/2,4.632592526983508,-4.632592526983508) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,6.203388853778405,-6.203388853778405) q[4];
u3(pi/2,3.051114785166407,-3.051114785166407) q[4];
u3(pi/2,6.006096835132967,-6.006096835132967) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,1.5048228810695108,-1.5048228810695108) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,6.246742832397945,-6.246742832397945) q[6];
u3(pi/2,1.2076282160399165,-1.2076282160399165) q[7];
u3(pi/2,1.2176813125314039,-1.2176813125314039) q[7];
rzz(-pi/2) q[6],q[7];
u3(pi/2,5.930070292916093,-5.930070292916093) q[7];
u3(pi/2,5.940751707938299,-5.940751707938299) q[7];
u3(pi/2,3.094468763785946,-3.094468763785946) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,1.5236724369910497,-1.5236724369910497) q[6];
u3(pi/2,1.5444069485047422,-1.5444069485047422) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,6.256795928889432,-6.256795928889432) q[5];
u3(pi/2,4.665265090580843,-4.665265090580843) q[6];
u3(pi/2,4.675946505603048,-4.675946505603048) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,6.246742832397945,-6.246742832397945) q[6];
u3(pi/2,6.25742424742015,-6.25742424742015) q[6];
u3(pi/2,6.2461145138672265,-6.2461145138672265) q[5];
rzz(-pi/2) q[6],q[5];
u3(pi/2,2.0282122171575705,-2.0282122171575705) q[0];
u3(pi/2,1.7731148936860792,-1.7731148936860792) q[1];
u3(pi/2,1.2937078547482768,-1.2937078547482768) q[3];
u3(pi/2,4.67531818707233,-4.67531818707233) q[5];
u3(pi,2.015017528012493,-2.015017528012493) q[6];
u3(pi,1.3144423662619695,-1.3144423662619695) q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
