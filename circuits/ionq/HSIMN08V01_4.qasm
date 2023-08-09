OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(pi/2,5.085610187631157,-5.085610187631157) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,0.3989822670059037,-0.3989822670059037) q[0];
u3(pi/2,0.1445132620651305,-0.1445132620651305) q[1];
u3(pi/2,3.275424500632718,-3.275424500632718) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,4.846220827427615,-4.846220827427615) q[1];
u3(pi/2,4.85690224244982,-4.85690224244982) q[1];
u3(pi/2,3.5512563356179023,-3.5512563356179023) q[0];
rzz(pi/2) q[1],q[0];
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
rzz(-pi/2) q[3],q[2];
u3(pi/2,1.5745662379792043,-1.5745662379792043) q[2];
u3(pi/2,1.7153095888600272,-1.7153095888600272) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,0.1445132620651305,-0.1445132620651305) q[1];
u3(pi/2,4.716158891568997,-4.716158891568997) q[2];
u3(pi/2,1.563884822956999,-1.563884822956999) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,3.134681149751896,-3.134681149751896) q[2];
u3(pi/2,3.1453625647741013,-3.1453625647741013) q[2];
u3(pi/2,3.296159012146411,-3.296159012146411) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,pi/2,-pi/2) q[5];
u3(pi/2,4.670919957357304,-4.670919957357304) q[5];
u3(pi/2,4.715530573038279,-4.715530573038279) q[4];
rzz(-pi/2) q[5],q[4];
u3(pi/2,6.235433098845021,-6.235433098845021) q[4];
u3(pi/2,6.235433098845021,-6.235433098845021) q[5];
u3(pi/2,3.083787348763741,-3.083787348763741) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,4.654583675558637,-4.654583675558637) q[5];
u3(pi/2,4.664636772050124,-4.664636772050124) q[5];
u3(pi/2,3.1045218602774334,-3.1045218602774334) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,1.533725533482537,-1.533725533482537) q[4];
u3(pi/2,1.4746635915950488,-1.4746635915950488) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,6.187052571979739,-6.187052571979739) q[3];
u3(pi/2,4.67531818707233,-4.67531818707233) q[4];
u3(pi/2,1.5230441184603318,-1.5230441184603318) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.0938404452552284,-3.0938404452552284) q[4];
u3(pi/2,3.1045218602774334,-3.1045218602774334) q[4];
u3(pi/2,3.0561413334121506,-3.0561413334121506) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,4.649557127312894,-4.649557127312894) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,3.108291771461741,-3.108291771461741) q[6];
u3(pi/2,4.297698750110837,-4.297698750110837) q[7];
u3(pi/2,1.1454246814988385,-1.1454246814988385) q[7];
rzz(pi/2) q[6],q[7];
u3(pi/2,2.716221008293735,-2.716221008293735) q[7];
u3(pi/2,2.7269024233159405,-2.7269024233159405) q[7];
u3(pi/2,6.260565840073739,-6.260565840073739) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,4.689769513278843,-4.689769513278843) q[6];
u3(pi/2,1.5230441184603318,-1.5230441184603318) q[5];
rzz(-pi/2) q[6],q[5];
u3(pi/2,3.0938404452552284,-3.0938404452552284) q[5];
u3(pi/2,4.689769513278843,-4.689769513278843) q[6];
u3(pi/2,1.5374954446668447,-1.5374954446668447) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,3.108291771461741,-3.108291771461741) q[6];
u3(pi/2,3.1189731864839465,-3.1189731864839465) q[6];
u3(pi/2,6.2461145138672265,-6.2461145138672265) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,1.7253626853515145,-1.7253626853515145) q[1];
u3(pi/2,5.121424343882081,-5.121424343882081) q[0];
u3(pi/2,0.41846014145816046,-0.41846014145816046) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,2.015017528012493,-2.015017528012493) q[0];
u3(pi/2,4.922247369644488,-4.922247369644488) q[1];
u3(pi/2,1.7699733010324894,-1.7699733010324894) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,0.19917697423759287,-0.19917697423759287) q[1];
u3(pi/2,3.3307165313358986,-3.3307165313358986) q[1];
u3(pi/2,5.145928766580081,-5.145928766580081) q[0];
rzz(pi/2) q[1],q[0];
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
rzz(-pi/2) q[3],q[2];
u3(pi/2,3.169866987472101,-3.169866987472101) q[2];
u3(pi/2,0.18912387774610553,-0.18912387774610553) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,4.901512858130795,-4.901512858130795) q[1];
u3(pi/2,0.028274333882308135,-0.028274333882308135) q[2];
u3(pi/2,3.1591855724498963,-3.1591855724498963) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,4.729981899244793,-4.729981899244793) q[2];
u3(pi/2,1.578336149163512,-1.578336149163512) q[2];
u3(pi/2,4.89083144310859,-4.89083144310859) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3.1045218602774334,-3.1045218602774334) q[5];
u3(pi/2,4.7167872100997155,-4.7167872100997155) q[5];
u3(pi/2,3.128397964444716,-3.128397964444716) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,1.5073361551923827,-1.5073361551923827) q[4];
u3(pi/2,6.28192867011815,-6.28192867011815) q[5];
u3(pi/2,3.129654601506152,-3.129654601506152) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,4.7004509283010485,-4.7004509283010485) q[5];
u3(pi/2,4.7111323433232535,-4.7111323433232535) q[5];
u3(pi/2,4.658981905273664,-4.658981905273664) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,3.0881855784787664,-3.0881855784787664) q[4];
u3(pi/2,4.531433243537918,-4.531433243537918) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,2.960636916743021,-2.960636916743021) q[3];
u3(pi/2,6.22977823206856,-6.22977823206856) q[4];
u3(pi/2,3.0781324819872795,-3.0781324819872795) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,4.648928808782176,-4.648928808782176) q[4];
u3(pi/2,4.658981905273664,-4.658981905273664) q[4];
u3(pi/2,6.112282666824301,-6.112282666824301) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,6.2190968170463545,-6.2190968170463545) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,4.677831461195202,-4.677831461195202) q[6];
u3(pi/2,1.171185741258275,-1.171185741258275) q[7];
u3(pi/2,1.1818671562804801,-1.1818671562804801) q[7];
rzz(pi/2) q[6],q[7];
u3(pi/2,2.7526634830753767,-2.7526634830753767) q[7];
u3(pi/2,2.763344898097582,-2.763344898097582) q[7];
u3(pi/2,4.6671500461729964,-4.6671500461729964) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,3.0963537193781003,-3.0963537193781003) q[6];
u3(pi/2,1.5695396897334606,-1.5695396897334606) q[5];
rzz(-pi/2) q[6],q[5];
u3(pi/2,3.1403360165283574,-3.1403360165283574) q[5];
u3(pi/2,3.0963537193781003,-3.0963537193781003) q[6];
u3(pi/2,3.1070351344003053,-3.1070351344003053) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,4.677831461195202,-4.677831461195202) q[6];
u3(pi/2,4.687884557686689,-4.687884557686689) q[6];
u3(pi/2,3.129654601506152,-3.129654601506152) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,3.320035116313693,-3.320035116313693) q[1];
u3(pi/2,0.43353978619539146,-0.43353978619539146) q[0];
u3(pi/2,2.0131325724203397,-2.0131325724203397) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.6096899589746725,-3.6096899589746725) q[0];
u3(pi/2,0.23373449342708058,-0.23373449342708058) q[1];
u3(pi/2,3.3646457319946683,-3.3646457319946683) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,4.935442058789564,-4.935442058789564) q[1];
u3(pi/2,1.7837963087082844,-1.7837963087082844) q[1];
u3(pi/2,3.599008543952467,-3.599008543952467) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,6.112282666824301,-6.112282666824301) q[3];
u3(pi/2,1.3615662560658164,-1.3615662560658164) q[3];
u3(pi/2,4.6677783647037145,-4.6677783647037145) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,3.0630528372500483,-3.0630528372500483) q[2];
u3(pi/2,2.8751855965653785,-2.8751855965653785) q[3];
u3(pi/2,6.006096835132967,-6.006096835132967) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,1.2937078547482768,-1.2937078547482768) q[3];
u3(pi/2,4.424619093315865,-4.424619093315865) q[3];
u3(pi/2,3.052999740758561,-3.052999740758561) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,1.4822034139636644,-1.4822034139636644) q[2];
u3(pi/2,1.7837963087082844,-1.7837963087082844) q[1];
rzz(-pi/2) q[2],q[1];
u3(pi/2,3.3545926355031814,-3.3545926355031814) q[1];
u3(pi/2,1.4822034139636644,-1.4822034139636644) q[2];
u3(pi/2,4.613114652531252,-4.613114652531252) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,6.183910979326148,-6.183910979326148) q[2];
u3(pi/2,3.03163691071415,-3.03163691071415) q[2];
u3(pi/2,3.343911220480976,-3.343911220480976) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,6.271247255095945,-6.271247255095945) q[5];
u3(pi/2,1.51738925168387,-1.51738925168387) q[5];
u3(pi/2,4.682858009440945,-4.682858009440945) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,3.0617962001886125,-3.0617962001886125) q[4];
u3(pi/2,3.0938404452552284,-3.0938404452552284) q[5];
u3(pi/2,3.1045218602774334,-3.1045218602774334) q[5];
rzz(-pi/2) q[4],q[5];
u3(pi/2,1.533725533482537,-1.533725533482537) q[5];
u3(pi/2,1.5444069485047422,-1.5444069485047422) q[5];
u3(pi/2,6.192707438756201,-6.192707438756201) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,4.621911111961304,-4.621911111961304) q[4];
u3(pi/2,1.2830264397260716,-1.2830264397260716) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,5.995415420110762,-5.995415420110762) q[3];
u3(pi/2,1.4803184583715105,-1.4803184583715105) q[4];
u3(pi/2,1.490999873393716,-1.490999873393716) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.0617962001886125,-3.0617962001886125) q[4];
u3(pi/2,6.192707438756201,-6.192707438756201) q[4];
u3(pi/2,2.8645041815431735,-2.8645041815431735) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,1.5048228810695108,-1.5048228810695108) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,6.246742832397945,-6.246742832397945) q[6];
u3(pi/2,1.2076282160399165,-1.2076282160399165) q[7];
u3(pi/2,1.2176813125314039,-1.2176813125314039) q[7];
rzz(pi/2) q[6],q[7];
u3(pi/2,2.7884776393263,-2.7884776393263) q[7];
u3(pi/2,2.7991590543485056,-2.7991590543485056) q[7];
u3(pi/2,6.23606141737574,-6.23606141737574) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,4.665265090580843,-4.665265090580843) q[6];
u3(pi/2,4.685999602094536,-4.685999602094536) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,3.115203275299639,-3.115203275299639) q[5];
u3(pi/2,1.5236724369910497,-1.5236724369910497) q[6];
u3(pi/2,1.534353852013255,-1.534353852013255) q[6];
rzz(-pi/2) q[5],q[6];
u3(pi/2,6.246742832397945,-6.246742832397945) q[6];
u3(pi/2,6.25742424742015,-6.25742424742015) q[6];
u3(pi/2,6.2461145138672265,-6.2461145138672265) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,5.169804870747363,-5.169804870747363) q[0];
u3(pi/2,1.7731148936860792,-1.7731148936860792) q[1];
u3(pi/2,1.2937078547482768,-1.2937078547482768) q[3];
u3(pi/2,1.533725533482537,-1.533725533482537) q[5];
u3(pi,0.44422120121759673,-0.44422120121759673) q[6];
u3(pi,4.456035019851763,-4.456035019851763) q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
