OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
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
u3(pi/2,4.795955344970178,-4.795955344970178) q[3];
u3(pi/2,1.6160352610065896,-1.6160352610065896) q[3];
u3(pi/2,4.743176588389869,-4.743176588389869) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,3.134681149751896,-3.134681149751896) q[2];
u3(pi/2,3.129654601506152,-3.129654601506152) q[3];
u3(pi/2,6.260565840073739,-6.260565840073739) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,1.5481768596890502,-1.5481768596890502) q[3];
u3(pi/2,1.5588582747112552,-1.5588582747112552) q[3];
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
u3(pi/2,1.5588582747112552,-1.5588582747112552) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,6.271247255095945,-6.271247255095945) q[3];
u3(pi/2,4.67531818707233,-4.67531818707233) q[4];
u3(pi/2,1.5230441184603318,-1.5230441184603318) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.0938404452552284,-3.0938404452552284) q[4];
u3(pi/2,3.1045218602774334,-3.1045218602774334) q[4];
u3(pi/2,3.1397076979976393,-3.1397076979976393) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,4.649557127312894,-4.649557127312894) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,3.108291771461741,-3.108291771461741) q[6];
u3(pi/2,3.1258846903218442,-3.1258846903218442) q[7];
u3(pi/2,6.256795928889432,-6.256795928889432) q[7];
rzz(pi/2) q[6],q[7];
u3(pi/2,1.5444069485047422,-1.5444069485047422) q[7];
u3(pi/2,1.5550883635269477,-1.5550883635269477) q[7];
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
u3(pi/2,1.928937889304133,-1.928937889304133) q[9];
u3(pi/2,5.029689838397259,-5.029689838397259) q[9];
u3(pi/2,4.674689868541612,-4.674689868541612) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,3.0724776152108175,-3.0724776152108175) q[8];
u3(pi/2,3.40108820677631,-3.40108820677631) q[9];
u3(pi/2,0.24881413816431164,-0.24881413816431164) q[9];
rzz(-pi/2) q[8],q[9];
u3(pi/2,4.961203118549001,-4.961203118549001) q[9];
u3(pi/2,4.971884533571207,-4.971884533571207) q[9];
u3(pi/2,3.083159030233023,-3.083159030233023) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,1.5123627034381264,-1.5123627034381264) q[8];
u3(pi/2,4.696681017116741,-4.696681017116741) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,3.1258846903218442,-3.1258846903218442) q[7];
u3(pi/2,4.65395535702792,-4.65395535702792) q[8];
u3(pi/2,1.501681288415921,-1.501681288415921) q[8];
rzz(pi/2) q[7],q[8];
u3(pi/2,3.0724776152108175,-3.0724776152108175) q[8];
u3(pi/2,3.083159030233023,-3.083159030233023) q[8];
u3(pi/2,6.278158758933842,-6.278158758933842) q[7];
rzz(-pi/2) q[8],q[7];
u3(pi/2,1.7253626853515145,-1.7253626853515145) q[1];
u3(pi/2,5.121424343882081,-5.121424343882081) q[0];
u3(pi/2,0.41846014145816046,-0.41846014145816046) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,2.015017528012493,-2.015017528012493) q[0];
u3(pi/2,4.922247369644488,-4.922247369644488) q[1];
u3(pi/2,1.7699733010324894,-1.7699733010324894) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,3.3407696278273855,-3.3407696278273855) q[1];
u3(pi/2,0.18912387774610553,-0.18912387774610553) q[1];
u3(pi/2,2.004336112990288,-2.004336112990288) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.1397076979976393,-3.1397076979976393) q[3];
u3(pi/2,4.67217659441874,-4.67217659441874) q[3];
u3(pi/2,0.05529203070318036,-0.05529203070318036) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,4.729981899244793,-4.729981899244793) q[2];
u3(pi/2,6.1857959349183025,-6.1857959349183025) q[3];
u3(pi/2,3.0335218663063044,-3.0335218663063044) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,4.604318193101201,-4.604318193101201) q[3];
u3(pi/2,4.614999608123406,-4.614999608123406) q[3];
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
rzz(-pi/2) q[2],q[1];
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
u3(pi/2,1.4734069545336128,-1.4734069545336128) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,3.0442032813285094,-3.0442032813285094) q[3];
u3(pi/2,3.0881855784787664,-3.0881855784787664) q[4];
u3(pi/2,6.219725135577073,-6.219725135577073) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,1.5073361551923827,-1.5073361551923827) q[4];
u3(pi/2,1.51738925168387,-1.51738925168387) q[4];
u3(pi/2,6.1964773499405075,-6.1964773499405075) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,1.565769778549153,-1.565769778549153) q[7];
u3(pi/2,6.2190968170463545,-6.2190968170463545) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,4.677831461195202,-4.677831461195202) q[6];
u3(pi/2,4.691654468870997,-4.691654468870997) q[7];
u3(pi/2,1.5400087187897167,-1.5400087187897167) q[7];
rzz(pi/2) q[6],q[7];
u3(pi/2,3.110805045584613,-3.110805045584613) q[7];
u3(pi/2,3.1208581420761003,-3.1208581420761003) q[7];
u3(pi/2,1.5462919040968963,-1.5462919040968963) q[6];
rzz(-pi/2) q[7],q[6];
u3(pi/2,3.1170882308917927,-3.1170882308917927) q[6];
u3(pi/2,1.5695396897334606,-1.5695396897334606) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,6.28192867011815,-6.28192867011815) q[5];
u3(pi/2,6.258680884481586,-6.258680884481586) q[6];
u3(pi/2,3.1070351344003053,-3.1070351344003053) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,4.677831461195202,-4.677831461195202) q[6];
u3(pi/2,4.687884557686689,-4.687884557686689) q[6];
u3(pi/2,3.1510174315505624,-3.1510174315505624) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,4.971884533571207,-4.971884533571207) q[9];
u3(pi/2,1.7894511754847462,-1.7894511754847462) q[9];
u3(pi/2,3.066194429903638,-3.066194429903638) q[8];
rzz(-pi/2) q[9],q[8];
u3(pi/2,4.605574830162636,-4.605574830162636) q[8];
u3(pi/2,0.27646015351590175,-0.27646015351590175) q[9];
u3(pi/2,0.2871415685381071,-0.2871415685381071) q[9];
rzz(pi/2) q[8],q[9];
u3(pi/2,1.8579378953330037,-1.8579378953330037) q[9];
u3(pi/2,1.867990991824491,-1.867990991824491) q[9];
u3(pi/2,4.595521733671149,-4.595521733671149) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,3.024725406876253,-3.024725406876253) q[8];
u3(pi/2,3.1208581420761003,-3.1208581420761003) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,1.550061815281204,-1.550061815281204) q[7];
u3(pi/2,6.166318060466046,-6.166318060466046) q[8];
u3(pi/2,6.176371156957533,-6.176371156957533) q[8];
rzz(-pi/2) q[7],q[8];
u3(pi/2,4.605574830162636,-4.605574830162636) q[8];
u3(pi/2,4.616256245184842,-4.616256245184842) q[8];
u3(pi/2,4.68160137237951,-4.68160137237951) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,3.320035116313693,-3.320035116313693) q[1];
u3(pi/2,3.5751324397851842,-3.5751324397851842) q[0];
u3(pi/2,5.154725226010132,-5.154725226010132) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,0.46809730538487915,-0.46809730538487915) q[0];
u3(pi/2,0.23373449342708058,-0.23373449342708058) q[1];
u3(pi/2,3.3646457319946683,-3.3646457319946683) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,4.935442058789564,-4.935442058789564) q[1];
u3(pi/2,1.7837963087082844,-1.7837963087082844) q[1];
u3(pi/2,0.4574158903626739,-0.4574158903626739) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,3.054884696350715,-3.054884696350715) q[3];
u3(pi/2,4.586725274241098,-4.586725274241098) q[3];
u3(pi/2,4.6677783647037145,-4.6677783647037145) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,3.0630528372500483,-3.0630528372500483) q[2];
u3(pi/2,6.10034461474066,-6.10034461474066) q[3];
u3(pi/2,2.948070546128662,-2.948070546128662) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,4.518866872923558,-4.518866872923558) q[3];
u3(pi/2,1.367221122842278,-1.367221122842278) q[3];
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
u3(pi/2,0.00942477796076938,-0.00942477796076938) q[5];
u3(pi/2,1.621690127783051,-1.621690127783051) q[5];
u3(pi/2,4.682858009440945,-4.682858009440945) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,3.0617962001886125,-3.0617962001886125) q[4];
u3(pi/2,3.1862032692707682,-3.1862032692707682) q[5];
u3(pi/2,0.03455751918948772,-0.03455751918948772) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,1.6053538459843844,-1.6053538459843844) q[5];
u3(pi/2,1.6154069424758717,-1.6154069424758717) q[5];
u3(pi/2,6.2134419502698925,-6.2134419502698925) q[4];
rzz(-pi/2) q[5],q[4];
u3(pi/2,1.501052969885203,-1.501052969885203) q[4];
u3(pi/2,4.5088137764320715,-4.5088137764320715) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,2.9380174496371745,-2.9380174496371745) q[3];
u3(pi/2,4.642645623474996,-4.642645623474996) q[4];
u3(pi/2,1.490999873393716,-1.490999873393716) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.0617962001886125,-3.0617962001886125) q[4];
u3(pi/2,6.192707438756201,-6.192707438756201) q[4];
u3(pi/2,2.927336034614969,-2.927336034614969) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3.110805045584613,-3.110805045584613) q[7];
u3(pi/2,1.5048228810695108,-1.5048228810695108) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,6.246742832397945,-6.246742832397945) q[6];
u3(pi/2,6.267477343911637,-6.267477343911637) q[7];
u3(pi/2,6.278158758933842,-6.278158758933842) q[7];
rzz(-pi/2) q[6],q[7];
u3(pi/2,4.707362432138946,-4.707362432138946) q[7];
u3(pi/2,4.717415528630434,-4.717415528630434) q[7];
u3(pi/2,3.094468763785946,-3.094468763785946) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,1.5236724369910497,-1.5236724369910497) q[6];
u3(pi/2,1.6154069424758717,-1.6154069424758717) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,0.04461061568097507,-0.04461061568097507) q[5];
u3(pi/2,4.665265090580843,-4.665265090580843) q[6];
u3(pi/2,4.675946505603048,-4.675946505603048) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,6.246742832397945,-6.246742832397945) q[6];
u3(pi/2,6.25742424742015,-6.25742424742015) q[6];
u3(pi/2,0.03455751918948772,-0.03455751918948772) q[5];
rzz(-pi/2) q[6],q[5];
u3(pi/2,1.867990991824491,-1.867990991824491) q[9];
u3(pi/2,4.9693712594483355,-4.9693712594483355) q[9];
u3(pi/2,1.4576989912656642,-1.4576989912656642) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,6.138672045114456,-6.138672045114456) q[8];
u3(pi/2,0.31478758388969724,-0.31478758388969724) q[9];
u3(pi/2,0.32484068038118463,-0.32484068038118463) q[9];
rzz(pi/2) q[8],q[9];
u3(pi/2,1.8956370071760813,-1.8956370071760813) q[9];
u3(pi/2,1.9063184221982865,-1.9063184221982865) q[9];
u3(pi/2,6.128618948622969,-6.128618948622969) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,4.557822621828072,-4.557822621828072) q[8];
u3(pi/2,1.5758228750406404,-1.5758228750406404) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,pi/625,-pi/625) q[7];
u3(pi/2,1.4162299682382786,-1.4162299682382786) q[8];
u3(pi/2,1.4262830647297662,-1.4262830647297662) q[8];
rzz(pi/2) q[7],q[8];
u3(pi/2,2.9970793915246623,-2.9970793915246623) q[8];
u3(pi/2,3.007760806546868,-3.007760806546868) q[8];
u3(pi/2,6.278158758933842,-6.278158758933842) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,2.0282122171575705,-2.0282122171575705) q[0];
u3(pi/2,1.7731148936860792,-1.7731148936860792) q[1];
u3(pi/2,1.3565397078200727,-1.3565397078200727) q[3];
u3(pi/2,4.746946499574177,-4.746946499574177) q[5];
u3(pi,4.074017353175243,-4.074017353175243) q[6];
u3(pi/2,1.565769778549153,-1.565769778549153) q[7];
u3(pi,3.276681137694154,-3.276681137694154) q[8];
u3(pi,2.384468824074653,-2.384468824074653) q[9];
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
