OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(pi/2,5.085610187631157,-5.085610187631157) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.540574920595697,-3.540574920595697) q[0];
u3(pi/2,3.2019112325387176,-3.2019112325387176) q[1];
u3(pi/2,0.05026548245743669,-0.05026548245743669) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,1.6210618092523332,-1.6210618092523332) q[1];
u3(pi/2,1.6311149057438206,-1.6311149057438206) q[1];
u3(pi/2,0.409663682028109,-0.409663682028109) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,1.532468896421101,-1.532468896421101) q[3];
u3(pi/2,4.743176588389869,-4.743176588389869) q[2];
rzz(-pi/2) q[3],q[2];
u3(pi/2,6.276273803341689,-6.276273803341689) q[2];
u3(pi/2,6.187052571979739,-6.187052571979739) q[3];
u3(pi/2,3.035406821898458,-3.035406821898458) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,4.606203148693354,-4.606203148693354) q[3];
u3(pi/2,4.616256245184842,-4.616256245184842) q[3];
u3(pi/2,3.1453625647741013,-3.1453625647741013) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,1.5745662379792043,-1.5745662379792043) q[2];
u3(pi/2,4.772707559333614,-4.772707559333614) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3.2019112325387176,-3.2019112325387176) q[1];
u3(pi/2,4.716158891568997,-4.716158891568997) q[2];
u3(pi/2,1.563884822956999,-1.563884822956999) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,3.134681149751896,-3.134681149751896) q[2];
u3(pi/2,3.1453625647741013,-3.1453625647741013) q[2];
u3(pi/2,0.07099999397112931,-0.07099999397112931) q[1];
rzz(-pi/2) q[2],q[1];
u3(pi/2,3.6850881826608273,-3.6850881826608273) q[5];
u3(pi/2,0.501398187512931,-0.501398187512931) q[5];
u3(pi/2,4.715530573038279,-4.715530573038279) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,3.0938404452552284,-3.0938404452552284) q[4];
u3(pi/2,5.208132301121159,-5.208132301121159) q[5];
u3(pi/2,2.0558582325091606,-2.0558582325091606) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,3.6266545593040576,-3.6266545593040576) q[5];
u3(pi/2,3.637335974326262,-3.637335974326262) q[5];
u3(pi/2,6.2461145138672265,-6.2461145138672265) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,4.67531818707233,-4.67531818707233) q[4];
u3(pi/2,1.4746635915950488,-1.4746635915950488) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,3.0454599183899456,-3.0454599183899456) q[3];
u3(pi/2,4.67531818707233,-4.67531818707233) q[4];
u3(pi/2,1.5230441184603318,-1.5230441184603318) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.0938404452552284,-3.0938404452552284) q[4];
u3(pi/2,3.1045218602774334,-3.1045218602774334) q[4];
u3(pi/2,6.1977339870019446,-6.1977339870019446) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,1.6417963207660258,-1.6417963207660258) q[1];
u3(pi/2,1.9798316902922877,-1.9798316902922877) q[0];
u3(pi/2,3.5600527950479535,-3.5600527950479535) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,5.156610181602287,-5.156610181602287) q[0];
u3(pi/2,4.838681005059,-4.838681005059) q[1];
u3(pi/2,1.6864069364470011,-1.6864069364470011) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,3.2572032632418972,-3.2572032632418972) q[1];
u3(pi/2,0.10555751316061704,-0.10555751316061704) q[1];
u3(pi/2,5.145928766580081,-5.145928766580081) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,3.0561413334121506,-3.0561413334121506) q[3];
u3(pi/2,4.588610229833251,-4.588610229833251) q[3];
u3(pi/2,3.1968846842929737,-3.1968846842929737) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,1.5883892456549995,-1.5883892456549995) q[2];
u3(pi/2,6.102229570332814,-6.102229570332814) q[3];
u3(pi/2,2.9499555017208157,-2.9499555017208157) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,4.520751828515713,-4.520751828515713) q[3];
u3(pi/2,4.531433243537918,-4.531433243537918) q[3];
u3(pi/2,4.740663314266998,-4.740663314266998) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,3.169866987472101,-3.169866987472101) q[2];
u3(pi/2,0.10555751316061704,-0.10555751316061704) q[1];
rzz(-pi/2) q[2],q[1];
u3(pi/2,1.6763538399555136,-1.6763538399555136) q[1];
u3(pi/2,3.169866987472101,-3.169866987472101) q[2];
u3(pi/2,0.017592918860102842,-0.017592918860102842) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,1.5883892456549995,-1.5883892456549995) q[2];
u3(pi/2,4.719928802753305,-4.719928802753305) q[2];
u3(pi/2,1.6656724249333084,-1.6656724249333084) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3.637335974326262,-3.637335974326262) q[5];
u3(pi/2,0.4542742977090841,-0.4542742977090841) q[5];
u3(pi/2,6.269990618034509,-6.269990618034509) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,4.648928808782176,-4.648928808782176) q[4];
u3(pi/2,2.030725491280442,-2.030725491280442) q[5];
u3(pi/2,2.041406906302648,-2.041406906302648) q[5];
rzz(-pi/2) q[4],q[5];
u3(pi/2,0.470610579507751,-0.470610579507751) q[5];
u3(pi/2,0.4812919945299563,-0.4812919945299563) q[5];
u3(pi/2,1.4966547401701773,-1.4966547401701773) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,6.209043720554867,-6.209043720554867) q[4];
u3(pi/2,1.3898405899481245,-1.3898405899481245) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,6.102229570332814,-6.102229570332814) q[3];
u3(pi/2,3.067451066965074,-3.067451066965074) q[4];
u3(pi/2,3.0781324819872795,-3.0781324819872795) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,4.648928808782176,-4.648928808782176) q[4];
u3(pi/2,4.658981905273664,-4.658981905273664) q[4];
u3(pi/2,6.091548155310609,-6.091548155310609) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,0.09487609813841176,-0.09487609813841176) q[1];
u3(pi/2,3.5751324397851842,-3.5751324397851842) q[0];
u3(pi/2,5.154725226010132,-5.154725226010132) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,0.46809730538487915,-0.46809730538487915) q[0];
u3(pi/2,3.2917607824313855,-3.2917607824313855) q[1];
u3(pi/2,0.13948671381938682,-0.13948671381938682) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,1.7102830406142833,-1.7102830406142833) q[1];
u3(pi/2,4.841822597712589,-4.841822597712589) q[1];
u3(pi/2,0.4574158903626739,-0.4574158903626739) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,6.091548155310609,-6.091548155310609) q[3];
u3(pi/2,1.4174866052997146,-1.4174866052997146) q[3];
u3(pi/2,1.5261857111139214,-1.5261857111139214) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,6.204645490839842,-6.204645490839842) q[2];
u3(pi/2,3.0454599183899456,-3.0454599183899456) q[3];
u3(pi/2,3.0561413334121506,-3.0561413334121506) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,4.626937660207048,-4.626937660207048) q[3];
u3(pi/2,1.4746635915950488,-1.4746635915950488) q[3];
u3(pi/2,3.0737342522722537,-3.0737342522722537) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,1.502937925477357,-1.502937925477357) q[2];
u3(pi/2,1.7002299441227962,-1.7002299441227962) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,0.1294336173278995,-0.1294336173278995) q[1];
u3(pi/2,4.64453057906715,-4.64453057906715) q[2];
u3(pi/2,4.655211994089355,-4.655211994089355) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,6.226008320884252,-6.226008320884252) q[2];
u3(pi/2,3.0737342522722537,-3.0737342522722537) q[2];
u3(pi/2,3.28107936740918,-3.28107936740918) q[1];
rzz(-pi/2) q[2],q[1];
u3(pi/2,0.4812919945299563,-0.4812919945299563) q[5];
u3(pi/2,3.580787306561646,-3.580787306561646) q[5];
u3(pi/2,4.682858009440945,-4.682858009440945) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,3.0617962001886125,-3.0617962001886125) q[4];
u3(pi/2,5.157866818663722,-5.157866818663722) q[5];
u3(pi/2,5.16791991515521,-5.16791991515521) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,0.45553093477052,-0.45553093477052) q[5];
u3(pi/2,0.4662123497927253,-0.4662123497927253) q[5];
u3(pi/2,3.051114785166407,-3.051114785166407) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,1.4803184583715105,-1.4803184583715105) q[4];
u3(pi/2,4.616256245184842,-4.616256245184842) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,6.187052571979739,-6.187052571979739) q[3];
u3(pi/2,1.4803184583715105,-1.4803184583715105) q[4];
u3(pi/2,1.490999873393716,-1.490999873393716) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.0617962001886125,-3.0617962001886125) q[4];
u3(pi/2,6.192707438756201,-6.192707438756201) q[4];
u3(pi/2,3.0561413334121506,-3.0561413334121506) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,5.169804870747363,-5.169804870747363) q[0];
u3(pi/2,1.7102830406142833,-1.7102830406142833) q[1];
u3(pi/2,4.626937660207048,-4.626937660207048) q[3];
u3(pi,2.005592750051724,-2.005592750051724) q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
