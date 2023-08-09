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
u3(pi/2,1.6311149057438206,-1.6311149057438206) q[1];
rzz(-pi/2) q[2],q[1];
u3(pi/2,3.2019112325387176,-3.2019112325387176) q[1];
u3(pi/2,4.716158891568997,-4.716158891568997) q[2];
u3(pi/2,1.563884822956999,-1.563884822956999) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,3.134681149751896,-3.134681149751896) q[2];
u3(pi/2,3.1453625647741013,-3.1453625647741013) q[2];
u3(pi/2,0.07099999397112931,-0.07099999397112931) q[1];
rzz(pi/2) q[2],q[1];
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
u3(pi/2,4.783388974355819,-4.783388974355819) q[1];
u3(pi/2,5.121424343882081,-5.121424343882081) q[0];
u3(pi/2,0.41846014145816046,-0.41846014145816046) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,5.156610181602287,-5.156610181602287) q[0];
u3(pi/2,4.838681005059,-4.838681005059) q[1];
u3(pi/2,1.6864069364470011,-1.6864069364470011) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,3.2572032632418972,-3.2572032632418972) q[1];
u3(pi/2,0.10555751316061704,-0.10555751316061704) q[1];
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
u3(pi/2,3.2471501667504103,-3.2471501667504103) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,1.6763538399555136,-1.6763538399555136) q[1];
u3(pi/2,0.028274333882308135,-0.028274333882308135) q[2];
u3(pi/2,3.1591855724498963,-3.1591855724498963) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,4.729981899244793,-4.729981899244793) q[2];
u3(pi/2,1.578336149163512,-1.578336149163512) q[2];
u3(pi/2,1.6656724249333084,-1.6656724249333084) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,0.4957433207364693,-0.4957433207364693) q[5];
u3(pi/2,3.5958669512988775,-3.5958669512988775) q[5];
u3(pi/2,6.269990618034509,-6.269990618034509) q[4];
rzz(-pi/2) q[5],q[4];
u3(pi/2,1.5073361551923827,-1.5073361551923827) q[4];
u3(pi/2,2.030725491280442,-2.030725491280442) q[5];
u3(pi/2,2.041406906302648,-2.041406906302648) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,3.612203233097544,-3.612203233097544) q[5];
u3(pi/2,3.6228846481197494,-3.6228846481197494) q[5];
u3(pi/2,1.4966547401701773,-1.4966547401701773) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,6.209043720554867,-6.209043720554867) q[4];
u3(pi/2,4.531433243537918,-4.531433243537918) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,2.960636916743021,-2.960636916743021) q[3];
u3(pi/2,3.067451066965074,-3.067451066965074) q[4];
u3(pi/2,3.0781324819872795,-3.0781324819872795) q[4];
rzz(-pi/2) q[3],q[4];
u3(pi/2,1.5073361551923827,-1.5073361551923827) q[4];
u3(pi/2,1.51738925168387,-1.51738925168387) q[4];
u3(pi/2,6.091548155310609,-6.091548155310609) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,0.09487609813841176,-0.09487609813841176) q[1];
u3(pi/2,0.43353978619539146,-0.43353978619539146) q[0];
u3(pi/2,2.0131325724203397,-2.0131325724203397) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.6096899589746725,-3.6096899589746725) q[0];
u3(pi/2,3.2917607824313855,-3.2917607824313855) q[1];
u3(pi/2,0.13948671381938682,-0.13948671381938682) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,1.7102830406142833,-1.7102830406142833) q[1];
u3(pi/2,4.841822597712589,-4.841822597712589) q[1];
u3(pi/2,3.599008543952467,-3.599008543952467) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,2.9499555017208157,-2.9499555017208157) q[3];
u3(pi/2,4.559079258889508,-4.559079258889508) q[3];
u3(pi/2,4.6677783647037145,-4.6677783647037145) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,3.0630528372500483,-3.0630528372500483) q[2];
u3(pi/2,6.187052571979739,-6.187052571979739) q[3];
u3(pi/2,6.1977339870019446,-6.1977339870019446) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,1.4853450066172542,-1.4853450066172542) q[3];
u3(pi/2,4.616256245184842,-4.616256245184842) q[3];
u3(pi/2,6.215326905862047,-6.215326905862047) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,4.64453057906715,-4.64453057906715) q[2];
u3(pi/2,4.841822597712589,-4.841822597712589) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3.2710262709176923,-3.2710262709176923) q[1];
u3(pi/2,1.502937925477357,-1.502937925477357) q[2];
u3(pi/2,1.5136193404995624,-1.5136193404995624) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,3.0844156672944587,-3.0844156672944587) q[2];
u3(pi/2,6.215326905862047,-6.215326905862047) q[2];
u3(pi/2,0.13948671381938682,-0.13948671381938682) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3.6228846481197494,-3.6228846481197494) q[5];
u3(pi/2,0.4391946529718531,-0.4391946529718531) q[5];
u3(pi/2,4.682858009440945,-4.682858009440945) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,3.0617962001886125,-3.0617962001886125) q[4];
u3(pi/2,2.0162741650739293,-2.0162741650739293) q[5];
u3(pi/2,2.0263272615654166,-2.0263272615654166) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,3.597123588360313,-3.597123588360313) q[5];
u3(pi/2,3.6078050033825186,-3.6078050033825186) q[5];
u3(pi/2,3.051114785166407,-3.051114785166407) q[4];
rzz(-pi/2) q[5],q[4];
u3(pi/2,4.621911111961304,-4.621911111961304) q[4];
u3(pi/2,1.4746635915950488,-1.4746635915950488) q[3];
rzz(pi/2) q[4],q[3];
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
u3(pi,0.4347964232568273,-0.4347964232568273) q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
