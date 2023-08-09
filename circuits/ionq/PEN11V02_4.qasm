OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[10];
rzz(-pi/2) q[9],q[10];
u3(pi/2,0.5887344632827273,-0.5887344632827273) q[10];
u3(pi/2,3.1170882308917927,-3.1170882308917927) q[10];
rzz(pi/2) q[9],q[10];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[8];
rzz(pi/2) q[8],q[10];
u3(pi/2,3.1170882308917927,-3.1170882308917927) q[10];
u3(pi/2,5.031574793989412,-5.031574793989412) q[10];
rzz(pi/2) q[8],q[10];
u3(pi/2,4.246176630591964,-4.246176630591964) q[7];
rzz(pi/2) q[7],q[10];
u3(pi/2,5.031574793989412,-5.031574793989412) q[10];
u3(pi/2,5.718955266594859,-5.718955266594859) q[10];
rzz(-pi/2) q[7],q[10];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
rzz(pi/2) q[6],q[10];
u3(pi/2,5.718955266594859,-5.718955266594859) q[10];
u3(pi/2,1.2026016677941727,-1.2026016677941727) q[10];
rzz(pi/2) q[6],q[10];
u3(pi/2,2.94555727200579,-2.94555727200579) q[5];
rzz(pi/2) q[5],q[10];
u3(pi/2,1.2026016677941727,-1.2026016677941727) q[10];
u3(pi/2,1.595300749492897,-1.595300749492897) q[10];
rzz(-pi/2) q[5],q[10];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[4],q[10];
u3(pi/2,1.595300749492897,-1.595300749492897) q[10];
u3(pi/2,3.951495239685242,-3.951495239685242) q[10];
rzz(pi/2) q[4],q[10];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,5.522291566480138,-5.522291566480138) q[10];
rzz(pi/2) q[3],q[10];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
u3(pi/2,7*pi/4,-7*pi/4) q[1];
rzz(-pi/2) q[0],q[1];
rzz(pi/2) q[0],q[2];
u3(pi/2,0,0) q[2];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
u3(pi/2,5.301751762198135,-5.301751762198135) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,2.944928953475072,-2.944928953475072) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[5];
u3(pi/2,4.466716434873968,-4.466716434873968) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,5.301123443667417,-5.301123443667417) q[6];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,4.246176630591964,-4.246176630591964) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,2.675380303797068,-2.675380303797068) q[7];
u3(pi/2,5.804406586772502,-5.804406586772502) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[8];
rzz(pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[8];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
rzz(-pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[9];
u3(pi/2,pi/4,-pi/4) q[1];
u3(pi/2,pi,-pi) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
u3(pi/2,13*pi/8,-13*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,5.301751762198135,-5.301751762198135) q[3];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[4];
u3(pi/2,5.203105752875415,-5.203105752875415) q[4];
rzz(-pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,4.466716434873968,-4.466716434873968) q[5];
u3(pi/2,1.2271060904921731,-1.2271060904921731) q[5];
rzz(pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[6];
u3(pi/2,5.227610175573416,-5.227610175573416) q[6];
rzz(-pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u3(pi/2,5.804406586772502,-5.804406586772502) q[7];
u3(pi/2,2.6383095104847083,-2.6383095104847083) q[7];
rzz(pi/2) q[1],q[7];
rzz(pi/2) q[1],q[8];
u3(pi/2,2.1048670779051615,-2.1048670779051615) q[8];
u3(pi/2,5.233893360880595,-5.233893360880595) q[8];
rzz(pi/2) q[1],q[8];
rzz(-pi/2) q[1],q[9];
rzz(pi/2) q[1],q[9];
u3(pi/2,pi/8,-pi/8) q[2];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[3];
u3(pi/2,4.123026198571244,-4.123026198571244) q[3];
rzz(pi/2) q[2],q[3];
rzz(-pi/2) q[2],q[4];
u3(pi/2,5.203105752875415,-5.203105752875415) q[4];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,1.2271060904921731,-1.2271060904921731) q[5];
u3(pi/2,4.172663362497963,-4.172663362497963) q[5];
rzz(pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[6];
u3(pi/2,5.227610175573416,-5.227610175573416) q[6];
u3(pi/2,1.9879998311916212,-1.9879998311916212) q[6];
rzz(pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,2.6383095104847083,-2.6383095104847083) q[7];
u3(pi/2,5.7308933186785005,-5.7308933186785005) q[7];
rzz(pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[8];
u3(pi/2,2.092300707290802,-2.092300707290802) q[8];
u3(pi/2,5.209388938182594,-5.209388938182594) q[8];
rzz(pi/2) q[2],q[8];
rzz(pi/2) q[2],q[9];
u3(pi/2,2.0954422999443922,-2.0954422999443922) q[9];
u3(pi/2,5.224468582919826,-5.224468582919826) q[9];
rzz(pi/2) q[2],q[9];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[4];
u3(pi/2,4.025008507779242,-4.025008507779242) q[4];
rzz(-pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,4.172663362497963,-4.172663362497963) q[5];
u3(pi/2,0.638371627209446,-0.638371627209446) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,1.9879998311916212,-1.9879998311916212) q[6];
u3(pi/2,4.933557103197411,-4.933557103197411) q[6];
rzz(-pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,5.7308933186785005,-5.7308933186785005) q[7];
u3(pi/2,2.491282974296706,-2.491282974296706) q[7];
rzz(pi/2) q[3],q[7];
rzz(pi/2) q[3],q[8];
u3(pi/2,5.209388938182594,-5.209388938182594) q[8];
u3(pi/2,2.018787439196801,-2.018787439196801) q[8];
rzz(pi/2) q[3],q[8];
rzz(-pi/2) q[3],q[9];
u3(pi/2,2.082875929330033,-2.082875929330033) q[9];
u3(pi/2,5.199964160221826,-5.199964160221826) q[9];
rzz(pi/2) q[3],q[9];
u3(pi/2,5.595804834574139,-5.595804834574139) q[4];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,0.638371627209446,-0.638371627209446) q[5];
u3(pi/2,2.994566117401791,-2.994566117401791) q[5];
rzz(pi/2) q[4],q[5];
rzz(-pi/2) q[4],q[6];
u3(pi/2,4.933557103197411,-4.933557103197411) q[6];
u3(pi/2,1.399265367908894,-1.399265367908894) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,2.491282974296706,-2.491282974296706) q[7];
u3(pi/2,5.436211927771778,-5.436211927771778) q[7];
rzz(pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,2.018787439196801,-2.018787439196801) q[8];
u3(pi/2,5.0623624019945925,-5.0623624019945925) q[8];
rzz(pi/2) q[4],q[8];
rzz(pi/2) q[4],q[9];
u3(pi/2,5.199964160221826,-5.199964160221826) q[9];
u3(pi/2,2.0093626612360316,-2.0093626612360316) q[9];
rzz(pi/2) q[4],q[9];
u3(pi/2,4.565362444196688,-4.565362444196688) q[5];
u3(pi/2,5.693822525366141,-5.693822525366141) q[5];
rzz(-pi/2) q[5],q[6];
u3(pi/2,4.540858021498687,-4.540858021498687) q[6];
u3(pi/2,0.6138672045114455,-0.6138672045114455) q[6];
rzz(pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,5.436211927771778,-5.436211927771778) q[7];
u3(pi/2,1.901920192483261,-1.901920192483261) q[7];
rzz(pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,5.0623624019945925,-5.0623624019945925) q[8];
u3(pi/2,1.7241060482900783,-1.7241060482900783) q[8];
rzz(pi/2) q[5],q[8];
rzz(pi/2) q[5],q[9];
u3(pi/2,2.0093626612360316,-2.0093626612360316) q[9];
u3(pi/2,5.052937624033824,-5.052937624033824) q[9];
rzz(pi/2) q[5],q[9];
u3(pi/2,5.326256184896136,-5.326256184896136) q[6];
u3(pi/2,5.546795989178139,-5.546795989178139) q[6];
rzz(pi/2) q[6],q[7];
u3(pi/2,1.901920192483261,-1.901920192483261) q[7];
u3(pi/2,4.258114682675606,-4.258114682675606) q[7];
rzz(-pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,1.7241060482900783,-1.7241060482900783) q[8];
u3(pi/2,4.472999620181147,-4.472999620181147) q[8];
rzz(pi/2) q[6],q[8];
rzz(pi/2) q[6],q[9];
u3(pi/2,5.052937624033824,-5.052937624033824) q[9];
u3(pi/2,1.7146812703293088,-1.7146812703293088) q[9];
rzz(-pi/2) q[6],q[9];
u3(pi/2,2.687318355880709,-2.687318355880709) q[7];
u3(pi/2,4.933557103197411,-4.933557103197411) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,4.472999620181147,-4.472999620181147) q[8];
u3(pi/2,0.5460088031939061,-0.5460088031939061) q[8];
rzz(pi/2) q[7],q[8];
rzz(pi/2) q[7],q[9];
u3(pi/2,4.856273923919102,-4.856273923919102) q[9];
u3(pi/2,1.321982188630585,-1.321982188630585) q[9];
rzz(pi/2) q[7],q[9];
u3(pi/2,5.258397783578595,-5.258397783578595) q[8];
u3(pi/2,1.7793980789932589,-1.7793980789932589) q[8];
rzz(-pi/2) q[8],q[9];
u3(pi/2,4.463574842220378,-4.463574842220378) q[9];
u3(pi/2,0.5365840252331366,-0.5365840252331366) q[9];
rzz(pi/2) q[8],q[9];
u3(pi,3*pi/2,-3*pi/2) q[1];
u3(pi,3*pi/4,-3*pi/4) q[2];
u3(pi,9*pi/8,-9*pi/8) q[3];
u3(pi,3*pi/2,-3*pi/2) q[4];
u3(pi,0.09801769079200154,-0.09801769079200154) q[5];
u3(pi,9*pi/8,-9*pi/8) q[6];
u3(pi,pi/4,-pi/4) q[7];
u3(pi,2.8469112626830704,-2.8469112626830704) q[8];
u3(pi/2,2.107380352028033,-2.107380352028033) q[9];
u3(pi/2,3.0617962001886125,-3.0617962001886125) q[9];
u3(pi,1.6198051721908973,-1.6198051721908973) q[10];
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
