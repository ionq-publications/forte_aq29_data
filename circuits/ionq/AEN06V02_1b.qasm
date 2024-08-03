OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u(pi/2,pi,-pi) q[5];
u(pi/2,0.7841415263360125,-0.7841415263360125) q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,-pi/8,pi/8) q[3];
u(pi/2,3.925734179925805,-3.925734179925805) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,-0.7866548004588843,0.7866548004588843) q[5];
rzz(pi/2) q[5],q[3];
u(pi/2,4.7111323433232535,-4.7111323433232535) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,-0.0012566370614359723,0.0012566370614359723) q[5];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[3],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[3],q[4];
u(pi/2,3.925734179925805,-3.925734179925805) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[3],q[4];
u(pi/2,-1.1787255636268905,1.1787255636268905) q[5];
rzz(-pi/2) q[3],q[5];
u(pi/2,0.39207076316800626,-0.39207076316800626) q[5];
u(pi/2,0,0) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,-pi/2,pi/2) q[5];
u(pi/2,3*pi/4,-3*pi/4) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,0,0) q[4];
u(pi/2,-1.4011503235010476,1.4011503235010476) q[5];
u(pi/2,3.8352563115024196,-3.8352563115024196) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,1.3087874994855078,-1.3087874994855078) q[5];
u(pi/2,-1.047406990706837,1.047406990706837) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,3.6649819896778526,-3.6649819896778526) q[5];
u(pi/2,1.3087874994855078,-1.3087874994855078) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[3],q[4];
u(pi/2,2.8795838262804043,-2.8795838262804043) q[5];
rzz(-pi/2) q[3],q[5];
u(pi/2,2.0941856628829565,-2.0941856628829565) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,2.0941856628829565,-2.0941856628829565) q[5];
u(pi/2,-0.654079590477395,0.654079590477395) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,4.058309389907295,-4.058309389907295) q[5];
u(pi/2,1.310044136546944,-1.310044136546944) q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[3],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[3],q[4];
u(pi/2,0.5246459731494957,-0.5246459731494957) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,-1.046150353645401,1.046150353645401) q[5];
u(pi/2,2.8808404633418405,-2.8808404633418405) q[5];
rzz(-pi/2) q[3],q[5];
u(pi/2,0,0) q[2];
rzz(pi/2) q[2],q[5];
u(pi/2,2.8808404633418405,-2.8808404633418405) q[5];
rzz(-pi/2) q[5],q[2];
u(pi/2,-1.046150353645401,1.046150353645401) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,0.5246459731494957,-0.5246459731494957) q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[2],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,1.310044136546944,-1.310044136546944) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,-0.652822953415959,0.652822953415959) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,4.059566026968731,-4.059566026968731) q[5];
u(pi/2,3.667495263800724,-3.667495263800724) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,-1.0448937165839651,1.0448937165839651) q[5];
u(pi/2,2.882097100403276,-2.882097100403276) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,-0.8746193947593984,0.8746193947593984) q[5];
u(pi/2,4.3611589217133515,-4.3611589217133515) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,-1.306274225362636,1.306274225362636) q[5];
u(pi/2,2.620716591624606,-2.620716591624606) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,1.049920264829709,-1.049920264829709) q[5];
u(pi/2,-1.306274225362636,1.306274225362636) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[2],q[4];
u(pi/2,3.4061147550220543,-3.4061147550220543) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,2.620716591624606,-2.620716591624606) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,2.620716591624606,-2.620716591624606) q[5];
u(pi/2,-0.12754866173574575,0.12754866173574575) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,1.4432476650591513,-1.4432476650591513) q[5];
u(pi/2,-1.3050175883012,1.3050175883012) q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[2],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[2],q[4];
rzz(-pi/2) q[2],q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[2],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,3.4073713920834896,-3.4073713920834896) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,1.444504302120587,-1.444504302120587) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,-0.12629202467430978,0.12629202467430978) q[5];
u(pi/2,-0.5183627878423158,0.5183627878423158) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,1.0524335389525805,-1.0524335389525805) q[5];
u(pi/2,-1.303760951239764,1.303760951239764) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,0,0) q[4];
u(pi/2,1.2220795422464295,-1.2220795422464295) q[5];
u(pi/2,0.17467255153959238,-0.17467255153959238) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,0.7904247116431922,-0.7904247116431922) q[5];
u(pi/2,-1.565769778549153,1.565769778549153) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,3.146619201835537,-3.146619201835537) q[5];
u(pi/2,0.7904247116431922,-0.7904247116431922) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[2],q[4];
u(pi/2,-0.7803716151517046,0.7803716151517046) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,-1.565769778549153,1.565769778549153) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,-1.565769778549153,1.565769778549153) q[5];
u(pi/2,1.9691502752700822,-1.9691502752700822) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,0.3983539484751859,-0.3983539484751859) q[5];
u(pi/2,3.9332740022944206,-3.9332740022944206) q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[2],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,0.0062831853071796395,-0.0062831853071796395) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,-1.564513141487717,1.564513141487717) q[5];
u(pi/2,2.3624776754995245,-2.3624776754995245) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,pi/2,-pi/2) q[1];
rzz(-pi/2) q[1],q[5];
u(pi/2,2.3624776754995245,-2.3624776754995245) q[5];
rzz(-pi/2) q[5],q[1];
u(pi/2,-1.564513141487717,1.564513141487717) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,0.0062831853071796395,-0.0062831853071796395) q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,3.9332740022944206,-3.9332740022944206) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,1.970406912331518,-1.970406912331518) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,0.39961058553662165,-0.39961058553662165) q[5];
u(pi/2,0.00753982236861539,-0.00753982236861539) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,-1.563256504426281,1.563256504426281) q[5];
u(pi/2,2.3637343125609602,-2.3637343125609602) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,0,0) q[4];
u(pi/2,1.748610470988079,-1.748610470988079) q[5];
u(pi/2,0.7012034802812415,-0.7012034802812415) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,1.3163273218541236,-1.3163273218541236) q[5];
u(pi/2,-1.0398671683382215,1.0398671683382215) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,3.672521812046468,-3.672521812046468) q[5];
u(pi/2,1.3163273218541236,-1.3163273218541236) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,-0.2544690049407734,0.2544690049407734) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,2.101725485251572,-2.101725485251572) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,2.101725485251572,-2.101725485251572) q[5];
u(pi/2,-0.6465397681087793,0.6465397681087793) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,4.065849212275911,-4.065849212275911) q[5];
u(pi/2,1.3175839589155594,-1.3175839589155594) q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[5];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,-0.25321236787933743,0.25321236787933743) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,0.9255131957475529,-0.9255131957475529) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,-0.6452831310473436,0.6452831310473436) q[5];
u(pi/2,-1.0373538942153497,1.0373538942153497) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,3.6750350861693395,-3.6750350861693395) q[5];
u(pi/2,1.318840595976995,-1.318840595976995) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,0.7037167544041134,-0.7037167544041134) q[5];
u(pi/2,-0.34369023630272344,0.34369023630272344) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,0.2720619238008761,-0.2720619238008761) q[5];
u(pi/2,4.199052740788117,-4.199052740788117) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,-0.5133362395965722,0.5133362395965722) q[5];
u(pi/2,3.413654577390669,-3.413654577390669) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,1.8428582505957727,-1.8428582505957727) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,4.199052740788117,-4.199052740788117) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,4.199052740788117,-4.199052740788117) q[5];
u(pi/2,1.4507874874277666,-1.4507874874277666) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,-0.12000883936713014,0.12000883936713014) q[5];
u(pi/2,3.414911214452105,-3.414911214452105) q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,1.8441148876572084,-1.8441148876572084) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,-0.11875220230569417,0.11875220230569417) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,4.593636778078995,-4.593636778078995) q[5];
u(pi/2,4.2015660149109895,-4.2015660149109895) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,2.630769688116093,-2.630769688116093) q[5];
u(pi/2,0.27457519792374807,-0.27457519792374807) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,0,0) q[4];
u(pi/2,-0.3411769621798515,0.3411769621798515) q[5];
u(pi/2,-1.3879556343559707,1.3879556343559707) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,-0.7728317927830891,0.7728317927830891) q[5];
u(pi/2,3.1541590242041524,-3.1541590242041524) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,1.5833626974092558,-1.5833626974092558) q[5];
u(pi/2,-0.7728317927830891,0.7728317927830891) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,3.9395571876016007,-3.9395571876016007) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,0.012566370614359279,-0.012566370614359279) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,0.012566370614359279,-0.012566370614359279) q[5];
u(pi/2,3.547486424433594,-3.547486424433594) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,1.9766900976386976,-1.9766900976386976) q[5];
u(pi/2,-0.7715751557216531,0.7715751557216531) q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[5];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0.7992211710732433,-0.7992211710732433) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,-1.1636459188896593,1.1636459188896593) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,3.5487430614950304,-3.5487430614950304) q[5];
u(pi/2,3.156672298327024,-3.156672298327024) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,1.5858759715321273,-1.5858759715321273) q[5];
u(pi/2,-0.7703185186602172,0.7703185186602172) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,-1.3854423602330987,1.3854423602330987) q[5];
u(pi/2,3.8503359562396504,-3.8503359562396504) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,4.465459797812532,-4.465459797812532) q[5];
u(pi/2,2.1092653076201873,-2.1092653076201873) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,0.5384689808252907,-0.5384689808252907) q[5];
u(pi/2,4.465459797812532,-4.465459797812532) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,2.894663471017635,-2.894663471017635) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,-1.032327345969606,1.032327345969606) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,2.1092653076201873,-2.1092653076201873) q[5];
u(pi/2,-0.6389999457401639,0.6389999457401639) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,0.9317963810547325,-0.9317963810547325) q[5];
u(pi/2,4.466716434873968,-4.466716434873968) q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,0.5397256178867265,-0.5397256178867265) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,-1.03107070890817,1.03107070890817) q[5];
u(pi/2,2.8959201080790713,-2.8959201080790713) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,pi,-pi) q[0];
rzz(pi/2) q[0],q[5];
u(pi/2,-0.24567254551072182,0.24567254551072182) q[5];
rzz(-pi/2) q[5],q[0];
u(pi/2,2.110521944681623,-2.110521944681623) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0.5397256178867265,-0.5397256178867265) q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,4.466716434873968,-4.466716434873968) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,2.503849344911065,-2.503849344911065) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0.9330530181161687,-0.9330530181161687) q[5];
u(pi/2,0.5409822549481627,-0.5409822549481627) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-1.0298140718467343,1.0298140718467343) q[5];
u(pi/2,2.8971767451405066,-2.8971767451405066) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,2.2820529035676254,-2.2820529035676254) q[5];
u(pi/2,1.2346459128607887,-1.2346459128607887) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,-1.291194580625405,1.291194580625405) q[5];
u(pi/2,2.635796236361836,-2.635796236361836) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,1.0649999095669398,-1.0649999095669398) q[5];
u(pi/2,-1.291194580625405,1.291194580625405) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,0.2796017461694915,-0.2796017461694915) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-0.5057964172279565,0.5057964172279565) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,2.635796236361836,-2.635796236361836) q[5];
u(pi/2,-0.11246901699851453,0.11246901699851453) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,1.4583273097963816,-1.4583273097963816) q[5];
u(pi/2,-1.289937943563969,1.289937943563969) q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[0],q[4];
rzz(pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0.2808583832309277,-0.2808583832309277) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,4.6011766004476105,-4.6011766004476105) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,3.030380273652714,-3.030380273652714) q[5];
u(pi/2,2.6383095104847083,-2.6383095104847083) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,4.209105837279605,-4.209105837279605) q[5];
u(pi/2,1.8529113470872605,-1.8529113470872605) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
u(pi/2,4.378751840573454,-4.378751840573454) q[5];
u(pi/2,3.331973168397335,-3.331973168397335) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,0.805504356380423,-0.805504356380423) q[5];
u(pi/2,-1.550690133811922,1.550690133811922) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0.020106192982974447,-0.020106192982974447) q[5];
u(pi/2,3.947097009970216,-3.947097009970216) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,2.3763006831753195,-2.3763006831753195) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,1.5909025197778712,-1.5909025197778712) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-1.550690133811922,1.550690133811922) q[5];
u(pi/2,1.984229920007313,-1.984229920007313) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0.4134335932124169,-0.4134335932124169) q[5];
u(pi/2,3.9483536470316514,-3.9483536470316514) q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(-pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,2.3775573202367553,-2.3775573202367553) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,0.41469023027385266,-0.41469023027385266) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-1.156106096521044,1.156106096521044) q[5];
u(pi/2,-1.54817685968905,1.54817685968905) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0.02261946710584639,-0.02261946710584639) q[5];
u(pi/2,3.9496102840930885,-3.9496102840930885) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,0.19289378893041342,-0.19289378893041342) q[5];
u(pi/2,-0.8545132017764236,0.8545132017764236) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,2.9028316119169686,-2.9028316119169686) q[5];
u(pi/2,0.5466371217246242,-0.5466371217246242) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-1.0241592050702726,1.0241592050702726) q[5];
u(pi/2,2.9028316119169686,-2.9028316119169686) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,1.3320352851220725,-1.3320352851220725) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0.5466371217246242,-0.5466371217246242) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,3.688229775314417,-3.688229775314417) q[5];
u(pi/2,0.939964521954066,-0.939964521954066) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-0.6308318048408303,0.6308318048408303) q[5];
u(pi/2,2.904088248978405,-2.904088248978405) q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,1.3332919221835082,-1.3332919221835082) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-0.6295751677793946,0.6295751677793946) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,4.082813812605296,-4.082813812605296) q[5];
u(pi/2,3.6907430494372893,-3.6907430494372893) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,2.1199467226423927,-2.1199467226423927) q[5];
u(pi/2,-0.23624776754995236,0.23624776754995236) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
u(pi/2,2.2895927259362407,-2.2895927259362407) q[5];
u(pi/2,1.2421857352294041,-1.2421857352294041) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,-1.2836547582567894,1.2836547582567894) q[5];
u(pi/2,2.6433360587304513,-2.6433360587304513) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,1.0725397319355556,-1.0725397319355556) q[5];
u(pi/2,-1.2836547582567894,1.2836547582567894) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3.4287342221278996,-3.4287342221278996) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,2.6433360587304513,-2.6433360587304513) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-0.49825659485934115,0.49825659485934115) q[5];
u(pi/2,3.036663458959894,-3.036663458959894) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,1.4658671321649974,-1.4658671321649974) q[5];
u(pi/2,-1.2823981211953535,1.2823981211953535) q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[5];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,3.4299908591893367,-3.4299908591893367) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,1.4671237692264332,-1.4671237692264332) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,3.03792009602133,-3.03792009602133) q[5];
u(pi/2,2.6458493328533237,-2.6458493328533237) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,1.0750530060584271,-1.0750530060584271) q[5];
u(pi/2,-1.2811414841339177,1.2811414841339177) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,1.244699009352276,-1.244699009352276) q[5];
u(pi/2,0.19792033717615687,-0.19792033717615687) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,3.9546368323388315,-3.9546368323388315) q[5];
u(pi/2,1.5984423421464866,-1.5984423421464866) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0.02764601535159028,-0.02764601535159028) q[5];
u(pi/2,3.9546368323388315,-3.9546368323388315) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,2.383840505543935,-2.383840505543935) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,1.5984423421464866,-1.5984423421464866) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-1.5431503114433063,1.5431503114433063) q[5];
u(pi/2,1.9917697423759284,-1.9917697423759284) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0.4209734155810323,-0.4209734155810323) q[5];
u(pi/2,3.9558934694002676,-3.9558934694002676) q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[0],q[4];
rzz(pi/2) q[0],q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-0.7564955109844222,0.7564955109844222) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,3.563822706232261,-3.563822706232261) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-1.1485662741524285,1.1485662741524285) q[5];
u(pi/2,-1.5406370373204346,1.5406370373204346) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,3.1717519430642556,-3.1717519430642556) q[5];
u(pi/2,0.8155574528719103,-0.8155574528719103) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0,0) q[4];
u(pi/2,3.3420262648888226,-3.3420262648888226) q[5];
u(pi/2,2.2946192741819846,-2.2946192741819846) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,-0.2312212193042087,0.2312212193042087) q[5];
u(pi/2,3.6957695976830323,-3.6957695976830323) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,2.1249732708881357,-2.1249732708881357) q[5];
u(pi/2,-0.2312212193042087,0.2312212193042087) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,1.3395751074906879,-1.3395751074906879) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0.5541769440932396,-0.5541769440932396) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,3.6957695976830323,-3.6957695976830323) q[5];
u(pi/2,0.9475043443226814,-0.9475043443226814) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-0.6232919824722151,0.6232919824722151) q[5];
u(pi/2,2.91162807134702,-2.91162807134702) q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(-pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,1.3408317445521236,-1.3408317445521236) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0,0) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-0.6220353454107791,0.6220353454107791) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,4.090353634973911,-4.090353634973911) q[5];
u(pi/2,3.6982828718059046,-3.6982828718059046) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,2.127486545011008,-2.127486545011008) q[5];
u(pi/2,-0.22870794518133697,0.22870794518133697) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,2.297132548304857,-2.297132548304857) q[5];
u(pi/2,1.24972555759802,-1.24972555759802) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,-1.276114935888174,1.276114935888174) q[5];
u(pi/2,2.6508758810990676,-2.6508758810990676) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,4.221672207893964,-4.221672207893964) q[5];
u(pi/2,1.8654777177016193,-1.8654777177016193) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,0.2946813909067225,-0.2946813909067225) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-0.49071677249072576,0.49071677249072576) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,2.6508758810990676,-2.6508758810990676) q[5];
u(pi/2,-0.09738937226128375,0.09738937226128375) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,4.614999608123406,-4.614999608123406) q[5];
u(pi/2,1.866734354763055,-1.866734354763055) q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[0],q[4];
rzz(pi/2) q[0],q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0.29593802796815827,-0.29593802796815827) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,4.616256245184842,-4.616256245184842) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,3.0454599183899456,-3.0454599183899456) q[5];
u(pi/2,2.653389155221939,-2.653389155221939) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,4.224185482016836,-4.224185482016836) q[5];
u(pi/2,1.8679909918244912,-1.8679909918244912) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0,0) q[4];
u(pi/2,4.393831485310685,-4.393831485310685) q[5];
u(pi/2,3.3470528131345656,-3.3470528131345656) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,0.8205840011176537,-0.8205840011176537) q[5];
u(pi/2,-1.5356104890746909,1.5356104890746909) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0.03518583772020567,-0.03518583772020567) q[5];
u(pi/2,3.9621766547074477,-3.9621766547074477) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,2.3913803279125507,-2.3913803279125507) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,1.6059821645151025,-1.6059821645151025) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-1.5356104890746909,1.5356104890746909) q[5];
u(pi/2,1.9993095647445447,-1.9993095647445447) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0.4285132379496477,-0.4285132379496477) q[5];
u(pi/2,3.963433291768883,-3.963433291768883) q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,3.1780351283714348,-3.1780351283714348) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,1.6072388015765382,-1.6072388015765382) q[5];
u(pi/2,-0.7489556886158066,0.7489556886158066) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[0];
u(pi/2,-pi/2,pi/2) q[1];
rzz(pi/2) q[0],q[1];
u(pi/2,0,0) q[1];
u(pi/2,5*pi/4,-5*pi/4) q[1];
rzz(pi/2) q[0],q[1];
u(pi/2,pi,-pi) q[2];
rzz(-pi/2) q[0],q[2];
u(pi/2,pi/2,-pi/2) q[2];
u(pi/2,-3*pi/8,3*pi/8) q[2];
rzz(pi/2) q[0],q[2];
u(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(pi/2) q[0],q[3];
u(pi/2,7*pi/8,-7*pi/8) q[3];
u(pi/2,-0.19666370011472112,0.19666370011472112) q[3];
rzz(pi/2) q[0],q[3];
u(pi/2,-pi/4,pi/4) q[1];
u(pi/2,-pi/2,pi/2) q[1];
rzz(pi/2) q[1],q[2];
u(pi/2,-3*pi/8,3*pi/8) q[2];
u(pi/2,7*pi/8,-7*pi/8) q[2];
rzz(-pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u(pi/2,-0.19603538158400302,0.19603538158400302) q[3];
u(pi/2,3.337628035173796,-3.337628035173796) q[3];
rzz(pi/2) q[1],q[3];
u(pi/2,11*pi/8,-11*pi/8) q[2];
u(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(pi/2) q[2],q[3];
u(pi/2,3.337628035173796,-3.337628035173796) q[3];
u(pi/2,0.9814335449814515,-0.9814335449814515) q[3];
rzz(pi/2) q[2],q[3];
u(pi,-pi/2,pi/2) q[0];
u(pi,pi/2,-pi/2) q[1];
u(pi,0,0) q[2];
u(pi/2,2.5522298717763476,-2.5522298717763476) q[3];
u(pi/2,pi/8,-pi/8) q[3];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,0.8218406381790899,-0.8218406381790899) q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
