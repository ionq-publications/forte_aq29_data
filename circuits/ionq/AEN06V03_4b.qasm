OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u(pi/2,pi,-pi) q[5];
u(pi/2,2.7495218904217866,-2.7495218904217866) q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/8,pi/8) q[3];
u(pi/2,2.7495218904217866,-2.7495218904217866) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,4.320318217216683,-4.320318217216683) q[5];
rzz(-pi/2) q[5],q[3];
u(pi/2,0.393327400229442,-0.393327400229442) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,1.9641237270243388,-1.9641237270243388) q[5];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[3],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[3],q[4];
u(pi/2,2.7495218904217866,-2.7495218904217866) q[5];
rzz(-pi/2) q[3],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[3],q[4];
u(pi/2,-0.19603538158400302,0.19603538158400302) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,4.516353598800687,-4.516353598800687) q[5];
u(pi/2,pi,-pi) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,pi/2,-pi/2) q[5];
u(pi/2,-pi/4,pi/4) q[5];
rzz(-pi/2) q[3],q[5];
u(pi/2,0,0) q[4];
u(pi/2,-1.4011503235010476,1.4011503235010476) q[5];
u(pi/2,3.8352563115024196,-3.8352563115024196) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,4.450380153075301,-4.450380153075301) q[5];
u(pi/2,2.0941856628829565,-2.0941856628829565) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,0.5233893360880595,-0.5233893360880595) q[5];
u(pi/2,4.450380153075301,-4.450380153075301) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[3],q[4];
u(pi/2,2.8795838262804043,-2.8795838262804043) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,-1.047406990706837,1.047406990706837) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,2.0941856628829565,-2.0941856628829565) q[5];
u(pi/2,-0.4580442088933918,0.4580442088933918) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,4.254344771491298,-4.254344771491298) q[5];
u(pi/2,1.70211489971495,-1.70211489971495) q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[3],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[3],q[4];
u(pi/2,4.058309389907295,-4.058309389907295) q[5];
rzz(-pi/2) q[3],q[5];
u(pi/2,-0.654079590477395,0.654079590477395) q[5];
u(pi/2,3.272911226509847,-3.272911226509847) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,0,0) q[2];
rzz(pi/2) q[2],q[5];
u(pi/2,0.13131857292005322,-0.13131857292005322) q[5];
rzz(-pi/2) q[5],q[2];
u(pi/2,2.4875130631123987,-2.4875130631123987) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,4.058309389907295,-4.058309389907295) q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-1.4394777538748431,1.4394777538748431) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,1.8981502812989532,-1.8981502812989532) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,0.32735395450405624,-0.32735395450405624) q[5];
u(pi/2,-1.047406990706837,1.047406990706837) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,3.6649819896778526,-3.6649819896778526) q[5];
u(pi/2,1.3087874994855078,-1.3087874994855078) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,3.8352563115024196,-3.8352563115024196) q[5];
u(pi/2,2.7878493207955826,-2.7878493207955826) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,0.2620088273093888,-0.2620088273093888) q[5];
u(pi/2,4.18899964429663,-4.18899964429663) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,2.6182033175017336,-2.6182033175017336) q[5];
u(pi/2,0.2620088273093888,-0.2620088273093888) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,-1.3087874994855078,1.3087874994855078) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,4.18899964429663,-4.18899964429663) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,4.18899964429663,-4.18899964429663) q[5];
u(pi/2,1.636141453989564,-1.636141453989564) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,3.206937780784461,-3.206937780784461) q[5];
u(pi/2,0.6547079090081129,-0.6547079090081129) q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[2],q[4];
rzz(-pi/2) q[2],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-0.9160884177867836,0.9160884177867836) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,2.4215396173870123,-2.4215396173870123) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,0.8507432905921162,-0.8507432905921162) q[5];
u(pi/2,-0.5233893360880597,0.5233893360880597) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,1.047406990706837,-1.047406990706837) q[5];
u(pi/2,-1.3087874994855078,1.3087874994855078) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,0,0) q[4];
u(pi/2,1.2170529940006856,-1.2170529940006856) q[5];
u(pi/2,0.16964600329384893,-0.16964600329384893) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[5];
u(pi/2,pi/2,-pi/2) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,0,0) q[5];
u(pi/2,5*pi/4,-5*pi/4) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,-pi/2,pi/2) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,-pi/2,pi/2) q[5];
u(pi/2,2.1601591086083416,-2.1601591086083416) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,0.589362781813445,-0.589362781813445) q[5];
u(pi/2,4.320318217216683,-4.320318217216683) q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[2],q[4];
u(pi/2,3.534920053819235,-3.534920053819235) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,1.9641237270243388,-1.9641237270243388) q[5];
u(pi/2,-0.39207076316800626,0.39207076316800626) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[1],q[5];
u(pi/2,-0.39207076316800626,0.39207076316800626) q[5];
rzz(pi/2) q[5],q[1];
u(pi/2,-1.1774689265654543,1.1774689265654543) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,3.534920053819235,-3.534920053819235) q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,1.1787255636268905,-1.1787255636268905) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,4.516353598800687,-4.516353598800687) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,2.9455572720057903,-2.9455572720057903) q[5];
u(pi/2,pi/2,-pi/2) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,pi,-pi) q[5];
u(pi/2,pi/4,-pi/4) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,pi,-pi) q[4];
u(pi/2,3.311238656883642,-3.311238656883642) q[5];
u(pi/2,2.264459984707523,-2.264459984707523) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,2.8795838262804043,-2.8795838262804043) q[5];
u(pi/2,0.5233893360880595,-0.5233893360880595) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,2.0941856628829565,-2.0941856628829565) q[5];
u(pi/2,-0.26200882730938857,0.26200882730938857) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,4.450380153075301,-4.450380153075301) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,3.6649819896778526,-3.6649819896778526) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,0.5233893360880595,-0.5233893360880595) q[5];
u(pi/2,4.254344771491298,-4.254344771491298) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,-0.4580442088933918,0.4580442088933918) q[5];
u(pi/2,3.272911226509847,-3.272911226509847) q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,1.70211489971495,-1.70211489971495) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,1.8981502812989532,-1.8981502812989532) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,0.32735395450405624,-0.32735395450405624) q[5];
u(pi/2,-1.047406990706837,1.047406990706837) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,3.6649819896778526,-3.6649819896778526) q[5];
u(pi/2,1.3087874994855078,-1.3087874994855078) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,0.6936636579126265,-0.6936636579126265) q[5];
u(pi/2,-0.35374333279421055,0.35374333279421055) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,0.2620088273093888,-0.2620088273093888) q[5];
u(pi/2,4.18899964429663,-4.18899964429663) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,-0.5233893360880597,0.5233893360880597) q[5];
u(pi/2,3.403601480899182,-3.403601480899182) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,1.832805154104285,-1.832805154104285) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,1.047406990706837,-1.047406990706837) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,4.18899964429663,-4.18899964429663) q[5];
u(pi/2,1.636141453989564,-1.636141453989564) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,3.206937780784461,-3.206937780784461) q[5];
u(pi/2,0.6547079090081129,-0.6547079090081129) q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-0.9160884177867836,0.9160884177867836) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,-0.7200530362027806,0.7200530362027806) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,3.9923359441819093,-3.9923359441819093) q[5];
u(pi/2,2.6182033175017336,-2.6182033175017336) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,1.047406990706837,-1.047406990706837) q[5];
u(pi/2,-1.3087874994855078,1.3087874994855078) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,0,0) q[4];
u(pi/2,1.2170529940006856,-1.2170529940006856) q[5];
u(pi/2,0.16964600329384893,-0.16964600329384893) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,pi/4,-pi/4) q[5];
u(pi/2,-pi/2,pi/2) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,0,0) q[5];
u(pi/2,5*pi/4,-5*pi/4) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,-pi/2,pi/2) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,pi/2,-pi/2) q[5];
u(pi/2,-0.9814335449814514,0.9814335449814514) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,0.589362781813445,-0.589362781813445) q[5];
u(pi/2,4.320318217216683,-4.320318217216683) q[5];
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
rzz(pi/2) q[1],q[5];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-0.39207076316800626,0.39207076316800626) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,-0.19603538158400302,0.19603538158400302) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,4.516353598800687,-4.516353598800687) q[5];
u(pi/2,pi,-pi) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,pi/2,-pi/2) q[5];
u(pi/2,-pi/4,pi/4) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,1.7404423300887455,-1.7404423300887455) q[5];
u(pi/2,0.6936636579126265,-0.6936636579126265) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,1.3087874994855078,-1.3087874994855078) q[5];
u(pi/2,-1.047406990706837,1.047406990706837) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,0.5233893360880595,-0.5233893360880595) q[5];
u(pi/2,4.450380153075301,-4.450380153075301) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,2.8795838262804043,-2.8795838262804043) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,2.0941856628829565,-2.0941856628829565) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,-1.047406990706837,1.047406990706837) q[5];
u(pi/2,2.6835484446964015,-2.6835484446964015) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,4.254344771491298,-4.254344771491298) q[5];
u(pi/2,1.70211489971495,-1.70211489971495) q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,0.9167167363175013,-0.9167167363175013) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,2.4875130631123987,-2.4875130631123987) q[5];
u(pi/2,0.13131857292005322,-0.13131857292005322) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,pi,-pi) q[0];
rzz(-pi/2) q[0],q[5];
u(pi/2,3.272911226509847,-3.272911226509847) q[5];
rzz(pi/2) q[5],q[0];
u(pi/2,2.4875130631123987,-2.4875130631123987) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,4.058309389907295,-4.058309389907295) q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,1.70211489971495,-1.70211489971495) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,1.8981502812989532,-1.8981502812989532) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,3.46894660809385,-3.46894660809385) q[5];
u(pi/2,2.0941856628829565,-2.0941856628829565) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0.5233893360880595,-0.5233893360880595) q[5];
u(pi/2,4.450380153075301,-4.450380153075301) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,0.6936636579126265,-0.6936636579126265) q[5];
u(pi/2,-0.35374333279421055,0.35374333279421055) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,3.403601480899182,-3.403601480899182) q[5];
u(pi/2,1.047406990706837,-1.047406990706837) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,2.6182033175017336,-2.6182033175017336) q[5];
u(pi/2,0.2620088273093888,-0.2620088273093888) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,1.832805154104285,-1.832805154104285) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,1.047406990706837,-1.047406990706837) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,1.047406990706837,-1.047406990706837) q[5];
u(pi/2,-1.5054511996002289,1.5054511996002289) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,3.206937780784461,-3.206937780784461) q[5];
u(pi/2,0.6547079090081129,-0.6547079090081129) q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[5],q[4];
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
u(pi/2,-0.9160884177867836,0.9160884177867836) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-0.7200530362027806,0.7200530362027806) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0.8507432905921162,-0.8507432905921162) q[5];
u(pi/2,-0.5233893360880597,0.5233893360880597) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,4.18899964429663,-4.18899964429663) q[5];
u(pi/2,1.832805154104285,-1.832805154104285) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
u(pi/2,4.358645647590479,-4.358645647590479) q[5];
u(pi/2,3.311238656883642,-3.311238656883642) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[5];
u(pi/2,pi/2,-pi/2) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0,0) q[5];
u(pi/2,5*pi/4,-5*pi/4) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-pi/4,pi/4) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-pi/2,pi/2) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[5];
u(pi/2,-0.9814335449814514,0.9814335449814514) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0.589362781813445,-0.589362781813445) q[5];
u(pi/2,4.320318217216683,-4.320318217216683) q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[5];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-0.39207076316800626,0.39207076316800626) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,-0.19603538158400302,0.19603538158400302) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,1.3747609452108933,-1.3747609452108933) q[5];
u(pi/2,0,0) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-pi/2,pi/2) q[5];
u(pi/2,3*pi/4,-3*pi/4) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,1.7404423300887455,-1.7404423300887455) q[5];
u(pi/2,0.6936636579126265,-0.6936636579126265) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,1.3087874994855078,-1.3087874994855078) q[5];
u(pi/2,-1.047406990706837,1.047406990706837) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,3.6649819896778526,-3.6649819896778526) q[5];
u(pi/2,1.3087874994855078,-1.3087874994855078) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,2.8795838262804043,-2.8795838262804043) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-1.047406990706837,1.047406990706837) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-1.047406990706837,1.047406990706837) q[5];
u(pi/2,2.6835484446964015,-2.6835484446964015) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,1.1127521179015045,-1.1127521179015045) q[5];
u(pi/2,-1.4394777538748431,1.4394777538748431) q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0.13131857292005322,-0.13131857292005322) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,0.32735395450405624,-0.32735395450405624) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,1.8981502812989532,-1.8981502812989532) q[5];
u(pi/2,0.5233893360880595,-0.5233893360880595) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,2.0941856628829565,-2.0941856628829565) q[5];
u(pi/2,-0.26200882730938857,0.26200882730938857) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
u(pi/2,2.264459984707523,-2.264459984707523) q[5];
u(pi/2,1.2170529940006856,-1.2170529940006856) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,1.832805154104285,-1.832805154104285) q[5];
u(pi/2,-0.5233893360880597,0.5233893360880597) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,4.18899964429663,-4.18899964429663) q[5];
u(pi/2,1.832805154104285,-1.832805154104285) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,3.403601480899182,-3.403601480899182) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-0.5233893360880597,0.5233893360880597) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-0.5233893360880597,0.5233893360880597) q[5];
u(pi/2,3.206937780784461,-3.206937780784461) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,1.636141453989564,-1.636141453989564) q[5];
u(pi/2,-0.9160884177867836,0.9160884177867836) q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[5];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0.6547079090081129,-0.6547079090081129) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,3.9923359441819093,-3.9923359441819093) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-0.7200530362027806,0.7200530362027806) q[5];
u(pi/2,4.18899964429663,-4.18899964429663) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,2.6182033175017336,-2.6182033175017336) q[5];
u(pi/2,0.2620088273093888,-0.2620088273093888) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,2.7878493207955826,-2.7878493207955826) q[5];
u(pi/2,1.7404423300887455,-1.7404423300887455) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[5];
u(pi/2,0,0) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[5];
u(pi/2,-pi/4,pi/4) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,pi/4,-pi/4) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0,0) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0,0) q[5];
u(pi/2,3.7309554354032386,-3.7309554354032386) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-0.9814335449814514,0.9814335449814514) q[5];
u(pi/2,2.7495218904217866,-2.7495218904217866) q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,4.320318217216683,-4.320318217216683) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,4.516353598800687,-4.516353598800687) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-0.19603538158400302,0.19603538158400302) q[5];
u(pi/2,-pi/2,pi/2) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[5];
u(pi/2,pi/4,-pi/4) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
u(pi/2,0.16964600329384893,-0.16964600329384893) q[5];
u(pi/2,-0.8771326688822703,0.8771326688822703) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,-0.26200882730938857,0.26200882730938857) q[5];
u(pi/2,3.6649819896778526,-3.6649819896778526) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,2.0941856628829565,-2.0941856628829565) q[5];
u(pi/2,-0.26200882730938857,0.26200882730938857) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,1.3087874994855078,-1.3087874994855078) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0.5233893360880595,-0.5233893360880595) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,3.6649819896778526,-3.6649819896778526) q[5];
u(pi/2,1.1127521179015045,-1.1127521179015045) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-0.4580442088933918,0.4580442088933918) q[5];
u(pi/2,3.272911226509847,-3.272911226509847) q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,1.70211489971495,-1.70211489971495) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,1.8981502812989532,-1.8981502812989532) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,3.46894660809385,-3.46894660809385) q[5];
u(pi/2,2.0941856628829565,-2.0941856628829565) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0.5233893360880595,-0.5233893360880595) q[5];
u(pi/2,4.450380153075301,-4.450380153075301) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,3.8352563115024196,-3.8352563115024196) q[5];
u(pi/2,2.7878493207955826,-2.7878493207955826) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,3.403601480899182,-3.403601480899182) q[5];
u(pi/2,1.047406990706837,-1.047406990706837) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,2.6182033175017336,-2.6182033175017336) q[5];
u(pi/2,0.2620088273093888,-0.2620088273093888) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,1.832805154104285,-1.832805154104285) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,1.047406990706837,-1.047406990706837) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,1.047406990706837,-1.047406990706837) q[5];
u(pi/2,-1.5054511996002289,1.5054511996002289) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,3.206937780784461,-3.206937780784461) q[5];
u(pi/2,0.6547079090081129,-0.6547079090081129) q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,2.225504235803009,-2.225504235803009) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,2.4215396173870123,-2.4215396173870123) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,3.9923359441819093,-3.9923359441819093) q[5];
u(pi/2,2.6182033175017336,-2.6182033175017336) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,1.047406990706837,-1.047406990706837) q[5];
u(pi/2,-1.3087874994855078,1.3087874994855078) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
u(pi/2,1.2170529940006856,-1.2170529940006856) q[5];
u(pi/2,0.16964600329384893,-0.16964600329384893) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[5];
u(pi/2,pi/2,-pi/2) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0,0) q[5];
u(pi/2,5*pi/4,-5*pi/4) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,-pi/4,pi/4) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-pi/2,pi/2) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-pi/2,pi/2) q[5];
u(pi/2,2.1601591086083416,-2.1601591086083416) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0.589362781813445,-0.589362781813445) q[5];
u(pi/2,4.320318217216683,-4.320318217216683) q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,0.393327400229442,-0.393327400229442) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-1.1774689265654543,1.1774689265654543) q[5];
u(pi/2,2.7495218904217866,-2.7495218904217866) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0,0) q[0];
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
u(pi/2,11*pi/8,-11*pi/8) q[3];
rzz(-pi/2) q[0],q[3];
u(pi/2,7*pi/8,-7*pi/8) q[3];
u(pi/2,-0.19666370011472112,0.19666370011472112) q[3];
rzz(-pi/2) q[0],q[3];
u(pi/2,-pi/4,pi/4) q[1];
u(pi/2,-pi/2,pi/2) q[1];
rzz(pi/2) q[1],q[2];
u(pi/2,-3*pi/8,3*pi/8) q[2];
u(pi/2,7*pi/8,-7*pi/8) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u(pi/2,2.9455572720057903,-2.9455572720057903) q[3];
u(pi/2,0.19603538158400302,-0.19603538158400302) q[3];
rzz(-pi/2) q[1],q[3];
u(pi/2,3*pi/8,-3*pi/8) q[2];
u(pi/2,pi/4,-pi/4) q[2];
rzz(pi/2) q[2],q[3];
u(pi/2,3.337628035173796,-3.337628035173796) q[3];
u(pi/2,0.9814335449814515,-0.9814335449814515) q[3];
rzz(pi/2) q[2],q[3];
u(pi,pi,-pi) q[0];
u(pi,pi,-pi) q[1];
u(pi,pi,-pi) q[2];
u(pi/2,2.5522298717763476,-2.5522298717763476) q[3];
u(pi/2,pi/8,-pi/8) q[3];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,1.1787255636268905,-1.1787255636268905) q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
