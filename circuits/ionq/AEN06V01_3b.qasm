OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u(pi/2,-pi/8,pi/8) q[3];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[3],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[3],q[4];
u(pi/2,-pi/4,pi/4) q[5];
rzz(-pi/2) q[3],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[3],q[4];
u(pi/2,1.7404423300887455,-1.7404423300887455) q[5];
u(pi/2,-0.35374333279421055,0.35374333279421055) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-0.01382300767579503,0.01382300767579503) q[5];
u(pi/2,-1.061229998382632,1.061229998382632) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,-0.44547783827903253,0.44547783827903253) q[5];
u(pi/2,3.481512978708209,-3.481512978708209) q[5];
rzz(-pi/2) q[3],q[5];
u(pi/2,-1.230876001676481,1.230876001676481) q[5];
u(pi/2,2.696114815310761,-2.696114815310761) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[3],q[4];
u(pi/2,4.266911142105657,-4.266911142105657) q[5];
rzz(-pi/2) q[3],q[5];
u(pi/2,1.9107166519133125,-1.9107166519133125) q[5];
u(pi/2,1.1253184885158638,-1.1253184885158638) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,-0.44547783827903253,0.44547783827903253) q[5];
u(pi/2,3.481512978708209,-3.481512978708209) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[3],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[3],q[4];
u(pi/2,-0.44547783827903253,0.44547783827903253) q[5];
rzz(pi/2) q[3],q[5];
u(pi/2,1.1253184885158638,-1.1253184885158638) q[5];
u(pi/2,-1.230876001676481,1.230876001676481) q[5];
rzz(-pi/2) q[3],q[5];
u(pi/2,0,0) q[2];
u(pi/2,3.481512978708209,-3.481512978708209) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,1.1253184885158638,-1.1253184885158638) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,0.5095663284122645,-0.5095663284122645) q[5];
u(pi/2,4.698565972708895,-4.698565972708895) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,0,0) q[4];
u(pi/2,1.8968936442375166,-1.8968936442375166) q[5];
u(pi/2,0.84948665353068,-0.84948665353068) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,1.4652388136342798,-1.4652388136342798) q[5];
u(pi/2,-0.8909556765580653,0.8909556765580653) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,3.821433303826624,-3.821433303826624) q[5];
u(pi/2,1.4652388136342798,-1.4652388136342798) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,3.036035140429176,-3.036035140429176) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,0.6798406502368315,-0.6798406502368315) q[5];
u(pi/2,-0.10555751316061701,0.10555751316061701) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,4.606831467224072,-4.606831467224072) q[5];
u(pi/2,2.2506369770317276,-2.2506369770317276) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,-0.10555751316061701,0.10555751316061701) q[5];
rzz(pi/2) q[5],q[2];
u(pi/2,4.606831467224072,-4.606831467224072) q[5];
u(pi/2,2.2506369770317276,-2.2506369770317276) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,3.821433303826624,-3.821433303826624) q[5];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,3.036035140429176,-3.036035140429176) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[2],q[4];
u(pi/2,-0.7213096732642166,0.7213096732642166) q[5];
u(pi/2,3.4676899710324136,-3.4676899710324136) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,3.8076102961508296,-3.8076102961508296) q[5];
u(pi/2,2.7602033054439925,-2.7602033054439925) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,3.375955465547592,-3.375955465547592) q[5];
u(pi/2,1.019760975355247,-1.019760975355247) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,-0.5510353514396498,0.5510353514396498) q[5];
u(pi/2,3.375955465547592,-3.375955465547592) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[2],q[4];
u(pi/2,1.8051591387526953,-1.8051591387526953) q[5];
rzz(-pi/2) q[2],q[5];
u(pi/2,-0.5510353514396498,0.5510353514396498) q[5];
u(pi/2,-1.336433514837098,1.336433514837098) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,3.375955465547592,-3.375955465547592) q[5];
u(pi/2,1.019760975355247,-1.019760975355247) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[2],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[2],q[4];
u(pi/2,0.2343628119577985,-0.2343628119577985) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,1.8051591387526953,-1.8051591387526953) q[5];
u(pi/2,-0.5510353514396498,0.5510353514396498) q[5];
rzz(pi/2) q[2],q[5];
u(pi/2,pi/2,-pi/2) q[1];
u(pi/2,1.019760975355247,-1.019760975355247) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,1.8051591387526953,-1.8051591387526953) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,4.330999632238889,-4.330999632238889) q[5];
u(pi/2,2.2368139693559326,-2.2368139693559326) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,-0.5648583591154448,0.5648583591154448) q[5];
u(pi/2,4.670919957357304,-4.670919957357304) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,-0.9971415082494003,0.9971415082494003) q[5];
u(pi/2,2.9298493087378414,-2.9298493087378414) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,1.3590529819429444,-1.3590529819429444) q[5];
u(pi/2,-0.9971415082494003,0.9971415082494003) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,3.7152474721352897,-3.7152474721352897) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,1.3590529819429444,-1.3590529819429444) q[5];
u(pi/2,0.5736548185454962,-0.5736548185454962) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,-0.9971415082494003,0.9971415082494003) q[5];
u(pi/2,2.9298493087378414,-2.9298493087378414) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,0.5736548185454962,-0.5736548185454962) q[5];
rzz(pi/2) q[5],q[1];
u(pi/2,-0.9971415082494003,0.9971415082494003) q[5];
u(pi/2,2.9298493087378414,-2.9298493087378414) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,1.3590529819429444,-1.3590529819429444) q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,3.7152474721352897,-3.7152474721352897) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,-0.04146902302738531,0.04146902302738531) q[5];
u(pi/2,4.147530621269245,-4.147530621269245) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,0,0) q[4];
u(pi/2,1.3452299742671494,-1.3452299742671494) q[5];
u(pi/2,0.2984513020910302,-0.2984513020910302) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,4.055167797253705,-4.055167797253705) q[5];
u(pi/2,1.69897330706136,-1.69897330706136) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,0.1281769802664634,-0.1281769802664634) q[5];
u(pi/2,4.055167797253705,-4.055167797253705) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,2.4843714704588082,-2.4843714704588082) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,3.2697696338562565,-3.2697696338562565) q[5];
u(pi/2,2.4843714704588082,-2.4843714704588082) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,4.055167797253705,-4.055167797253705) q[5];
u(pi/2,1.69897330706136,-1.69897330706136) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,-0.6572211831309848,0.6572211831309848) q[5];
rzz(pi/2) q[5],q[1];
u(pi/2,4.055167797253705,-4.055167797253705) q[5];
u(pi/2,1.69897330706136,-1.69897330706136) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,3.2697696338562565,-3.2697696338562565) q[5];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-0.6572211831309848,0.6572211831309848) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,0,0) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,1.869247628885927,-1.869247628885927) q[5];
u(pi/2,-0.22556635252774715,0.22556635252774715) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,0.11435397259066837,-0.11435397259066837) q[5];
u(pi/2,-0.9324246995854506,0.9324246995854506) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,2.824291795577224,-2.824291795577224) q[5];
u(pi/2,0.46809730538487937,-0.46809730538487937) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,-1.1026990214100174,1.1026990214100174) q[5];
u(pi/2,2.824291795577224,-2.824291795577224) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,4.395088122372121,-4.395088122372121) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,-1.1026990214100174,1.1026990214100174) q[5];
u(pi/2,4.395088122372121,-4.395088122372121) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,-0.31730085801256913,0.31730085801256913) q[5];
u(pi/2,3.6096899589746725,-3.6096899589746725) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[1],q[4];
u(pi/2,4.395088122372121,-4.395088122372121) q[5];
rzz(-pi/2) q[5],q[1];
u(pi/2,-0.31730085801256913,0.31730085801256913) q[5];
u(pi/2,3.6096899589746725,-3.6096899589746725) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,2.038893632179776,-2.038893632179776) q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,4.395088122372121,-4.395088122372121) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,0.638371627209446,-0.638371627209446) q[5];
u(pi/2,-1.4564423542042282,1.4564423542042282) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,pi,-pi) q[4];
u(pi/2,-1.1165220290858124,1.1165220290858124) q[5];
u(pi/2,4.119884605917655,-4.119884605917655) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,-1.54817685968905,1.54817685968905) q[5];
u(pi/2,2.3788139572981915,-2.3788139572981915) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,3.9496102840930885,-3.9496102840930885) q[5];
u(pi/2,1.5934157939007432,-1.5934157939007432) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,0.02261946710584639,-0.02261946710584639) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,3.9496102840930885,-3.9496102840930885) q[5];
u(pi/2,3.16421212069564,-3.16421212069564) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,-1.54817685968905,1.54817685968905) q[5];
u(pi/2,2.3788139572981915,-2.3788139572981915) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[1],q[4];
u(pi/2,-1.54817685968905,1.54817685968905) q[5];
rzz(-pi/2) q[1],q[5];
u(pi/2,3.16421212069564,-3.16421212069564) q[5];
u(pi/2,0.8080176305032949,-0.8080176305032949) q[5];
rzz(pi/2) q[1],q[5];
u(pi/2,pi,-pi) q[0];
u(pi/2,2.3788139572981915,-2.3788139572981915) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,3.16421212069564,-3.16421212069564) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-0.592504374467035,0.592504374467035) q[5];
u(pi/2,3.595866951298877,-3.595866951298877) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
u(pi/2,0.7941946228274999,-0.7941946228274999) q[5];
u(pi/2,-0.25321236787933743,0.25321236787933743) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,3.5041324458140553,-3.5041324458140553) q[5];
u(pi/2,1.1479379556217104,-1.1479379556217104) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,2.718734282416607,-2.718734282416607) q[5];
u(pi/2,0.3625397922242619,-0.3625397922242619) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-1.2082565345706344,1.2082565345706344) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,2.718734282416607,-2.718734282416607) q[5];
u(pi/2,1.9333361190191587,-1.9333361190191587) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0.3625397922242619,-0.3625397922242619) q[5];
u(pi/2,4.2895306092115035,-4.2895306092115035) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,1.9333361190191587,-1.9333361190191587) q[5];
rzz(-pi/2) q[5],q[0];
u(pi/2,3.5041324458140553,-3.5041324458140553) q[5];
u(pi/2,1.1479379556217104,-1.1479379556217104) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-0.42285837117318614,0.42285837117318614) q[5];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-1.2082565345706344,1.2082565345706344) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,4.4591766125053525,-4.4591766125053525) q[5];
u(pi/2,2.364990949622396,-2.364990949622396) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,2.7049112747408115,-2.7049112747408115) q[5];
u(pi/2,1.6575042840339749,-1.6575042840339749) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,2.273256444137574,-2.273256444137574) q[5];
u(pi/2,-0.0829380460547704,0.0829380460547704) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,4.6294509343299195,-4.6294509343299195) q[5];
u(pi/2,2.273256444137574,-2.273256444137574) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3.844052770932471,-3.844052770932471) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,1.487858280740126,-1.487858280740126) q[5];
u(pi/2,0.7024601173426777,-0.7024601173426777) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-0.8683362094522188,0.8683362094522188) q[5];
u(pi/2,3.058654607535023,-3.058654607535023) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3.844052770932471,-3.844052770932471) q[5];
rzz(pi/2) q[5],q[0];
u(pi/2,2.273256444137574,-2.273256444137574) q[5];
u(pi/2,-0.0829380460547704,0.0829380460547704) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,1.487858280740126,-1.487858280740126) q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,3.844052770932471,-3.844052770932471) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,0.08670795723907809,-0.08670795723907809) q[5];
u(pi/2,4.275707601535708,-4.275707601535708) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0,0) q[4];
u(pi/2,4.615627926654124,-4.615627926654124) q[5];
u(pi/2,3.568220935947287,-3.568220935947287) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,1.0423804424610932,-1.0423804424610932) q[5];
u(pi/2,-1.3138140477312514,1.3138140477312514) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,3.398574932653439,-3.398574932653439) q[5];
u(pi/2,1.0423804424610932,-1.0423804424610932) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,-0.5284158843338032,0.5284158843338032) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,3.398574932653439,-3.398574932653439) q[5];
u(pi/2,2.6131767692559906,-2.6131767692559906) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,1.0423804424610932,-1.0423804424610932) q[5];
u(pi/2,-1.3138140477312514,1.3138140477312514) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,pi/4,-pi/4) q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,2.6131767692559906,-2.6131767692559906) q[5];
rzz(pi/2) q[5],q[0];
u(pi/2,1.0423804424610932,-1.0423804424610932) q[5];
u(pi/2,-1.3138140477312514,1.3138140477312514) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,3.398574932653439,-3.398574932653439) q[5];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-0.5284158843338032,0.5284158843338032) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0,0) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,1.9974246091523904,-1.9974246091523904) q[5];
u(pi/2,-0.09676105373056565,0.09676105373056565) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,3.3847519249776425,-3.3847519249776425) q[5];
u(pi/2,2.337344934270806,-2.337344934270806) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,2.9530970943744057,-2.9530970943744057) q[5];
u(pi/2,0.5969026041820604,-0.5969026041820604) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,2.167698930976957,-2.167698930976957) q[5];
u(pi/2,-0.18849555921538763,0.18849555921538763) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,4.523893421169302,-4.523893421169302) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-0.9738937226128359,0.9738937226128359) q[5];
u(pi/2,4.523893421169302,-4.523893421169302) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,2.9530970943744057,-2.9530970943744057) q[5];
u(pi/2,0.5969026041820604,-0.5969026041820604) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,1.3823007675795087,-1.3823007675795087) q[5];
rzz(pi/2) q[5],q[0];
u(pi/2,-0.18849555921538763,0.18849555921538763) q[5];
u(pi/2,3.738495257771854,-3.738495257771854) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,2.167698930976957,-2.167698930976957) q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,1.3823007675795087,-1.3823007675795087) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,0.7665486074759094,-0.7665486074759094) q[5];
u(pi/2,-1.3276370554070467,1.3276370554070467) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0,0) q[4];
u(pi/2,-0.987716730288631,0.987716730288631) q[5];
u(pi/2,4.248061586184119,-4.248061586184119) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,1.7215927741672066,-1.7215927741672066) q[5];
u(pi/2,-0.6346017160251383,0.6346017160251383) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,4.077787264359552,-4.077787264359552) q[5];
u(pi/2,1.7215927741672066,-1.7215927741672066) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3.2923891009621036,-3.2923891009621036) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0.9361946107697583,-0.9361946107697583) q[5];
u(pi/2,0.15079644737231024,-0.15079644737231024) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,1.7215927741672066,-1.7215927741672066) q[5];
u(pi/2,-0.6346017160251383,0.6346017160251383) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,0.15079644737231024,-0.15079644737231024) q[5];
rzz(pi/2) q[5],q[0];
u(pi/2,-1.4199998794225865,1.4199998794225865) q[5];
u(pi/2,2.5069909375646553,-2.5069909375646553) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0.9361946107697583,-0.9361946107697583) q[5];
u(pi/2,0,0) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0.15079644737231024,-0.15079644737231024) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,2.6772652593892223,-2.6772652593892223) q[5];
u(pi/2,0.5830795965062654,-0.5830795965062654) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,0.9223716030939633,-0.9223716030939633) q[5];
u(pi/2,-0.12440706908215593,0.12440706908215593) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,0.49071677249072554,-0.49071677249072554) q[5];
u(pi/2,4.417707589477967,-4.417707589477967) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-0.2946813909067225,0.2946813909067225) q[5];
u(pi/2,3.6323094260805187,-3.6323094260805187) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,2.061513099285622,-2.061513099285622) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-0.2946813909067225,0.2946813909067225) q[5];
u(pi/2,-1.0800795543041708,1.0800795543041708) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,3.6323094260805187,-3.6323094260805187) q[5];
u(pi/2,1.2761149358881738,-1.2761149358881738) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-1.0800795543041708,1.0800795543041708) q[5];
rzz(-pi/2) q[5],q[0];
u(pi/2,0.49071677249072554,-0.49071677249072554) q[5];
u(pi/2,4.417707589477967,-4.417707589477967) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,2.8469112626830704,-2.8469112626830704) q[5];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,2.061513099285622,-2.061513099285622) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-pi/2,pi/2) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,1.4463892577127409,-1.4463892577127409) q[5];
u(pi/2,-0.6484247237009332,0.6484247237009332) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0,0) q[4];
u(pi/2,-0.30850439858251777,0.30850439858251777) q[5];
u(pi/2,-1.3552830707586367,1.3552830707586367) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,-0.7401592291857552,0.7401592291857552) q[5];
u(pi/2,3.1868315878014863,-3.1868315878014863) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,1.6160352610065893,-1.6160352610065893) q[5];
u(pi/2,-0.7401592291857552,0.7401592291857552) q[5];
rzz(-pi/2) q[4],q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,0.8306370976091411,-0.8306370976091411) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-1.5255573925832036,1.5255573925832036) q[5];
u(pi/2,3.9722297511989346,-3.9722297511989346) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,2.4014334244040376,-2.4014334244040376) q[5];
u(pi/2,0.045238934211693005,-0.045238934211693005) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3.9722297511989346,-3.9722297511989346) q[5];
rzz(pi/2) q[5],q[0];
u(pi/2,2.4014334244040376,-2.4014334244040376) q[5];
u(pi/2,0.045238934211693005,-0.045238934211693005) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,1.6160352610065893,-1.6160352610065893) q[5];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,0,0) q[4];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
u(pi/2,0,0) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,3.9722297511989346,-3.9722297511989346) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi,-pi) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,0.21551325603625981,-0.21551325603625981) q[5];
u(pi/2,4.4038845818021715,-4.4038845818021715) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,-1.5393804002589986,1.5393804002589986) q[5];
u(pi/2,3.6963979162137512,-3.6963979162137512) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,4.312150076317351,-4.312150076317351) q[5];
u(pi/2,1.9559555861250053,-1.9559555861250053) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0.38515925933010875,-0.38515925933010875) q[5];
u(pi/2,4.312150076317351,-4.312150076317351) q[5];
rzz(pi/2) q[4],q[5];
u(pi/2,5*pi/4,-5*pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,2.741353749522454,-2.741353749522454) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,0.38515925933010875,-0.38515925933010875) q[5];
u(pi/2,-0.40023890406733975,0.40023890406733975) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,1.170557422727557,-1.170557422727557) q[5];
u(pi/2,-1.1856370674647878,1.1856370674647878) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-pi/4,pi/4) q[4];
rzz(-pi/2) q[5],q[4];
u(pi/2,pi/2,-pi/2) q[4];
u(pi/2,pi/4,-pi/4) q[4];
rzz(-pi/2) q[0],q[4];
u(pi/2,-pi/4,pi/4) q[4];
u(pi/2,pi,-pi) q[4];
rzz(pi/2) q[5],q[4];
u(pi/2,-pi/2,pi/2) q[4];
u(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,4.312150076317351,-4.312150076317351) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,-0.40023890406733975,0.40023890406733975) q[5];
u(pi/2,3.5267519129199023,-3.5267519129199023) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,0,0) q[0];
u(pi/2,-pi/2,pi/2) q[1];
rzz(-pi/2) q[0],q[1];
u(pi/2,pi,-pi) q[1];
u(pi/2,pi/4,-pi/4) q[1];
rzz(-pi/2) q[0],q[1];
u(pi/2,0,0) q[2];
rzz(pi/2) q[0],q[2];
u(pi/2,pi/2,-pi/2) q[2];
u(pi/2,-3*pi/8,3*pi/8) q[2];
rzz(pi/2) q[0],q[2];
u(pi/2,-pi/8,pi/8) q[3];
rzz(pi/2) q[0],q[3];
u(pi/2,3*pi/8,-3*pi/8) q[3];
u(pi/2,4.516353598800687,-4.516353598800687) q[3];
rzz(-pi/2) q[0],q[3];
u(pi/2,-pi/4,pi/4) q[1];
u(pi/2,-pi/2,pi/2) q[1];
rzz(pi/2) q[1],q[2];
u(pi/2,-3*pi/8,3*pi/8) q[2];
u(pi/2,7*pi/8,-7*pi/8) q[2];
rzz(-pi/2) q[1],q[2];
rzz(-pi/2) q[1],q[3];
u(pi/2,4.516353598800687,-4.516353598800687) q[3];
u(pi/2,1.7668317083788998,-1.7668317083788998) q[3];
rzz(-pi/2) q[1],q[3];
u(pi/2,11*pi/8,-11*pi/8) q[2];
u(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(-pi/2) q[2],q[3];
u(pi/2,1.7668317083788998,-1.7668317083788998) q[3];
u(pi/2,-0.5893627818134451,0.5893627818134451) q[3];
rzz(pi/2) q[2],q[3];
u(pi,pi,-pi) q[0];
u(pi,pi,-pi) q[1];
u(pi,pi,-pi) q[2];
u(pi/2,0.9814335449814515,-0.9814335449814515) q[3];
u(pi/2,-3*pi/8,3*pi/8) q[3];
u(pi/2,5*pi/4,-5*pi/4) q[4];
u(pi/2,pi,-pi) q[4];
u(pi/2,-1.1856370674647878,1.1856370674647878) q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];