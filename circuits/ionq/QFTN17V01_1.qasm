OPENQASM 2.0;
include "qelib1.inc";
qreg q[17];
creg c[17];
u3(pi/2,3.9357872764172925,-3.9357872764172925) q[16];
u3(pi/2,4.721185439814741,-4.721185439814741) q[16];
rzz(pi/2) q[15],q[16];
rzz(-pi/2) q[14],q[16];
u3(pi/2,4.721185439814741,-4.721185439814741) q[16];
u3(pi/2,1.186893704526224,-1.186893704526224) q[16];
rzz(-pi/2) q[14],q[16];
rzz(pi/2) q[14],q[15];
u3(pi/2,0.4209734155810323,-0.4209734155810323) q[15];
u3(pi/2,1.2063715789784806,-1.2063715789784806) q[15];
rzz(-pi/2) q[14],q[15];
rzz(pi/2) q[13],q[16];
u3(pi/2,4.328486358116017,-4.328486358116017) q[16];
u3(pi/2,0.9902300044115028,-0.9902300044115028) q[16];
rzz(-pi/2) q[13],q[16];
rzz(-pi/2) q[13],q[15];
u3(pi/2,4.347964232568273,-4.347964232568273) q[15];
u3(pi/2,0.8136724972797564,-0.8136724972797564) q[15];
rzz(-pi/2) q[13],q[15];
rzz(-pi/2) q[13],q[14];
u3(pi/2,1.6047255274536665,-1.6047255274536665) q[14];
u3(pi/2,2.3901236908511145,-2.3901236908511145) q[14];
rzz(-pi/2) q[13],q[14];
rzz(-pi/2) q[12],q[16];
u3(pi/2,4.131822658001296,-4.131822658001296) q[16];
u3(pi/2,0.8922123136195012,-0.8922123136195012) q[16];
rzz(-pi/2) q[12],q[16];
rzz(-pi/2) q[12],q[15];
u3(pi/2,3.955265150869549,-3.955265150869549) q[15];
u3(pi/2,0.6170087971650353,-0.6170087971650353) q[15];
rzz(pi/2) q[12],q[15];
rzz(-pi/2) q[12],q[14];
u3(pi/2,2.3901236908511145,-2.3901236908511145) q[14];
u3(pi/2,5.1390172627421835,-5.1390172627421835) q[14];
rzz(-pi/2) q[12],q[14];
rzz(pi/2) q[12],q[13];
u3(pi/2,0.1124690169985146,-0.1124690169985146) q[13];
u3(pi/2,0.8978671803959629,-0.8978671803959629) q[13];
rzz(pi/2) q[12],q[13];
rzz(-pi/2) q[11],q[16];
u3(pi/2,0.8922123136195012,-0.8922123136195012) q[16];
u3(pi/2,3.9847961218132935,-3.9847961218132935) q[16];
rzz(-pi/2) q[11],q[16];
rzz(-pi/2) q[11],q[15];
u3(pi/2,3.758601450754828,-3.758601450754828) q[15];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[15];
rzz(-pi/2) q[11],q[15];
rzz(pi/2) q[11],q[14];
u3(pi/2,1.9974246091523906,-1.9974246091523906) q[14];
u3(pi/2,4.94298188115818,-4.94298188115818) q[14];
rzz(-pi/2) q[11],q[14];
rzz(-pi/2) q[11],q[13];
u3(pi/2,0.8978671803959629,-0.8978671803959629) q[13];
u3(pi/2,3.6467607522870322,-3.6467607522870322) q[13];
rzz(-pi/2) q[11],q[13];
rzz(-pi/2) q[11],q[12];
u3(pi/2,1.7071414479606937,-1.7071414479606937) q[12];
u3(pi/2,2.4925396113581417,-2.4925396113581417) q[12];
rzz(-pi/2) q[11],q[12];
rzz(-pi/2) q[10],q[16];
u3(pi/2,0.8432034682235006,-0.8432034682235006) q[16];
u3(pi/2,3.960291699115293,-3.960291699115293) q[16];
rzz(-pi/2) q[10],q[16];
rzz(pi/2) q[10],q[15];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[15];
u3(pi/2,3.611574914566826,-3.611574914566826) q[15];
rzz(-pi/2) q[10],q[15];
rzz(-pi/2) q[10],q[14];
u3(pi/2,1.8013892275683874,-1.8013892275683874) q[14];
u3(pi/2,4.844964190366179,-4.844964190366179) q[14];
rzz(-pi/2) q[10],q[14];
rzz(pi/2) q[10],q[13];
u3(pi/2,3.6467607522870322,-3.6467607522870322) q[13];
u3(pi/2,0.30913271711323564,-0.30913271711323564) q[13];
rzz(-pi/2) q[10],q[13];
rzz(-pi/2) q[10],q[12];
u3(pi/2,2.4925396113581417,-2.4925396113581417) q[12];
u3(pi/2,5.241433183249211,-5.241433183249211) q[12];
rzz(-pi/2) q[10],q[12];
rzz(pi/2) q[10],q[11];
u3(pi/2,2.021929031850391,-2.021929031850391) q[11];
u3(pi/2,2.807327195247839,-2.807327195247839) q[11];
rzz(-pi/2) q[10],q[11];
rzz(-pi/2) q[9],q[16];
u3(pi/2,3.960291699115293,-3.960291699115293) q[16];
u3(pi/2,0.8061326749111409,-0.8061326749111409) q[16];
rzz(-pi/2) q[9],q[16];
rzz(-pi/2) q[9],q[15];
u3(pi/2,3.611574914566826,-3.611574914566826) q[15];
u3(pi/2,0.4454778382790327,-0.4454778382790327) q[15];
rzz(pi/2) q[9],q[15];
rzz(-pi/2) q[9],q[14];
u3(pi/2,4.844964190366179,-4.844964190366179) q[14];
u3(pi/2,1.653734372849667,-1.653734372849667) q[14];
rzz(-pi/2) q[9],q[14];
rzz(-pi/2) q[9],q[13];
u3(pi/2,0.30913271711323564,-0.30913271711323564) q[13];
u3(pi/2,3.352707679911027,-3.352707679911027) q[13];
rzz(pi/2) q[9],q[13];
rzz(-pi/2) q[9],q[12];
u3(pi/2,5.241433183249211,-5.241433183249211) q[12];
u3(pi/2,1.9038051480754146,-1.9038051480754146) q[12];
rzz(-pi/2) q[9],q[12];
rzz(-pi/2) q[9],q[11];
u3(pi/2,5.948919848837632,-5.948919848837632) q[11];
u3(pi/2,2.414628113549115,-2.414628113549115) q[11];
rzz(-pi/2) q[9],q[11];
rzz(-pi/2) q[9],q[10];
u3(pi/2,5.258397783578595,-5.258397783578595) q[10];
u3(pi/2,6.043795946976044,-6.043795946976044) q[10];
rzz(-pi/2) q[9],q[10];
rzz(-pi/2) q[8],q[16];
rzz(-pi/2) q[8],q[16];
rzz(pi/2) q[8],q[15];
u3(pi/2,0.4454778382790327,-0.4454778382790327) q[15];
u3(pi/2,3.574504121254466,-3.574504121254466) q[15];
rzz(-pi/2) q[8],q[15];
rzz(-pi/2) q[8],q[14];
u3(pi/2,1.653734372849667,-1.653734372849667) q[14];
u3(pi/2,4.770822603741459,-4.770822603741459) q[14];
rzz(-pi/2) q[8],q[14];
rzz(pi/2) q[8],q[13];
u3(pi/2,3.352707679911027,-3.352707679911027) q[13];
u3(pi/2,0.16210618092523332,-0.16210618092523332) q[13];
rzz(-pi/2) q[8],q[13];
rzz(-pi/2) q[8],q[12];
u3(pi/2,1.9038051480754146,-1.9038051480754146) q[12];
u3(pi/2,4.9473801108732065,-4.9473801108732065) q[12];
rzz(-pi/2) q[8],q[12];
rzz(pi/2) q[8],q[11];
u3(pi/2,5.556220767138908,-5.556220767138908) q[11];
u3(pi/2,2.2179644134343937,-2.2179644134343937) q[11];
rzz(-pi/2) q[8],q[11];
rzz(-pi/2) q[8],q[10];
u3(pi/2,2.902203293386251,-2.902203293386251) q[10];
u3(pi/2,5.65109686527732,-5.65109686527732) q[10];
rzz(-pi/2) q[8],q[10];
rzz(pi/2) q[8],q[9];
u3(pi/2,3.730955435403238,-3.730955435403238) q[9];
u3(pi/2,4.516353598800687,-4.516353598800687) q[9];
rzz(pi/2) q[8],q[9];
rzz(-pi/2) q[7],q[16];
rzz(-pi/2) q[7],q[16];
rzz(-pi/2) q[7],q[15];
rzz(pi/2) q[7],q[15];
rzz(-pi/2) q[7],q[14];
u3(pi/2,4.770822603741459,-4.770822603741459) q[14];
u3(pi/2,1.6172918980680255,-1.6172918980680255) q[14];
rzz(-pi/2) q[7],q[14];
rzz(-pi/2) q[7],q[13];
u3(pi/2,0.16210618092523332,-0.16210618092523332) q[13];
u3(pi/2,3.279194411817026,-3.279194411817026) q[13];
rzz(-pi/2) q[7],q[13];
rzz(pi/2) q[7],q[12];
u3(pi/2,1.805787457283413,-1.805787457283413) q[12];
u3(pi/2,4.8977429469464875,-4.8977429469464875) q[12];
rzz(-pi/2) q[7],q[12];
rzz(-pi/2) q[7],q[11];
u3(pi/2,2.2179644134343937,-2.2179644134343937) q[11];
u3(pi/2,5.261539376232186,-5.261539376232186) q[11];
rzz(-pi/2) q[7],q[11];
rzz(pi/2) q[7],q[10];
u3(pi/2,2.509504211687527,-2.509504211687527) q[10];
u3(pi/2,5.455061483693316,-5.455061483693316) q[10];
rzz(-pi/2) q[7],q[10];
rzz(-pi/2) q[7],q[9];
u3(pi/2,4.516353598800687,-4.516353598800687) q[9];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[9];
rzz(-pi/2) q[7],q[9];
rzz(pi/2) q[7],q[8];
u3(pi/2,3.7554598581012386,-3.7554598581012386) q[8];
u3(pi/2,4.540858021498687,-4.540858021498687) q[8];
rzz(-pi/2) q[7],q[8];
rzz(-pi/2) q[6],q[16];
rzz(-pi/2) q[6],q[16];
rzz(pi/2) q[6],q[15];
rzz(-pi/2) q[6],q[15];
rzz(-pi/2) q[6],q[14];
rzz(-pi/2) q[6],q[14];
rzz(-pi/2) q[6],q[13];
u3(pi/2,0.13760175822723295,-0.13760175822723295) q[13];
u3(pi/2,3.266628041202667,-3.266628041202667) q[13];
rzz(-pi/2) q[6],q[13];
rzz(-pi/2) q[6],q[12];
u3(pi/2,1.7561502933566946,-1.7561502933566946) q[12];
u3(pi/2,4.873238524248487,-4.873238524248487) q[12];
rzz(-pi/2) q[6],q[12];
rzz(-pi/2) q[6],q[11];
u3(pi/2,2.1199467226423923,-2.1199467226423923) q[11];
u3(pi/2,5.212530530836185,-5.212530530836185) q[11];
rzz(pi/2) q[6],q[11];
rzz(-pi/2) q[6],q[10];
u3(pi/2,2.313468830103524,-2.313468830103524) q[10];
u3(pi/2,5.356415474370597,-5.356415474370597) q[10];
rzz(-pi/2) q[6],q[10];
rzz(-pi/2) q[6],q[9];
u3(pi/2,4.123026198571244,-4.123026198571244) q[9];
u3(pi/2,pi/4,-pi/4) q[9];
rzz(pi/2) q[6],q[9];
rzz(-pi/2) q[6],q[8];
u3(pi/2,4.540858021498687,-4.540858021498687) q[8];
u3(pi/2,1.0065662862101699,-1.0065662862101699) q[8];
rzz(-pi/2) q[6],q[8];
rzz(pi/2) q[6],q[7];
u3(pi/2,2.503221026380347,-2.503221026380347) q[7];
u3(pi/2,3.2886191897777954,-3.2886191897777954) q[7];
rzz(-pi/2) q[6],q[7];
rzz(pi/2) q[5],q[16];
rzz(-pi/2) q[5],q[16];
rzz(-pi/2) q[5],q[15];
rzz(-pi/2) q[5],q[15];
rzz(pi/2) q[5],q[14];
rzz(-pi/2) q[5],q[14];
rzz(-pi/2) q[5],q[13];
rzz(-pi/2) q[5],q[13];
rzz(-pi/2) q[5],q[12];
u3(pi/2,1.7316458706586941,-1.7316458706586941) q[12];
u3(pi/2,4.861300472164846,-4.861300472164846) q[12];
rzz(pi/2) q[5],q[12];
rzz(-pi/2) q[5],q[11];
u3(pi/2,5.212530530836185,-5.212530530836185) q[11];
u3(pi/2,2.0464334545483913,-2.0464334545483913) q[11];
rzz(-pi/2) q[5],q[11];
rzz(pi/2) q[5],q[10];
u3(pi/2,5.356415474370597,-5.356415474370597) q[10];
u3(pi/2,2.1658139753848036,-2.1658139753848036) q[10];
rzz(-pi/2) q[5],q[10];
rzz(-pi/2) q[5],q[9];
u3(pi/2,pi/4,-pi/4) q[9];
u3(pi/2,3.82897312619524,-3.82897312619524) q[9];
rzz(-pi/2) q[5],q[9];
rzz(-pi/2) q[5],q[8];
u3(pi/2,4.148158939799963,-4.148158939799963) q[8];
u3(pi/2,0.8099025860954486,-0.8099025860954486) q[8];
rzz(pi/2) q[5],q[8];
rzz(-pi/2) q[5],q[7];
u3(pi/2,3.2886191897777954,-3.2886191897777954) q[7];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
rzz(-pi/2) q[5],q[7];
rzz(pi/2) q[5],q[6];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[6];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[6];
rzz(pi/2) q[5],q[6];
rzz(-pi/2) q[4],q[16];
rzz(-pi/2) q[4],q[16];
rzz(-pi/2) q[4],q[15];
rzz(-pi/2) q[4],q[15];
rzz(pi/2) q[4],q[14];
rzz(-pi/2) q[4],q[14];
rzz(-pi/2) q[4],q[13];
rzz(-pi/2) q[4],q[13];
rzz(pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[11];
u3(pi/2,2.0464334545483913,-2.0464334545483913) q[11];
u3(pi/2,5.175459737523825,-5.175459737523825) q[11];
rzz(-pi/2) q[4],q[11];
rzz(-pi/2) q[4],q[10];
u3(pi/2,2.1658139753848036,-2.1658139753848036) q[10];
u3(pi/2,5.282902206276596,-5.282902206276596) q[10];
rzz(-pi/2) q[4],q[10];
rzz(-pi/2) q[4],q[9];
u3(pi/2,3.82897312619524,-3.82897312619524) q[9];
u3(pi/2,0.638371627209446,-0.638371627209446) q[9];
rzz(-pi/2) q[4],q[9];
rzz(-pi/2) q[4],q[8];
u3(pi/2,3.951495239685242,-3.951495239685242) q[8];
u3(pi/2,0.7118848953034471,-0.7118848953034471) q[8];
rzz(pi/2) q[4],q[8];
rzz(-pi/2) q[4],q[7];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
u3(pi/2,2.6998847264950685,-2.6998847264950685) q[7];
rzz(-pi/2) q[4],q[7];
rzz(-pi/2) q[4],q[6];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[6];
u3(pi/2,4.417707589477967,-4.417707589477967) q[6];
rzz(-pi/2) q[4],q[6];
rzz(pi/2) q[4],q[5];
u3(pi/2,5.301751762198135,-5.301751762198135) q[5];
u3(pi/2,6.087149925595583,-6.087149925595583) q[5];
rzz(-pi/2) q[4],q[5];
rzz(-pi/2) q[3],q[16];
rzz(pi/2) q[3],q[16];
rzz(-pi/2) q[3],q[15];
rzz(-pi/2) q[3],q[15];
rzz(-pi/2) q[3],q[14];
rzz(pi/2) q[3],q[14];
rzz(-pi/2) q[3],q[13];
rzz(-pi/2) q[3],q[13];
rzz(-pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[12];
rzz(pi/2) q[3],q[11];
rzz(-pi/2) q[3],q[11];
rzz(-pi/2) q[3],q[10];
u3(pi/2,2.141309552686803,-2.141309552686803) q[10];
u3(pi/2,5.270964154192955,-5.270964154192955) q[10];
rzz(-pi/2) q[3],q[10];
rzz(pi/2) q[3],q[9];
u3(pi/2,0.638371627209446,-0.638371627209446) q[9];
u3(pi/2,3.7554598581012386,-3.7554598581012386) q[9];
rzz(-pi/2) q[3],q[9];
rzz(-pi/2) q[3],q[8];
u3(pi/2,0.7118848953034471,-0.7118848953034471) q[8];
u3(pi/2,3.8044687034972395,-3.8044687034972395) q[8];
rzz(-pi/2) q[3],q[8];
rzz(-pi/2) q[3],q[7];
u3(pi/2,5.841477380084861,-5.841477380084861) q[7];
u3(pi/2,2.6018670357030667,-2.6018670357030667) q[7];
rzz(pi/2) q[3],q[7];
rzz(-pi/2) q[3],q[6];
u3(pi/2,1.276114935888174,-1.276114935888174) q[6];
u3(pi/2,4.221672207893964,-4.221672207893964) q[6];
rzz(-pi/2) q[3],q[6];
rzz(-pi/2) q[3],q[5];
u3(pi/2,6.087149925595583,-6.087149925595583) q[5];
u3(pi/2,2.552229871776348,-2.552229871776348) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[4];
u3(pi/2,13*pi/8,-13*pi/8) q[4];
u3(pi/2,15*pi/8,-15*pi/8) q[4];
rzz(-pi/2) q[3],q[4];
rzz(-pi/2) q[2],q[16];
rzz(pi/2) q[2],q[16];
rzz(-pi/2) q[2],q[15];
rzz(-pi/2) q[2],q[15];
rzz(-pi/2) q[2],q[14];
rzz(pi/2) q[2],q[14];
rzz(-pi/2) q[2],q[13];
rzz(-pi/2) q[2],q[13];
rzz(-pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[12];
rzz(pi/2) q[2],q[11];
rzz(-pi/2) q[2],q[11];
rzz(-pi/2) q[2],q[10];
rzz(-pi/2) q[2],q[10];
rzz(pi/2) q[2],q[9];
u3(pi/2,0.6138672045114455,-0.6138672045114455) q[9];
u3(pi/2,3.7428934874868798,-3.7428934874868798) q[9];
rzz(-pi/2) q[2],q[9];
rzz(-pi/2) q[2],q[8];
u3(pi/2,3.8044687034972395,-3.8044687034972395) q[8];
u3(pi/2,0.638371627209446,-0.638371627209446) q[8];
rzz(-pi/2) q[2],q[8];
rzz(-pi/2) q[2],q[7];
u3(pi/2,5.74345968929286,-5.74345968929286) q[7];
u3(pi/2,2.5528581903070657,-2.5528581903070657) q[7];
rzz(-pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[6];
u3(pi/2,4.221672207893964,-4.221672207893964) q[6];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[6];
rzz(-pi/2) q[2],q[6];
rzz(-pi/2) q[2],q[5];
u3(pi/2,5.693822525366141,-5.693822525366141) q[5];
u3(pi/2,3*pi/4,-3*pi/4) q[5];
rzz(pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[4];
u3(pi/2,7*pi/8,-7*pi/8) q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
rzz(-pi/2) q[2],q[4];
rzz(pi/2) q[2],q[3];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[2],q[3];
rzz(-pi/2) q[1],q[16];
rzz(-pi/2) q[1],q[16];
rzz(-pi/2) q[1],q[15];
rzz(pi/2) q[1],q[15];
rzz(-pi/2) q[1],q[14];
rzz(-pi/2) q[1],q[14];
rzz(-pi/2) q[1],q[13];
rzz(pi/2) q[1],q[13];
rzz(-pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[11];
rzz(pi/2) q[1],q[10];
rzz(-pi/2) q[1],q[10];
rzz(-pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[9];
rzz(pi/2) q[1],q[8];
u3(pi/2,3.779964280799239,-3.779964280799239) q[8];
u3(pi/2,0.6258052565950868,-0.6258052565950868) q[8];
rzz(-pi/2) q[1],q[8];
rzz(-pi/2) q[1],q[7];
u3(pi/2,2.5528581903070657,-2.5528581903070657) q[7];
u3(pi/2,5.66931810266814,-5.66931810266814) q[7];
rzz(-pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[6];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[6];
u3(pi/2,4.074017353175243,-4.074017353175243) q[6];
rzz(pi/2) q[1],q[6];
rzz(-pi/2) q[1],q[5];
u3(pi/2,7*pi/4,-7*pi/4) q[5];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[5];
rzz(-pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[4];
rzz(pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[3];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
rzz(-pi/2) q[1],q[3];
rzz(pi/2) q[1],q[2];
u3(pi/2,0,0) q[2];
u3(pi/2,pi/4,-pi/4) q[2];
rzz(-pi/2) q[1],q[2];
rzz(pi/2) q[0],q[16];
rzz(-pi/2) q[0],q[15];
rzz(-pi/2) q[0],q[14];
rzz(pi/2) q[0],q[13];
rzz(-pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[9];
rzz(pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[7];
rzz(-pi/2) q[0],q[6];
rzz(-pi/2) q[0],q[5];
rzz(pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[3];
rzz(-pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi/2,0,0) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,5.301751762198135,-5.301751762198135) q[4];
u3(pi/2,5*pi/8,-5*pi/8) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,5.399769452990137,-5.399769452990137) q[5];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,4.074017353175243,-4.074017353175243) q[6];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,2.5277254490783476,-2.5277254490783476) q[7];
u3(pi/2,5.64481367997014,-5.64481367997014) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,0.6258052565950868,-0.6258052565950868) q[8];
u3(pi/2,3.7554598581012386,-3.7554598581012386) q[8];
rzz(pi/2) q[0],q[8];
rzz(pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[10];
rzz(pi/2) q[0],q[11];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[14];
rzz(-pi/2) q[0],q[15];
rzz(pi/2) q[0],q[16];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,0,0) q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(-pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,9*pi/8,-9*pi/8) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,13*pi/8,-13*pi/8) q[4];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[4];
rzz(pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[5];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[5];
u3(pi/2,5.203105752875415,-5.203105752875415) q[5];
rzz(pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[6];
u3(pi/2,3.9759996623832423,-3.9759996623832423) q[6];
rzz(pi/2) q[1],q[6];
rzz(-pi/2) q[1],q[7];
u3(pi/2,5.64481367997014,-5.64481367997014) q[7];
u3(pi/2,2.4787166036823467,-2.4787166036823467) q[7];
rzz(pi/2) q[1],q[7];
rzz(pi/2) q[1],q[8];
u3(pi/2,0.6138672045114455,-0.6138672045114455) q[8];
u3(pi/2,3.7428934874868798,-3.7428934874868798) q[8];
rzz(pi/2) q[1],q[8];
rzz(pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[9];
rzz(pi/2) q[1],q[10];
rzz(pi/2) q[1],q[10];
rzz(pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[11];
rzz(pi/2) q[1],q[12];
rzz(pi/2) q[1],q[12];
rzz(pi/2) q[1],q[13];
rzz(pi/2) q[1],q[13];
rzz(-pi/2) q[1],q[14];
rzz(pi/2) q[1],q[14];
rzz(pi/2) q[1],q[15];
rzz(pi/2) q[1],q[15];
rzz(-pi/2) q[1],q[16];
rzz(pi/2) q[1],q[16];
u3(pi/2,pi/4,-pi/4) q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,9*pi/8,-9*pi/8) q[3];
u3(pi/2,15*pi/8,-15*pi/8) q[3];
rzz(pi/2) q[2],q[3];
rzz(-pi/2) q[2],q[4];
u3(pi/2,4.909052680499411,-4.909052680499411) q[4];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,5.203105752875415,-5.203105752875415) q[5];
u3(pi/2,1.865477717701619,-1.865477717701619) q[5];
rzz(pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[6];
u3(pi/2,0.8344070087934491,-0.8344070087934491) q[6];
u3(pi/2,3.8779819715912405,-3.8779819715912405) q[6];
rzz(pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,2.4787166036823467,-2.4787166036823467) q[7];
u3(pi/2,5.571300411876139,-5.571300411876139) q[7];
rzz(pi/2) q[2],q[7];
rzz(pi/2) q[2],q[8];
u3(pi/2,3.7428934874868798,-3.7428934874868798) q[8];
u3(pi/2,0.5767964111990861,-0.5767964111990861) q[8];
rzz(pi/2) q[2],q[8];
rzz(pi/2) q[2],q[9];
u3(pi/2,3.7246722500960585,-3.7246722500960585) q[9];
u3(pi/2,0.5705132258919065,-0.5705132258919065) q[9];
rzz(pi/2) q[2],q[9];
rzz(pi/2) q[2],q[10];
rzz(-pi/2) q[2],q[10];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[12];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[14];
rzz(pi/2) q[2],q[14];
rzz(-pi/2) q[2],q[15];
rzz(pi/2) q[2],q[15];
rzz(pi/2) q[2],q[16];
rzz(pi/2) q[2],q[16];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(-pi/2) q[3],q[4];
u3(pi/2,4.516353598800687,-4.516353598800687) q[4];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,1.865477717701619,-1.865477717701619) q[5];
u3(pi/2,4.614371289592689,-4.614371289592689) q[5];
rzz(pi/2) q[3],q[5];
rzz(-pi/2) q[3],q[6];
u3(pi/2,0.7363893180014475,-0.7363893180014475) q[6];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[6];
rzz(pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,5.571300411876139,-5.571300411876139) q[7];
u3(pi/2,2.3316900674943444,-2.3316900674943444) q[7];
rzz(pi/2) q[3],q[7];
rzz(-pi/2) q[3],q[8];
u3(pi/2,3.718389064788879,-3.718389064788879) q[8];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[8];
rzz(pi/2) q[3],q[8];
rzz(pi/2) q[3],q[9];
u3(pi/2,0.5705132258919065,-0.5705132258919065) q[9];
u3(pi/2,3.687601456783699,-3.687601456783699) q[9];
rzz(pi/2) q[3],q[9];
rzz(pi/2) q[3],q[10];
u3(pi/2,5.2489730056178265,-5.2489730056178265) q[10];
u3(pi/2,2.0954422999443922,-2.0954422999443922) q[10];
rzz(-pi/2) q[3],q[10];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[12];
rzz(pi/2) q[3],q[13];
rzz(pi/2) q[3],q[13];
rzz(pi/2) q[3],q[14];
rzz(pi/2) q[3],q[14];
rzz(-pi/2) q[3],q[15];
rzz(pi/2) q[3],q[15];
rzz(pi/2) q[3],q[16];
rzz(pi/2) q[3],q[16];
u3(pi/2,5.301751762198135,-5.301751762198135) q[4];
u3(pi/2,7*pi/4,-7*pi/4) q[4];
rzz(-pi/2) q[4],q[5];
u3(pi/2,1.472778636002895,-1.472778636002895) q[5];
u3(pi/2,3.82897312619524,-3.82897312619524) q[5];
rzz(pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[6];
u3(pi/2,0.14702653618800232,-0.14702653618800232) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,2.3316900674943444,-2.3316900674943444) q[7];
u3(pi/2,5.276619020969417,-5.276619020969417) q[7];
rzz(pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[8];
u3(pi/2,3.571362528600877,-3.571362528600877) q[8];
rzz(pi/2) q[4],q[8];
rzz(-pi/2) q[4],q[9];
u3(pi/2,0.5460088031939061,-0.5460088031939061) q[9];
u3(pi/2,3.638592611387698,-3.638592611387698) q[9];
rzz(pi/2) q[4],q[9];
rzz(pi/2) q[4],q[10];
u3(pi/2,5.237034953534185,-5.237034953534185) q[10];
u3(pi/2,2.070937877246392,-2.070937877246392) q[10];
rzz(pi/2) q[4],q[10];
rzz(pi/2) q[4],q[11];
u3(pi/2,2.020044076258237,-2.020044076258237) q[11];
u3(pi/2,5.1496986777643885,-5.1496986777643885) q[11];
rzz(pi/2) q[4],q[11];
rzz(pi/2) q[4],q[12];
rzz(pi/2) q[4],q[12];
rzz(pi/2) q[4],q[13];
rzz(-pi/2) q[4],q[13];
rzz(pi/2) q[4],q[14];
rzz(pi/2) q[4],q[14];
rzz(pi/2) q[4],q[15];
rzz(-pi/2) q[4],q[15];
rzz(pi/2) q[4],q[16];
rzz(pi/2) q[4],q[16];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[5];
u3(pi/2,3*pi/4,-3*pi/4) q[5];
rzz(pi/2) q[5],q[6];
u3(pi/2,0.14702653618800232,-0.14702653618800232) q[6];
u3(pi/2,2.503221026380347,-2.503221026380347) q[6];
rzz(pi/2) q[5],q[6];
rzz(-pi/2) q[5],q[7];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[7];
u3(pi/2,4.883919939270692,-4.883919939270692) q[7];
rzz(pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,3.571362528600877,-3.571362528600877) q[8];
u3(pi/2,0.23310617489636265,-0.23310617489636265) q[8];
rzz(pi/2) q[5],q[8];
rzz(-pi/2) q[5],q[9];
u3(pi/2,0.4969999577979053,-0.4969999577979053) q[9];
u3(pi/2,3.540574920595697,-3.540574920595697) q[9];
rzz(pi/2) q[5],q[9];
rzz(pi/2) q[5],q[10];
u3(pi/2,2.070937877246392,-2.070937877246392) q[10];
u3(pi/2,5.163521685440184,-5.163521685440184) q[10];
rzz(pi/2) q[5],q[10];
rzz(pi/2) q[5],q[11];
u3(pi/2,5.1496986777643885,-5.1496986777643885) q[11];
u3(pi/2,1.9836016014765951,-1.9836016014765951) q[11];
rzz(pi/2) q[5],q[11];
rzz(pi/2) q[5],q[12];
u3(pi/2,4.843707553304744,-4.843707553304744) q[12];
u3(pi/2,1.6895485291005905,-1.6895485291005905) q[12];
rzz(pi/2) q[5],q[12];
rzz(-pi/2) q[5],q[13];
rzz(pi/2) q[5],q[13];
rzz(pi/2) q[5],q[14];
rzz(pi/2) q[5],q[14];
rzz(pi/2) q[5],q[15];
rzz(-pi/2) q[5],q[15];
rzz(pi/2) q[5],q[16];
rzz(pi/2) q[5],q[16];
u3(pi/2,4.074017353175243,-4.074017353175243) q[6];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[6];
rzz(pi/2) q[6],q[7];
u3(pi/2,4.883919939270692,-4.883919939270692) q[7];
u3(pi/2,0.956929122283451,-0.956929122283451) q[7];
rzz(-pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,0.23310617489636265,-0.23310617489636265) q[8];
u3(pi/2,2.981999746787432,-2.981999746787432) q[8];
rzz(pi/2) q[6],q[8];
rzz(pi/2) q[6],q[9];
u3(pi/2,3.540574920595697,-3.540574920595697) q[9];
u3(pi/2,0.20231856689118266,-0.20231856689118266) q[9];
rzz(pi/2) q[6],q[9];
rzz(-pi/2) q[6],q[10];
u3(pi/2,2.021929031850391,-2.021929031850391) q[10];
u3(pi/2,5.065503994648183,-5.065503994648183) q[10];
rzz(pi/2) q[6],q[10];
rzz(pi/2) q[6],q[11];
u3(pi/2,1.9836016014765951,-1.9836016014765951) q[11];
u3(pi/2,5.076185409670387,-5.076185409670387) q[11];
rzz(pi/2) q[6],q[11];
rzz(-pi/2) q[6],q[12];
u3(pi/2,4.831141182690384,-4.831141182690384) q[12];
u3(pi/2,1.6650441064025905,-1.6650441064025905) q[12];
rzz(pi/2) q[6],q[12];
rzz(pi/2) q[6],q[13];
u3(pi/2,0.11435397259066847,-0.11435397259066847) q[13];
u3(pi/2,3.2440085740968203,-3.2440085740968203) q[13];
rzz(pi/2) q[6],q[13];
rzz(pi/2) q[6],q[14];
rzz(pi/2) q[6],q[14];
rzz(pi/2) q[6],q[15];
rzz(pi/2) q[6],q[15];
rzz(pi/2) q[6],q[16];
rzz(pi/2) q[6],q[16];
u3(pi/2,2.5277254490783476,-2.5277254490783476) q[7];
u3(pi/2,2.5528581903070657,-2.5528581903070657) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,2.981999746787432,-2.981999746787432) q[8];
u3(pi/2,5.338194236979777,-5.338194236979777) q[8];
rzz(pi/2) q[7],q[8];
rzz(pi/2) q[7],q[9];
u3(pi/2,0.20231856689118266,-0.20231856689118266) q[9];
u3(pi/2,2.951212138782252,-2.951212138782252) q[9];
rzz(pi/2) q[7],q[9];
rzz(pi/2) q[7],q[10];
u3(pi/2,5.065503994648183,-5.065503994648183) q[10];
u3(pi/2,1.7272476409436681,-1.7272476409436681) q[10];
rzz(pi/2) q[7],q[10];
rzz(pi/2) q[7],q[11];
u3(pi/2,5.076185409670387,-5.076185409670387) q[11];
u3(pi/2,1.8359467467578752,-1.8359467467578752) q[11];
rzz(-pi/2) q[7],q[11];
rzz(pi/2) q[7],q[12];
u3(pi/2,1.6650441064025905,-1.6650441064025905) q[12];
u3(pi/2,4.757627914596383,-4.757627914596383) q[12];
rzz(pi/2) q[7],q[12];
rzz(pi/2) q[7],q[13];
u3(pi/2,3.2440085740968203,-3.2440085740968203) q[13];
u3(pi/2,0.07791149780902687,-0.07791149780902687) q[13];
rzz(-pi/2) q[7],q[13];
rzz(pi/2) q[7],q[14];
u3(pi/2,4.7387783586748435,-4.7387783586748435) q[14];
u3(pi/2,1.5852476530014097,-1.5852476530014097) q[14];
rzz(pi/2) q[7],q[14];
rzz(pi/2) q[7],q[15];
rzz(pi/2) q[7],q[15];
rzz(-pi/2) q[7],q[16];
rzz(pi/2) q[7],q[16];
u3(pi/2,3.7673979101848802,-3.7673979101848802) q[8];
u3(pi/2,3.779964280799239,-3.779964280799239) q[8];
rzz(pi/2) q[8],q[9];
u3(pi/2,2.951212138782252,-2.951212138782252) q[9];
u3(pi/2,5.307406628974596,-5.307406628974596) q[9];
rzz(pi/2) q[8],q[9];
rzz(-pi/2) q[8],q[10];
u3(pi/2,4.868840294533461,-4.868840294533461) q[10];
u3(pi/2,1.3345485592449442,-1.3345485592449442) q[10];
rzz(pi/2) q[8],q[10];
rzz(pi/2) q[8],q[11];
u3(pi/2,4.977539400347668,-4.977539400347668) q[11];
u3(pi/2,1.6399113651738721,-1.6399113651738721) q[11];
rzz(pi/2) q[8],q[11];
rzz(pi/2) q[8],q[12];
u3(pi/2,4.757627914596383,-4.757627914596383) q[12];
u3(pi/2,1.5180175702145882,-1.5180175702145882) q[12];
rzz(-pi/2) q[8],q[12];
rzz(pi/2) q[8],q[13];
u3(pi/2,3.21950415139882,-3.21950415139882) q[13];
u3(pi/2,0.028902652413026097,-0.028902652413026097) q[13];
rzz(pi/2) q[8],q[13];
rzz(pi/2) q[8],q[14];
u3(pi/2,1.5852476530014097,-1.5852476530014097) q[14];
u3(pi/2,4.702335883893202,-4.702335883893202) q[14];
rzz(pi/2) q[8],q[14];
rzz(pi/2) q[8],q[15];
u3(pi/2,3.5562828838636453,-3.5562828838636453) q[15];
u3(pi/2,0.40212385965949354,-0.40212385965949354) q[15];
rzz(pi/2) q[8],q[15];
rzz(pi/2) q[8],q[16];
rzz(-pi/2) q[8],q[16];
u3(pi,3.7366103021797,-3.7366103021797) q[9];
rzz(pi/2) q[9],q[10];
u3(pi/2,1.3345485592449442,-1.3345485592449442) q[10];
u3(pi/2,3.6907430494372893,-3.6907430494372893) q[10];
rzz(pi/2) q[9],q[10];
rzz(pi/2) q[9],q[11];
u3(pi/2,1.6399113651738721,-1.6399113651738721) q[11];
u3(pi/2,4.388804937064941,-4.388804937064941) q[11];
rzz(-pi/2) q[9],q[11];
rzz(pi/2) q[9],q[12];
u3(pi/2,4.659610223804381,-4.659610223804381) q[12];
u3(pi/2,1.3213538700998668,-1.3213538700998668) q[12];
rzz(pi/2) q[9],q[12];
rzz(pi/2) q[9],q[13];
u3(pi/2,0.028902652413026097,-0.028902652413026097) q[13];
u3(pi/2,3.0724776152108175,-3.0724776152108175) q[13];
rzz(pi/2) q[9],q[13];
rzz(-pi/2) q[9],q[14];
u3(pi/2,1.5607432303034092,-1.5607432303034092) q[14];
u3(pi/2,4.653327038497202,-4.653327038497202) q[14];
rzz(pi/2) q[9],q[14];
rzz(pi/2) q[9],q[15];
u3(pi/2,0.40212385965949354,-0.40212385965949354) q[15];
u3(pi/2,3.5192120905512865,-3.5192120905512865) q[15];
rzz(pi/2) q[9],q[15];
rzz(-pi/2) q[9],q[16];
u3(pi/2,0.7954512598889355,-0.7954512598889355) q[16];
u3(pi/2,3.92447754286437,-3.92447754286437) q[16];
rzz(pi/2) q[9],q[16];
u3(pi,2.1199467226423923,-2.1199467226423923) q[10];
rzz(pi/2) q[10],q[11];
u3(pi/2,1.247212283475148,-1.247212283475148) q[11];
u3(pi/2,3.603406773667493,-3.603406773667493) q[11];
rzz(pi/2) q[10],q[11];
rzz(pi/2) q[10],q[12];
u3(pi/2,1.3213538700998668,-1.3213538700998668) q[12];
u3(pi/2,4.0702474419909365,-4.0702474419909365) q[12];
rzz(pi/2) q[10],q[12];
rzz(pi/2) q[10],q[13];
u3(pi/2,3.0724776152108175,-3.0724776152108175) q[13];
u3(pi/2,6.01740656868589,-6.01740656868589) q[13];
rzz(pi/2) q[10],q[13];
rzz(pi/2) q[10],q[14];
u3(pi/2,4.653327038497202,-4.653327038497202) q[14];
u3(pi/2,1.4130883755846888,-1.4130883755846888) q[14];
rzz(-pi/2) q[10],q[14];
rzz(pi/2) q[10],q[15];
u3(pi/2,3.5192120905512865,-3.5192120905512865) q[15];
u3(pi/2,0.3286105915654924,-0.3286105915654924) q[15];
rzz(pi/2) q[10],q[15];
rzz(pi/2) q[10],q[16];
u3(pi/2,3.92447754286437,-3.92447754286437) q[16];
u3(pi/2,0.7583804665765761,-0.7583804665765761) q[16];
rzz(pi/2) q[10],q[16];
rzz(pi/2) q[11],q[12];
u3(pi/2,4.0702474419909365,-4.0702474419909365) q[12];
u3(pi/2,0.14325662500369457,-0.14325662500369457) q[12];
rzz(pi/2) q[11],q[12];
rzz(pi/2) q[11],q[13];
u3(pi/2,6.01740656868589,-6.01740656868589) q[13];
u3(pi/2,2.4831148333973725,-2.4831148333973725) q[13];
rzz(-pi/2) q[11],q[13];
rzz(pi/2) q[11],q[14];
u3(pi/2,4.554681029174482,-4.554681029174482) q[14];
u3(pi/2,1.217052994000686,-1.217052994000686) q[14];
rzz(pi/2) q[11],q[14];
rzz(pi/2) q[11],q[15];
u3(pi/2,0.3286105915654924,-0.3286105915654924) q[15];
u3(pi/2,3.3721855543632837,-3.3721855543632837) q[15];
rzz(-pi/2) q[11],q[15];
rzz(pi/2) q[11],q[16];
u3(pi/2,0.7583804665765761,-0.7583804665765761) q[16];
u3(pi/2,3.8509642747703685,-3.8509642747703685) q[16];
rzz(pi/2) q[11],q[16];
rzz(pi/2) q[12],q[13];
u3(pi/2,5.624707486987165,-5.624707486987165) q[13];
u3(pi/2,1.6977166699999242,-1.6977166699999242) q[13];
rzz(pi/2) q[12],q[13];
rzz(-pi/2) q[12],q[14];
u3(pi/2,4.358645647590479,-4.358645647590479) q[14];
u3(pi/2,0.8243539123019618,-0.8243539123019618) q[14];
rzz(pi/2) q[12],q[14];
rzz(pi/2) q[12],q[15];
u3(pi/2,0.23059290077349084,-0.23059290077349084) q[15];
u3(pi/2,3.175521854248563,-3.175521854248563) q[15];
rzz(pi/2) q[12],q[15];
rzz(-pi/2) q[12],q[16];
u3(pi/2,0.7093716211805753,-0.7093716211805753) q[16];
u3(pi/2,3.752946583978367,-3.752946583978367) q[16];
rzz(pi/2) q[12],q[16];
u3(pi,0.12692034320502762,-0.12692034320502762) q[13];
rzz(pi/2) q[13],q[14];
u3(pi/2,0.8243539123019618,-0.8243539123019618) q[14];
u3(pi/2,3.1805484024943063,-3.1805484024943063) q[14];
rzz(pi/2) q[13],q[14];
rzz(pi/2) q[13],q[15];
u3(pi/2,3.175521854248563,-3.175521854248563) q[15];
u3(pi/2,5.924415426139632,-5.924415426139632) q[15];
rzz(-pi/2) q[13],q[15];
rzz(pi/2) q[13],q[16];
u3(pi/2,3.752946583978367,-3.752946583978367) q[16];
u3(pi/2,0.4146902302738527,-0.4146902302738527) q[16];
rzz(pi/2) q[13],q[16];
u3(pi,1.60975207569941,-1.60975207569941) q[14];
rzz(pi/2) q[14],q[15];
u3(pi/2,2.7828227725498387,-2.7828227725498387) q[15];
u3(pi/2,5.1390172627421835,-5.1390172627421835) q[15];
rzz(-pi/2) q[14],q[15];
rzz(pi/2) q[14],q[16];
u3(pi/2,0.4146902302738527,-0.4146902302738527) q[16];
u3(pi/2,3.163583802164921,-3.163583802164921) q[16];
rzz(pi/2) q[14],q[16];
u3(pi,0.4266282823574939,-0.4266282823574939) q[15];
rzz(pi/2) q[15],q[16];
u3(pi/2,3.163583802164921,-3.163583802164921) q[16];
u3(pi/2,5.519778292357266,-5.519778292357266) q[16];
rzz(-pi/2) q[15],q[16];
u3(pi,pi,-pi) q[1];
u3(pi,0,0) q[2];
u3(pi,3*pi/4,-3*pi/4) q[3];
u3(pi,0,0) q[4];
u3(pi,7*pi/8,-7*pi/8) q[5];
u3(pi,4.221672207893964,-4.221672207893964) q[6];
u3(pi,4.221672207893964,-4.221672207893964) q[7];
u3(pi,1.472778636002895,-1.472778636002895) q[8];
u3(pi,4.945495155281052,-4.945495155281052) q[9];
u3(pi,1.1535928223981722,-1.1535928223981722) q[10];
u3(pi,3.1170882308917927,-3.1170882308917927) q[11];
u3(pi,0.6873804726054468,-0.6873804726054468) q[12];
u3(pi,5.896141092257324,-5.896141092257324) q[13];
u3(pi,2.9235661234306614,-2.9235661234306614) q[14];
u3(pi,1.0838494654884785,-1.0838494654884785) q[15];
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
measure q[11] -> c[11];
measure q[12] -> c[12];
measure q[13] -> c[13];
measure q[14] -> c[14];
measure q[15] -> c[15];
measure q[16] -> c[16];
