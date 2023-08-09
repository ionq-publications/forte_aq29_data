OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(pi/2,4.829884545628948,-4.829884545628948) q[0];
u3(pi/2,3.335114761050925,-3.335114761050925) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi,-pi) q[2];
u3(pi/2,1.4953981031087415,-1.4953981031087415) q[0];
rzz(pi/2) q[2],q[0];
u3(pi/2,4.636990756698535,-4.636990756698535) q[0];
u3(pi/2,4.905911087845821,-4.905911087845821) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[1];
u3(pi/2,3.4658050154402598,-3.4658050154402598) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
u3(pi/2,1.895008688645363,-1.895008688645363) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.4658050154402598,-3.4658050154402598) q[0];
u3(pi/2,6.2077870834934314,-6.2077870834934314) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,5.691309251243269,-5.691309251243269) q[3];
u3(pi/2,2.549716597653476,-2.549716597653476) q[4];
rzz(-pi/2) q[3],q[4];
u3(pi/2,0,0) q[5];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[5],q[3];
u3(pi/2,pi,-pi) q[3];
u3(pi/2,0.9789202708585795,-0.9789202708585795) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,1.957212223186441,-1.957212223186441) q[4];
u3(pi/2,3.812008525865855,-3.812008525865855) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[5];
rzz(-pi/2) q[4],q[5];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
u3(pi/2,5.382804852660752,-5.382804852660752) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,0.6704158722760619,-0.6704158722760619) q[3];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[5];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[1];
rzz(-pi/2) q[3],q[1];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
u3(pi/2,1.4953981031087415,-1.4953981031087415) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.066194429903638,-3.066194429903638) q[0];
u3(pi/2,5.74031809663927,-5.74031809663927) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[3],q[1];
u3(pi/2,5.098804876776234,-5.098804876776234) q[4];
u3(pi/2,pi,-pi) q[3];
rzz(-pi/2) q[4],q[3];
rzz(pi/2) q[4],q[1];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,2.598725443049477,-2.598725443049477) q[0];
u3(pi/2,5.273477428315827,-5.273477428315827) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,5.098804876776234,-5.098804876776234) q[4];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(-pi/2) q[4],q[1];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,2.1318847747260334,-2.1318847747260334) q[0];
u3(pi/2,4.806636759992384,-4.806636759992384) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,5.098804876776234,-5.098804876776234) q[4];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[1];
rzz(pi/2) q[4],q[1];
u3(pi/2,pi,-pi) q[3];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
rzz(-pi/2) q[3],q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,5.098804876776234,-5.098804876776234) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,1.957212223186441,-1.957212223186441) q[4];
rzz(pi/2) q[4],q[1];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,4.806636759992384,-4.806636759992384) q[0];
u3(pi/2,1.198203438079147,-1.198203438079147) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[1];
rzz(pi/2) q[4],q[1];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[3],q[4];
rzz(-pi/2) q[3],q[1];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
u3(pi/2,2.7689997648740436,-2.7689997648740436) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.198203438079147,-1.198203438079147) q[0];
u3(pi/2,3.8729554233454966,-3.8729554233454966) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(-pi/2) q[3],q[1];
u3(pi/2,pi,-pi) q[3];
rzz(-pi/2) q[4],q[3];
rzz(pi/2) q[4],q[1];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(-pi/2) q[1],q[0];
u3(pi/2,3.8729554233454966,-3.8729554233454966) q[0];
u3(pi/2,0.26452210143226057,-0.26452210143226057) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[4],q[1];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.406114755022054,-3.406114755022054) q[0];
u3(pi/2,6.080866740288403,-6.080866740288403) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[1];
rzz(-pi/2) q[4],q[1];
u3(pi/2,pi,-pi) q[3];
u3(pi/2,5.098804876776234,-5.098804876776234) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
rzz(pi/2) q[4],q[1];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[1];
rzz(-pi/2) q[1],q[0];
u3(pi/2,6.080866740288403,-6.080866740288403) q[0];
u3(pi/2,2.4724334183751675,-2.4724334183751675) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[1];
rzz(pi/2) q[4],q[1];
u3(pi/2,pi/2,-pi/2) q[5];
rzz(pi/2) q[5],q[3];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[1];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
u3(pi/2,0.9016370915802705,-0.9016370915802705) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,5.61402607196496,-5.61402607196496) q[0];
u3(pi/2,2.2280175099258814,-2.2280175099258814) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(-pi/2) q[3],q[1];
u3(pi/2,1.957212223186441,-1.957212223186441) q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[3],q[1];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(-pi/2) q[1],q[0];
u3(pi/2,2.2280175099258814,-2.2280175099258814) q[0];
u3(pi/2,5.125194255066388,-5.125194255066388) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(pi/2) q[3],q[1];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.9836016014765951,-1.9836016014765951) q[0];
u3(pi/2,4.8814066651478205,-4.8814066651478205) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[1];
rzz(pi/2) q[3],q[1];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(-pi/2) q[3],q[1];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,4.8814066651478205,-4.8814066651478205) q[0];
u3(pi/2,1.4953981031087415,-1.4953981031087415) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[1];
rzz(pi/2) q[3],q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[3],q[1];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
u3(pi/2,3.066194429903638,-3.066194429903638) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.4953981031087415,-1.4953981031087415) q[0];
u3(pi/2,4.392574848249249,-4.392574848249249) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[3],q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[3],q[1];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.2509821946594557,-1.2509821946594557) q[0];
u3(pi/2,4.148158939799963,-4.148158939799963) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[3],q[1];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.0065662862101699,-1.0065662862101699) q[0];
u3(pi/2,3.9037430313506767,-3.9037430313506767) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[1];
rzz(pi/2) q[3],q[1];
u3(pi/2,pi,-pi) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[3],q[1];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[1];
rzz(-pi/2) q[1],q[0];
u3(pi/2,3.9037430313506767,-3.9037430313506767) q[0];
u3(pi/2,0.517734469311598,-0.517734469311598) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[1];
rzz(pi/2) q[3],q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,pi/2,-pi/2) q[5];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
u3(pi/2,0,0) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,0,0) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
u3(pi/2,2.0885307961064945,-2.0885307961064945) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.659327122901391,-3.659327122901391) q[0];
u3(pi/2,0.28023006470020956,-0.28023006470020956) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,1.957212223186441,-1.957212223186441) q[4];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[4],q[3];
rzz(pi/2) q[4],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(-pi/2) q[2],q[1];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,0.28023006470020956,-0.28023006470020956) q[0];
u3(pi/2,3.1843183136786144,-3.1843183136786144) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(-pi/2) q[2],q[1];
u3(pi/2,5.098804876776234,-5.098804876776234) q[4];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[4],q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[2],q[1];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,0.042725660088821185,-0.042725660088821185) q[0];
u3(pi/2,2.9468139090672256,-2.9468139090672256) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(-pi/2) q[2],q[1];
u3(pi/2,1.957212223186441,-1.957212223186441) q[4];
u3(pi/2,0,0) q[2];
rzz(pi/2) q[4],q[2];
u3(pi/2,pi,-pi) q[3];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,5.098804876776234,-5.098804876776234) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,1.957212223186441,-1.957212223186441) q[4];
rzz(-pi/2) q[4],q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[2],q[1];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,2.9468139090672256,-2.9468139090672256) q[0];
u3(pi/2,5.850902158045631,-5.850902158045631) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,pi,-pi) q[2];
rzz(-pi/2) q[4],q[2];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[3],q[4];
rzz(-pi/2) q[3],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[2],q[1];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
u3(pi/2,1.138513177660941,-1.138513177660941) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,2.7093095044558377,-2.7093095044558377) q[0];
u3(pi/2,5.613397753434242,-5.613397753434242) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(-pi/2) q[3],q[2];
u3(pi/2,0,0) q[3];
rzz(-pi/2) q[4],q[3];
rzz(pi/2) q[4],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[2],q[1];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,5.613397753434242,-5.613397753434242) q[0];
u3(pi/2,2.234300695233061,-2.234300695233061) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(-pi/2) q[2],q[1];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/2) q[4],q[2];
u3(pi/2,0,0) q[2];
rzz(pi/2) q[2],q[1];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,5.375893348822854,-5.375893348822854) q[0];
u3(pi/2,1.9967962906216727,-1.9967962906216727) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
u3(pi/2,0,0) q[2];
rzz(pi/2) q[4],q[2];
u3(pi/2,0,0) q[3];
u3(pi/2,1.957212223186441,-1.957212223186441) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
rzz(-pi/2) q[3],q[4];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
rzz(pi/2) q[4],q[2];
u3(pi/2,0,0) q[2];
rzz(pi/2) q[2],q[1];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,5.138388944211465,-5.138388944211465) q[0];
u3(pi/2,1.7592918860102844,-1.7592918860102844) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,0,0) q[2];
rzz(pi/2) q[4],q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[5];
rzz(pi/2) q[5],q[3];
u3(pi/2,pi,-pi) q[3];
rzz(-pi/2) q[3],q[4];
rzz(-pi/2) q[3],q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/2) q[2],q[1];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
u3(pi/2,3.330088212805181,-3.330088212805181) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.7592918860102844,-1.7592918860102844) q[0];
u3(pi/2,4.863185427757,-4.863185427757) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(-pi/2) q[2],q[1];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,5.098804876776234,-5.098804876776234) q[4];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
u3(pi/2,pi/2,-pi/2) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,0,0) q[3];
rzz(-pi/2) q[3],q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/2) q[2],q[1];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.7215927741672068,-1.7215927741672068) q[0];
u3(pi/2,4.825486315913922,-4.825486315913922) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(-pi/2) q[2],q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[2],q[1];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.6838936623241292,-1.6838936623241292) q[0];
u3(pi/2,4.787158885540127,-4.787158885540127) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
rzz(-pi/2) q[5],q[4];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[3],q[2];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[2],q[1];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,4.787158885540127,-4.787158885540127) q[0];
u3(pi/2,1.6078671201072563,-1.6078671201072563) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,pi,-pi) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
rzz(-pi/2) q[5],q[4];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[3],q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/2) q[2],q[1];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
u3(pi/2,1.6078671201072563,-1.6078671201072563) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,0.03707079331235956,-0.03707079331235956) q[0];
u3(pi/2,3.140964335059075,-3.140964335059075) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,pi,-pi) q[3];
rzz(pi/2) q[3],q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/2) q[2],q[1];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(-pi/2) q[1],q[0];
u3(pi/2,3.140964335059075,-3.140964335059075) q[0];
u3(pi/2,6.244857876805791,-6.244857876805791) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,0,0) q[2];
rzz(-pi/2) q[2],q[1];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.1032652232159976,-3.1032652232159976) q[0];
u3(pi/2,6.207158764962713,-6.207158764962713) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,5.098804876776234,-5.098804876776234) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,pi,-pi) q[2];
rzz(-pi/2) q[3],q[2];
u3(pi/2,0,0) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
rzz(pi/2) q[4],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[3];
rzz(pi/2) q[3],q[2];
u3(pi/2,0,0) q[2];
rzz(-pi/2) q[2],q[1];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.0655661113729202,-3.0655661113729202) q[0];
u3(pi/2,6.168831334588917,-6.168831334588917) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.957212223186441,-1.957212223186441) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,pi,-pi) q[2];
rzz(-pi/2) q[3],q[2];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3*pi/2,-3*pi/2) q[5];
u3(pi/2,3.5280085499813376,-3.5280085499813376) q[4];
rzz(pi/2) q[5],q[4];
u3(pi,1.229619364615045,-1.229619364615045) q[0];
u3(pi/2,0.38641589639154456,-0.38641589639154456) q[1];
u3(pi,3*pi/2,-3*pi/2) q[2];
u3(pi,2.262575029115369,-2.262575029115369) q[3];
u3(pi/2,1.957212223186441,-1.957212223186441) q[4];
u3(pi/2,pi,-pi) q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
