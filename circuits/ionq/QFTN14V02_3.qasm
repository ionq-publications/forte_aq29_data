OPENQASM 2.0;
include "qelib1.inc";
qreg q[14];
creg c[14];
u3(pi/2,1.7065131294299756,-1.7065131294299756) q[13];
u3(pi/2,2.491911292827424,-2.491911292827424) q[13];
rzz(-pi/2) q[12],q[13];
rzz(-pi/2) q[11],q[13];
u3(pi/2,5.633503946417217,-5.633503946417217) q[13];
u3(pi/2,2.0992122111287,-2.0992122111287) q[13];
rzz(-pi/2) q[11],q[13];
rzz(-pi/2) q[11],q[12];
u3(pi/2,1.8007609090376695,-1.8007609090376695) q[12];
u3(pi/2,2.586159072435118,-2.586159072435118) q[12];
rzz(-pi/2) q[11],q[12];
rzz(-pi/2) q[10],q[13];
u3(pi/2,5.240804864718492,-5.240804864718492) q[13];
u3(pi/2,1.9031768295446967,-1.9031768295446967) q[13];
rzz(-pi/2) q[10],q[13];
rzz(-pi/2) q[10],q[12];
u3(pi/2,2.586159072435118,-2.586159072435118) q[12];
u3(pi/2,5.335052644326186,-5.335052644326186) q[12];
rzz(pi/2) q[10],q[12];
rzz(pi/2) q[10],q[11];
u3(pi/2,5.252114598271416,-5.252114598271416) q[11];
u3(pi/2,6.037512761668864,-6.037512761668864) q[11];
rzz(-pi/2) q[10],q[11];
rzz(-pi/2) q[9],q[13];
u3(pi/2,1.9031768295446967,-1.9031768295446967) q[13];
u3(pi/2,4.94612347381177,-4.94612347381177) q[13];
rzz(pi/2) q[9],q[13];
rzz(-pi/2) q[9],q[12];
u3(pi/2,2.1934599907363936,-2.1934599907363936) q[12];
u3(pi/2,5.1390172627421835,-5.1390172627421835) q[12];
rzz(-pi/2) q[9],q[12];
rzz(-pi/2) q[9],q[11];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[11];
u3(pi/2,5.64481367997014,-5.64481367997014) q[11];
rzz(-pi/2) q[9],q[11];
rzz(-pi/2) q[9],q[10];
u3(pi/2,2.491282974296706,-2.491282974296706) q[10];
u3(pi/2,3.276681137694154,-3.276681137694154) q[10];
rzz(-pi/2) q[9],q[10];
rzz(-pi/2) q[8],q[13];
u3(pi/2,1.8045308202219772,-1.8045308202219772) q[13];
u3(pi/2,4.897114628415769,-4.897114628415769) q[13];
rzz(-pi/2) q[8],q[13];
rzz(pi/2) q[8],q[12];
u3(pi/2,1.9974246091523906,-1.9974246091523906) q[12];
u3(pi/2,5.0403712534194645,-5.0403712534194645) q[12];
rzz(-pi/2) q[8],q[12];
rzz(-pi/2) q[8],q[11];
u3(pi/2,5.64481367997014,-5.64481367997014) q[11];
u3(pi/2,2.3071856447963444,-2.3071856447963444) q[11];
rzz(-pi/2) q[8],q[11];
rzz(pi/2) q[8],q[10];
u3(pi/2,3.276681137694154,-3.276681137694154) q[10];
u3(pi/2,6.025574709585223,-6.025574709585223) q[10];
rzz(-pi/2) q[8],q[10];
rzz(pi/2) q[8],q[9];
u3(pi/2,3.730955435403238,-3.730955435403238) q[9];
u3(pi/2,4.516353598800687,-4.516353598800687) q[9];
rzz(-pi/2) q[8],q[9];
rzz(pi/2) q[7],q[13];
u3(pi/2,1.7555219748259763,-1.7555219748259763) q[13];
u3(pi/2,4.872610205717769,-4.872610205717769) q[13];
rzz(-pi/2) q[7],q[13];
rzz(-pi/2) q[7],q[12];
u3(pi/2,5.0403712534194645,-5.0403712534194645) q[12];
u3(pi/2,1.8497697544336702,-1.8497697544336702) q[12];
rzz(-pi/2) q[7],q[12];
rzz(-pi/2) q[7],q[11];
u3(pi/2,2.3071856447963444,-2.3071856447963444) q[11];
u3(pi/2,5.350760607594136,-5.350760607594136) q[11];
rzz(-pi/2) q[7],q[11];
rzz(-pi/2) q[7],q[10];
u3(pi/2,6.025574709585223,-6.025574709585223) q[10];
u3(pi/2,2.687318355880709,-2.687318355880709) q[10];
rzz(-pi/2) q[7],q[10];
rzz(-pi/2) q[7],q[9];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[9];
u3(pi/2,4.123654517101962,-4.123654517101962) q[9];
rzz(pi/2) q[7],q[9];
rzz(pi/2) q[7],q[8];
u3(pi/2,3.82897312619524,-3.82897312619524) q[8];
u3(pi/2,4.614371289592689,-4.614371289592689) q[8];
rzz(-pi/2) q[7],q[8];
rzz(-pi/2) q[6],q[13];
u3(pi/2,1.7310175521279763,-1.7310175521279763) q[13];
u3(pi/2,4.860672153634128,-4.860672153634128) q[13];
rzz(pi/2) q[6],q[13];
rzz(-pi/2) q[6],q[12];
u3(pi/2,4.9913624080234635,-4.9913624080234635) q[12];
u3(pi/2,1.8252653317356697,-1.8252653317356697) q[12];
rzz(-pi/2) q[6],q[12];
rzz(-pi/2) q[6],q[11];
u3(pi/2,2.2091679540043425,-2.2091679540043425) q[11];
u3(pi/2,5.301123443667417,-5.301123443667417) q[11];
rzz(-pi/2) q[6],q[11];
rzz(pi/2) q[6],q[10];
u3(pi/2,2.687318355880709,-2.687318355880709) q[10];
u3(pi/2,5.7308933186785005,-5.7308933186785005) q[10];
rzz(-pi/2) q[6],q[10];
rzz(-pi/2) q[6],q[9];
u3(pi/2,4.123654517101962,-4.123654517101962) q[9];
u3(pi/2,pi/4,-pi/4) q[9];
rzz(-pi/2) q[6],q[9];
rzz(pi/2) q[6],q[8];
u3(pi/2,1.472778636002895,-1.472778636002895) q[8];
u3(pi/2,4.221672207893964,-4.221672207893964) q[8];
rzz(-pi/2) q[6],q[8];
rzz(pi/2) q[6],q[7];
u3(pi/2,pi/4,-pi/4) q[7];
u3(pi/2,pi/2,-pi/2) q[7];
rzz(-pi/2) q[6],q[7];
rzz(pi/2) q[5],q[13];
rzz(-pi/2) q[5],q[13];
rzz(-pi/2) q[5],q[12];
u3(pi/2,1.8252653317356697,-1.8252653317356697) q[12];
u3(pi/2,4.954919933241822,-4.954919933241822) q[12];
rzz(-pi/2) q[5],q[12];
rzz(pi/2) q[5],q[11];
u3(pi/2,2.159530790077624,-2.159530790077624) q[11];
u3(pi/2,5.276619020969417,-5.276619020969417) q[11];
rzz(-pi/2) q[5],q[11];
rzz(-pi/2) q[5],q[10];
u3(pi/2,5.7308933186785005,-5.7308933186785005) q[10];
u3(pi/2,2.5402918196927065,-2.5402918196927065) q[10];
rzz(-pi/2) q[5],q[10];
rzz(-pi/2) q[5],q[9];
u3(pi/2,pi/4,-pi/4) q[9];
u3(pi/2,3.82897312619524,-3.82897312619524) q[9];
rzz(pi/2) q[5],q[9];
rzz(-pi/2) q[5],q[8];
u3(pi/2,4.221672207893964,-4.221672207893964) q[8];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[8];
rzz(-pi/2) q[5],q[8];
rzz(-pi/2) q[5],q[7];
u3(pi/2,3*pi/2,-3*pi/2) q[7];
u3(pi/2,3*pi/8,-3*pi/8) q[7];
rzz(pi/2) q[5],q[7];
rzz(pi/2) q[5],q[6];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
u3(pi/2,4.516353598800687,-4.516353598800687) q[6];
rzz(-pi/2) q[5],q[6];
rzz(-pi/2) q[4],q[13];
rzz(-pi/2) q[4],q[13];
rzz(pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[11];
u3(pi/2,2.1350263673796235,-2.1350263673796235) q[11];
u3(pi/2,5.264680968885775,-5.264680968885775) q[11];
rzz(-pi/2) q[4],q[11];
rzz(pi/2) q[4],q[10];
u3(pi/2,2.5402918196927065,-2.5402918196927065) q[10];
u3(pi/2,5.657380050584499,-5.657380050584499) q[10];
rzz(-pi/2) q[4],q[10];
rzz(-pi/2) q[4],q[9];
u3(pi/2,3.82897312619524,-3.82897312619524) q[9];
u3(pi/2,0.638371627209446,-0.638371627209446) q[9];
rzz(-pi/2) q[4],q[9];
rzz(-pi/2) q[4],q[8];
u3(pi/2,4.025008507779242,-4.025008507779242) q[8];
u3(pi/2,pi/4,-pi/4) q[8];
rzz(-pi/2) q[4],q[8];
rzz(-pi/2) q[4],q[7];
u3(pi/2,3*pi/8,-3*pi/8) q[7];
u3(pi/2,4.123654517101962,-4.123654517101962) q[7];
rzz(-pi/2) q[4],q[7];
rzz(pi/2) q[4],q[6];
u3(pi/2,1.3747609452108935,-1.3747609452108935) q[6];
u3(pi/2,4.123026198571244,-4.123026198571244) q[6];
rzz(-pi/2) q[4],q[6];
rzz(pi/2) q[4],q[5];
u3(pi/2,3*pi/2,-3*pi/2) q[5];
u3(pi/2,7*pi/4,-7*pi/4) q[5];
rzz(-pi/2) q[4],q[5];
rzz(-pi/2) q[3],q[13];
rzz(-pi/2) q[3],q[13];
rzz(-pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[11];
rzz(pi/2) q[3],q[11];
rzz(-pi/2) q[3],q[10];
u3(pi/2,2.515787396994706,-2.515787396994706) q[10];
u3(pi/2,5.64481367997014,-5.64481367997014) q[10];
rzz(-pi/2) q[3],q[10];
rzz(-pi/2) q[3],q[9];
u3(pi/2,3.779964280799239,-3.779964280799239) q[9];
u3(pi/2,0.6138672045114455,-0.6138672045114455) q[9];
rzz(pi/2) q[3],q[9];
rzz(-pi/2) q[3],q[8];
u3(pi/2,5*pi/4,-5*pi/4) q[8];
u3(pi/2,0.7363893180014475,-0.7363893180014475) q[8];
rzz(-pi/2) q[3],q[8];
rzz(-pi/2) q[3],q[7];
u3(pi/2,0.9820618635121693,-0.9820618635121693) q[7];
u3(pi/2,4.025008507779242,-4.025008507779242) q[7];
rzz(-pi/2) q[3],q[7];
rzz(pi/2) q[3],q[6];
u3(pi/2,4.123026198571244,-4.123026198571244) q[6];
u3(pi/2,pi/4,-pi/4) q[6];
rzz(-pi/2) q[3],q[6];
rzz(-pi/2) q[3],q[5];
u3(pi/2,7*pi/4,-7*pi/4) q[5];
u3(pi/2,5*pi/8,-5*pi/8) q[5];
rzz(-pi/2) q[3],q[5];
rzz(-pi/2) q[3],q[4];
u3(pi/2,pi/4,-pi/4) q[4];
u3(pi/2,pi/2,-pi/2) q[4];
rzz(-pi/2) q[3],q[4];
rzz(-pi/2) q[2],q[13];
rzz(-pi/2) q[2],q[13];
rzz(-pi/2) q[2],q[12];
rzz(pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[11];
rzz(-pi/2) q[2],q[11];
rzz(pi/2) q[2],q[10];
rzz(-pi/2) q[2],q[10];
rzz(-pi/2) q[2],q[9];
u3(pi/2,3.7554598581012386,-3.7554598581012386) q[9];
u3(pi/2,0.6013008338970863,-0.6013008338970863) q[9];
rzz(-pi/2) q[2],q[9];
rzz(-pi/2) q[2],q[8];
u3(pi/2,0.7363893180014475,-0.7363893180014475) q[8];
u3(pi/2,3.85347754889324,-3.85347754889324) q[8];
rzz(pi/2) q[2],q[8];
rzz(-pi/2) q[2],q[7];
u3(pi/2,4.025008507779242,-4.025008507779242) q[7];
u3(pi/2,0.8344070087934491,-0.8344070087934491) q[7];
rzz(-pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[6];
u3(pi/2,pi/4,-pi/4) q[6];
u3(pi/2,3.82897312619524,-3.82897312619524) q[6];
rzz(pi/2) q[2],q[6];
rzz(-pi/2) q[2],q[5];
u3(pi/2,5*pi/8,-5*pi/8) q[5];
u3(pi/2,4.9084243619686925,-4.9084243619686925) q[5];
rzz(-pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[4];
u3(pi/2,3*pi/2,-3*pi/2) q[4];
u3(pi/2,3*pi/8,-3*pi/8) q[4];
rzz(-pi/2) q[2],q[4];
rzz(-pi/2) q[2],q[3];
u3(pi/2,pi/4,-pi/4) q[3];
u3(pi/2,pi/2,-pi/2) q[3];
rzz(-pi/2) q[2],q[3];
rzz(-pi/2) q[1],q[13];
rzz(-pi/2) q[1],q[13];
rzz(pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[11];
rzz(-pi/2) q[1],q[10];
rzz(-pi/2) q[1],q[10];
rzz(-pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[9];
rzz(-pi/2) q[1],q[8];
u3(pi/2,3.85347754889324,-3.85347754889324) q[8];
u3(pi/2,0.699318524689088,-0.699318524689088) q[8];
rzz(pi/2) q[1],q[8];
rzz(-pi/2) q[1],q[7];
u3(pi/2,3.9759996623832423,-3.9759996623832423) q[7];
u3(pi/2,0.8099025860954486,-0.8099025860954486) q[7];
rzz(-pi/2) q[1],q[7];
rzz(-pi/2) q[1],q[6];
u3(pi/2,3.82897312619524,-3.82897312619524) q[6];
u3(pi/2,0.638371627209446,-0.638371627209446) q[6];
rzz(-pi/2) q[1],q[6];
rzz(-pi/2) q[1],q[5];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[5];
u3(pi/2,4.810406671176691,-4.810406671176691) q[5];
rzz(-pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[4];
u3(pi/2,11*pi/8,-11*pi/8) q[4];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[4];
rzz(pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[3];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
rzz(-pi/2) q[1],q[3];
rzz(pi/2) q[1],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
rzz(pi/2) q[1],q[2];
rzz(-pi/2) q[0],q[13];
rzz(-pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[10];
rzz(pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[7];
rzz(-pi/2) q[0],q[6];
rzz(pi/2) q[0],q[5];
rzz(-pi/2) q[0],q[4];
rzz(-pi/2) q[0],q[3];
rzz(-pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi/2,pi,-pi) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,3*pi/4,-3*pi/4) q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/2) q[0],q[2];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
u3(pi/2,pi/4,-pi/4) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,4.123026198571244,-4.123026198571244) q[4];
u3(pi/2,pi/4,-pi/4) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[5];
u3(pi/2,3*pi/2,-3*pi/2) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,0.638371627209446,-0.638371627209446) q[6];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,0.8099025860954486,-0.8099025860954486) q[7];
u3(pi/2,5*pi/4,-5*pi/4) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,3.8409111782788807,-3.8409111782788807) q[8];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[8];
rzz(pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[11];
rzz(pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[13];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,pi/2,-pi/2) q[2];
u3(pi/2,5*pi/4,-5*pi/4) q[2];
rzz(pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
u3(pi/2,pi/8,-pi/8) q[3];
rzz(pi/2) q[1],q[3];
rzz(-pi/2) q[1],q[4];
u3(pi/2,pi/4,-pi/4) q[4];
u3(pi/2,3.730955435403238,-3.730955435403238) q[4];
rzz(pi/2) q[1],q[4];
rzz(pi/2) q[1],q[5];
u3(pi/2,3*pi/2,-3*pi/2) q[5];
u3(pi/2,1.472778636002895,-1.472778636002895) q[5];
rzz(pi/2) q[1],q[5];
rzz(-pi/2) q[1],q[6];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[6];
rzz(pi/2) q[1],q[6];
rzz(pi/2) q[1],q[7];
u3(pi/2,pi/4,-pi/4) q[7];
u3(pi/2,3.902486394289241,-3.902486394289241) q[7];
rzz(pi/2) q[1],q[7];
rzz(pi/2) q[1],q[8];
u3(pi/2,3.82897312619524,-3.82897312619524) q[8];
u3(pi/2,0.6748141019910875,-0.6748141019910875) q[8];
rzz(pi/2) q[1],q[8];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[9];
rzz(pi/2) q[1],q[10];
rzz(-pi/2) q[1],q[10];
rzz(pi/2) q[1],q[11];
rzz(pi/2) q[1],q[11];
rzz(pi/2) q[1],q[12];
rzz(-pi/2) q[1],q[12];
rzz(pi/2) q[1],q[13];
rzz(pi/2) q[1],q[13];
u3(pi/2,7*pi/4,-7*pi/4) q[2];
u3(pi/2,pi/2,-pi/2) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,pi/8,-pi/8) q[3];
u3(pi/2,7*pi/8,-7*pi/8) q[3];
rzz(-pi/2) q[2],q[3];
rzz(pi/2) q[2],q[4];
u3(pi/2,3.730955435403238,-3.730955435403238) q[4];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,1.472778636002895,-1.472778636002895) q[5];
u3(pi/2,4.417707589477967,-4.417707589477967) q[5];
rzz(-pi/2) q[2],q[5];
rzz(pi/2) q[2],q[6];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[6];
u3(pi/2,3.583300580684518,-3.583300580684518) q[6];
rzz(pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,3.902486394289241,-3.902486394289241) q[7];
u3(pi/2,0.7118848953034471,-0.7118848953034471) q[7];
rzz(pi/2) q[2],q[7];
rzz(-pi/2) q[2],q[8];
u3(pi/2,3.816406755580881,-3.816406755580881) q[8];
u3(pi/2,0.6503096792930871,-0.6503096792930871) q[8];
rzz(pi/2) q[2],q[8];
rzz(pi/2) q[2],q[9];
u3(pi/2,3.7246722500960585,-3.7246722500960585) q[9];
u3(pi/2,0.5705132258919065,-0.5705132258919065) q[9];
rzz(pi/2) q[2],q[9];
rzz(-pi/2) q[2],q[10];
rzz(pi/2) q[2],q[10];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[12];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[13];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
u3(pi/2,5*pi/4,-5*pi/4) q[3];
rzz(pi/2) q[3],q[4];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[4];
u3(pi/2,2.552229871776348,-2.552229871776348) q[4];
rzz(-pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,1.276114935888174,-1.276114935888174) q[5];
u3(pi/2,4.025008507779242,-4.025008507779242) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,3.583300580684518,-3.583300580684518) q[6];
u3(pi/2,0.24567254551072185,-0.24567254551072185) q[6];
rzz(pi/2) q[3],q[6];
rzz(-pi/2) q[3],q[7];
u3(pi/2,3.85347754889324,-3.85347754889324) q[7];
u3(pi/2,0.6138672045114455,-0.6138672045114455) q[7];
rzz(pi/2) q[3],q[7];
rzz(pi/2) q[3],q[8];
u3(pi/2,0.6503096792930871,-0.6503096792930871) q[8];
u3(pi/2,3.7428934874868798,-3.7428934874868798) q[8];
rzz(-pi/2) q[3],q[8];
rzz(pi/2) q[3],q[9];
u3(pi/2,0.5705132258919065,-0.5705132258919065) q[9];
u3(pi/2,3.687601456783699,-3.687601456783699) q[9];
rzz(pi/2) q[3],q[9];
rzz(pi/2) q[3],q[10];
u3(pi/2,5.629734035232909,-5.629734035232909) q[10];
u3(pi/2,2.475575011028757,-2.475575011028757) q[10];
rzz(pi/2) q[3],q[10];
rzz(-pi/2) q[3],q[11];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[12];
rzz(pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[13];
rzz(pi/2) q[3],q[13];
u3(pi/2,4.123026198571244,-4.123026198571244) q[4];
u3(pi/2,11*pi/8,-11*pi/8) q[4];
rzz(pi/2) q[4],q[5];
u3(pi/2,4.025008507779242,-4.025008507779242) q[5];
u3(pi/2,0.09801769079200154,-0.09801769079200154) q[5];
rzz(pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,0.24567254551072185,-0.24567254551072185) q[6];
u3(pi/2,2.994566117401791,-2.994566117401791) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,0.6138672045114455,-0.6138672045114455) q[7];
u3(pi/2,3.5587961579865177,-3.5587961579865177) q[7];
rzz(pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,0.6013008338970863,-0.6013008338970863) q[8];
u3(pi/2,3.6448757966948775,-3.6448757966948775) q[8];
rzz(-pi/2) q[4],q[8];
rzz(pi/2) q[4],q[9];
u3(pi/2,3.687601456783699,-3.687601456783699) q[9];
u3(pi/2,0.4969999577979053,-0.4969999577979053) q[9];
rzz(pi/2) q[4],q[9];
rzz(pi/2) q[4],q[10];
u3(pi/2,2.475575011028757,-2.475575011028757) q[10];
u3(pi/2,5.592663241920549,-5.592663241920549) q[10];
rzz(-pi/2) q[4],q[10];
rzz(pi/2) q[4],q[11];
u3(pi/2,5.24771636855639,-5.24771636855639) q[11];
u3(pi/2,2.094185662882956,-2.094185662882956) q[11];
rzz(pi/2) q[4],q[11];
rzz(pi/2) q[4],q[12];
rzz(pi/2) q[4],q[12];
rzz(-pi/2) q[4],q[13];
rzz(pi/2) q[4],q[13];
u3(pi/2,1.6688140175868982,-1.6688140175868982) q[5];
u3(pi/2,3*pi/2,-3*pi/2) q[5];
rzz(pi/2) q[5],q[6];
u3(pi/2,2.994566117401791,-2.994566117401791) q[6];
u3(pi/2,5.350760607594136,-5.350760607594136) q[6];
rzz(-pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,3.5587961579865177,-3.5587961579865177) q[7];
u3(pi/2,0.024504422698000385,-0.024504422698000385) q[7];
rzz(pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,0.5032831431050849,-0.5032831431050849) q[8];
u3(pi/2,3.4482120965801566,-3.4482120965801566) q[8];
rzz(pi/2) q[5],q[8];
rzz(-pi/2) q[5],q[9];
u3(pi/2,3.638592611387698,-3.638592611387698) q[9];
u3(pi/2,0.3989822670059037,-0.3989822670059037) q[9];
rzz(pi/2) q[5],q[9];
rzz(pi/2) q[5],q[10];
u3(pi/2,2.4510705883307566,-2.4510705883307566) q[10];
u3(pi/2,5.543654396524548,-5.543654396524548) q[10];
rzz(pi/2) q[5],q[10];
rzz(-pi/2) q[5],q[11];
u3(pi/2,5.235778316472749,-5.235778316472749) q[11];
u3(pi/2,2.0690529216542375,-2.0690529216542375) q[11];
rzz(pi/2) q[5],q[11];
rzz(pi/2) q[5],q[12];
u3(pi/2,4.93732701438172,-4.93732701438172) q[12];
u3(pi/2,1.7831679901775666,-1.7831679901775666) q[12];
rzz(pi/2) q[5],q[12];
rzz(pi/2) q[5],q[13];
rzz(-pi/2) q[5],q[13];
u3(pi/2,3.779964280799239,-3.779964280799239) q[6];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[6];
rzz(pi/2) q[6],q[7];
u3(pi/2,0.024504422698000385,-0.024504422698000385) q[7];
u3(pi/2,2.3806989128903453,-2.3806989128903453) q[7];
rzz(pi/2) q[6],q[7];
rzz(pi/2) q[6],q[8];
u3(pi/2,3.4482120965801566,-3.4482120965801566) q[8];
u3(pi/2,6.197105668471226,-6.197105668471226) q[8];
rzz(-pi/2) q[6],q[8];
rzz(pi/2) q[6],q[9];
u3(pi/2,0.3989822670059037,-0.3989822670059037) q[9];
u3(pi/2,3.343911220480976,-3.343911220480976) q[9];
rzz(pi/2) q[6],q[9];
rzz(pi/2) q[6],q[10];
u3(pi/2,5.543654396524548,-5.543654396524548) q[10];
u3(pi/2,2.3040440521427543,-2.3040440521427543) q[10];
rzz(pi/2) q[6],q[10];
rzz(-pi/2) q[6],q[11];
u3(pi/2,5.2106455752440315,-5.2106455752440315) q[11];
u3(pi/2,2.020044076258237,-2.020044076258237) q[11];
rzz(pi/2) q[6],q[11];
rzz(pi/2) q[6],q[12];
u3(pi/2,1.7831679901775666,-1.7831679901775666) q[12];
u3(pi/2,4.90025622106936,-4.90025622106936) q[12];
rzz(pi/2) q[6],q[12];
rzz(-pi/2) q[6],q[13];
u3(pi/2,4.846849145958332,-4.846849145958332) q[13];
u3(pi/2,1.6933184402848986,-1.6933184402848986) q[13];
rzz(pi/2) q[6],q[13];
u3(pi/2,0.8099025860954486,-0.8099025860954486) q[7];
u3(pi/2,0.8344070087934491,-0.8344070087934491) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,3.055513014881433,-3.055513014881433) q[8];
u3(pi/2,5.411707505073777,-5.411707505073777) q[8];
rzz(pi/2) q[7],q[8];
rzz(-pi/2) q[7],q[9];
u3(pi/2,0.20231856689118266,-0.20231856689118266) q[9];
u3(pi/2,2.951212138782252,-2.951212138782252) q[9];
rzz(pi/2) q[7],q[9];
rzz(pi/2) q[7],q[10];
u3(pi/2,2.3040440521427543,-2.3040440521427543) q[10];
u3(pi/2,5.2489730056178265,-5.2489730056178265) q[10];
rzz(pi/2) q[7],q[10];
rzz(-pi/2) q[7],q[11];
u3(pi/2,5.161636729848031,-5.161636729848031) q[11];
u3(pi/2,1.9220263854662354,-1.9220263854662354) q[11];
rzz(pi/2) q[7],q[11];
rzz(pi/2) q[7],q[12];
u3(pi/2,4.90025622106936,-4.90025622106936) q[12];
u3(pi/2,1.7096547220835654,-1.7096547220835654) q[12];
rzz(pi/2) q[7],q[12];
rzz(pi/2) q[7],q[13];
u3(pi/2,1.6933184402848986,-1.6933184402848986) q[13];
u3(pi/2,4.810406671176691,-4.810406671176691) q[13];
rzz(pi/2) q[7],q[13];
u3(pi/2,3.8409111782788807,-3.8409111782788807) q[8];
u3(pi/2,3.85347754889324,-3.85347754889324) q[8];
rzz(pi/2) q[8],q[9];
u3(pi/2,2.951212138782252,-2.951212138782252) q[9];
u3(pi/2,5.307406628974596,-5.307406628974596) q[9];
rzz(pi/2) q[8],q[9];
rzz(pi/2) q[8],q[10];
u3(pi/2,5.2489730056178265,-5.2489730056178265) q[10];
u3(pi/2,1.7146812703293088,-1.7146812703293088) q[10];
rzz(-pi/2) q[8],q[10];
rzz(pi/2) q[8],q[11];
u3(pi/2,1.9220263854662354,-1.9220263854662354) q[11];
u3(pi/2,4.867583657472026,-4.867583657472026) q[11];
rzz(pi/2) q[8],q[11];
rzz(pi/2) q[8],q[12];
u3(pi/2,1.7096547220835654,-1.7096547220835654) q[12];
u3(pi/2,4.753229684881357,-4.753229684881357) q[12];
rzz(-pi/2) q[8],q[12];
rzz(pi/2) q[8],q[13];
u3(pi/2,4.810406671176691,-4.810406671176691) q[13];
u3(pi/2,1.6198051721908973,-1.6198051721908973) q[13];
rzz(pi/2) q[8],q[13];
u3(pi,3.7366103021797,-3.7366103021797) q[9];
rzz(pi/2) q[9],q[10];
u3(pi/2,4.856273923919102,-4.856273923919102) q[10];
u3(pi/2,0.9292831069318608,-0.9292831069318608) q[10];
rzz(pi/2) q[9],q[10];
rzz(-pi/2) q[9],q[11];
u3(pi/2,1.7259910038822324,-1.7259910038822324) q[11];
u3(pi/2,4.474884575773301,-4.474884575773301) q[11];
rzz(pi/2) q[9],q[11];
rzz(pi/2) q[9],q[12];
u3(pi/2,1.611637031291564,-1.611637031291564) q[12];
u3(pi/2,4.556565984766635,-4.556565984766635) q[12];
rzz(pi/2) q[9],q[12];
rzz(-pi/2) q[9],q[13];
u3(pi/2,4.761397825780691,-4.761397825780691) q[13];
u3(pi/2,1.521159162868178,-1.521159162868178) q[13];
rzz(pi/2) q[9],q[13];
u3(pi,5.64167208731655,-5.64167208731655) q[10];
rzz(pi/2) q[10],q[11];
u3(pi/2,4.474884575773301,-4.474884575773301) q[11];
u3(pi/2,0.54789375878606,-0.54789375878606) q[11];
rzz(pi/2) q[10],q[11];
rzz(-pi/2) q[10],q[12];
u3(pi/2,1.4149733311768429,-1.4149733311768429) q[12];
u3(pi/2,4.163866903067912,-4.163866903067912) q[12];
rzz(pi/2) q[10],q[12];
rzz(pi/2) q[10],q[13];
u3(pi/2,1.521159162868178,-1.521159162868178) q[13];
u3(pi/2,4.466716434873968,-4.466716434873968) q[13];
rzz(pi/2) q[10],q[13];
rzz(-pi/2) q[11],q[12];
u3(pi/2,1.0222742494781187,-1.0222742494781187) q[12];
u3(pi/2,3.3784687396704633,-3.3784687396704633) q[12];
rzz(pi/2) q[11],q[12];
rzz(pi/2) q[11],q[13];
u3(pi/2,4.466716434873968,-4.466716434873968) q[13];
u3(pi/2,0.9324246995854506,-0.9324246995854506) q[13];
rzz(pi/2) q[11],q[13];
rzz(pi/2) q[12],q[13];
u3(pi/2,0.9324246995854506,-0.9324246995854506) q[13];
u3(pi/2,3.2886191897777954,-3.2886191897777954) q[13];
rzz(-pi/2) q[12],q[13];
u3(pi,pi,-pi) q[1];
u3(pi,pi,-pi) q[2];
u3(pi,pi,-pi) q[3];
u3(pi,3*pi/2,-3*pi/2) q[4];
u3(pi,1.7668317083788996,-1.7668317083788996) q[5];
u3(pi,4.614371289592689,-4.614371289592689) q[6];
u3(pi,5.9885039162728635,-5.9885039162728635) q[7];
u3(pi,3.2886191897777954,-3.2886191897777954) q[8];
u3(pi,3.460778467194516,-3.460778467194516) q[9];
u3(pi,7*pi/4,-7*pi/4) q[10];
u3(pi,3.620371373996878,-3.620371373996878) q[11];
u3(pi,5.698849073611885,-5.698849073611885) q[12];
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
