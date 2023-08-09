OPENQASM 2.0;
include "qelib1.inc";
qreg q[24];
creg c[24];
u3(pi,1.6543626913803848,-1.6543626913803848) q[1];
u3(pi/2,1.9440175340413641,-1.9440175340413641) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.489681119607542,-3.489681119607542) q[0];
u3(pi/2,3.169866987472101,-3.169866987472101) q[1];
u3(pi/2,0.018221237390820797,-0.018221237390820797) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,4.730610217775511,-4.730610217775511) q[1];
u3(pi/2,4.740663314266998,-4.740663314266998) q[1];
u3(pi/2,0.35814156250923646,-0.35814156250923646) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi,1.487229962209408,-1.487229962209408) q[11];
u3(pi/2,1.5902742012471531,-1.5902742012471531) q[10];
rzz(pi/2) q[11],q[10];
u3(pi/2,3.1196015050146646,-3.1196015050146646) q[10];
u3(pi/2,3.0856723043558945,-3.0856723043558945) q[11];
u3(pi/2,6.217211861454201,-6.217211861454201) q[11];
rzz(-pi/2) q[10],q[11];
u3(pi/2,4.646415534659305,-4.646415534659305) q[11];
u3(pi/2,4.6564686311507915,-4.6564686311507915) q[11];
u3(pi/2,6.271875573626663,-6.271875573626663) q[10];
rzz(pi/2) q[11],q[10];
u3(pi/2,1.487229962209408,-1.487229962209408) q[13];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[13];
u3(pi/2,1.6286016316209486,-1.6286016316209486) q[12];
rzz(-pi/2) q[13],q[12];
u3(pi/2,0.03015928947446201,-0.03015928947446201) q[12];
u3(pi/2,6.256795928889432,-6.256795928889432) q[13];
u3(pi/2,3.1051501788081515,-3.1051501788081515) q[13];
rzz(-pi/2) q[12],q[13];
u3(pi/2,1.534353852013255,-1.534353852013255) q[13];
u3(pi/2,1.5444069485047422,-1.5444069485047422) q[13];
u3(pi/2,3.18243335808646,-3.18243335808646) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,1.611637031291564,-1.611637031291564) q[12];
u3(pi/2,4.6564686311507915,-4.6564686311507915) q[11];
rzz(pi/2) q[12],q[11];
u3(pi/2,6.227264957945688,-6.227264957945688) q[11];
u3(pi/2,1.611637031291564,-1.611637031291564) q[12];
u3(pi/2,4.742548269859152,-4.742548269859152) q[12];
rzz(pi/2) q[11],q[12];
u3(pi/2,3.171751943064255,-3.171751943064255) q[12];
u3(pi/2,3.18243335808646,-3.18243335808646) q[12];
u3(pi/2,3.0756192078644076,-3.0756192078644076) q[11];
rzz(-pi/2) q[12],q[11];
u3(pi/2,4.628822615799201,-4.628822615799201) q[15];
u3(pi/2,4.6564686311507915,-4.6564686311507915) q[15];
u3(pi/2,1.5349821705439728,-1.5349821705439728) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,3.074362570802972,-3.074362570802972) q[14];
u3(pi/2,3.052999740758561,-3.052999740758561) q[15];
u3(pi/2,6.183910979326148,-6.183910979326148) q[15];
rzz(-pi/2) q[14],q[15];
u3(pi/2,4.613114652531252,-4.613114652531252) q[15];
u3(pi/2,4.6237960675534575,-4.6237960675534575) q[15];
u3(pi/2,6.226008320884252,-6.226008320884252) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,4.655211994089355,-4.655211994089355) q[14];
u3(pi/2,1.5444069485047422,-1.5444069485047422) q[13];
rzz(-pi/2) q[14],q[13];
u3(pi/2,6.256795928889432,-6.256795928889432) q[13];
u3(pi/2,1.5136193404995624,-1.5136193404995624) q[14];
u3(pi/2,4.6451588975978675,-4.6451588975978675) q[14];
rzz(pi/2) q[13],q[14];
u3(pi/2,3.074362570802972,-3.074362570802972) q[14];
u3(pi/2,3.0844156672944587,-3.0844156672944587) q[14];
u3(pi/2,3.1051501788081515,-3.1051501788081515) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,1.487229962209408,-1.487229962209408) q[17];
u3(pi/2,1.5418936743818705,-1.5418936743818705) q[17];
u3(pi/2,1.534353852013255,-1.534353852013255) q[16];
rzz(pi/2) q[17],q[16];
u3(pi/2,3.168610350410665,-3.168610350410665) q[16];
u3(pi/2,3.0567696519428686,-3.0567696519428686) q[17];
u3(pi/2,6.188309209041175,-6.188309209041175) q[17];
rzz(-pi/2) q[16],q[17];
u3(pi/2,4.617512882246278,-4.617512882246278) q[17];
u3(pi/2,4.627565978737765,-4.627565978737765) q[17];
u3(pi/2,0.03769911184307752,-0.03769911184307752) q[16];
rzz(-pi/2) q[17],q[16];
u3(pi/2,1.6084954386379742,-1.6084954386379742) q[16];
u3(pi/2,4.6237960675534575,-4.6237960675534575) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,6.194592394348354,-6.194592394348354) q[15];
u3(pi/2,1.6084954386379742,-1.6084954386379742) q[16];
u3(pi/2,4.739406677205562,-4.739406677205562) q[16];
rzz(pi/2) q[15],q[16];
u3(pi/2,3.168610350410665,-3.168610350410665) q[16];
u3(pi/2,3.1792917654328705,-3.1792917654328705) q[16];
u3(pi/2,3.042318325736356,-3.042318325736356) q[15];
rzz(pi/2) q[16],q[15];
u3(pi,3*pi/2,-3*pi/2) q[19];
u3(pi/2,1.5286989852367932,-1.5286989852367932) q[18];
rzz(-pi/2) q[19],q[18];
u3(pi/2,6.250512743582252,-6.250512743582252) q[18];
u3(pi/2,0.028274333882308135,-0.028274333882308135) q[19];
u3(pi/2,3.1591855724498963,-3.1591855724498963) q[19];
rzz(-pi/2) q[18],q[19];
u3(pi/2,1.5883892456549995,-1.5883892456549995) q[19];
u3(pi/2,1.5990706606772047,-1.5990706606772047) q[19];
u3(pi/2,3.1196015050146646,-3.1196015050146646) q[18];
rzz(pi/2) q[19],q[18];
u3(pi/2,1.548805178219768,-1.548805178219768) q[18];
u3(pi/2,1.485973325147972,-1.485973325147972) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,3.0567696519428686,-3.0567696519428686) q[17];
u3(pi/2,1.548805178219768,-1.548805178219768) q[18];
u3(pi/2,4.679716416787356,-4.679716416787356) q[18];
rzz(-pi/2) q[17],q[18];
u3(pi/2,6.250512743582252,-6.250512743582252) q[18];
u3(pi/2,6.261194158604458,-6.261194158604458) q[18];
u3(pi/2,3.0467165554513813,-3.0467165554513813) q[17];
rzz(pi/2) q[18],q[17];
u3(pi,pi/2,-pi/2) q[21];
u3(pi/2,1.6166635795373074,-1.6166635795373074) q[20];
rzz(pi/2) q[21],q[20];
u3(pi/2,3.206309462253743,-3.206309462253743) q[20];
u3(pi/2,3.0988669935009723,-3.0988669935009723) q[21];
u3(pi/2,6.22977823206856,-6.22977823206856) q[21];
rzz(-pi/2) q[20],q[21];
u3(pi/2,4.658981905273664,-4.658981905273664) q[21];
u3(pi/2,4.669663320295868,-4.669663320295868) q[21];
u3(pi/2,0.07539822368615504,-0.07539822368615504) q[20];
rzz(-pi/2) q[21],q[20];
u3(pi/2,1.6461945504810516,-1.6461945504810516) q[20];
u3(pi/2,1.5990706606772047,-1.5990706606772047) q[19];
rzz(pi/2) q[20],q[19];
u3(pi/2,3.169866987472101,-3.169866987472101) q[19];
u3(pi/2,1.6461945504810516,-1.6461945504810516) q[20];
u3(pi/2,4.7771057890486395,-4.7771057890486395) q[20];
rzz(pi/2) q[19],q[20];
u3(pi/2,3.206309462253743,-3.206309462253743) q[20];
u3(pi/2,3.2169908772759483,-3.2169908772759483) q[20];
u3(pi/2,0.017592918860102842,-0.017592918860102842) q[19];
rzz(pi/2) q[20],q[19];
u3(pi/2,4.031291693086422,-4.031291693086422) q[23];
u3(pi/2,4.079043901420987,-4.079043901420987) q[23];
u3(pi/2,1.5129910219688443,-1.5129910219688443) q[22];
rzz(pi/2) q[23],q[22];
u3(pi/2,3.151645750081281,-3.151645750081281) q[22];
u3(pi/2,2.4787166036823467,-2.4787166036823467) q[23];
u3(pi/2,5.609627842249935,-5.609627842249935) q[23];
rzz(-pi/2) q[22],q[23];
u3(pi/2,4.038831515455039,-4.038831515455039) q[23];
u3(pi/2,4.049512930477243,-4.049512930477243) q[23];
u3(pi/2,0.020734511513692634,-0.020734511513692634) q[22];
rzz(pi/2) q[23],q[22];
u3(pi/2,4.733123491898382,-4.733123491898382) q[22];
u3(pi/2,1.5280706667060753,-1.5280706667060753) q[21];
rzz(pi/2) q[22],q[21];
u3(pi/2,3.0988669935009723,-3.0988669935009723) q[21];
u3(pi/2,4.733123491898382,-4.733123491898382) q[22];
u3(pi/2,1.580849423286384,-1.580849423286384) q[22];
rzz(pi/2) q[21],q[22];
u3(pi/2,0.010053096491487338,-0.010053096491487338) q[22];
u3(pi/2,0.020734511513692634,-0.020734511513692634) q[22];
u3(pi/2,6.22977823206856,-6.22977823206856) q[21];
rzz(-pi/2) q[22],q[21];
u3(pi/2,1.6543626913803848,-1.6543626913803848) q[3];
u3(pi/2,1.6933184402848986,-1.6933184402848986) q[3];
u3(pi/2,1.6015839348000767,-1.6015839348000767) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,3.206309462253743,-3.206309462253743) q[2];
u3(pi/2,3.321291753375129,-3.321291753375129) q[3];
u3(pi/2,0.16901768476313087,-0.16901768476313087) q[3];
rzz(-pi/2) q[2],q[3];
u3(pi/2,4.8814066651478205,-4.8814066651478205) q[3];
u3(pi/2,4.8920880801700255,-4.8920880801700255) q[3];
u3(pi/2,0.07539822368615504,-0.07539822368615504) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,4.7877872040708445,-4.7877872040708445) q[2];
u3(pi/2,1.5990706606772047,-1.5990706606772047) q[1];
rzz(-pi/2) q[2],q[1];
u3(pi/2,0.028274333882308135,-0.028274333882308135) q[1];
u3(pi/2,1.6461945504810516,-1.6461945504810516) q[2];
u3(pi/2,4.7771057890486395,-4.7771057890486395) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,3.206309462253743,-3.206309462253743) q[2];
u3(pi/2,3.2169908772759483,-3.2169908772759483) q[2];
u3(pi/2,3.159813890980614,-3.159813890980614) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,4.628822615799201,-4.628822615799201) q[5];
u3(pi/2,4.670291638826586,-4.670291638826586) q[5];
u3(pi/2,1.5739379194484864,-1.5739379194484864) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,3.1931147731086655,-3.1931147731086655) q[4];
u3(pi/2,3.1051501788081515,-3.1051501788081515) q[5];
u3(pi/2,6.236689735906458,-6.236689735906458) q[5];
rzz(-pi/2) q[4],q[5];
u3(pi/2,4.665893409111561,-4.665893409111561) q[5];
u3(pi/2,4.675946505603048,-4.675946505603048) q[5];
u3(pi/2,0.061575216010359944,-0.061575216010359944) q[4];
rzz(-pi/2) q[5],q[4];
u3(pi/2,1.6323715428052563,-1.6323715428052563) q[4];
u3(pi/2,4.8920880801700255,-4.8920880801700255) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,0.17969909978533616,-0.17969909978533616) q[3];
u3(pi/2,1.6323715428052563,-1.6323715428052563) q[4];
u3(pi/2,4.763911099903562,-4.763911099903562) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.1931147731086655,-3.1931147731086655) q[4];
u3(pi/2,3.2031678696001533,-3.2031678696001533) q[4];
u3(pi/2,3.3106103383529244,-3.3106103383529244) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,1.5293273037675112,-1.5293273037675112) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,6.209672039085585,-6.209672039085585) q[6];
u3(pi/2,3.042318325736356,-3.042318325736356) q[7];
u3(pi/2,3.052371422227843,-3.052371422227843) q[7];
rzz(-pi/2) q[6],q[7];
u3(pi/2,4.623167749022739,-4.623167749022739) q[7];
u3(pi/2,4.633849164044945,-4.633849164044945) q[7];
u3(pi/2,3.0781324819872795,-3.0781324819872795) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,1.5073361551923827,-1.5073361551923827) q[6];
u3(pi/2,1.534353852013255,-1.534353852013255) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,3.1051501788081515,-3.1051501788081515) q[5];
u3(pi/2,1.5073361551923827,-1.5073361551923827) q[6];
u3(pi/2,4.638875712290688,-4.638875712290688) q[6];
rzz(-pi/2) q[5],q[6];
u3(pi/2,6.209672039085585,-6.209672039085585) q[6];
u3(pi/2,6.219725135577073,-6.219725135577073) q[6];
u3(pi/2,3.095097082316664,-3.095097082316664) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,4.628822615799201,-4.628822615799201) q[9];
u3(pi/2,4.669035001765151,-4.669035001765151) q[9];
u3(pi/2,1.5538317264655117,-1.5538317264655117) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,3.1573006168577415,-3.1573006168577415) q[8];
u3(pi/2,3.156043979796306,-3.156043979796306) q[9];
u3(pi/2,0.0043982297150257105,-0.0043982297150257105) q[9];
rzz(-pi/2) q[8],q[9];
u3(pi/2,4.7167872100997155,-4.7167872100997155) q[9];
u3(pi/2,4.726840306591202,-4.726840306591202) q[9];
u3(pi/2,0.02638937829015426,-0.02638937829015426) q[8];
rzz(-pi/2) q[9],q[8];
u3(pi/2,4.701079246831767,-4.701079246831767) q[10];
u3(pi/2,1.5852476530014097,-1.5852476530014097) q[9];
rzz(pi/2) q[10],q[9];
u3(pi/2,1.5971857050850506,-1.5971857050850506) q[8];
u3(pi/2,4.633849164044945,-4.633849164044945) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,6.204645490839842,-6.204645490839842) q[7];
u3(pi/2,1.5971857050850506,-1.5971857050850506) q[8];
u3(pi/2,4.7280969436526386,-4.7280969436526386) q[8];
rzz(pi/2) q[7],q[8];
u3(pi/2,3.1573006168577415,-3.1573006168577415) q[8];
u3(pi/2,3.1679820318799474,-3.1679820318799474) q[8];
u3(pi/2,3.052371422227843,-3.052371422227843) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,3.156043979796306,-3.156043979796306) q[9];
u3(pi/2,4.701079246831767,-4.701079246831767) q[10];
u3(pi/2,1.548805178219768,-1.548805178219768) q[10];
rzz(-pi/2) q[9],q[10];
u3(pi/2,3.1196015050146646,-3.1196015050146646) q[10];
u3(pi/2,3.1302829200368696,-3.1302829200368696) q[10];
u3(pi/2,3.145990883304819,-3.145990883304819) q[9];
rzz(pi/2) q[10],q[9];
u3(pi/2,4.730610217775511,-4.730610217775511) q[1];
u3(pi/2,1.9295662078348508,-1.9295662078348508) q[0];
u3(pi/2,3.5091589940597987,-3.5091589940597987) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.913229926036184,-1.913229926036184) q[0];
u3(pi/2,4.67531818707233,-4.67531818707233) q[1];
u3(pi/2,1.5230441184603318,-1.5230441184603318) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,6.235433098845021,-6.235433098845021) q[1];
u3(pi/2,3.083159030233023,-3.083159030233023) q[1];
u3(pi/2,1.9025485110139788,-1.9025485110139788) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,1.5048228810695108,-1.5048228810695108) q[11];
u3(pi/2,3.1101767270538954,-3.1101767270538954) q[10];
rzz(pi/2) q[11],q[10];
u3(pi/2,1.4979113772316133,-1.4979113772316133) q[10];
u3(pi/2,1.5224157999296137,-1.5224157999296137) q[11];
u3(pi/2,4.653327038497202,-4.653327038497202) q[11];
rzz(-pi/2) q[10],q[11];
u3(pi/2,3.082530711702305,-3.082530711702305) q[11];
u3(pi/2,3.0932121267245103,-3.0932121267245103) q[11];
u3(pi/2,4.649557127312894,-4.649557127312894) q[10];
rzz(pi/2) q[11],q[10];
u3(pi/2,3.1051501788081515,-3.1051501788081515) q[13];
u3(pi/2,4.640760667882843,-4.640760667882843) q[13];
u3(pi/2,6.266220706850201,-6.266220706850201) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,4.668406683234433,-4.668406683234433) q[12];
u3(pi/2,3.0925838081937926,-3.0925838081937926) q[13];
u3(pi/2,6.224123365292098,-6.224123365292098) q[13];
rzz(pi/2) q[12],q[13];
u3(pi/2,1.5117343849074085,-1.5117343849074085) q[13];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[13];
u3(pi/2,4.679088098256638,-4.679088098256638) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,3.108291771461741,-3.108291771461741) q[12];
u3(pi/2,3.0932121267245103,-3.0932121267245103) q[11];
rzz(pi/2) q[12],q[11];
u3(pi/2,4.664008453519407,-4.664008453519407) q[11];
u3(pi/2,3.108291771461741,-3.108291771461741) q[12];
u3(pi/2,6.239203010029329,-6.239203010029329) q[12];
rzz(pi/2) q[11],q[12];
u3(pi/2,4.668406683234433,-4.668406683234433) q[12];
u3(pi/2,4.679088098256638,-4.679088098256638) q[12];
u3(pi/2,1.5117343849074085,-1.5117343849074085) q[11];
rzz(-pi/2) q[12],q[11];
u3(pi/2,3.042318325736356,-3.042318325736356) q[15];
u3(pi/2,4.640760667882843,-4.640760667882843) q[15];
u3(pi/2,3.1202298235453823,-3.1202298235453823) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,1.51738925168387,-1.51738925168387) q[14];
u3(pi/2,3.037291777490612,-3.037291777490612) q[15];
u3(pi/2,6.1682030160582,-6.1682030160582) q[15];
rzz(-pi/2) q[14],q[15];
u3(pi/2,4.597406689263304,-4.597406689263304) q[15];
u3(pi/2,4.608088104285509,-4.608088104285509) q[15];
u3(pi/2,4.669663320295868,-4.669663320295868) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,3.0988669935009723,-3.0988669935009723) q[14];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,3.0925838081937926,-3.0925838081937926) q[13];
u3(pi/2,3.0988669935009723,-3.0988669935009723) q[14];
u3(pi/2,6.22977823206856,-6.22977823206856) q[14];
rzz(pi/2) q[13],q[14];
u3(pi/2,4.658981905273664,-4.658981905273664) q[14];
u3(pi/2,4.669663320295868,-4.669663320295868) q[14];
u3(pi/2,6.224123365292098,-6.224123365292098) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,3.0467165554513813,-3.0467165554513813) q[17];
u3(pi/2,4.562220851543097,-4.562220851543097) q[17];
u3(pi/2,3.2157342402145126,-3.2157342402145126) q[16];
rzz(pi/2) q[17],q[16];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[16];
u3(pi/2,2.9355041755143025,-2.9355041755143025) q[17];
u3(pi/2,6.066415414081891,-6.066415414081891) q[17];
rzz(pi/2) q[16],q[17];
u3(pi/2,1.3540264336972008,-1.3540264336972008) q[17];
u3(pi/2,1.3647078487194062,-1.3647078487194062) q[17];
u3(pi/2,1.7190795000443349,-1.7190795000443349) q[16];
rzz(pi/2) q[17],q[16];
u3(pi/2,0.14828317324943824,-0.14828317324943824) q[16];
u3(pi/2,4.608088104285509,-4.608088104285509) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,6.178884431080405,-6.178884431080405) q[15];
u3(pi/2,0.14828317324943824,-0.14828317324943824) q[16];
u3(pi/2,3.2798227303477443,-3.2798227303477443) q[16];
rzz(pi/2) q[15],q[16];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[16];
u3(pi/2,1.7190795000443349,-1.7190795000443349) q[16];
u3(pi/2,3.0266103624684066,-3.0266103624684066) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,1.5883892456549995,-1.5883892456549995) q[19];
u3(pi/2,0.020106192982974676,-0.020106192982974676) q[18];
rzz(pi/2) q[19],q[18];
u3(pi/2,4.741919951328434,-4.741919951328434) q[18];
u3(pi/2,1.6166635795373074,-1.6166635795373074) q[19];
u3(pi/2,4.747574818104896,-4.747574818104896) q[19];
rzz(-pi/2) q[18],q[19];
u3(pi/2,3.176778491309999,-3.176778491309999) q[19];
u3(pi/2,3.187459906332204,-3.187459906332204) q[19];
u3(pi/2,1.6103803942301278,-1.6103803942301278) q[18];
rzz(pi/2) q[19],q[18];
u3(pi/2,0.03958406743523139,-0.03958406743523139) q[18];
u3(pi/2,1.3647078487194062,-1.3647078487194062) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,2.9355041755143025,-2.9355041755143025) q[17];
u3(pi/2,0.03958406743523139,-0.03958406743523139) q[18];
u3(pi/2,3.1711236245335375,-3.1711236245335375) q[18];
rzz(-pi/2) q[17],q[18];
u3(pi/2,4.741919951328434,-4.741919951328434) q[18];
u3(pi/2,4.751973047819921,-4.751973047819921) q[18];
u3(pi/2,2.9248227604920976,-2.9248227604920976) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,4.658981905273664,-4.658981905273664) q[21];
u3(pi/2,3.1503891130198443,-3.1503891130198443) q[20];
rzz(pi/2) q[21],q[20];
u3(pi/2,1.5984423421464868,-1.5984423421464868) q[20];
u3(pi/2,4.606831467224072,-4.606831467224072) q[21];
u3(pi/2,1.4545573986120743,-1.4545573986120743) q[21];
rzz(-pi/2) q[20],q[21];
u3(pi/2,6.166946378996764,-6.166946378996764) q[21];
u3(pi/2,6.177627794018969,-6.177627794018969) q[21];
u3(pi/2,4.750716410758485,-4.750716410758485) q[20];
rzz(-pi/2) q[21],q[20];
u3(pi/2,0.038327430373795476,-0.038327430373795476) q[20];
u3(pi/2,3.187459906332204,-3.187459906332204) q[19];
rzz(pi/2) q[20],q[19];
u3(pi/2,4.7582562331271,-4.7582562331271) q[19];
u3(pi/2,0.038327430373795476,-0.038327430373795476) q[20];
u3(pi/2,3.169238668941383,-3.169238668941383) q[20];
rzz(pi/2) q[19],q[20];
u3(pi/2,1.5984423421464868,-1.5984423421464868) q[20];
u3(pi/2,1.609123757168692,-1.609123757168692) q[20];
u3(pi/2,1.6059821645151022,-1.6059821645151022) q[19];
rzz(pi/2) q[20],q[19];
u3(pi/2,4.049512930477243,-4.049512930477243) q[23];
u3(pi/2,0.8601680685528853,-0.8601680685528853) q[23];
u3(pi/2,3.199397958415845,-3.199397958415845) q[22];
rzz(pi/2) q[23],q[22];
u3(pi/2,1.6970883514692063,-1.6970883514692063) q[22];
u3(pi/2,5.602088019881319,-5.602088019881319) q[23];
u3(pi/2,5.612769434903524,-5.612769434903524) q[23];
rzz(pi/2) q[22],q[23];
u3(pi/2,0.9003804545188347,-0.9003804545188347) q[23];
u3(pi/2,0.9104335510103221,-0.9104335510103221) q[23];
u3(pi/2,4.828627908567512,-4.828627908567512) q[22];
rzz(pi/2) q[23],q[22];
u3(pi/2,3.2578315817726153,-3.2578315817726153) q[22];
u3(pi/2,3.0360351404291763,-3.0360351404291763) q[21];
rzz(pi/2) q[22],q[21];
u3(pi/2,4.606831467224072,-4.606831467224072) q[21];
u3(pi/2,3.2578315817726153,-3.2578315817726153) q[22];
u3(pi/2,3.2678846782641027,-3.2678846782641027) q[22];
rzz(pi/2) q[21],q[22];
u3(pi/2,1.6970883514692063,-1.6970883514692063) q[22];
u3(pi/2,1.7077697664914115,-1.7077697664914115) q[22];
u3(pi/2,4.617512882246278,-4.617512882246278) q[21];
rzz(-pi/2) q[22],q[21];
u3(pi/2,3.3106103383529244,-3.3106103383529244) q[3];
u3(pi/2,4.8430792347740255,-4.8430792347740255) q[3];
u3(pi/2,3.164840439226358,-3.164840439226358) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,1.6279733130902307,-1.6279733130902307) q[2];
u3(pi/2,3.329459894274463,-3.329459894274463) q[3];
u3(pi/2,0.17718582566246432,-0.17718582566246432) q[3];
rzz(-pi/2) q[2],q[3];
u3(pi/2,4.889574806047154,-4.889574806047154) q[3];
u3(pi/2,4.90025622106936,-4.90025622106936) q[3];
u3(pi/2,4.780247381702229,-4.780247381702229) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,3.209451054907333,-3.209451054907333) q[2];
u3(pi/2,6.2247516838228165,-6.2247516838228165) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,1.5123627034381264,-1.5123627034381264) q[1];
u3(pi/2,3.209451054907333,-3.209451054907333) q[2];
u3(pi/2,0.057176986295334235,-0.057176986295334235) q[2];
rzz(-pi/2) q[1],q[2];
u3(pi/2,1.6279733130902307,-1.6279733130902307) q[2];
u3(pi/2,4.759512870188536,-4.759512870188536) q[2];
u3(pi/2,4.664636772050124,-4.664636772050124) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3.095097082316664,-3.095097082316664) q[5];
u3(pi/2,4.707362432138946,-4.707362432138946) q[5];
u3(pi/2,3.1792917654328705,-3.1792917654328705) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,1.6568759655032568,-1.6568759655032568) q[4];
u3(pi/2,3.142220972120511,-3.142220972120511) q[5];
u3(pi/2,6.273760529218817,-6.273760529218817) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,1.561371548834127,-1.561371548834127) q[5];
u3(pi/2,1.5714246453256144,-1.5714246453256144) q[5];
u3(pi/2,1.6675573805254624,-1.6675573805254624) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,0.09676105373056564,-0.09676105373056564) q[4];
u3(pi/2,4.90025622106936,-4.90025622106936) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,0.18786724068466962,-0.18786724068466962) q[3];
u3(pi/2,0.09676105373056564,-0.09676105373056564) q[4];
u3(pi/2,3.2276722922981538,-3.2276722922981538) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,1.6568759655032568,-1.6568759655032568) q[4];
u3(pi/2,1.6675573805254624,-1.6675573805254624) q[4];
u3(pi/2,3.3187784792522574,-3.3187784792522574) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,4.623167749022739,-4.623167749022739) q[7];
u3(pi/2,6.261822477135176,-6.261822477135176) q[6];
rzz(-pi/2) q[7],q[6];
u3(pi/2,1.51738925168387,-1.51738925168387) q[6];
u3(pi/2,1.4972830587008954,-1.4972830587008954) q[7];
u3(pi/2,4.628194297268483,-4.628194297268483) q[7];
rzz(-pi/2) q[6],q[7];
u3(pi/2,3.0573979704735863,-3.0573979704735863) q[7];
u3(pi/2,3.068079385495792,-3.068079385495792) q[7];
u3(pi/2,4.669035001765151,-4.669035001765151) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,3.0982386749702537,-3.0982386749702537) q[6];
u3(pi/2,1.5714246453256144,-1.5714246453256144) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,3.142220972120511,-3.142220972120511) q[5];
u3(pi/2,3.0982386749702537,-3.0982386749702537) q[6];
u3(pi/2,6.22977823206856,-6.22977823206856) q[6];
rzz(-pi/2) q[5],q[6];
u3(pi/2,1.51738925168387,-1.51738925168387) q[6];
u3(pi/2,1.5274423481753574,-1.5274423481753574) q[6];
u3(pi/2,3.132167875629024,-3.132167875629024) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,3.145990883304819,-3.145990883304819) q[9];
u3(pi/2,4.757627914596383,-4.757627914596383) q[9];
u3(pi/2,3.1843183136786144,-3.1843183136786144) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,1.6461945504810516,-1.6461945504810516) q[8];
u3(pi/2,3.244636892627538,-3.244636892627538) q[9];
u3(pi/2,0.09236282401553991,-0.09236282401553991) q[9];
rzz(-pi/2) q[8],q[9];
u3(pi/2,4.8047518044002295,-4.8047518044002295) q[9];
u3(pi/2,4.8154332194224345,-4.8154332194224345) q[9];
u3(pi/2,4.797840300562331,-4.797840300562331) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,3.078760800517997,-3.078760800517997) q[10];
u3(pi/2,4.8154332194224345,-4.8154332194224345) q[9];
rzz(-pi/2) q[10],q[9];
u3(pi/2,3.2270439737674352,-3.2270439737674352) q[8];
u3(pi/2,3.068079385495792,-3.068079385495792) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,4.638875712290688,-4.638875712290688) q[7];
u3(pi/2,3.2270439737674352,-3.2270439737674352) q[8];
u3(pi/2,0.07539822368615504,-0.07539822368615504) q[8];
rzz(pi/2) q[7],q[8];
u3(pi/2,4.7877872040708445,-4.7877872040708445) q[8];
u3(pi/2,4.797840300562331,-4.797840300562331) q[8];
u3(pi/2,1.4866016436786902,-1.4866016436786902) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,3.244636892627538,-3.244636892627538) q[9];
u3(pi/2,6.220353454107791,-6.220353454107791) q[10];
u3(pi/2,3.06870770402651,-3.06870770402651) q[10];
rzz(-pi/2) q[9],q[10];
u3(pi/2,4.639504030821406,-4.639504030821406) q[10];
u3(pi/2,4.649557127312894,-4.649557127312894) q[10];
u3(pi/2,3.2339554776053334,-3.2339554776053334) q[9];
rzz(pi/2) q[10],q[9];
u3(pi/2,6.235433098845021,-6.235433098845021) q[1];
u3(pi/2,3.473344837808875,-3.473344837808875) q[0];
u3(pi/2,5.053565942564541,-5.053565942564541) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.4570085560102086,-3.4570085560102086) q[0];
u3(pi/2,6.180141068141841,-6.180141068141841) q[1];
u3(pi/2,3.0278669995298424,-3.0278669995298424) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,1.457070672734946,-1.457070672734946) q[1];
u3(pi/2,4.588610229833251,-4.588610229833251) q[1];
u3(pi/2,3.446327140988003,-3.446327140988003) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,6.224123365292098,-6.224123365292098) q[11];
u3(pi/2,4.630079252860637,-4.630079252860637) q[10];
rzz(pi/2) q[11],q[10];
u3(pi/2,3.0171855845076374,-3.0171855845076374) q[10];
u3(pi/2,6.242344602682919,-6.242344602682919) q[11];
u3(pi/2,3.0900705340709207,-3.0900705340709207) q[11];
rzz(-pi/2) q[10],q[11];
u3(pi/2,1.519274207276024,-1.519274207276024) q[11];
u3(pi/2,1.5299556222982291,-1.5299556222982291) q[11];
u3(pi/2,6.1694596531196355,-6.1694596531196355) q[10];
rzz(pi/2) q[11],q[10];
u3(pi/2,6.224123365292098,-6.224123365292098) q[13];
u3(pi/2,1.4765485471872026,-1.4765485471872026) q[13];
u3(pi/2,1.4796901398407925,-1.4796901398407925) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,6.16506142340461,-6.16506142340461) q[12];
u3(pi/2,6.211556994677739,-6.211556994677739) q[13];
u3(pi/2,3.0599112445964582,-3.0599112445964582) q[13];
rzz(pi/2) q[12],q[13];
u3(pi/2,4.630707571391355,-4.630707571391355) q[13];
u3(pi/2,4.640760667882843,-4.640760667882843) q[13];
u3(pi/2,6.1751145198960975,-6.1751145198960975) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,4.604318193101201,-4.604318193101201) q[12];
u3(pi/2,1.5299556222982291,-1.5299556222982291) q[11];
rzz(pi/2) q[12],q[11];
u3(pi/2,3.1007519490931257,-3.1007519490931257) q[11];
u3(pi/2,4.604318193101201,-4.604318193101201) q[12];
u3(pi/2,1.4526724430199203,-1.4526724430199203) q[12];
rzz(pi/2) q[11],q[12];
u3(pi/2,6.16506142340461,-6.16506142340461) q[12];
u3(pi/2,6.1751145198960975,-6.1751145198960975) q[12];
u3(pi/2,6.231663187660714,-6.231663187660714) q[11];
rzz(pi/2) q[12],q[11];
u3(pi/2,3.0266103624684066,-3.0266103624684066) q[15];
u3(pi/2,4.625052704614894,-4.625052704614894) q[15];
u3(pi/2,4.705477476546792,-4.705477476546792) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,3.10263690468528,-3.10263690468528) q[14];
u3(pi/2,3.021583814222663,-3.021583814222663) q[15];
u3(pi/2,6.1531233713209685,-6.1531233713209685) q[15];
rzz(-pi/2) q[14],q[15];
u3(pi/2,4.582327044526072,-4.582327044526072) q[15];
u3(pi/2,4.59238014101756,-4.59238014101756) q[15];
u3(pi/2,6.254910973297278,-6.254910973297278) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,4.6841146465023815,-4.6841146465023815) q[14];
u3(pi/2,4.640760667882843,-4.640760667882843) q[13];
rzz(-pi/2) q[14],q[13];
u3(pi/2,3.0699643410879456,-3.0699643410879456) q[13];
u3(pi/2,1.5425219929125884,-1.5425219929125884) q[14];
u3(pi/2,4.6734332314801765,-4.6734332314801765) q[14];
rzz(pi/2) q[13],q[14];
u3(pi/2,3.10263690468528,-3.10263690468528) q[14];
u3(pi/2,3.113318319707485,-3.113318319707485) q[14];
u3(pi/2,6.201503898186251,-6.201503898186251) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,2.9248227604920976,-2.9248227604920976) q[17];
u3(pi/2,4.440955375114531,-4.440955375114531) q[17];
u3(pi/2,1.7561502933566946,-1.7561502933566946) q[16];
rzz(pi/2) q[17],q[16];
u3(pi/2,0.24881413816431164,-0.24881413816431164) q[16];
u3(pi/2,2.814238699085737,-2.814238699085737) q[17];
u3(pi/2,5.945149937653325,-5.945149937653325) q[17];
rzz(pi/2) q[16],q[17];
u3(pi/2,1.2327609572686349,-1.2327609572686349) q[17];
u3(pi/2,1.24344237229084,-1.24344237229084) q[17];
u3(pi/2,0.25949555318651696,-0.25949555318651696) q[16];
rzz(pi/2) q[17],q[16];
u3(pi/2,4.971884533571207,-4.971884533571207) q[16];
u3(pi/2,4.59238014101756,-4.59238014101756) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,6.163176467812456,-6.163176467812456) q[15];
u3(pi/2,4.971884533571207,-4.971884533571207) q[16];
u3(pi/2,1.8196104649592084,-1.8196104649592084) q[16];
rzz(pi/2) q[15],q[16];
u3(pi/2,0.24881413816431164,-0.24881413816431164) q[16];
u3(pi/2,0.25949555318651696,-0.25949555318651696) q[16];
u3(pi/2,3.011530717731176,-3.011530717731176) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,3.176778491309999,-3.176778491309999) q[19];
u3(pi/2,4.7940703893780245,-4.7940703893780245) q[18];
rzz(-pi/2) q[19],q[18];
u3(pi/2,0.09173450548482195,-0.09173450548482195) q[18];
u3(pi/2,0.06346017160251381,-0.06346017160251381) q[19];
u3(pi/2,3.1943714101701013,-3.1943714101701013) q[19];
rzz(-pi/2) q[18],q[19];
u3(pi/2,1.6235750833752052,-1.6235750833752052) q[19];
u3(pi/2,1.6342564983974104,-1.6342564983974104) q[19];
u3(pi/2,3.243380255566102,-3.243380255566102) q[18];
rzz(pi/2) q[19],q[18];
u3(pi/2,1.6725839287712059,-1.6725839287712059) q[18];
u3(pi/2,1.24344237229084,-1.24344237229084) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,2.814238699085737,-2.814238699085737) q[17];
u3(pi/2,1.6725839287712059,-1.6725839287712059) q[18];
u3(pi/2,4.804123485869511,-4.804123485869511) q[18];
rzz(-pi/2) q[17],q[18];
u3(pi/2,0.09173450548482195,-0.09173450548482195) q[18];
u3(pi/2,0.10178760197630929,-0.10178760197630929) q[18];
u3(pi/2,2.8035572840635314,-2.8035572840635314) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,3.0467165554513813,-3.0467165554513813) q[21];
u3(pi/2,1.5425219929125884,-1.5425219929125884) q[20];
rzz(pi/2) q[21],q[20];
u3(pi/2,6.273760529218817,-6.273760529218817) q[20];
u3(pi/2,3.0988669935009723,-3.0988669935009723) q[21];
u3(pi/2,3.1095484085231773,-3.1095484085231773) q[21];
rzz(-pi/2) q[20],q[21];
u3(pi/2,1.5387520817282807,-1.5387520817282807) q[21];
u3(pi/2,1.548805178219768,-1.548805178219768) q[21];
u3(pi/2,6.26370743272733,-6.26370743272733) q[20];
rzz(pi/2) q[21],q[20];
u3(pi/2,4.692911105932433,-4.692911105932433) q[20];
u3(pi/2,1.6342564983974104,-1.6342564983974104) q[19];
rzz(-pi/2) q[20],q[19];
u3(pi/2,0.06346017160251381,-0.06346017160251381) q[19];
u3(pi/2,1.55131845234264,-1.55131845234264) q[20];
u3(pi/2,1.561371548834127,-1.561371548834127) q[20];
rzz(pi/2) q[19],q[20];
u3(pi/2,6.273760529218817,-6.273760529218817) q[20];
u3(pi/2,pi/2500,-pi/2500) q[20];
u3(pi/2,0.07351326809400116,-0.07351326809400116) q[19];
rzz(pi/2) q[20],q[19];
u3(pi/2,0.9104335510103221,-0.9104335510103221) q[23];
u3(pi/2,4.0042739962655505,-4.0042739962655505) q[23];
u3(pi/2,4.8864332133935635,-4.8864332133935635) q[22];
rzz(pi/2) q[23],q[22];
u3(pi/2,3.384123606446925,-3.384123606446925) q[22];
u3(pi/2,2.463008640414398,-2.463008640414398) q[23];
u3(pi/2,2.4736900554366033,-2.4736900554366033) q[23];
rzz(pi/2) q[22],q[23];
u3(pi/2,4.0444863822315,-4.0444863822315) q[23];
u3(pi/2,4.055167797253705,-4.055167797253705) q[23];
u3(pi/2,0.23184953783492673,-0.23184953783492673) q[22];
rzz(pi/2) q[23],q[22];
u3(pi/2,4.944238518219617,-4.944238518219617) q[22];
u3(pi/2,1.548805178219768,-1.548805178219768) q[21];
rzz(pi/2) q[22],q[21];
u3(pi/2,3.1196015050146646,-3.1196015050146646) q[21];
u3(pi/2,4.944238518219617,-4.944238518219617) q[22];
u3(pi/2,4.954919933241822,-4.954919933241822) q[22];
rzz(pi/2) q[21],q[22];
u3(pi/2,3.384123606446925,-3.384123606446925) q[22];
u3(pi/2,3.3948050214691303,-3.3948050214691303) q[22];
u3(pi/2,3.1302829200368696,-3.1302829200368696) q[21];
rzz(-pi/2) q[22],q[21];
u3(pi/2,3.3187784792522574,-3.3187784792522574) q[3];
u3(pi/2,4.851247375673359,-4.851247375673359) q[3];
u3(pi/2,4.81103498970741,-4.81103498970741) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,3.2779377747555904,-3.2779377747555904) q[2];
u3(pi/2,3.3376280351737964,-3.3376280351737964) q[3];
u3(pi/2,0.18535396656179778,-0.18535396656179778) q[3];
rzz(-pi/2) q[2],q[3];
u3(pi/2,4.8977429469464875,-4.8977429469464875) q[3];
u3(pi/2,1.7460971968652068,-1.7460971968652068) q[3];
u3(pi/2,3.267256359733385,-3.267256359733385) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,1.6964600329384885,-1.6964600329384885) q[2];
u3(pi/2,4.588610229833251,-4.588610229833251) q[1];
rzz(-pi/2) q[2],q[1];
u3(pi/2,3.0178139030383555,-3.0178139030383555) q[1];
u3(pi/2,4.838052686528282,-4.838052686528282) q[2];
u3(pi/2,1.6857786179162828,-1.6857786179162828) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,0.11498229112138643,-0.11498229112138643) q[2];
u3(pi/2,3.2465218482196927,-3.2465218482196927) q[2];
u3(pi/2,3.0278669995298424,-3.0278669995298424) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3.132167875629024,-3.132167875629024) q[5];
u3(pi/2,4.7444332254513055,-4.7444332254513055) q[5];
u3(pi/2,1.643052957827462,-1.643052957827462) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,0.12063715789784804,-0.12063715789784804) q[4];
u3(pi/2,3.1792917654328705,-3.1792917654328705) q[5];
u3(pi/2,0.02764601535159018,-0.02764601535159018) q[5];
rzz(-pi/2) q[4],q[5];
u3(pi/2,4.74003499573628,-4.74003499573628) q[5];
u3(pi/2,4.7500880922277675,-4.7500880922277675) q[5];
u3(pi/2,3.2729112265098466,-3.2729112265098466) q[4];
rzz(-pi/2) q[5],q[4];
u3(pi/2,4.843707553304744,-4.843707553304744) q[4];
u3(pi/2,1.7460971968652068,-1.7460971968652068) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3.3168935236601036,-3.3168935236601036) q[3];
u3(pi/2,4.843707553304744,-4.843707553304744) q[4];
u3(pi/2,1.6914334846927446,-1.6914334846927446) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,0.12063715789784804,-0.12063715789784804) q[4];
u3(pi/2,3.2521767149961534,-3.2521767149961534) q[4];
u3(pi/2,3.326946620151591,-3.326946620151591) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3.0573979704735863,-3.0573979704735863) q[7];
u3(pi/2,1.5695396897334606,-1.5695396897334606) q[6];
rzz(-pi/2) q[7],q[6];
u3(pi/2,3.108291771461741,-3.108291771461741) q[6];
u3(pi/2,6.214070268800611,-6.214070268800611) q[7];
u3(pi/2,3.06242451871933,-3.06242451871933) q[7];
rzz(-pi/2) q[6],q[7];
u3(pi/2,1.4916281919244339,-1.4916281919244339) q[7];
u3(pi/2,1.501681288415921,-1.501681288415921) q[7];
u3(pi/2,6.259937521543022,-6.259937521543022) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,4.689141194748125,-4.689141194748125) q[6];
u3(pi/2,1.6084954386379742,-1.6084954386379742) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,3.1792917654328705,-3.1792917654328705) q[5];
u3(pi/2,4.689141194748125,-4.689141194748125) q[6];
u3(pi/2,1.5374954446668447,-1.5374954446668447) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,6.249884425051534,-6.249884425051534) q[6];
u3(pi/2,6.259937521543022,-6.259937521543022) q[6];
u3(pi/2,0.02764601535159018,-0.02764601535159018) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,3.2339554776053334,-3.2339554776053334) q[9];
u3(pi/2,4.845592508896897,-4.845592508896897) q[9];
u3(pi/2,4.814804900891716,-4.814804900891716) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,3.276681137694154,-3.276681137694154) q[8];
u3(pi/2,3.3326014869280525,-3.3326014869280525) q[9];
u3(pi/2,0.18032741831605412,-0.18032741831605412) q[9];
rzz(-pi/2) q[8],q[9];
u3(pi/2,4.892716398700744,-4.892716398700744) q[9];
u3(pi/2,4.903397813722949,-4.903397813722949) q[9];
u3(pi/2,0.14514158059584845,-0.14514158059584845) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,4.598663326324739,-4.598663326324739) q[10];
u3(pi/2,4.903397813722949,-4.903397813722949) q[9];
rzz(-pi/2) q[10],q[9];
u3(pi/2,4.857530560980538,-4.857530560980538) q[8];
u3(pi/2,1.501681288415921,-1.501681288415921) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,3.0724776152108175,-3.0724776152108175) q[7];
u3(pi/2,4.857530560980538,-4.857530560980538) q[8];
u3(pi/2,1.7058848108992577,-1.7058848108992577) q[8];
rzz(pi/2) q[7],q[8];
u3(pi/2,0.13508848410436108,-0.13508848410436108) q[8];
u3(pi/2,0.14514158059584845,-0.14514158059584845) q[8];
u3(pi/2,6.204017172309124,-6.204017172309124) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,3.3326014869280525,-3.3326014869280525) q[9];
u3(pi/2,1.457070672734946,-1.457070672734946) q[10];
u3(pi/2,4.587981911302534,-4.587981911302534) q[10];
rzz(-pi/2) q[9],q[10];
u3(pi/2,6.1587782380974305,-6.1587782380974305) q[10];
u3(pi/2,6.1694596531196355,-6.1694596531196355) q[10];
u3(pi/2,3.321920071905847,-3.321920071905847) q[9];
rzz(pi/2) q[10],q[9];
u3(pi/2,5.017123467782899,-5.017123467782899) q[0];
u3(pi/2,1.457070672734946,-1.457070672734946) q[1];
u3(pi,2.4039466985269096,-2.4039466985269096) q[2];
u3(pi/2,1.7561502933566946,-1.7561502933566946) q[3];
u3(pi,3.151645750081281,-3.151645750081281) q[4];
u3(pi/2,4.74003499573628,-4.74003499573628) q[5];
u3(pi/2,4.633220845514227,-4.633220845514227) q[7];
u3(pi/2,1.7511237451109507,-1.7511237451109507) q[9];
u3(pi/2,4.660866860865817,-4.660866860865817) q[11];
u3(pi/2,4.630707571391355,-4.630707571391355) q[13];
u3(pi/2,1.440734390936279,-1.440734390936279) q[15];
u3(pi/2,1.2327609572686349,-1.2327609572686349) q[17];
u3(pi/2,1.6443095948888977,-1.6443095948888977) q[19];
u3(pi/2,1.5594865932419733,-1.5594865932419733) q[21];
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
measure q[17] -> c[17];
measure q[18] -> c[18];
measure q[19] -> c[19];
measure q[20] -> c[20];
measure q[21] -> c[21];
measure q[22] -> c[22];
measure q[23] -> c[23];
