OPENQASM 2.0;
include "qelib1.inc";
qreg q[25];
creg c[25];
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
rzz(pi/2) q[16],q[17];
u3(pi/2,1.4759202286564848,-1.4759202286564848) q[17];
u3(pi/2,1.485973325147972,-1.485973325147972) q[17];
u3(pi/2,3.1792917654328705,-3.1792917654328705) q[16];
rzz(pi/2) q[17],q[16];
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
rzz(-pi/2) q[16],q[15];
u3(pi,3*pi/2,-3*pi/2) q[19];
u3(pi/2,1.5286989852367932,-1.5286989852367932) q[18];
rzz(pi/2) q[19],q[18];
u3(pi/2,3.108920089992459,-3.108920089992459) q[18];
u3(pi/2,3.169866987472101,-3.169866987472101) q[19];
u3(pi/2,0.017592918860102842,-0.017592918860102842) q[19];
rzz(-pi/2) q[18],q[19];
u3(pi/2,4.729981899244793,-4.729981899244793) q[19];
u3(pi/2,4.740663314266998,-4.740663314266998) q[19];
u3(pi/2,6.261194158604458,-6.261194158604458) q[18];
rzz(pi/2) q[19],q[18];
u3(pi/2,4.690397831809562,-4.690397831809562) q[18];
u3(pi/2,1.485973325147972,-1.485973325147972) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,3.0567696519428686,-3.0567696519428686) q[17];
u3(pi/2,4.690397831809562,-4.690397831809562) q[18];
u3(pi/2,1.5381237631975626,-1.5381237631975626) q[18];
rzz(pi/2) q[17],q[18];
u3(pi/2,6.250512743582252,-6.250512743582252) q[18];
u3(pi/2,6.261194158604458,-6.261194158604458) q[18];
u3(pi/2,6.188309209041175,-6.188309209041175) q[17];
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
u3(pi/2,4.740663314266998,-4.740663314266998) q[19];
rzz(pi/2) q[20],q[19];
u3(pi/2,0.028274333882308135,-0.028274333882308135) q[19];
u3(pi/2,1.6461945504810516,-1.6461945504810516) q[20];
u3(pi/2,4.7771057890486395,-4.7771057890486395) q[20];
rzz(pi/2) q[19],q[20];
u3(pi/2,3.206309462253743,-3.206309462253743) q[20];
u3(pi/2,3.2169908772759483,-3.2169908772759483) q[20];
u3(pi/2,3.1591855724498963,-3.1591855724498963) q[19];
rzz(pi/2) q[20],q[19];
u3(pi/2,3*pi/2,-3*pi/2) q[23];
u3(pi/2,4.760141188719255,-4.760141188719255) q[23];
u3(pi/2,1.5129910219688443,-1.5129910219688443) q[22];
rzz(-pi/2) q[23],q[22];
u3(pi/2,0.010053096491487338,-0.010053096491487338) q[22];
u3(pi/2,0.018221237390820797,-0.018221237390820797) q[23];
u3(pi/2,3.1491324759584085,-3.1491324759584085) q[23];
rzz(-pi/2) q[22],q[23];
u3(pi/2,1.578336149163512,-1.578336149163512) q[23];
u3(pi/2,1.5890175641857174,-1.5890175641857174) q[23];
u3(pi/2,3.1623271651034854,-3.1623271651034854) q[22];
rzz(pi/2) q[23],q[22];
u3(pi/2,1.5915308383085893,-1.5915308383085893) q[22];
u3(pi/2,1.5280706667060753,-1.5280706667060753) q[21];
rzz(pi/2) q[22],q[21];
u3(pi/2,3.0988669935009723,-3.0988669935009723) q[21];
u3(pi/2,1.5915308383085893,-1.5915308383085893) q[22];
u3(pi/2,4.7224420768761775,-4.7224420768761775) q[22];
rzz(pi/2) q[21],q[22];
u3(pi/2,3.151645750081281,-3.151645750081281) q[22];
u3(pi/2,3.1623271651034854,-3.1623271651034854) q[22];
u3(pi/2,6.22977823206856,-6.22977823206856) q[21];
rzz(-pi/2) q[22],q[21];
u3(pi/2,5.306778310443879,-5.306778310443879) q[24];
u3(pi/2,5.370238482046393,-5.370238482046393) q[24];
u3(pi/2,1.5890175641857174,-1.5890175641857174) q[23];
rzz(pi/2) q[24],q[23];
u3(pi/2,3.159813890980614,-3.159813890980614) q[23];
u3(pi/2,3.766141273123444,-3.766141273123444) q[24];
u3(pi/2,3.776822688145649,-3.776822688145649) q[24];
rzz(pi/2) q[23],q[24];
u3(pi/2,2.206026361350753,-2.206026361350753) q[24];
u3(pi/2,2.2160794578422403,-2.2160794578422403) q[24];
u3(pi/2,3.1704953060028194,-3.1704953060028194) q[23];
rzz(pi/2) q[24],q[23];
u3(pi/2,1.6543626913803848,-1.6543626913803848) q[3];
u3(pi/2,1.6933184402848986,-1.6933184402848986) q[3];
u3(pi/2,1.6015839348000767,-1.6015839348000767) q[2];
rzz(-pi/2) q[3],q[2];
u3(pi/2,0.06471680866394974,-0.06471680866394974) q[2];
u3(pi/2,0.17969909978533616,-0.17969909978533616) q[3];
u3(pi/2,3.3106103383529244,-3.3106103383529244) q[3];
rzz(-pi/2) q[2],q[3];
u3(pi/2,1.7398140115580274,-1.7398140115580274) q[3];
u3(pi/2,1.7504954265802328,-1.7504954265802328) q[3];
u3(pi/2,3.2169908772759483,-3.2169908772759483) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,1.6461945504810516,-1.6461945504810516) q[2];
u3(pi/2,1.5990706606772047,-1.5990706606772047) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3.169866987472101,-3.169866987472101) q[1];
u3(pi/2,1.6461945504810516,-1.6461945504810516) q[2];
u3(pi/2,4.7771057890486395,-4.7771057890486395) q[2];
rzz(-pi/2) q[1],q[2];
u3(pi/2,0.06471680866394974,-0.06471680866394974) q[2];
u3(pi/2,0.07539822368615504,-0.07539822368615504) q[2];
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
u3(pi/2,1.7504954265802328,-1.7504954265802328) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3.321291753375129,-3.321291753375129) q[3];
u3(pi/2,1.6323715428052563,-1.6323715428052563) q[4];
u3(pi/2,4.763911099903562,-4.763911099903562) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.1931147731086655,-3.1931147731086655) q[4];
u3(pi/2,3.2031678696001533,-3.2031678696001533) q[4];
u3(pi/2,0.16901768476313087,-0.16901768476313087) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,1.5293273037675112,-1.5293273037675112) q[6];
rzz(-pi/2) q[7],q[6];
u3(pi/2,3.068079385495792,-3.068079385495792) q[6];
u3(pi/2,6.183910979326148,-6.183910979326148) q[7];
u3(pi/2,6.193964075817636,-6.193964075817636) q[7];
rzz(pi/2) q[6],q[7];
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
rzz(pi/2) q[5],q[6];
u3(pi/2,3.068079385495792,-3.068079385495792) q[6];
u3(pi/2,3.0781324819872795,-3.0781324819872795) q[6];
u3(pi/2,6.236689735906458,-6.236689735906458) q[5];
rzz(-pi/2) q[6],q[5];
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
rzz(pi/2) q[9],q[8];
u3(pi/2,4.701079246831767,-4.701079246831767) q[10];
u3(pi/2,4.726840306591202,-4.726840306591202) q[9];
rzz(pi/2) q[10],q[9];
u3(pi/2,4.7387783586748435,-4.7387783586748435) q[8];
u3(pi/2,4.633849164044945,-4.633849164044945) q[7];
rzz(-pi/2) q[8],q[7];
u3(pi/2,3.0630528372500483,-3.0630528372500483) q[7];
u3(pi/2,1.5971857050850506,-1.5971857050850506) q[8];
u3(pi/2,4.7280969436526386,-4.7280969436526386) q[8];
rzz(pi/2) q[7],q[8];
u3(pi/2,3.1573006168577415,-3.1573006168577415) q[8];
u3(pi/2,3.1679820318799474,-3.1679820318799474) q[8];
u3(pi/2,6.193964075817636,-6.193964075817636) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,0.014451326206513048,-0.014451326206513048) q[9];
u3(pi/2,4.701079246831767,-4.701079246831767) q[10];
u3(pi/2,1.548805178219768,-1.548805178219768) q[10];
rzz(pi/2) q[9],q[10];
u3(pi/2,6.261194158604458,-6.261194158604458) q[10];
u3(pi/2,6.271875573626663,-6.271875573626663) q[10];
u3(pi/2,3.145990883304819,-3.145990883304819) q[9];
rzz(-pi/2) q[10],q[9];
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
rzz(pi/2) q[1],q[0];
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
rzz(pi/2) q[12],q[11];
u3(pi/2,6.183910979326148,-6.183910979326148) q[15];
u3(pi/2,1.4991680142930492,-1.4991680142930492) q[15];
u3(pi/2,3.1202298235453823,-3.1202298235453823) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,1.51738925168387,-1.51738925168387) q[14];
u3(pi/2,6.178884431080405,-6.178884431080405) q[15];
u3(pi/2,3.0266103624684066,-3.0266103624684066) q[15];
rzz(-pi/2) q[14],q[15];
u3(pi/2,1.45581403567351,-1.45581403567351) q[15];
u3(pi/2,1.4664954506957153,-1.4664954506957153) q[15];
u3(pi/2,4.669663320295868,-4.669663320295868) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,3.0988669935009723,-3.0988669935009723) q[14];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,3.0925838081937926,-3.0925838081937926) q[13];
u3(pi/2,3.0988669935009723,-3.0988669935009723) q[14];
u3(pi/2,6.22977823206856,-6.22977823206856) q[14];
rzz(-pi/2) q[13],q[14];
u3(pi/2,1.51738925168387,-1.51738925168387) q[14];
u3(pi/2,1.5280706667060753,-1.5280706667060753) q[14];
u3(pi/2,3.082530711702305,-3.082530711702305) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,6.188309209041175,-6.188309209041175) q[17];
u3(pi/2,1.4206281979533044,-1.4206281979533044) q[17];
u3(pi/2,0.07414158662471912,-0.07414158662471912) q[16];
rzz(pi/2) q[17],q[16];
u3(pi/2,4.850619057142641,-4.850619057142641) q[16];
u3(pi/2,6.077096829104096,-6.077096829104096) q[17];
u3(pi/2,2.9248227604920976,-2.9248227604920976) q[17];
rzz(-pi/2) q[16],q[17];
u3(pi/2,1.3540264336972008,-1.3540264336972008) q[17];
u3(pi/2,1.3647078487194062,-1.3647078487194062) q[17];
u3(pi/2,1.7190795000443349,-1.7190795000443349) q[16];
rzz(-pi/2) q[17],q[16];
u3(pi/2,3.289875826839231,-3.289875826839231) q[16];
u3(pi/2,1.4664954506957153,-1.4664954506957153) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,3.037291777490612,-3.037291777490612) q[15];
u3(pi/2,3.289875826839231,-3.289875826839231) q[16];
u3(pi/2,0.13823007675795088,-0.13823007675795088) q[16];
rzz(pi/2) q[15],q[16];
u3(pi/2,4.850619057142641,-4.850619057142641) q[16];
u3(pi/2,4.860672153634128,-4.860672153634128) q[16];
u3(pi/2,6.1682030160582,-6.1682030160582) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,4.729981899244793,-4.729981899244793) q[19];
u3(pi/2,0.020106192982974676,-0.020106192982974676) q[18];
rzz(pi/2) q[19],q[18];
u3(pi/2,4.741919951328434,-4.741919951328434) q[18];
u3(pi/2,4.7582562331271,-4.7582562331271) q[19];
u3(pi/2,1.6059821645151022,-1.6059821645151022) q[19];
rzz(pi/2) q[18],q[19];
u3(pi/2,3.176778491309999,-3.176778491309999) q[19];
u3(pi/2,3.187459906332204,-3.187459906332204) q[19];
u3(pi/2,4.751973047819921,-4.751973047819921) q[18];
rzz(pi/2) q[19],q[18];
u3(pi/2,3.1811767210250244,-3.1811767210250244) q[18];
u3(pi/2,4.506300502309199,-4.506300502309199) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,6.077096829104096,-6.077096829104096) q[17];
u3(pi/2,3.1811767210250244,-3.1811767210250244) q[18];
u3(pi/2,0.029530970943744055,-0.029530970943744055) q[18];
rzz(pi/2) q[17],q[18];
u3(pi/2,4.741919951328434,-4.741919951328434) q[18];
u3(pi/2,4.751973047819921,-4.751973047819921) q[18];
u3(pi/2,2.9248227604920976,-2.9248227604920976) q[17];
rzz(-pi/2) q[18],q[17];
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
rzz(pi/2) q[21],q[20];
u3(pi/2,3.1799200839635886,-3.1799200839635886) q[20];
u3(pi/2,3.187459906332204,-3.187459906332204) q[19];
rzz(pi/2) q[20],q[19];
u3(pi/2,4.7582562331271,-4.7582562331271) q[19];
u3(pi/2,3.1799200839635886,-3.1799200839635886) q[20];
u3(pi/2,0.02764601535159018,-0.02764601535159018) q[20];
rzz(-pi/2) q[19],q[20];
u3(pi/2,1.5984423421464868,-1.5984423421464868) q[20];
u3(pi/2,1.609123757168692,-1.609123757168692) q[20];
u3(pi/2,4.747574818104896,-4.747574818104896) q[19];
rzz(pi/2) q[20],q[19];
u3(pi/2,3.1704953060028194,-3.1704953060028194) q[23];
u3(pi/2,4.693539424463151,-4.693539424463151) q[23];
u3(pi/2,0.057805304826052194,-0.057805304826052194) q[22];
rzz(pi/2) q[23],q[22];
u3(pi/2,4.838681005059,-4.838681005059) q[22];
u3(pi/2,3.1522740686119985,-3.1522740686119985) q[23];
u3(pi/2,3.1623271651034854,-3.1623271651034854) q[23];
rzz(pi/2) q[22],q[23];
u3(pi/2,4.733123491898382,-4.733123491898382) q[23];
u3(pi/2,4.743804906920587,-4.743804906920587) q[23];
u3(pi/2,1.687035254977719,-1.687035254977719) q[22];
rzz(pi/2) q[23],q[22];
u3(pi/2,0.11623892818282235,-0.11623892818282235) q[22];
u3(pi/2,6.177627794018969,-6.177627794018969) q[21];
rzz(pi/2) q[22],q[21];
u3(pi/2,1.4652388136342795,-1.4652388136342795) q[21];
u3(pi/2,0.11623892818282235,-0.11623892818282235) q[22];
u3(pi/2,0.1262920246743097,-0.1262920246743097) q[22];
rzz(pi/2) q[21],q[22];
u3(pi/2,4.838681005059,-4.838681005059) q[22];
u3(pi/2,4.849362420081205,-4.849362420081205) q[22];
u3(pi/2,1.4759202286564848,-1.4759202286564848) q[21];
rzz(pi/2) q[22],q[21];
u3(pi/2,5.357672111432033,-5.357672111432033) q[24];
u3(pi/2,2.1532476047704443,-2.1532476047704443) q[24];
u3(pi/2,4.743804906920587,-4.743804906920587) q[23];
rzz(-pi/2) q[24],q[23];
u3(pi/2,3.173008580125691,-3.173008580125691) q[23];
u3(pi/2,0.5491503958474959,-0.5491503958474959) q[24];
u3(pi/2,0.5598318108697011,-0.5598318108697011) q[24];
rzz(pi/2) q[23],q[24];
u3(pi/2,5.272220791254391,-5.272220791254391) q[24];
u3(pi/2,5.282273887745879,-5.282273887745879) q[24];
u3(pi/2,3.1836899951478967,-3.1836899951478967) q[23];
rzz(pi/2) q[24],q[23];
u3(pi/2,0.16901768476313087,-0.16901768476313087) q[3];
u3(pi/2,1.701486581184232,-1.701486581184232) q[3];
u3(pi/2,0.02324778563656447,-0.02324778563656447) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,4.769565966680024,-4.769565966680024) q[2];
u3(pi/2,0.18786724068466962,-0.18786724068466962) q[3];
u3(pi/2,3.3187784792522574,-3.3187784792522574) q[3];
rzz(pi/2) q[2],q[3];
u3(pi/2,4.889574806047154,-4.889574806047154) q[3];
u3(pi/2,4.90025622106936,-4.90025622106936) q[3];
u3(pi/2,4.780247381702229,-4.780247381702229) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,3.209451054907333,-3.209451054907333) q[2];
u3(pi/2,3.083159030233023,-3.083159030233023) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,4.65395535702792,-4.65395535702792) q[1];
u3(pi/2,3.209451054907333,-3.209451054907333) q[2];
u3(pi/2,0.057176986295334235,-0.057176986295334235) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,4.769565966680024,-4.769565966680024) q[2];
u3(pi/2,1.6179202165987434,-1.6179202165987434) q[2];
u3(pi/2,4.664636772050124,-4.664636772050124) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3.095097082316664,-3.095097082316664) q[5];
u3(pi/2,4.707362432138946,-4.707362432138946) q[5];
u3(pi/2,3.1792917654328705,-3.1792917654328705) q[4];
rzz(-pi/2) q[5],q[4];
u3(pi/2,4.79846861909305,-4.79846861909305) q[4];
u3(pi/2,pi/5000,-pi/5000) q[5];
u3(pi/2,3.132167875629024,-3.132167875629024) q[5];
rzz(-pi/2) q[4],q[5];
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
rzz(-pi/2) q[3],q[4];
u3(pi/2,4.79846861909305,-4.79846861909305) q[4];
u3(pi/2,4.809150034115255,-4.809150034115255) q[4];
u3(pi/2,0.17718582566246432,-0.17718582566246432) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,1.4815750954329465,-1.4815750954329465) q[7];
u3(pi/2,6.261822477135176,-6.261822477135176) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,4.658981905273664,-4.658981905273664) q[6];
u3(pi/2,1.4972830587008954,-1.4972830587008954) q[7];
u3(pi/2,4.628194297268483,-4.628194297268483) q[7];
rzz(-pi/2) q[6],q[7];
u3(pi/2,3.0573979704735863,-3.0573979704735863) q[7];
u3(pi/2,3.068079385495792,-3.068079385495792) q[7];
u3(pi/2,1.5274423481753574,-1.5274423481753574) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,6.239831328560047,-6.239831328560047) q[6];
u3(pi/2,1.5714246453256144,-1.5714246453256144) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,3.142220972120511,-3.142220972120511) q[5];
u3(pi/2,6.239831328560047,-6.239831328560047) q[6];
u3(pi/2,3.0881855784787664,-3.0881855784787664) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,1.51738925168387,-1.51738925168387) q[6];
u3(pi/2,1.5274423481753574,-1.5274423481753574) q[6];
u3(pi/2,6.273760529218817,-6.273760529218817) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,0.0043982297150257105,-0.0043982297150257105) q[9];
u3(pi/2,1.6160352610065896,-1.6160352610065896) q[9];
u3(pi/2,3.1843183136786144,-3.1843183136786144) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,1.6461945504810516,-1.6461945504810516) q[8];
u3(pi/2,0.10304423903774522,-0.10304423903774522) q[9];
u3(pi/2,3.2339554776053334,-3.2339554776053334) q[9];
rzz(pi/2) q[8],q[9];
u3(pi/2,4.8047518044002295,-4.8047518044002295) q[9];
u3(pi/2,4.8154332194224345,-4.8154332194224345) q[9];
u3(pi/2,1.656247646972539,-1.656247646972539) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,3.078760800517997,-3.078760800517997) q[10];
u3(pi/2,4.8154332194224345,-4.8154332194224345) q[9];
rzz(pi/2) q[10],q[9];
u3(pi/2,0.08545132017764237,-0.08545132017764237) q[8];
u3(pi/2,3.068079385495792,-3.068079385495792) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,4.638875712290688,-4.638875712290688) q[7];
u3(pi/2,0.08545132017764237,-0.08545132017764237) q[8];
u3(pi/2,3.2169908772759483,-3.2169908772759483) q[8];
rzz(pi/2) q[7],q[8];
u3(pi/2,1.6461945504810516,-1.6461945504810516) q[8];
u3(pi/2,1.656247646972539,-1.656247646972539) q[8];
u3(pi/2,1.4866016436786902,-1.4866016436786902) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,0.10304423903774522,-0.10304423903774522) q[9];
u3(pi/2,3.078760800517997,-3.078760800517997) q[10];
u3(pi/2,6.210300357616303,-6.210300357616303) q[10];
rzz(pi/2) q[9],q[10];
u3(pi/2,4.639504030821406,-4.639504030821406) q[10];
u3(pi/2,4.649557127312894,-4.649557127312894) q[10];
u3(pi/2,3.2339554776053334,-3.2339554776053334) q[9];
rzz(pi/2) q[10],q[9];
u3(pi/2,6.235433098845021,-6.235433098845021) q[1];
u3(pi/2,0.33175218421908215,-0.33175218421908215) q[0];
u3(pi/2,1.9119732889747483,-1.9119732889747483) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,3.4563802374794905,-3.4563802374794905) q[0];
u3(pi/2,3.038548414552048,-3.038548414552048) q[1];
u3(pi/2,6.1694596531196355,-6.1694596531196355) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,4.598663326324739,-4.598663326324739) q[1];
u3(pi/2,1.4470175762434587,-1.4470175762434587) q[1];
u3(pi/2,3.446327140988003,-3.446327140988003) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.082530711702305,-3.082530711702305) q[11];
u3(pi/2,4.630079252860637,-4.630079252860637) q[10];
rzz(pi/2) q[11],q[10];
u3(pi/2,3.0171855845076374,-3.0171855845076374) q[10];
u3(pi/2,3.1007519490931257,-3.1007519490931257) q[11];
u3(pi/2,6.231663187660714,-6.231663187660714) q[11];
rzz(pi/2) q[10],q[11];
u3(pi/2,1.519274207276024,-1.519274207276024) q[11];
u3(pi/2,1.5299556222982291,-1.5299556222982291) q[11];
u3(pi/2,3.0278669995298424,-3.0278669995298424) q[10];
rzz(pi/2) q[11],q[10];
u3(pi/2,3.082530711702305,-3.082530711702305) q[13];
u3(pi/2,4.6181412007769955,-4.6181412007769955) q[13];
u3(pi/2,4.621282793430586,-4.621282793430586) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,3.023468769814817,-3.023468769814817) q[12];
u3(pi/2,3.0699643410879456,-3.0699643410879456) q[13];
u3(pi/2,6.201503898186251,-6.201503898186251) q[13];
rzz(-pi/2) q[12],q[13];
u3(pi/2,4.630707571391355,-4.630707571391355) q[13];
u3(pi/2,4.640760667882843,-4.640760667882843) q[13];
u3(pi/2,6.1751145198960975,-6.1751145198960975) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,4.604318193101201,-4.604318193101201) q[12];
u3(pi/2,1.5299556222982291,-1.5299556222982291) q[11];
rzz(-pi/2) q[12],q[11];
u3(pi/2,6.242344602682919,-6.242344602682919) q[11];
u3(pi/2,1.4627255395114078,-1.4627255395114078) q[12];
u3(pi/2,4.594265096609713,-4.594265096609713) q[12];
rzz(pi/2) q[11],q[12];
u3(pi/2,3.023468769814817,-3.023468769814817) q[12];
u3(pi/2,3.0335218663063044,-3.0335218663063044) q[12];
u3(pi/2,3.0900705340709207,-3.0900705340709207) q[11];
rzz(pi/2) q[12],q[11];
u3(pi/2,6.1682030160582,-6.1682030160582) q[15];
u3(pi/2,1.4834600510251004,-1.4834600510251004) q[15];
u3(pi/2,1.563884822956999,-1.563884822956999) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,6.244229558275073,-6.244229558275073) q[14];
u3(pi/2,6.163176467812456,-6.163176467812456) q[15];
u3(pi/2,3.011530717731176,-3.011530717731176) q[15];
rzz(pi/2) q[14],q[15];
u3(pi/2,4.582327044526072,-4.582327044526072) q[15];
u3(pi/2,4.59238014101756,-4.59238014101756) q[15];
u3(pi/2,6.254910973297278,-6.254910973297278) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,4.6841146465023815,-4.6841146465023815) q[14];
u3(pi/2,4.640760667882843,-4.640760667882843) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,6.211556994677739,-6.211556994677739) q[13];
u3(pi/2,4.6841146465023815,-4.6841146465023815) q[14];
u3(pi/2,1.531840577890383,-1.531840577890383) q[14];
rzz(pi/2) q[13],q[14];
u3(pi/2,6.244229558275073,-6.244229558275073) q[14];
u3(pi/2,6.254910973297278,-6.254910973297278) q[14];
u3(pi/2,3.0599112445964582,-3.0599112445964582) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,6.066415414081891,-6.066415414081891) q[17];
u3(pi/2,1.2993627215247385,-1.2993627215247385) q[17];
u3(pi/2,4.8977429469464875,-4.8977429469464875) q[16];
rzz(-pi/2) q[17],q[16];
u3(pi/2,0.24881413816431164,-0.24881413816431164) q[16];
u3(pi/2,2.814238699085737,-2.814238699085737) q[17];
u3(pi/2,5.945149937653325,-5.945149937653325) q[17];
rzz(-pi/2) q[16],q[17];
u3(pi/2,4.374353610858428,-4.374353610858428) q[17];
u3(pi/2,4.385035025880633,-4.385035025880633) q[17];
u3(pi/2,3.40108820677631,-3.40108820677631) q[16];
rzz(pi/2) q[17],q[16];
u3(pi/2,1.8302918799814134,-1.8302918799814134) q[16];
u3(pi/2,4.59238014101756,-4.59238014101756) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,6.163176467812456,-6.163176467812456) q[15];
u3(pi/2,1.8302918799814134,-1.8302918799814134) q[16];
u3(pi/2,4.961203118549001,-4.961203118549001) q[16];
rzz(-pi/2) q[15],q[16];
u3(pi/2,0.24881413816431164,-0.24881413816431164) q[16];
u3(pi/2,0.25949555318651696,-0.25949555318651696) q[16];
u3(pi/2,6.1531233713209685,-6.1531233713209685) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,0.035185837720205684,-0.035185837720205684) q[19];
u3(pi/2,1.6524777357882312,-1.6524777357882312) q[18];
rzz(pi/2) q[19],q[18];
u3(pi/2,0.09173450548482195,-0.09173450548482195) q[18];
u3(pi/2,0.06346017160251381,-0.06346017160251381) q[19];
u3(pi/2,3.1943714101701013,-3.1943714101701013) q[19];
rzz(-pi/2) q[18],q[19];
u3(pi/2,1.6235750833752052,-1.6235750833752052) q[19];
u3(pi/2,1.6342564983974104,-1.6342564983974104) q[19];
u3(pi/2,3.243380255566102,-3.243380255566102) q[18];
rzz(-pi/2) q[19],q[18];
u3(pi/2,4.814176582360999,-4.814176582360999) q[18];
u3(pi/2,4.385035025880633,-4.385035025880633) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,5.9558313526755295,-5.9558313526755295) q[17];
u3(pi/2,4.814176582360999,-4.814176582360999) q[18];
u3(pi/2,1.6625308322797185,-1.6625308322797185) q[18];
rzz(pi/2) q[17],q[18];
u3(pi/2,0.09173450548482195,-0.09173450548482195) q[18];
u3(pi/2,0.10178760197630929,-0.10178760197630929) q[18];
u3(pi/2,2.8035572840635314,-2.8035572840635314) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,3.0467165554513813,-3.0467165554513813) q[21];
u3(pi/2,1.5425219929125884,-1.5425219929125884) q[20];
rzz(-pi/2) q[21],q[20];
u3(pi/2,3.132167875629024,-3.132167875629024) q[20];
u3(pi/2,6.2404596470907645,-6.2404596470907645) q[21];
u3(pi/2,6.25114106211297,-6.25114106211297) q[21];
rzz(-pi/2) q[20],q[21];
u3(pi/2,4.680344735318074,-4.680344735318074) q[21];
u3(pi/2,4.690397831809562,-4.690397831809562) q[21];
u3(pi/2,3.1221147791375365,-3.1221147791375365) q[20];
rzz(pi/2) q[21],q[20];
u3(pi/2,1.55131845234264,-1.55131845234264) q[20];
u3(pi/2,4.775849151987203,-4.775849151987203) q[19];
rzz(pi/2) q[20],q[19];
u3(pi/2,0.06346017160251381,-0.06346017160251381) q[19];
u3(pi/2,1.55131845234264,-1.55131845234264) q[20];
u3(pi/2,1.561371548834127,-1.561371548834127) q[20];
rzz(pi/2) q[19],q[20];
u3(pi/2,6.273760529218817,-6.273760529218817) q[20];
u3(pi/2,pi/2500,-pi/2500) q[20];
u3(pi/2,0.07351326809400116,-0.07351326809400116) q[19];
rzz(-pi/2) q[20],q[19];
u3(pi/2,3.1836899951478967,-3.1836899951478967) q[23];
u3(pi/2,4.706734113608228,-4.706734113608228) q[23];
u3(pi/2,4.8864332133935635,-4.8864332133935635) q[22];
rzz(pi/2) q[23],q[22];
u3(pi/2,3.384123606446925,-3.384123606446925) q[22];
u3(pi/2,3.165468757757076,-3.165468757757076) q[23];
u3(pi/2,3.175521854248563,-3.175521854248563) q[23];
rzz(-pi/2) q[22],q[23];
u3(pi/2,1.6047255274536665,-1.6047255274536665) q[23];
u3(pi/2,1.6154069424758717,-1.6154069424758717) q[23];
u3(pi/2,3.37344219142472,-3.37344219142472) q[22];
rzz(pi/2) q[23],q[22];
u3(pi/2,1.8026458646298231,-1.8026458646298231) q[22];
u3(pi/2,4.690397831809562,-4.690397831809562) q[21];
rzz(-pi/2) q[22],q[21];
u3(pi/2,3.1196015050146646,-3.1196015050146646) q[21];
u3(pi/2,4.944238518219617,-4.944238518219617) q[22];
u3(pi/2,4.954919933241822,-4.954919933241822) q[22];
rzz(pi/2) q[21],q[22];
u3(pi/2,3.384123606446925,-3.384123606446925) q[22];
u3(pi/2,3.3948050214691303,-3.3948050214691303) q[22];
u3(pi/2,3.1302829200368696,-3.1302829200368696) q[21];
rzz(pi/2) q[22],q[21];
u3(pi/2,2.140681234156085,-2.140681234156085) q[24];
u3(pi/2,5.219442034674082,-5.219442034674082) q[24];
u3(pi/2,1.6154069424758717,-1.6154069424758717) q[23];
rzz(pi/2) q[24],q[23];
u3(pi/2,3.1862032692707682,-3.1862032692707682) q[23];
u3(pi/2,0.47375217216134075,-0.47375217216134075) q[24];
u3(pi/2,0.4844335871835461,-0.4844335871835461) q[24];
rzz(pi/2) q[23],q[24];
u3(pi/2,5.196822567568235,-5.196822567568235) q[24];
u3(pi/2,5.206875664059723,-5.206875664059723) q[24];
u3(pi/2,3.1968846842929737,-3.1968846842929737) q[23];
rzz(pi/2) q[24],q[23];
u3(pi/2,0.17718582566246432,-0.17718582566246432) q[3];
u3(pi/2,1.7096547220835654,-1.7096547220835654) q[3];
u3(pi/2,1.669442336117616,-1.669442336117616) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,0.13634512116579703,-0.13634512116579703) q[2];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[3];
u3(pi/2,3.326946620151591,-3.326946620151591) q[3];
rzz(-pi/2) q[2],q[3];
u3(pi/2,1.7561502933566946,-1.7561502933566946) q[3];
u3(pi/2,4.887689850455001,-4.887689850455001) q[3];
u3(pi/2,pi/25,-pi/25) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,4.838052686528282,-4.838052686528282) q[2];
u3(pi/2,1.4470175762434587,-1.4470175762434587) q[1];
rzz(-pi/2) q[2],q[1];
u3(pi/2,6.159406556628148,-6.159406556628148) q[1];
u3(pi/2,1.6964600329384885,-1.6964600329384885) q[2];
u3(pi/2,4.827371271506076,-4.827371271506076) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,3.2565749447111796,-3.2565749447111796) q[2];
u3(pi/2,0.10492919462989908,-0.10492919462989908) q[2];
u3(pi/2,6.1694596531196355,-6.1694596531196355) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,6.273760529218817,-6.273760529218817) q[5];
u3(pi/2,1.6028405718615124,-1.6028405718615124) q[5];
u3(pi/2,4.784645611417255,-4.784645611417255) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,3.262229811487641,-3.262229811487641) q[4];
u3(pi/2,0.03769911184307752,-0.03769911184307752) q[5];
u3(pi/2,3.169238668941383,-3.169238668941383) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,4.74003499573628,-4.74003499573628) q[5];
u3(pi/2,4.7500880922277675,-4.7500880922277675) q[5];
u3(pi/2,3.2729112265098466,-3.2729112265098466) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,1.7021148997149498,-1.7021148997149498) q[4];
u3(pi/2,4.887689850455001,-4.887689850455001) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,0.17530087007031048,-0.17530087007031048) q[3];
u3(pi/2,1.7021148997149498,-1.7021148997149498) q[4];
u3(pi/2,4.833026138282538,-4.833026138282538) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,3.262229811487641,-3.262229811487641) q[4];
u3(pi/2,0.11058406140636072,-0.11058406140636072) q[4];
u3(pi/2,0.18535396656179778,-0.18535396656179778) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,3.0573979704735863,-3.0573979704735863) q[7];
u3(pi/2,1.5695396897334606,-1.5695396897334606) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,6.249884425051534,-6.249884425051534) q[6];
u3(pi/2,3.0724776152108175,-3.0724776152108175) q[7];
u3(pi/2,6.204017172309124,-6.204017172309124) q[7];
rzz(-pi/2) q[6],q[7];
u3(pi/2,4.633220845514227,-4.633220845514227) q[7];
u3(pi/2,4.643273942005714,-4.643273942005714) q[7];
u3(pi/2,3.118344867953229,-3.118344867953229) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,1.547548541158332,-1.547548541158332) q[6];
u3(pi/2,4.7500880922277675,-4.7500880922277675) q[5];
rzz(-pi/2) q[6],q[5];
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
u3(pi/2,1.6732122473019237,-1.6732122473019237) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,0.13508848410436108,-0.13508848410436108) q[8];
u3(pi/2,3.3326014869280525,-3.3326014869280525) q[9];
u3(pi/2,0.18032741831605412,-0.18032741831605412) q[9];
rzz(-pi/2) q[8],q[9];
u3(pi/2,4.892716398700744,-4.892716398700744) q[9];
u3(pi/2,4.903397813722949,-4.903397813722949) q[9];
u3(pi/2,3.2867342341856416,-3.2867342341856416) q[8];
rzz(-pi/2) q[9],q[8];
u3(pi/2,1.457070672734946,-1.457070672734946) q[10];
u3(pi/2,1.761805160133156,-1.761805160133156) q[9];
rzz(pi/2) q[10],q[9];
u3(pi/2,4.857530560980538,-4.857530560980538) q[8];
u3(pi/2,4.643273942005714,-4.643273942005714) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,6.214070268800611,-6.214070268800611) q[7];
u3(pi/2,4.857530560980538,-4.857530560980538) q[8];
u3(pi/2,1.7058848108992577,-1.7058848108992577) q[8];
rzz(pi/2) q[7],q[8];
u3(pi/2,0.13508848410436108,-0.13508848410436108) q[8];
u3(pi/2,0.14514158059584845,-0.14514158059584845) q[8];
u3(pi/2,3.06242451871933,-3.06242451871933) q[7];
rzz(-pi/2) q[8],q[7];
u3(pi/2,3.3326014869280525,-3.3326014869280525) q[9];
u3(pi/2,1.457070672734946,-1.457070672734946) q[10];
u3(pi/2,4.587981911302534,-4.587981911302534) q[10];
rzz(pi/2) q[9],q[10];
u3(pi/2,3.0171855845076374,-3.0171855845076374) q[10];
u3(pi/2,3.0278669995298424,-3.0278669995298424) q[10];
u3(pi/2,0.18032741831605412,-0.18032741831605412) q[9];
rzz(pi/2) q[10],q[9];
u3(pi/2,5.017123467782899,-5.017123467782899) q[0];
u3(pi/2,4.598663326324739,-4.598663326324739) q[1];
u3(pi,5.545539352116704,-5.545539352116704) q[2];
u3(pi/2,1.7561502933566946,-1.7561502933566946) q[3];
u3(pi,1.580849423286384,-1.580849423286384) q[4];
u3(pi/2,4.74003499573628,-4.74003499573628) q[5];
u3(pi/2,4.633220845514227,-4.633220845514227) q[7];
u3(pi/2,4.892716398700744,-4.892716398700744) q[9];
u3(pi/2,1.519274207276024,-1.519274207276024) q[11];
u3(pi/2,1.489114917801562,-1.489114917801562) q[13];
u3(pi/2,4.582327044526072,-4.582327044526072) q[15];
u3(pi/2,1.2327609572686349,-1.2327609572686349) q[17];
u3(pi/2,4.785902248478691,-4.785902248478691) q[19];
u3(pi/2,4.701079246831767,-4.701079246831767) q[21];
u3(pi/2,4.76768101108787,-4.76768101108787) q[23];
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
measure q[24] -> c[24];
