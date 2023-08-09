OPENQASM 2.0;
include "qelib1.inc";
qreg q[23];
creg c[23];
u3(pi,1.6543626913803848,-1.6543626913803848) q[1];
u3(pi/2,1.9440175340413641,-1.9440175340413641) q[0];
rzz(-pi/2) q[1],q[0];
u3(pi/2,0.34808846601774907,-0.34808846601774907) q[0];
u3(pi/2,0.028274333882308135,-0.028274333882308135) q[1];
u3(pi/2,3.159813890980614,-3.159813890980614) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,1.5890175641857174,-1.5890175641857174) q[1];
u3(pi/2,1.5990706606772047,-1.5990706606772047) q[1];
u3(pi/2,3.49973421609903,-3.49973421609903) q[0];
rzz(pi/2) q[1],q[0];
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
rzz(-pi/2) q[11],q[10];
u3(pi/2,1.487229962209408,-1.487229962209408) q[13];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[13];
u3(pi/2,1.6286016316209486,-1.6286016316209486) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,3.171751943064255,-3.171751943064255) q[12];
u3(pi/2,3.115203275299639,-3.115203275299639) q[13];
u3(pi/2,6.246742832397945,-6.246742832397945) q[13];
rzz(-pi/2) q[12],q[13];
u3(pi/2,4.675946505603048,-4.675946505603048) q[13];
u3(pi/2,4.685999602094536,-4.685999602094536) q[13];
u3(pi/2,0.04084070449666731,-0.04084070449666731) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,4.753229684881357,-4.753229684881357) q[12];
u3(pi/2,1.5148759775609983,-1.5148759775609983) q[11];
rzz(-pi/2) q[12],q[11];
u3(pi/2,6.227264957945688,-6.227264957945688) q[11];
u3(pi/2,1.611637031291564,-1.611637031291564) q[12];
u3(pi/2,4.742548269859152,-4.742548269859152) q[12];
rzz(pi/2) q[11],q[12];
u3(pi/2,3.171751943064255,-3.171751943064255) q[12];
u3(pi/2,3.18243335808646,-3.18243335808646) q[12];
u3(pi/2,3.0756192078644076,-3.0756192078644076) q[11];
rzz(pi/2) q[12],q[11];
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
rzz(-pi/2) q[15],q[14];
u3(pi/2,1.5136193404995624,-1.5136193404995624) q[14];
u3(pi/2,4.685999602094536,-4.685999602094536) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,6.256795928889432,-6.256795928889432) q[13];
u3(pi/2,1.5136193404995624,-1.5136193404995624) q[14];
u3(pi/2,4.6451588975978675,-4.6451588975978675) q[14];
rzz(pi/2) q[13],q[14];
u3(pi/2,3.074362570802972,-3.074362570802972) q[14];
u3(pi/2,3.0844156672944587,-3.0844156672944587) q[14];
u3(pi/2,3.1051501788081515,-3.1051501788081515) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,pi/2,-pi/2) q[17];
u3(pi/2,1.6260883574980767,-1.6260883574980767) q[17];
u3(pi/2,1.534353852013255,-1.534353852013255) q[16];
rzz(-pi/2) q[17],q[16];
u3(pi/2,0.02701769682087222,-0.02701769682087222) q[16];
u3(pi/2,6.282556988648868,-6.282556988648868) q[17];
u3(pi/2,3.1302829200368696,-3.1302829200368696) q[17];
rzz(-pi/2) q[16],q[17];
u3(pi/2,1.5594865932419733,-1.5594865932419733) q[17];
u3(pi/2,1.5701680082641787,-1.5701680082641787) q[17];
u3(pi/2,3.1792917654328705,-3.1792917654328705) q[16];
rzz(pi/2) q[17],q[16];
u3(pi/2,1.6084954386379742,-1.6084954386379742) q[16];
u3(pi/2,1.4822034139636644,-1.4822034139636644) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,3.052999740758561,-3.052999740758561) q[15];
u3(pi/2,1.6084954386379742,-1.6084954386379742) q[16];
u3(pi/2,4.739406677205562,-4.739406677205562) q[16];
rzz(pi/2) q[15],q[16];
u3(pi/2,3.168610350410665,-3.168610350410665) q[16];
u3(pi/2,3.1792917654328705,-3.1792917654328705) q[16];
u3(pi/2,6.183910979326148,-6.183910979326148) q[15];
rzz(-pi/2) q[16],q[15];
u3(pi,3*pi/2,-3*pi/2) q[19];
u3(pi/2,1.5079644737231006,-1.5079644737231006) q[18];
rzz(pi/2) q[19],q[18];
u3(pi/2,3.0881855784787664,-3.0881855784787664) q[18];
u3(pi/2,3.169866987472101,-3.169866987472101) q[19];
u3(pi/2,0.017592918860102842,-0.017592918860102842) q[19];
rzz(-pi/2) q[18],q[19];
u3(pi/2,4.729981899244793,-4.729981899244793) q[19];
u3(pi/2,4.740663314266998,-4.740663314266998) q[19];
u3(pi/2,6.239831328560047,-6.239831328560047) q[18];
rzz(pi/2) q[19],q[18];
u3(pi/2,4.669035001765151,-4.669035001765151) q[18];
u3(pi/2,1.5701680082641787,-1.5701680082641787) q[17];
rzz(-pi/2) q[18],q[17];
u3(pi/2,6.282556988648868,-6.282556988648868) q[17];
u3(pi/2,1.5274423481753574,-1.5274423481753574) q[18];
u3(pi/2,4.658981905273664,-4.658981905273664) q[18];
rzz(pi/2) q[17],q[18];
u3(pi/2,3.0881855784787664,-3.0881855784787664) q[18];
u3(pi/2,3.0982386749702537,-3.0982386749702537) q[18];
u3(pi/2,3.1302829200368696,-3.1302829200368696) q[17];
rzz(pi/2) q[18],q[17];
u3(pi,pi/2,-pi/2) q[21];
u3(pi/2,1.6166635795373074,-1.6166635795373074) q[20];
rzz(pi/2) q[21],q[20];
u3(pi/2,3.206309462253743,-3.206309462253743) q[20];
u3(pi/2,3.0988669935009723,-3.0988669935009723) q[21];
u3(pi/2,6.22977823206856,-6.22977823206856) q[21];
rzz(pi/2) q[20],q[21];
u3(pi/2,1.51738925168387,-1.51738925168387) q[21];
u3(pi/2,1.5280706667060753,-1.5280706667060753) q[21];
u3(pi/2,3.2169908772759483,-3.2169908772759483) q[20];
rzz(pi/2) q[21],q[20];
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
rzz(-pi/2) q[20],q[19];
u3(pi/2,4.143760710084937,-4.143760710084937) q[22];
u3(pi/2,4.209734155810323,-4.209734155810323) q[22];
u3(pi/2,1.5280706667060753,-1.5280706667060753) q[21];
rzz(pi/2) q[22],q[21];
u3(pi/2,3.0988669935009723,-3.0988669935009723) q[21];
u3(pi/2,2.6760086223277857,-2.6760086223277857) q[22];
u3(pi/2,2.686690037349991,-2.686690037349991) q[22];
rzz(pi/2) q[21],q[22];
u3(pi/2,1.1158937105550946,-1.1158937105550946) q[22];
u3(pi/2,1.125946807046582,-1.125946807046582) q[22];
u3(pi/2,3.1095484085231773,-3.1095484085231773) q[21];
rzz(pi/2) q[22],q[21];
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
rzz(pi/2) q[2],q[1];
u3(pi/2,3.169866987472101,-3.169866987472101) q[1];
u3(pi/2,4.7877872040708445,-4.7877872040708445) q[2];
u3(pi/2,1.6355131354588461,-1.6355131354588461) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,0.06471680866394974,-0.06471680866394974) q[2];
u3(pi/2,0.07539822368615504,-0.07539822368615504) q[2];
u3(pi/2,0.018221237390820797,-0.018221237390820797) q[1];
rzz(-pi/2) q[2],q[1];
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
rzz(pi/2) q[5],q[4];
u3(pi/2,4.77396419639505,-4.77396419639505) q[4];
u3(pi/2,4.8920880801700255,-4.8920880801700255) q[3];
rzz(-pi/2) q[4],q[3];
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
rzz(-pi/2) q[6],q[7];
u3(pi/2,1.4815750954329465,-1.4815750954329465) q[7];
u3(pi/2,1.4922565104551517,-1.4922565104551517) q[7];
u3(pi/2,6.219725135577073,-6.219725135577073) q[6];
rzz(-pi/2) q[7],q[6];
u3(pi/2,1.5073361551923827,-1.5073361551923827) q[6];
u3(pi/2,4.675946505603048,-4.675946505603048) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,6.246742832397945,-6.246742832397945) q[5];
u3(pi/2,1.5073361551923827,-1.5073361551923827) q[6];
u3(pi/2,4.638875712290688,-4.638875712290688) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,3.068079385495792,-3.068079385495792) q[6];
u3(pi/2,3.0781324819872795,-3.0781324819872795) q[6];
u3(pi/2,3.095097082316664,-3.095097082316664) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,4.628822615799201,-4.628822615799201) q[9];
u3(pi/2,4.669035001765151,-4.669035001765151) q[9];
u3(pi/2,1.5538317264655117,-1.5538317264655117) q[8];
rzz(-pi/2) q[9],q[8];
u3(pi/2,pi/200,-pi/200) q[8];
u3(pi/2,0.014451326206513048,-0.014451326206513048) q[9];
u3(pi/2,3.145990883304819,-3.145990883304819) q[9];
rzz(-pi/2) q[8],q[9];
u3(pi/2,1.5751945565099221,-1.5751945565099221) q[9];
u3(pi/2,1.5852476530014097,-1.5852476530014097) q[9];
u3(pi/2,3.1679820318799474,-3.1679820318799474) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,1.5594865932419733,-1.5594865932419733) q[10];
u3(pi/2,1.5852476530014097,-1.5852476530014097) q[9];
rzz(pi/2) q[10],q[9];
u3(pi/2,1.5971857050850506,-1.5971857050850506) q[8];
u3(pi/2,4.633849164044945,-4.633849164044945) q[7];
rzz(-pi/2) q[8],q[7];
u3(pi/2,3.0630528372500483,-3.0630528372500483) q[7];
u3(pi/2,4.7387783586748435,-4.7387783586748435) q[8];
u3(pi/2,1.5865042900628454,-1.5865042900628454) q[8];
rzz(pi/2) q[7],q[8];
u3(pi/2,pi/200,-pi/200) q[8];
u3(pi/2,0.02638937829015426,-0.02638937829015426) q[8];
u3(pi/2,6.193964075817636,-6.193964075817636) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,3.156043979796306,-3.156043979796306) q[9];
u3(pi/2,1.5594865932419733,-1.5594865932419733) q[10];
u3(pi/2,4.690397831809562,-4.690397831809562) q[10];
rzz(pi/2) q[9],q[10];
u3(pi/2,3.1196015050146646,-3.1196015050146646) q[10];
u3(pi/2,3.1302829200368696,-3.1302829200368696) q[10];
u3(pi/2,0.0043982297150257105,-0.0043982297150257105) q[9];
rzz(-pi/2) q[10],q[9];
u3(pi/2,4.730610217775511,-4.730610217775511) q[1];
u3(pi/2,1.9295662078348508,-1.9295662078348508) q[0];
u3(pi/2,3.5091589940597987,-3.5091589940597987) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,1.9126016075054662,-1.9126016075054662) q[0];
u3(pi/2,4.669663320295868,-4.669663320295868) q[1];
u3(pi/2,1.5180175702145882,-1.5180175702145882) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,6.230406550599278,-6.230406550599278) q[1];
u3(pi/2,3.0781324819872795,-3.0781324819872795) q[1];
u3(pi/2,1.901920192483261,-1.901920192483261) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,4.646415534659305,-4.646415534659305) q[11];
u3(pi/2,6.2517693806436885,-6.2517693806436885) q[10];
rzz(pi/2) q[11],q[10];
u3(pi/2,4.639504030821406,-4.639504030821406) q[10];
u3(pi/2,4.664008453519407,-4.664008453519407) q[11];
u3(pi/2,1.5117343849074085,-1.5117343849074085) q[11];
rzz(pi/2) q[10],q[11];
u3(pi/2,3.082530711702305,-3.082530711702305) q[11];
u3(pi/2,3.0932121267245103,-3.0932121267245103) q[11];
u3(pi/2,4.649557127312894,-4.649557127312894) q[10];
rzz(pi/2) q[11],q[10];
u3(pi/2,3.1051501788081515,-3.1051501788081515) q[13];
u3(pi/2,4.640760667882843,-4.640760667882843) q[13];
u3(pi/2,3.1246280532604085,-3.1246280532604085) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,1.5268140296446395,-1.5268140296446395) q[12];
u3(pi/2,3.0925838081937926,-3.0925838081937926) q[13];
u3(pi/2,6.224123365292098,-6.224123365292098) q[13];
rzz(-pi/2) q[12],q[13];
u3(pi/2,4.653327038497202,-4.653327038497202) q[13];
u3(pi/2,4.663380134988689,-4.663380134988689) q[13];
u3(pi/2,4.679088098256638,-4.679088098256638) q[12];
rzz(-pi/2) q[13],q[12];
u3(pi/2,6.249884425051534,-6.249884425051534) q[12];
u3(pi/2,3.0932121267245103,-3.0932121267245103) q[11];
rzz(pi/2) q[12],q[11];
u3(pi/2,4.664008453519407,-4.664008453519407) q[11];
u3(pi/2,6.249884425051534,-6.249884425051534) q[12];
u3(pi/2,3.097610356439536,-3.097610356439536) q[12];
rzz(pi/2) q[11],q[12];
u3(pi/2,1.5268140296446395,-1.5268140296446395) q[12];
u3(pi/2,1.5374954446668447,-1.5374954446668447) q[12];
u3(pi/2,1.5117343849074085,-1.5117343849074085) q[11];
rzz(pi/2) q[12],q[11];
u3(pi/2,3.042318325736356,-3.042318325736356) q[15];
u3(pi/2,4.640760667882843,-4.640760667882843) q[15];
u3(pi/2,3.1202298235453823,-3.1202298235453823) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,1.51738925168387,-1.51738925168387) q[14];
u3(pi/2,3.037291777490612,-3.037291777490612) q[15];
u3(pi/2,6.1682030160582,-6.1682030160582) q[15];
rzz(pi/2) q[14],q[15];
u3(pi/2,1.45581403567351,-1.45581403567351) q[15];
u3(pi/2,1.4664954506957153,-1.4664954506957153) q[15];
u3(pi/2,1.5280706667060753,-1.5280706667060753) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,6.2404596470907645,-6.2404596470907645) q[14];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,3.0925838081937926,-3.0925838081937926) q[13];
u3(pi/2,6.2404596470907645,-6.2404596470907645) q[14];
u3(pi/2,3.0881855784787664,-3.0881855784787664) q[14];
rzz(pi/2) q[13],q[14];
u3(pi/2,1.51738925168387,-1.51738925168387) q[14];
u3(pi/2,1.5280706667060753,-1.5280706667060753) q[14];
u3(pi/2,6.224123365292098,-6.224123365292098) q[13];
rzz(-pi/2) q[14],q[13];
u3(pi/2,3.1302829200368696,-3.1302829200368696) q[17];
u3(pi/2,4.645787216128586,-4.645787216128586) q[17];
u3(pi/2,0.07414158662471912,-0.07414158662471912) q[16];
rzz(pi/2) q[17],q[16];
u3(pi/2,4.850619057142641,-4.850619057142641) q[16];
u3(pi/2,3.0190705400997913,-3.0190705400997913) q[17];
u3(pi/2,6.150610097198097,-6.150610097198097) q[17];
rzz(-pi/2) q[16],q[17];
u3(pi/2,4.5798137704032005,-4.5798137704032005) q[17];
u3(pi/2,4.589866866894688,-4.589866866894688) q[17];
u3(pi/2,1.7190795000443349,-1.7190795000443349) q[16];
rzz(pi/2) q[17],q[16];
u3(pi/2,0.14828317324943824,-0.14828317324943824) q[16];
u3(pi/2,1.4664954506957153,-1.4664954506957153) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,3.037291777490612,-3.037291777490612) q[15];
u3(pi/2,0.14828317324943824,-0.14828317324943824) q[16];
u3(pi/2,3.2798227303477443,-3.2798227303477443) q[16];
rzz(pi/2) q[15],q[16];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[16];
u3(pi/2,1.7190795000443349,-1.7190795000443349) q[16];
u3(pi/2,6.1682030160582,-6.1682030160582) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,1.5883892456549995,-1.5883892456549995) q[19];
u3(pi/2,3.1403360165283574,-3.1403360165283574) q[18];
rzz(pi/2) q[19],q[18];
u3(pi/2,1.5789644676942303,-1.5789644676942303) q[18];
u3(pi/2,1.6204334907216154,-1.6204334907216154) q[19];
u3(pi/2,4.751344729289203,-4.751344729289203) q[19];
rzz(pi/2) q[18],q[19];
u3(pi/2,0.038955748904513435,-0.038955748904513435) q[19];
u3(pi/2,0.049637163926718735,-0.049637163926718735) q[19];
u3(pi/2,1.5896458827164353,-1.5896458827164353) q[18];
rzz(pi/2) q[19],q[18];
u3(pi/2,0.01884955592153876,-0.01884955592153876) q[18];
u3(pi/2,4.589866866894688,-4.589866866894688) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,6.160663193689585,-6.160663193689585) q[17];
u3(pi/2,0.01884955592153876,-0.01884955592153876) q[18];
u3(pi/2,3.149760794489126,-3.149760794489126) q[18];
rzz(pi/2) q[17],q[18];
u3(pi/2,1.5789644676942303,-1.5789644676942303) q[18];
u3(pi/2,1.5896458827164353,-1.5896458827164353) q[18];
u3(pi/2,3.009017443608304,-3.009017443608304) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,4.680344735318074,-4.680344735318074) q[21];
u3(pi/2,0.008796459430051421,-0.008796459430051421) q[20];
rzz(pi/2) q[21],q[20];
u3(pi/2,4.74003499573628,-4.74003499573628) q[20];
u3(pi/2,4.732495173367664,-4.732495173367664) q[21];
u3(pi/2,4.743176588389869,-4.743176588389869) q[21];
rzz(-pi/2) q[20],q[21];
u3(pi/2,3.172380261594973,-3.172380261594973) q[21];
u3(pi/2,3.1830616766171786,-3.1830616766171786) q[21];
u3(pi/2,4.729353580714075,-4.729353580714075) q[20];
rzz(pi/2) q[21],q[20];
u3(pi/2,3.158557253919178,-3.158557253919178) q[20];
u3(pi/2,0.049637163926718735,-0.049637163926718735) q[19];
rzz(pi/2) q[20],q[19];
u3(pi/2,1.6204334907216154,-1.6204334907216154) q[19];
u3(pi/2,3.158557253919178,-3.158557253919178) q[20];
u3(pi/2,3.169238668941383,-3.169238668941383) q[20];
rzz(-pi/2) q[19],q[20];
u3(pi/2,4.74003499573628,-4.74003499573628) q[20];
u3(pi/2,4.750716410758485,-4.750716410758485) q[20];
u3(pi/2,4.772079240802896,-4.772079240802896) q[19];
rzz(pi/2) q[20],q[19];
u3(pi/2,4.2675394606363755,-4.2675394606363755) q[22];
u3(pi/2,1.0599733613211961,-1.0599733613211961) q[22];
u3(pi/2,3.1830616766171786,-3.1830616766171786) q[21];
rzz(pi/2) q[22],q[21];
u3(pi/2,4.753858003412075,-4.753858003412075) q[21];
u3(pi/2,2.6684687999591703,-2.6684687999591703) q[22];
u3(pi/2,2.6785218964506576,-2.6785218964506576) q[22];
rzz(pi/2) q[21],q[22];
u3(pi/2,1.107725569655761,-1.107725569655761) q[22];
u3(pi/2,1.1184069846779663,-1.1184069846779663) q[22];
u3(pi/2,4.763911099903562,-4.763911099903562) q[21];
rzz(-pi/2) q[22],q[21];
u3(pi/2,0.16901768476313087,-0.16901768476313087) q[3];
u3(pi/2,1.701486581184232,-1.701486581184232) q[3];
u3(pi/2,3.164840439226358,-3.164840439226358) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,1.6279733130902307,-1.6279733130902307) q[2];
u3(pi/2,0.18786724068466962,-0.18786724068466962) q[3];
u3(pi/2,3.3187784792522574,-3.3187784792522574) q[3];
rzz(-pi/2) q[2],q[3];
u3(pi/2,1.7479821524573609,-1.7479821524573609) q[3];
u3(pi/2,1.758663567479566,-1.758663567479566) q[3];
u3(pi/2,4.780247381702229,-4.780247381702229) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,3.209451054907333,-3.209451054907333) q[2];
u3(pi/2,3.0781324819872795,-3.0781324819872795) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,4.648928808782176,-4.648928808782176) q[1];
u3(pi/2,3.209451054907333,-3.209451054907333) q[2];
u3(pi/2,0.057176986295334235,-0.057176986295334235) q[2];
rzz(-pi/2) q[1],q[2];
u3(pi/2,1.6279733130902307,-1.6279733130902307) q[2];
u3(pi/2,4.759512870188536,-4.759512870188536) q[2];
u3(pi/2,1.5180175702145882,-1.5180175702145882) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3.095097082316664,-3.095097082316664) q[5];
u3(pi/2,4.707362432138946,-4.707362432138946) q[5];
u3(pi/2,3.1792917654328705,-3.1792917654328705) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,1.6568759655032568,-1.6568759655032568) q[4];
u3(pi/2,3.142220972120511,-3.142220972120511) q[5];
u3(pi/2,6.273760529218817,-6.273760529218817) q[5];
rzz(-pi/2) q[4],q[5];
u3(pi/2,4.702964202423921,-4.702964202423921) q[5];
u3(pi/2,4.713017298915408,-4.713017298915408) q[5];
u3(pi/2,4.809150034115255,-4.809150034115255) q[4];
rzz(-pi/2) q[5],q[4];
u3(pi/2,0.09676105373056564,-0.09676105373056564) q[4];
u3(pi/2,1.758663567479566,-1.758663567479566) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,3.329459894274463,-3.329459894274463) q[3];
u3(pi/2,0.09676105373056564,-0.09676105373056564) q[4];
u3(pi/2,3.2276722922981538,-3.2276722922981538) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,1.6568759655032568,-1.6568759655032568) q[4];
u3(pi/2,1.6675573805254624,-1.6675573805254624) q[4];
u3(pi/2,0.17718582566246432,-0.17718582566246432) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,1.4815750954329465,-1.4815750954329465) q[7];
u3(pi/2,3.1202298235453823,-3.1202298235453823) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,1.51738925168387,-1.51738925168387) q[6];
u3(pi/2,1.4972830587008954,-1.4972830587008954) q[7];
u3(pi/2,4.628194297268483,-4.628194297268483) q[7];
rzz(pi/2) q[6],q[7];
u3(pi/2,6.19899062406338,-6.19899062406338) q[7];
u3(pi/2,6.209672039085585,-6.209672039085585) q[7];
u3(pi/2,1.5274423481753574,-1.5274423481753574) q[6];
rzz(pi/2) q[7],q[6];
u3(pi/2,6.239831328560047,-6.239831328560047) q[6];
u3(pi/2,1.5714246453256144,-1.5714246453256144) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,3.142220972120511,-3.142220972120511) q[5];
u3(pi/2,6.239831328560047,-6.239831328560047) q[6];
u3(pi/2,3.0881855784787664,-3.0881855784787664) q[6];
rzz(-pi/2) q[5],q[6];
u3(pi/2,4.658981905273664,-4.658981905273664) q[6];
u3(pi/2,4.669035001765151,-4.669035001765151) q[6];
u3(pi/2,3.132167875629024,-3.132167875629024) q[5];
rzz(pi/2) q[6],q[5];
u3(pi/2,3.145990883304819,-3.145990883304819) q[9];
u3(pi/2,4.757627914596383,-4.757627914596383) q[9];
u3(pi/2,0.042725660088821185,-0.042725660088821185) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,4.7877872040708445,-4.7877872040708445) q[8];
u3(pi/2,3.244636892627538,-3.244636892627538) q[9];
u3(pi/2,0.09236282401553991,-0.09236282401553991) q[9];
rzz(-pi/2) q[8],q[9];
u3(pi/2,4.8047518044002295,-4.8047518044002295) q[9];
u3(pi/2,4.8154332194224345,-4.8154332194224345) q[9];
u3(pi/2,1.656247646972539,-1.656247646972539) q[8];
rzz(pi/2) q[9],q[8];
u3(pi/2,3.078760800517997,-3.078760800517997) q[10];
u3(pi/2,4.8154332194224345,-4.8154332194224345) q[9];
rzz(-pi/2) q[10],q[9];
u3(pi/2,0.08545132017764237,-0.08545132017764237) q[8];
u3(pi/2,6.209672039085585,-6.209672039085585) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,1.4972830587008954,-1.4972830587008954) q[7];
u3(pi/2,0.08545132017764237,-0.08545132017764237) q[8];
u3(pi/2,3.2169908772759483,-3.2169908772759483) q[8];
rzz(pi/2) q[7],q[8];
u3(pi/2,1.6461945504810516,-1.6461945504810516) q[8];
u3(pi/2,1.656247646972539,-1.656247646972539) q[8];
u3(pi/2,4.628194297268483,-4.628194297268483) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,3.244636892627538,-3.244636892627538) q[9];
u3(pi/2,6.220353454107791,-6.220353454107791) q[10];
u3(pi/2,3.06870770402651,-3.06870770402651) q[10];
rzz(-pi/2) q[9],q[10];
u3(pi/2,4.639504030821406,-4.639504030821406) q[10];
u3(pi/2,4.649557127312894,-4.649557127312894) q[10];
u3(pi/2,3.2339554776053334,-3.2339554776053334) q[9];
rzz(pi/2) q[10],q[9];
u3(pi/2,3.0888138970094845,-3.0888138970094845) q[1];
u3(pi/2,0.33112386568836416,-0.33112386568836416) q[0];
u3(pi/2,1.9113449704440304,-1.9113449704440304) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,0.31478758388969724,-0.31478758388969724) q[0];
u3(pi/2,3.0335218663063044,-3.0335218663063044) q[1];
u3(pi/2,6.164433104873892,-6.164433104873892) q[1];
rzz(-pi/2) q[0],q[1];
u3(pi/2,4.593636778078995,-4.593636778078995) q[1];
u3(pi/2,1.441362709466997,-1.441362709466997) q[1];
u3(pi/2,0.304106168867492,-0.304106168867492) q[0];
rzz(pi/2) q[1],q[0];
u3(pi/2,3.082530711702305,-3.082530711702305) q[11];
u3(pi/2,4.630079252860637,-4.630079252860637) q[10];
rzz(-pi/2) q[11],q[10];
u3(pi/2,6.1587782380974305,-6.1587782380974305) q[10];
u3(pi/2,6.242344602682919,-6.242344602682919) q[11];
u3(pi/2,3.0900705340709207,-3.0900705340709207) q[11];
rzz(-pi/2) q[10],q[11];
u3(pi/2,1.519274207276024,-1.519274207276024) q[11];
u3(pi/2,1.5299556222982291,-1.5299556222982291) q[11];
u3(pi/2,3.0278669995298424,-3.0278669995298424) q[10];
rzz(pi/2) q[11],q[10];
u3(pi/2,3.082530711702305,-3.082530711702305) q[13];
u3(pi/2,4.6181412007769955,-4.6181412007769955) q[13];
u3(pi/2,1.4796901398407925,-1.4796901398407925) q[12];
rzz(pi/2) q[13],q[12];
u3(pi/2,6.16506142340461,-6.16506142340461) q[12];
u3(pi/2,3.0699643410879456,-3.0699643410879456) q[13];
u3(pi/2,6.201503898186251,-6.201503898186251) q[13];
rzz(pi/2) q[12],q[13];
u3(pi/2,1.489114917801562,-1.489114917801562) q[13];
u3(pi/2,1.4991680142930492,-1.4991680142930492) q[13];
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
u3(pi/2,6.1682030160582,-6.1682030160582) q[15];
u3(pi/2,1.4834600510251004,-1.4834600510251004) q[15];
u3(pi/2,4.705477476546792,-4.705477476546792) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,3.10263690468528,-3.10263690468528) q[14];
u3(pi/2,6.163176467812456,-6.163176467812456) q[15];
u3(pi/2,3.011530717731176,-3.011530717731176) q[15];
rzz(-pi/2) q[14],q[15];
u3(pi/2,1.440734390936279,-1.440734390936279) q[15];
u3(pi/2,1.4507874874277664,-1.4507874874277664) q[15];
u3(pi/2,6.254910973297278,-6.254910973297278) q[14];
rzz(pi/2) q[15],q[14];
u3(pi/2,4.6841146465023815,-4.6841146465023815) q[14];
u3(pi/2,1.4991680142930492,-1.4991680142930492) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,3.0699643410879456,-3.0699643410879456) q[13];
u3(pi/2,4.6841146465023815,-4.6841146465023815) q[14];
u3(pi/2,1.531840577890383,-1.531840577890383) q[14];
rzz(-pi/2) q[13],q[14];
u3(pi/2,3.10263690468528,-3.10263690468528) q[14];
u3(pi/2,3.113318319707485,-3.113318319707485) q[14];
u3(pi/2,3.0599112445964582,-3.0599112445964582) q[13];
rzz(pi/2) q[14],q[13];
u3(pi/2,3.009017443608304,-3.009017443608304) q[17];
u3(pi/2,4.5245217397000195,-4.5245217397000195) q[17];
u3(pi/2,1.7561502933566946,-1.7561502933566946) q[16];
rzz(pi/2) q[17],q[16];
u3(pi/2,0.24881413816431164,-0.24881413816431164) q[16];
u3(pi/2,2.897805063671225,-2.897805063671225) q[17];
u3(pi/2,6.028716302238813,-6.028716302238813) q[17];
rzz(-pi/2) q[16],q[17];
u3(pi/2,4.457919975443916,-4.457919975443916) q[17];
u3(pi/2,4.468601390466122,-4.468601390466122) q[17];
u3(pi/2,3.40108820677631,-3.40108820677631) q[16];
rzz(pi/2) q[17],q[16];
u3(pi/2,1.8302918799814134,-1.8302918799814134) q[16];
u3(pi/2,1.4507874874277664,-1.4507874874277664) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,3.021583814222663,-3.021583814222663) q[15];
u3(pi/2,1.8302918799814134,-1.8302918799814134) q[16];
u3(pi/2,4.961203118549001,-4.961203118549001) q[16];
rzz(pi/2) q[15],q[16];
u3(pi/2,3.3904067917541045,-3.3904067917541045) q[16];
u3(pi/2,3.40108820677631,-3.40108820677631) q[16];
u3(pi/2,6.1531233713209685,-6.1531233713209685) q[15];
rzz(pi/2) q[16],q[15];
u3(pi/2,0.05969026041820607,-0.05969026041820607) q[19];
u3(pi/2,1.6317432242745384,-1.6317432242745384) q[18];
rzz(-pi/2) q[19],q[18];
u3(pi/2,3.2119643290302045,-3.2119643290302045) q[18];
u3(pi/2,3.169866987472101,-3.169866987472101) q[19];
u3(pi/2,3.1799200839635886,-3.1799200839635886) q[19];
rzz(-pi/2) q[18],q[19];
u3(pi/2,1.609123757168692,-1.609123757168692) q[19];
u3(pi/2,1.6198051721908973,-1.6198051721908973) q[19];
u3(pi/2,3.2019112325387176,-3.2019112325387176) q[18];
rzz(pi/2) q[19],q[18];
u3(pi/2,1.6311149057438206,-1.6311149057438206) q[18];
u3(pi/2,4.468601390466122,-4.468601390466122) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,6.039397717261019,-6.039397717261019) q[17];
u3(pi/2,1.6311149057438206,-1.6311149057438206) q[18];
u3(pi/2,1.641168002235308,-1.641168002235308) q[18];
rzz(-pi/2) q[17],q[18];
u3(pi/2,3.2119643290302045,-3.2119643290302045) q[18];
u3(pi/2,3.22264574405241,-3.22264574405241) q[18];
u3(pi/2,2.90848647869343,-2.90848647869343) q[17];
rzz(pi/2) q[18],q[17];
u3(pi/2,3.1931147731086655,-3.1931147731086655) q[21];
u3(pi/2,4.6841146465023815,-4.6841146465023815) q[20];
rzz(pi/2) q[21],q[20];
u3(pi/2,3.132167875629024,-3.132167875629024) q[20];
u3(pi/2,3.2458935296889737,-3.2458935296889737) q[21];
u3(pi/2,3.2559466261804615,-3.2559466261804615) q[21];
rzz(-pi/2) q[20],q[21];
u3(pi/2,1.685150299385565,-1.685150299385565) q[21];
u3(pi/2,1.6958317144077701,-1.6958317144077701) q[21];
u3(pi/2,3.1221147791375365,-3.1221147791375365) q[20];
rzz(pi/2) q[21],q[20];
u3(pi/2,1.55131845234264,-1.55131845234264) q[20];
u3(pi/2,1.6198051721908973,-1.6198051721908973) q[19];
rzz(-pi/2) q[20],q[19];
u3(pi/2,0.04900884539600077,-0.04900884539600077) q[19];
u3(pi/2,4.692911105932433,-4.692911105932433) q[20];
u3(pi/2,4.702964202423921,-4.702964202423921) q[20];
rzz(pi/2) q[19],q[20];
u3(pi/2,3.132167875629024,-3.132167875629024) q[20];
u3(pi/2,3.142849290651229,-3.142849290651229) q[20];
u3(pi/2,0.05969026041820607,-0.05969026041820607) q[19];
rzz(pi/2) q[20],q[19];
u3(pi/2,1.1184069846779663,-1.1184069846779663) q[22];
u3(pi/2,4.194026192542374,-4.194026192542374) q[22];
u3(pi/2,1.6958317144077701,-1.6958317144077701) q[21];
rzz(pi/2) q[22],q[21];
u3(pi/2,3.266628041202667,-3.266628041202667) q[21];
u3(pi/2,5.8018933126496295,-5.8018933126496295) q[22];
u3(pi/2,5.812574727671835,-5.812574727671835) q[22];
rzz(-pi/2) q[21],q[22];
u3(pi/2,1.1001857472871457,-1.1001857472871457) q[22];
u3(pi/2,1.1102388437786328,-1.1102388437786328) q[22];
u3(pi/2,0.13571680263507907,-0.13571680263507907) q[21];
rzz(pi/2) q[22],q[21];
u3(pi/2,0.17718582566246432,-0.17718582566246432) q[3];
u3(pi/2,1.7096547220835654,-1.7096547220835654) q[3];
u3(pi/2,4.81103498970741,-4.81103498970741) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,3.2779377747555904,-3.2779377747555904) q[2];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[3];
u3(pi/2,3.326946620151591,-3.326946620151591) q[3];
rzz(-pi/2) q[2],q[3];
u3(pi/2,1.7561502933566946,-1.7561502933566946) q[3];
u3(pi/2,4.887689850455001,-4.887689850455001) q[3];
u3(pi/2,3.267256359733385,-3.267256359733385) q[2];
rzz(pi/2) q[3],q[2];
u3(pi/2,1.6964600329384885,-1.6964600329384885) q[2];
u3(pi/2,1.441362709466997,-1.441362709466997) q[1];
rzz(-pi/2) q[2],q[1];
u3(pi/2,6.153751689851687,-6.153751689851687) q[1];
u3(pi/2,4.838052686528282,-4.838052686528282) q[2];
u3(pi/2,1.6857786179162828,-1.6857786179162828) q[2];
rzz(pi/2) q[1],q[2];
u3(pi/2,0.11498229112138643,-0.11498229112138643) q[2];
u3(pi/2,3.2465218482196927,-3.2465218482196927) q[2];
u3(pi/2,6.164433104873892,-6.164433104873892) q[1];
rzz(pi/2) q[2],q[1];
u3(pi/2,3.132167875629024,-3.132167875629024) q[5];
u3(pi/2,4.7444332254513055,-4.7444332254513055) q[5];
u3(pi/2,1.643052957827462,-1.643052957827462) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,0.12063715789784804,-0.12063715789784804) q[4];
u3(pi/2,3.1792917654328705,-3.1792917654328705) q[5];
u3(pi/2,0.02764601535159018,-0.02764601535159018) q[5];
rzz(pi/2) q[4],q[5];
u3(pi/2,1.5984423421464868,-1.5984423421464868) q[5];
u3(pi/2,1.6084954386379742,-1.6084954386379742) q[5];
u3(pi/2,0.13131857292005333,-0.13131857292005333) q[4];
rzz(pi/2) q[5],q[4];
u3(pi/2,4.843707553304744,-4.843707553304744) q[4];
u3(pi/2,4.887689850455001,-4.887689850455001) q[3];
rzz(pi/2) q[4],q[3];
u3(pi/2,0.17530087007031048,-0.17530087007031048) q[3];
u3(pi/2,4.843707553304744,-4.843707553304744) q[4];
u3(pi/2,1.6914334846927446,-1.6914334846927446) q[4];
rzz(pi/2) q[3],q[4];
u3(pi/2,0.12063715789784804,-0.12063715789784804) q[4];
u3(pi/2,3.2521767149961534,-3.2521767149961534) q[4];
u3(pi/2,0.18535396656179778,-0.18535396656179778) q[3];
rzz(-pi/2) q[4],q[3];
u3(pi/2,6.19899062406338,-6.19899062406338) q[7];
u3(pi/2,4.7111323433232535,-4.7111323433232535) q[6];
rzz(pi/2) q[7],q[6];
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
rzz(-pi/2) q[6],q[5];
u3(pi/2,0.03769911184307752,-0.03769911184307752) q[5];
u3(pi/2,1.547548541158332,-1.547548541158332) q[6];
u3(pi/2,4.679088098256638,-4.679088098256638) q[6];
rzz(pi/2) q[5],q[6];
u3(pi/2,3.108291771461741,-3.108291771461741) q[6];
u3(pi/2,3.118344867953229,-3.118344867953229) q[6];
u3(pi/2,3.169238668941383,-3.169238668941383) q[5];
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
u3(pi/2,1.501681288415921,-1.501681288415921) q[7];
rzz(pi/2) q[8],q[7];
u3(pi/2,3.0724776152108175,-3.0724776152108175) q[7];
u3(pi/2,4.857530560980538,-4.857530560980538) q[8];
u3(pi/2,1.7058848108992577,-1.7058848108992577) q[8];
rzz(pi/2) q[7],q[8];
u3(pi/2,0.13508848410436108,-0.13508848410436108) q[8];
u3(pi/2,0.14514158059584845,-0.14514158059584845) q[8];
u3(pi/2,6.204017172309124,-6.204017172309124) q[7];
rzz(-pi/2) q[8],q[7];
u3(pi/2,3.3326014869280525,-3.3326014869280525) q[9];
u3(pi/2,1.457070672734946,-1.457070672734946) q[10];
u3(pi/2,4.587981911302534,-4.587981911302534) q[10];
rzz(pi/2) q[9],q[10];
u3(pi/2,3.0171855845076374,-3.0171855845076374) q[10];
u3(pi/2,3.0278669995298424,-3.0278669995298424) q[10];
u3(pi/2,0.18032741831605412,-0.18032741831605412) q[9];
rzz(pi/2) q[10],q[9];
u3(pi/2,1.8749024956623885,-1.8749024956623885) q[0];
u3(pi/2,4.593636778078995,-4.593636778078995) q[1];
u3(pi,2.4039466985269096,-2.4039466985269096) q[2];
u3(pi/2,1.7561502933566946,-1.7561502933566946) q[3];
u3(pi,4.7224420768761775,-4.7224420768761775) q[4];
u3(pi/2,1.5984423421464868,-1.5984423421464868) q[5];
u3(pi/2,1.4916281919244339,-1.4916281919244339) q[7];
u3(pi/2,4.892716398700744,-4.892716398700744) q[9];
u3(pi/2,4.660866860865817,-4.660866860865817) q[11];
u3(pi/2,1.489114917801562,-1.489114917801562) q[13];
u3(pi/2,4.582327044526072,-4.582327044526072) q[15];
u3(pi/2,4.479282805488327,-4.479282805488327) q[17];
u3(pi/2,1.6304865872131027,-1.6304865872131027) q[19];
u3(pi/2,1.7065131294299756,-1.7065131294299756) q[21];
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
