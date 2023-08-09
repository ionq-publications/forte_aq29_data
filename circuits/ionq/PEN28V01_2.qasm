OPENQASM 2.0;
include "qelib1.inc";
qreg q[28];
creg c[28];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[26];
u3(pi,1.458327309796382,-1.458327309796382) q[27];
rzz(0.696317738218232) q[26],q[27];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[25];
rzz(1.3926846109455662) q[25],q[27];
u3(pi/2,2.9543537314358415,-2.9543537314358415) q[24];
u3(pi,4.101035049996115,-4.101035049996115) q[27];
rzz(0.3559047485398804) q[24],q[27];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[23];
rzz(0.711882130701912) q[23],q[27];
u3(pi/2,0.8186990455255001,-0.8186990455255001) q[22];
rzz(1.4236189941595216) q[22],q[27];
u3(pi/2,3.669380219392878,-3.669380219392878) q[21];
u3(pi,1.6461945504810516,-1.6461945504810516) q[27];
rzz(0.29409908100041643) q[21],q[27];
u3(pi/2,4.842450916243307,-4.842450916243307) q[20];
rzz(0.5883323394230677) q[20],q[27];
u3(pi/2,3.669380219392878,-3.669380219392878) q[19];
rzz(1.1763963240016657) q[19],q[27];
u3(pi/2,2.0891591146372126,-2.0891591146372126) q[18];
u3(pi,1.6895485291005905,-1.6895485291005905) q[27];
rzz(0.7884961664529695) q[18],q[27];
u3(pi/2,3.669380219392878,-3.669380219392878) q[17];
u3(pi,1.0536901760140165,-1.0536901760140165) q[27];
rzz(1.564519424673024) q[17],q[27];
u3(pi/2,3.6411058855105702,-3.6411058855105702) q[16];
u3(pi,1.8931237330532094,-1.8931237330532094) q[27];
rzz(0.01266696441112712) q[16],q[27];
u3(pi/2,3.669380219392878,-3.669380219392878) q[15];
rzz(0.02529549143613757) q[15],q[27];
u3(pi/2,3.5650793432936974,-3.5650793432936974) q[14];
rzz(0.050611450348413266) q[14],q[27];
u3(pi/2,3.669380219392878,-3.669380219392878) q[13];
rzz(0.10126171507406165) q[13],q[27];
u3(pi/2,3.2603448558954873,-3.2603448558954873) q[12];
rzz(0.20244580139365306) q[12],q[27];
u3(pi/2,3.669380219392878,-3.669380219392878) q[11];
rzz(0.4050468602962466) q[11],q[27];
u3(pi/2,5.188026108138184,-5.188026108138184) q[10];
rzz(0.8100087876351032) q[10],q[27];
u3(pi/2,3.6756634047000576,-3.6756634047000576) q[9];
u3(pi,3.7133625165431354,-3.7133625165431354) q[27];
rzz(1.5218605077201284) q[9],q[27];
u3(pi/2,0.33112386568836416,-0.33112386568836416) q[8];
u3(pi,1.0807078728348887,-1.0807078728348887) q[27];
rzz(pi/32) q[8],q[27];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[7];
rzz(0.19640131429629326) q[7],q[27];
u3(pi/2,6.037512761668864,-6.037512761668864) q[6];
rzz(pi/8) q[6],q[27];
u3(pi/2,3.82897312619524,-3.82897312619524) q[5];
rzz(0.7853830994606742) q[5],q[27];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
u3(pi,2.3675042237452684,-2.3675042237452684) q[27];
rzz(pi/2) q[4],q[27];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
rzz(0.7853830994606742) q[0],q[1];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/8) q[0],q[2];
u3(pi/2,11*pi/8,-11*pi/8) q[3];
rzz(0.19640131429629326) q[0],q[3];
rzz(pi/32) q[0],q[4];
rzz(0.049100328574073315) q[0],q[5];
rzz(pi/128) q[0],q[6];
rzz(0.012275082143518329) q[0],q[7];
u3(pi/2,3*pi/2,-3*pi/2) q[0];
u3(pi/2,0.33112386568836416,-0.33112386568836416) q[8];
rzz(0) q[0],q[8];
u3(pi,1.9974246091523906,-1.9974246091523906) q[0];
u3(pi/2,3.6756634047000576,-3.6756634047000576) q[9];
rzz(0) q[0],q[9];
u3(pi,5.992902145987889,-5.992902145987889) q[0];
u3(pi/2,2.0464334545483913,-2.0464334545483913) q[10];
rzz(0) q[0],q[10];
u3(pi,3.705194375643802,-3.705194375643802) q[0];
u3(pi/2,0.5284158843338032,-0.5284158843338032) q[11];
rzz(0) q[0],q[11];
u3(pi,1.4174866052997146,-1.4174866052997146) q[0];
u3(pi/2,0.11938052083641214,-0.11938052083641214) q[12];
rzz(0) q[0],q[12];
u3(pi,5.4129641421352135,-5.4129641421352135) q[0];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[13];
rzz(0) q[0],q[13];
u3(pi,3.125256371791126,-3.125256371791126) q[0];
u3(pi/2,0.42348668970390413,-0.42348668970390413) q[14];
rzz(0) q[0],q[14];
u3(pi,0.8375486014470389,-0.8375486014470389) q[0];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[15];
rzz(0) q[0],q[15];
u3(pi,4.833026138282538,-4.833026138282538) q[0];
u3(pi/2,3.6411058855105702,-3.6411058855105702) q[16];
rzz(0) q[0],q[16];
u3(pi,2.5453183679384503,-2.5453183679384503) q[0];
u3(pi/2,3.669380219392878,-3.669380219392878) q[17];
rzz(0) q[0],q[17];
u3(pi,0.25761059759436306,-0.25761059759436306) q[0];
u3(pi/2,2.0891591146372126,-2.0891591146372126) q[18];
rzz(0) q[0],q[18];
u3(pi,4.253088134429862,-4.253088134429862) q[0];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[19];
rzz(0) q[0],q[19];
u3(pi,1.9653803640857748,-1.9653803640857748) q[0];
u3(pi/2,1.700858262653514,-1.700858262653514) q[20];
rzz(0) q[0],q[20];
u3(pi,5.960857900921273,-5.960857900921273) q[0];
u3(pi/2,3.669380219392878,-3.669380219392878) q[21];
rzz(0) q[0],q[21];
u3(pi,3.673150130577186,-3.673150130577186) q[0];
u3(pi/2,3.960291699115293,-3.960291699115293) q[22];
rzz(0) q[0],q[22];
u3(pi/2,0.9581857593448869,-0.9581857593448869) q[0];
u3(pi/2,3.666238626739289,-3.666238626739289) q[23];
rzz(pi/2) q[0],q[23];
rzz(-pi/2) q[0],q[23];
u3(pi/2,4.09977841293468,-4.09977841293468) q[0];
u3(pi/2,2.9543537314358415,-2.9543537314358415) q[24];
rzz(0) q[0],q[24];
u3(pi/2,0.9581857593448869,-0.9581857593448869) q[0];
u3(pi/2,3.666238626739289,-3.666238626739289) q[25];
rzz(pi/2) q[0],q[25];
rzz(pi/2) q[0],q[25];
u3(pi/2,3*pi/4,-3*pi/4) q[1];
rzz(0.7853830994606742) q[1],q[2];
rzz(pi/8) q[1],q[3];
rzz(0.19640131429629326) q[1],q[4];
rzz(pi/32) q[1],q[5];
rzz(0.049100328574073315) q[1],q[6];
rzz(pi/128) q[1],q[7];
u3(pi/2,3.472716519278157,-3.472716519278157) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi/2,pi/4,-pi/4) q[1];
rzz(0) q[1],q[9];
u3(pi,3.387265199100515,-3.387265199100515) q[1];
rzz(0) q[1],q[10];
u3(pi,5.450034935447572,-5.450034935447572) q[1];
rzz(0) q[1],q[11];
u3(pi,1.229619364615045,-1.229619364615045) q[1];
rzz(0) q[1],q[12];
u3(pi,3.292389100962103,-3.292389100962103) q[1];
rzz(0) q[1],q[13];
u3(pi,5.354530518778443,-5.354530518778443) q[1];
rzz(0) q[1],q[14];
u3(pi,1.1341149479459154,-1.1341149479459154) q[1];
rzz(0) q[1],q[15];
u3(pi,3.1968846842929737,-3.1968846842929737) q[1];
rzz(0) q[1],q[16];
u3(pi,5.2596544206400315,-5.2596544206400315) q[1];
rzz(0) q[1],q[17];
u3(pi,1.0386105312767857,-1.0386105312767857) q[1];
rzz(0) q[1],q[18];
u3(pi,3.101380267623844,-3.101380267623844) q[1];
rzz(0) q[1],q[19];
u3(pi,5.164150003970902,-5.164150003970902) q[1];
rzz(0) q[1],q[20];
u3(pi,0.9437344331383739,-0.9437344331383739) q[1];
rzz(0) q[1],q[21];
u3(pi,3.005875850954714,-3.005875850954714) q[1];
rzz(0) q[1],q[22];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,0.5309291584566751,-0.5309291584566751) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi/2,3*pi/8,-3*pi/8) q[2];
rzz(0) q[2],q[10];
u3(pi,4.746318181043459,-4.746318181043459) q[2];
rzz(0) q[2],q[11];
u3(pi,2.458610410699372,-2.458610410699372) q[2];
rzz(0) q[2],q[12];
u3(pi/2,6.0274596651773775,-6.0274596651773775) q[2];
rzz(pi/2) q[2],q[13];
rzz(-pi/2) q[2],q[13];
rzz(pi/2) q[2],q[14];
rzz(pi/2) q[2],q[14];
rzz(pi/2) q[2],q[15];
rzz(pi/2) q[2],q[15];
rzz(-pi/2) q[2],q[16];
rzz(pi/2) q[2],q[16];
rzz(pi/2) q[2],q[17];
rzz(pi/2) q[2],q[17];
rzz(-pi/2) q[2],q[18];
rzz(pi/2) q[2],q[18];
rzz(pi/2) q[2],q[19];
rzz(pi/2) q[2],q[19];
rzz(pi/2) q[2],q[20];
rzz(-pi/2) q[2],q[20];
rzz(pi/2) q[2],q[21];
rzz(pi/2) q[2],q[21];
rzz(pi/2) q[2],q[22];
rzz(-pi/2) q[2],q[22];
u3(pi/2,2.885867011587584,-2.885867011587584) q[2];
rzz(0) q[2],q[23];
u3(pi/2,2.552229871776348,-2.552229871776348) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,5.183627878423159,-5.183627878423159) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,0.9814335449814514,-0.9814335449814514) q[3];
rzz(0) q[3],q[11];
u3(pi,3.3206634348444113,-3.3206634348444113) q[3];
rzz(0) q[3],q[12];
u3(pi,4.85690224244982,-4.85690224244982) q[3];
rzz(0) q[3],q[13];
u3(pi,0.10995574287564278,-0.10995574287564278) q[3];
rzz(0) q[3],q[14];
u3(pi,1.6461945504810516,-1.6461945504810516) q[3];
rzz(0) q[3],q[15];
u3(pi,3.18243335808646,-3.18243335808646) q[3];
rzz(0) q[3],q[16];
u3(pi,4.718672165691869,-4.718672165691869) q[3];
rzz(0) q[3],q[17];
u3(pi,6.254910973297278,-6.254910973297278) q[3];
rzz(0) q[3],q[18];
u3(pi,1.5079644737231006,-1.5079644737231006) q[3];
rzz(0) q[3],q[19];
u3(pi,3.0442032813285094,-3.0442032813285094) q[3];
rzz(0) q[3],q[20];
u3(pi,4.580442088933919,-4.580442088933919) q[3];
rzz(0) q[3],q[21];
u3(pi,6.116680896539328,-6.116680896539328) q[3];
rzz(0) q[3],q[22];
u3(pi/2,2.172097160691983,-2.172097160691983) q[3];
rzz(pi/2) q[3],q[23];
rzz(pi/2) q[3],q[23];
u3(pi/2,2.172097160691983,-2.172097160691983) q[3];
rzz(0) q[3],q[24];
u3(pi/2,3.82897312619524,-3.82897312619524) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
rzz(pi/2) q[4],q[7];
u3(pi/2,4.466716434873968,-4.466716434873968) q[7];
u3(pi/2,1.1290883997001717,-1.1290883997001717) q[7];
rzz(-pi/2) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,3.6630970340856988,-3.6630970340856988) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,2.2581767994003434,-2.2581767994003434) q[4];
rzz(0) q[4],q[12];
u3(pi,4.783388974355819,-4.783388974355819) q[4];
rzz(0) q[4],q[13];
u3(pi,0.409663682028109,-0.409663682028109) q[4];
rzz(0) q[4],q[14];
u3(pi/2,2.9355041755143025,-2.9355041755143025) q[4];
rzz(pi/2) q[4],q[15];
rzz(pi/2) q[4],q[15];
rzz(pi/2) q[4],q[16];
rzz(-pi/2) q[4],q[16];
rzz(pi/2) q[4],q[17];
rzz(pi/2) q[4],q[17];
rzz(pi/2) q[4],q[18];
rzz(pi/2) q[4],q[18];
rzz(pi/2) q[4],q[19];
rzz(pi/2) q[4],q[19];
rzz(pi/2) q[4],q[20];
rzz(pi/2) q[4],q[20];
rzz(-pi/2) q[4],q[21];
rzz(pi/2) q[4],q[21];
rzz(pi/2) q[4],q[22];
rzz(pi/2) q[4],q[22];
u3(pi/2,2.9355041755143025,-2.9355041755143025) q[4];
rzz(0) q[4],q[23];
u3(pi/2,6.077096829104096,-6.077096829104096) q[4];
rzz(-pi/2) q[4],q[24];
rzz(pi/2) q[4],q[24];
u3(pi/2,2.9355041755143025,-2.9355041755143025) q[4];
rzz(0) q[4],q[25];
u3(pi/2,1.4237697906068942,-1.4237697906068942) q[5];
rzz(0.7853830994606742) q[5],q[6];
u3(pi/2,2.6998847264950685,-2.6998847264950685) q[7];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,3.255318307649744,-3.255318307649744) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,6.136158770991584,-6.136158770991584) q[5];
rzz(0) q[5],q[13];
u3(pi,2.3009024594891647,-2.3009024594891647) q[5];
rzz(0) q[5],q[14];
u3(pi,4.055796115784423,-4.055796115784423) q[5];
rzz(0) q[5],q[15];
u3(pi,5.810061453548963,-5.810061453548963) q[5];
rzz(0) q[5],q[16];
u3(pi,1.2817698026646356,-1.2817698026646356) q[5];
rzz(0) q[5],q[17];
u3(pi,3.0360351404291763,-3.0360351404291763) q[5];
rzz(0) q[5],q[18];
u3(pi,4.790928796724434,-4.790928796724434) q[5];
rzz(0) q[5],q[19];
u3(pi,0.26263714584010667,-0.26263714584010667) q[5];
rzz(0) q[5],q[20];
u3(pi,2.0169024836046474,-2.0169024836046474) q[5];
rzz(0) q[5],q[21];
u3(pi,3.771796139899905,-3.771796139899905) q[5];
rzz(0) q[5],q[22];
u3(pi/2,6.219725135577073,-6.219725135577073) q[5];
rzz(pi/2) q[5],q[23];
rzz(pi/2) q[5],q[23];
u3(pi/2,6.219725135577073,-6.219725135577073) q[5];
rzz(0) q[5],q[24];
u3(pi/2,3.0781324819872795,-3.0781324819872795) q[5];
rzz(pi/2) q[5],q[25];
rzz(-pi/2) q[5],q[25];
u3(pi/2,6.219725135577073,-6.219725135577073) q[5];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[26];
rzz(0) q[5],q[26];
u3(pi/2,4.049512930477243,-4.049512930477243) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,3.6630970340856988,-3.6630970340856988) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,5.620309257272139,-5.620309257272139) q[6];
rzz(0) q[6],q[14];
u3(pi,1.9396193043263381,-1.9396193043263381) q[6];
rzz(0) q[6],q[15];
u3(pi,4.001760722142679,-4.001760722142679) q[6];
rzz(0) q[6],q[16];
u3(pi/2,0.3210707691968768,-0.3210707691968768) q[6];
rzz(pi/2) q[6],q[17];
rzz(pi/2) q[6],q[17];
rzz(pi/2) q[6],q[18];
rzz(-pi/2) q[6],q[18];
rzz(pi/2) q[6],q[19];
rzz(pi/2) q[6],q[19];
rzz(pi/2) q[6],q[20];
rzz(-pi/2) q[6],q[20];
rzz(pi/2) q[6],q[21];
rzz(pi/2) q[6],q[21];
rzz(pi/2) q[6],q[22];
rzz(-pi/2) q[6],q[22];
u3(pi/2,3.46266342278667,-3.46266342278667) q[6];
rzz(0) q[6],q[23];
u3(pi/2,0.3210707691968768,-0.3210707691968768) q[6];
rzz(pi/2) q[6],q[24];
rzz(pi/2) q[6],q[24];
u3(pi/2,0.3210707691968768,-0.3210707691968768) q[6];
rzz(0) q[6],q[25];
u3(pi/2,3.46266342278667,-3.46266342278667) q[6];
rzz(pi/2) q[6],q[26];
rzz(pi/2) q[6],q[26];
u3(pi/2,0.9204866475018093,-0.9204866475018093) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,0.41783182292744253,-0.41783182292744253) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,5.632875627886499,-5.632875627886499) q[7];
rzz(0) q[7],q[15];
u3(pi,3.476486430462465,-3.476486430462465) q[7];
rzz(0) q[7],q[16];
u3(pi,2.305929007734908,-2.305929007734908) q[7];
rzz(0) q[7],q[17];
u3(pi,1.1353715850073511,-1.1353715850073511) q[7];
rzz(0) q[7],q[18];
u3(pi,6.247371150928663,-6.247371150928663) q[7];
rzz(0) q[7],q[19];
u3(pi,5.076813728201106,-5.076813728201106) q[7];
rzz(0) q[7],q[20];
u3(pi,3.906256305473549,-3.906256305473549) q[7];
rzz(0) q[7],q[21];
u3(pi,2.735070564215274,-2.735070564215274) q[7];
rzz(0) q[7],q[22];
u3(pi/2,0.5793096853219579,-0.5793096853219579) q[7];
rzz(-pi/2) q[7],q[23];
rzz(pi/2) q[7],q[23];
u3(pi/2,3.720902338911751,-3.720902338911751) q[7];
rzz(0) q[7],q[24];
u3(pi/2,0.5793096853219579,-0.5793096853219579) q[7];
rzz(pi/2) q[7],q[25];
rzz(pi/2) q[7],q[25];
u3(pi/2,0.5793096853219579,-0.5793096853219579) q[7];
rzz(0) q[7],q[26];
u3(pi/2,4.933557103197411,-4.933557103197411) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,3.6630970340856988,-3.6630970340856988) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,0.22116812281272144,-0.22116812281272144) q[8];
rzz(0) q[8],q[16];
u3(pi,3.960291699115293,-3.960291699115293) q[8];
rzz(0) q[8],q[17];
u3(pi,2.015017528012493,-2.015017528012493) q[8];
rzz(0) q[8],q[18];
u3(pi/2,5.754141104315065,-5.754141104315065) q[8];
rzz(-pi/2) q[8],q[19];
rzz(pi/2) q[8],q[19];
rzz(pi/2) q[8],q[20];
rzz(pi/2) q[8],q[20];
rzz(pi/2) q[8],q[21];
rzz(-pi/2) q[8],q[21];
rzz(pi/2) q[8],q[22];
rzz(pi/2) q[8],q[22];
u3(pi/2,5.754141104315065,-5.754141104315065) q[8];
rzz(0) q[8],q[23];
u3(pi/2,2.612548450725272,-2.612548450725272) q[8];
rzz(pi/2) q[8],q[24];
rzz(-pi/2) q[8],q[24];
u3(pi/2,5.754141104315065,-5.754141104315065) q[8];
rzz(0) q[8],q[25];
u3(pi/2,2.612548450725272,-2.612548450725272) q[8];
rzz(pi/2) q[8],q[26];
rzz(pi/2) q[8],q[26];
u3(pi/2,3.6140881886896983,-3.6140881886896983) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,0.4938583651443155,-0.4938583651443155) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,2.0432918618948013,-2.0432918618948013) q[9];
rzz(0) q[9],q[17];
u3(pi,5.611512797842089,-5.611512797842089) q[9];
rzz(0) q[9],q[18];
u3(pi,3.3238050274980013,-3.3238050274980013) q[9];
rzz(0) q[9],q[19];
u3(pi,1.0360972571539138,-1.0360972571539138) q[9];
rzz(0) q[9],q[20];
u3(pi,5.031574793989412,-5.031574793989412) q[9];
rzz(0) q[9],q[21];
u3(pi,2.743867023645325,-2.743867023645325) q[9];
rzz(0) q[9],q[22];
u3(pi/2,0.029530970943744055,-0.029530970943744055) q[9];
rzz(pi/2) q[9],q[23];
rzz(pi/2) q[9],q[23];
u3(pi/2,0.029530970943744055,-0.029530970943744055) q[9];
rzz(0) q[9],q[24];
u3(pi/2,3.1711236245335375,-3.1711236245335375) q[9];
rzz(pi/2) q[9],q[25];
rzz(pi/2) q[9],q[25];
u3(pi/2,3.1711236245335375,-3.1711236245335375) q[9];
rzz(0) q[9],q[26];
u3(pi/2,2.79476082463348,-2.79476082463348) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,3.6630970340856988,-3.6630970340856988) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi/2,4.365557151428376,-4.365557151428376) q[10];
rzz(0) q[10],q[18];
u3(pi,1.6512210987267952,-1.6512210987267952) q[10];
rzz(0) q[10],q[19];
u3(pi,5.646698635562294,-5.646698635562294) q[10];
rzz(0) q[10],q[20];
u3(pi/2,2.931734264329995,-2.931734264329995) q[10];
rzz(pi/2) q[10],q[21];
rzz(-pi/2) q[10],q[21];
rzz(pi/2) q[10],q[22];
rzz(pi/2) q[10],q[22];
u3(pi/2,6.073326917919788,-6.073326917919788) q[10];
rzz(0) q[10],q[23];
u3(pi/2,2.931734264329995,-2.931734264329995) q[10];
rzz(pi/2) q[10],q[24];
rzz(pi/2) q[10],q[24];
u3(pi/2,2.931734264329995,-2.931734264329995) q[10];
rzz(0) q[10],q[25];
u3(pi/2,6.073326917919788,-6.073326917919788) q[10];
rzz(pi/2) q[10],q[26];
rzz(pi/2) q[10],q[26];
u3(pi/2,1.6807520696705394,-1.6807520696705394) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,2.0835042478607506,-2.0835042478607506) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,0.10995574287564278,-0.10995574287564278) q[11];
rzz(0) q[11],q[19];
u3(pi,3.5858138548073897,-3.5858138548073897) q[11];
rzz(0) q[11],q[20];
u3(pi,1.1133804364322226,-1.1133804364322226) q[11];
rzz(0) q[11],q[21];
u3(pi,4.924132325236641,-4.924132325236641) q[11];
rzz(0) q[11],q[22];
u3(pi/2,2.1168051299888027,-2.1168051299888027) q[11];
rzz(pi/2) q[11],q[23];
rzz(pi/2) q[11],q[23];
u3(pi/2,2.1168051299888027,-2.1168051299888027) q[11];
rzz(0) q[11],q[24];
u3(pi/2,5.258397783578595,-5.258397783578595) q[11];
rzz(-pi/2) q[11],q[25];
rzz(pi/2) q[11],q[25];
u3(pi/2,2.1168051299888027,-2.1168051299888027) q[11];
rzz(0) q[11],q[26];
u3(pi/2,1.4752919101257669,-1.4752919101257669) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/128) q[12],q[18];
u3(pi/2,0.5215043804959056,-0.5215043804959056) q[19];
rzz(0.012275082143518329) q[12],q[19];
u3(pi/2,3.0460882369206637,-3.0460882369206637) q[12];
rzz(0) q[12],q[20];
u3(pi,5.385318126783623,-5.385318126783623) q[12];
rzz(0) q[12],q[21];
u3(pi,0.638371627209446,-0.638371627209446) q[12];
rzz(0) q[12],q[22];
u3(pi,2.1746104348148547,-2.1746104348148547) q[12];
rzz(0) q[12],q[23];
u3(pi/2,4.513212006147097,-4.513212006147097) q[12];
rzz(pi/2) q[12],q[24];
rzz(pi/2) q[12],q[24];
u3(pi/2,4.513212006147097,-4.513212006147097) q[12];
rzz(0) q[12],q[25];
u3(pi/2,1.3716193525573037,-1.3716193525573037) q[12];
rzz(-pi/2) q[12],q[26];
rzz(pi/2) q[12],q[26];
u3(pi/2,1.9848582385380313,-1.9848582385380313) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/128) q[13],q[19];
u3(pi/2,4.8374243679975635,-4.8374243679975635) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi/2,0.4140619117431347,-0.4140619117431347) q[13];
rzz(0) q[13],q[21];
u3(pi,3.015928947446201,-3.015928947446201) q[13];
rzz(0) q[13],q[22];
u3(pi,5.0786986837932595,-5.0786986837932595) q[13];
rzz(0) q[13],q[23];
u3(pi,0.8582831129607315,-0.8582831129607315) q[13];
rzz(0) q[13],q[24];
u3(pi/2,3.4601501486637978,-3.4601501486637978) q[13];
rzz(pi/2) q[13],q[25];
rzz(pi/2) q[13],q[25];
u3(pi/2,3.4601501486637978,-3.4601501486637978) q[13];
rzz(0) q[13],q[26];
u3(pi/2,5.073672135547516,-5.073672135547516) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/128) q[14],q[20];
u3(pi/2,0.5215043804959056,-0.5215043804959056) q[21];
rzz(0.012275082143518329) q[14],q[21];
u3(pi/2,0.3612831551628262,-0.3612831551628262) q[14];
rzz(0) q[14],q[22];
u3(pi,2.5415484567541426,-2.5415484567541426) q[14];
rzz(0) q[14],q[23];
u3(pi,3.7611147248777006,-3.7611147248777006) q[14];
rzz(0) q[14],q[24];
u3(pi,4.980680993001258,-4.980680993001258) q[14];
rzz(0) q[14],q[25];
u3(pi/2,0.8777609874129881,-0.8777609874129881) q[14];
rzz(pi/2) q[14],q[26];
rzz(-pi/2) q[14],q[26];
u3(pi/2,2.0608847807549044,-2.0608847807549044) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/32) q[15],q[19];
rzz(0.049100328574073315) q[15],q[20];
rzz(pi/128) q[15],q[21];
u3(pi/2,3.955265150869549,-3.955265150869549) q[22];
rzz(0.012275082143518329) q[15],q[22];
u3(pi/2,0.4900884539600077,-0.4900884539600077) q[15];
rzz(0) q[15],q[23];
u3(pi,4.0589377084380125,-4.0589377084380125) q[15];
rzz(0) q[15],q[24];
u3(pi,1.7712299380939251,-1.7712299380939251) q[15];
rzz(0) q[15],q[25];
u3(pi,5.766079156398706,-5.766079156398706) q[15];
rzz(0) q[15],q[26];
u3(pi/2,2.045805136017673,-2.045805136017673) q[16];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/8) q[16],q[18];
rzz(0.19640131429629326) q[16],q[19];
rzz(pi/32) q[16],q[20];
rzz(0.049100328574073315) q[16],q[21];
rzz(pi/128) q[16],q[22];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[23];
rzz(0.012275082143518329) q[16],q[23];
u3(pi/2,3.61660146281257,-3.61660146281257) q[16];
rzz(0) q[16],q[24];
u3(pi,0.9789202708585795,-0.9789202708585795) q[16];
rzz(0) q[16],q[25];
u3(pi,5.12770752918926,-5.12770752918926) q[16];
rzz(0) q[16],q[26];
u3(pi/2,0.5089380098815465,-0.5089380098815465) q[17];
rzz(0.7853830994606742) q[17],q[18];
rzz(pi/8) q[17],q[19];
rzz(0.19640131429629326) q[17],q[20];
rzz(pi/32) q[17],q[21];
rzz(0.049100328574073315) q[17],q[22];
rzz(pi/128) q[17],q[23];
u3(pi/2,2.9493271831900976,-2.9493271831900976) q[24];
rzz(0.012275082143518329) q[17],q[24];
u3(pi/2,5.221326990266236,-5.221326990266236) q[17];
rzz(0) q[17],q[25];
u3(pi,2.506990937564655,-2.506990937564655) q[17];
rzz(0) q[17],q[26];
u3(pi/2,2.8601059518281478,-2.8601059518281478) q[18];
rzz(0.7853830994606742) q[18],q[19];
rzz(pi/8) q[18],q[20];
rzz(0.19640131429629326) q[18],q[21];
rzz(pi/32) q[18],q[22];
rzz(0.049100328574073315) q[18],q[23];
rzz(pi/128) q[18],q[24];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[25];
rzz(0.012275082143518329) q[18],q[25];
u3(pi/2,4.430902278623044,-4.430902278623044) q[18];
rzz(0) q[18],q[26];
u3(pi/2,0.12126547642856603,-0.12126547642856603) q[19];
rzz(0.7853830994606742) q[19],q[20];
rzz(pi/8) q[19],q[21];
rzz(0.19640131429629326) q[19],q[22];
rzz(pi/32) q[19],q[23];
rzz(0.049100328574073315) q[19],q[24];
rzz(pi/128) q[19],q[25];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[26];
rzz(0.012275082143518329) q[19],q[26];
u3(pi/2,3.848451000647497,-3.848451000647497) q[20];
rzz(0.7853830994606742) q[20],q[21];
rzz(pi/8) q[20],q[22];
rzz(0.19640131429629326) q[20],q[23];
rzz(pi/32) q[20],q[24];
rzz(0.049100328574073315) q[20],q[25];
rzz(pi/128) q[20],q[26];
u3(pi/2,2.3806989128903453,-2.3806989128903453) q[21];
rzz(0.7853830994606742) q[21],q[22];
rzz(pi/8) q[21],q[23];
rzz(0.19640131429629326) q[21],q[24];
rzz(pi/32) q[21],q[25];
rzz(0.049100328574073315) q[21],q[26];
u3(pi/2,0.9544158481605792,-0.9544158481605792) q[22];
rzz(0.7853830994606742) q[22],q[23];
rzz(pi/8) q[22],q[24];
rzz(0.19640131429629326) q[22],q[25];
rzz(pi/32) q[22],q[26];
u3(pi/2,4.513840324677815,-4.513840324677815) q[23];
rzz(0.7853830994606742) q[23],q[24];
rzz(pi/8) q[23],q[25];
rzz(0.19640131429629326) q[23],q[26];
u3(pi/2,4.1575837177607315,-4.1575837177607315) q[24];
rzz(0.7853830994606742) q[24],q[25];
rzz(pi/8) q[24],q[26];
u3(pi/2,0.33489377687267197,-0.33489377687267197) q[25];
rzz(0.7853830994606742) q[25],q[26];
u3(pi/2,5.608371205188498,-5.608371205188498) q[1];
u3(pi/2,6.0274596651773775,-6.0274596651773775) q[2];
u3(pi/2,5.313689814281776,-5.313689814281776) q[3];
u3(pi/2,6.077096829104096,-6.077096829104096) q[4];
u3(pi/2,3.0781324819872795,-3.0781324819872795) q[5];
u3(pi/2,3.720902338911751,-3.720902338911751) q[7];
u3(pi/2,0.029530970943744055,-0.029530970943744055) q[9];
u3(pi/2,5.258397783578595,-5.258397783578595) q[11];
u3(pi/2,0.31855749507400505,-0.31855749507400505) q[13];
u3(pi/2,3.051743103697125,-3.051743103697125) q[15];
u3(pi/2,2.49002633723527,-2.49002633723527) q[16];
u3(pi/2,6.075211873511942,-6.075211873511942) q[17];
u3(pi/2,1.289309625033251,-1.289309625033251) q[18];
u3(pi/2,2.7803094984269667,-2.7803094984269667) q[26];
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
measure q[25] -> c[25];
measure q[26] -> c[26];
measure q[27] -> c[27];
