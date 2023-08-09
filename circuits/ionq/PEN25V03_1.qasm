OPENQASM 2.0;
include "qelib1.inc";
qreg q[25];
creg c[25];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[23];
u3(pi,1.288052987971815,-1.288052987971815) q[24];
rzz(1.208145008031432) q[23],q[24];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[22];
u3(pi,3.1931147731086655,-3.1931147731086655) q[24];
rzz(0.7250550171939733) q[22],q[24];
u3(pi/2,3.910654535188574,-3.910654535188574) q[21];
rzz(1.450361503171903) q[21],q[24];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[20];
u3(pi,5.908079144340965,-5.908079144340965) q[24];
rzz(0.2411286087284226) q[20],q[24];
u3(pi/2,1.4922565104551517,-1.4922565104551517) q[19];
rzz(0.4822765853755547) q[19],q[24];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[18];
rzz(0.9645144349136904) q[18],q[24];
u3(pi/2,4.386291662942069,-4.386291662942069) q[17];
u3(pi,3.007760806546868,-3.007760806546868) q[24];
rzz(1.2122813074589645) q[17],q[24];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[16];
u3(pi,3.790017377290726,-3.790017377290726) q[24];
rzz(0.716754756615593) q[16],q[24];
u3(pi/2,3.3948050214691303,-3.3948050214691303) q[15];
rzz(1.4335057276120384) q[15],q[24];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[14];
u3(pi,3.426849266535746,-3.426849266535746) q[24];
rzz(0.2744558174029115) q[14],q[24];
u3(pi/2,5.71267208128768,-5.71267208128768) q[13];
rzz(0.5488349956530386) q[13],q[24];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[12];
rzz(1.0975606481737692) q[12],q[24];
u3(pi/2,2.419026343264141,-2.419026343264141) q[11];
u3(pi,4.035689922801448,-4.035689922801448) q[24];
rzz(0.9463834868957338) q[11],q[24];
u3(pi/2,3.669380219392878,-3.669380219392878) q[10];
u3(pi,3.2597165373647696,-3.2597165373647696) q[24];
rzz(1.248705419631546) q[10],q[24];
u3(pi/2,1.8164688723056186,-1.8164688723056186) q[9];
u3(pi,1.5557166820576656,-1.5557166820576656) q[24];
rzz(0.6442026745019208) q[9],q[24];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[8];
rzz(1.2885124301894393) q[8],q[24];
u3(pi/2,4.883919939270692,-4.883919939270692) q[7];
u3(pi,2.382583868482499,-2.382583868482499) q[24];
rzz(0.5645445925491676) q[7],q[24];
u3(pi/2,3.730955435403238,-3.730955435403238) q[6];
rzz(1.1290331705013215) q[6],q[24];
u3(pi/2,3*pi/4,-3*pi/4) q[5];
u3(pi,4.691654468870997,-4.691654468870997) q[24];
rzz(0.8835729338221293) q[5],q[24];
u3(pi/2,5*pi/4,-5*pi/4) q[4];
u3(pi,4.207849200218169,-4.207849200218169) q[24];
rzz(1.3744046257721234) q[4],q[24];
u3(pi/2,pi/2,-pi/2) q[3];
u3(pi,2.017530802135365,-2.017530802135365) q[24];
rzz(pi/8) q[3],q[24];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(0.7853830994606742) q[2],q[24];
u3(pi/2,3*pi/2,-3*pi/2) q[1];
u3(pi,2.704911274740812,-2.704911274740812) q[24];
rzz(-pi/2) q[1],q[24];
rzz(0.7853830994606742) q[0],q[1];
rzz(pi/8) q[0],q[2];
rzz(0.19640131429629326) q[0],q[3];
rzz(pi/32) q[0],q[4];
rzz(0.049100328574073315) q[0],q[5];
rzz(pi/128) q[0],q[6];
rzz(0.012275082143518329) q[0],q[7];
u3(pi/2,pi,-pi) q[0];
u3(pi/2,0.5397256178867265,-0.5397256178867265) q[8];
rzz(0) q[0],q[8];
u3(pi,0.4266282823574939,-0.4266282823574939) q[0];
u3(pi/2,1.8164688723056186,-1.8164688723056186) q[9];
rzz(0) q[0],q[9];
u3(pi,4.422105819192993,-4.422105819192993) q[0];
u3(pi/2,3.669380219392878,-3.669380219392878) q[10];
rzz(0) q[0],q[10];
u3(pi,2.1343980488489054,-2.1343980488489054) q[0];
u3(pi/2,2.4196546617948584,-2.4196546617948584) q[11];
rzz(0) q[0],q[11];
u3(pi,6.129875585684404,-6.129875585684404) q[0];
u3(pi/2,3.6668669452700065,-3.6668669452700065) q[12];
rzz(0) q[0],q[12];
u3(pi,3.8421678153403174,-3.8421678153403174) q[0];
u3(pi/2,2.5717077462286047,-2.5717077462286047) q[13];
rzz(0) q[0],q[13];
u3(pi,1.5544600449962296,-1.5544600449962296) q[0];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[14];
rzz(0) q[0],q[14];
u3(pi,5.5499375818317285,-5.5499375818317285) q[0];
u3(pi/2,0.2532123678793373,-0.2532123678793373) q[15];
rzz(0) q[0],q[15];
u3(pi,3.262229811487641,-3.262229811487641) q[0];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[16];
rzz(0) q[0],q[16];
u3(pi,0.9745220411435538,-0.9745220411435538) q[0];
u3(pi/2,4.386291662942069,-4.386291662942069) q[17];
rzz(0) q[0],q[17];
u3(pi,4.969999577979053,-4.969999577979053) q[0];
u3(pi/2,3.666238626739289,-3.666238626739289) q[18];
rzz(0) q[0],q[18];
u3(pi,2.6822918076349653,-2.6822918076349653) q[0];
u3(pi/2,4.633849164044945,-4.633849164044945) q[19];
rzz(0) q[0],q[19];
u3(pi,0.394584037290878,-0.394584037290878) q[0];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[20];
rzz(0) q[0],q[20];
u3(pi,4.390061574126377,-4.390061574126377) q[0];
u3(pi/2,0.7690618815987813,-0.7690618815987813) q[21];
rzz(0) q[0],q[21];
u3(pi,2.1023538037822895,-2.1023538037822895) q[0];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[22];
rzz(0) q[0],q[22];
u3(pi/2,5.6705747397295765,-5.6705747397295765) q[0];
u3(pi/2,0.5246459731494955,-0.5246459731494955) q[23];
rzz(pi/2) q[0],q[23];
rzz(pi/2) q[0],q[23];
u3(pi/2,5*pi/4,-5*pi/4) q[1];
rzz(0.7853830994606742) q[1],q[2];
rzz(pi/8) q[1],q[3];
rzz(0.19640131429629326) q[1],q[4];
rzz(pi/32) q[1],q[5];
rzz(0.049100328574073315) q[1],q[6];
rzz(pi/128) q[1],q[7];
u3(pi/2,3.6813182714765196,-3.6813182714765196) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi/2,3*pi/4,-3*pi/4) q[1];
rzz(0) q[1],q[9];
u3(pi,5.924415426139632,-5.924415426139632) q[1];
rzz(0) q[1],q[10];
u3(pi,3.6367076557955444,-3.6367076557955444) q[1];
rzz(0) q[1],q[11];
u3(pi,1.348999885451457,-1.348999885451457) q[1];
rzz(0) q[1],q[12];
u3(pi,5.344477422286956,-5.344477422286956) q[1];
rzz(0) q[1],q[13];
u3(pi,3.0567696519428686,-3.0567696519428686) q[1];
rzz(0) q[1],q[14];
u3(pi,0.7690618815987813,-0.7690618815987813) q[1];
rzz(0) q[1],q[15];
u3(pi,4.76453941843428,-4.76453941843428) q[1];
rzz(0) q[1],q[16];
u3(pi,2.476831648090193,-2.476831648090193) q[1];
rzz(0) q[1],q[17];
u3(pi,0.18912387774610553,-0.18912387774610553) q[1];
rzz(0) q[1],q[18];
u3(pi,4.184601414581604,-4.184601414581604) q[1];
rzz(0) q[1],q[19];
u3(pi,1.896893644237517,-1.896893644237517) q[1];
rzz(0) q[1],q[20];
u3(pi,5.892371181073016,-5.892371181073016) q[1];
rzz(0) q[1],q[21];
u3(pi,3.6046634107289286,-3.6046634107289286) q[1];
rzz(0) q[1],q[22];
u3(pi/2,9*pi/8,-9*pi/8) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,4.954919933241822,-4.954919933241822) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi/2,5*pi/8,-5*pi/8) q[2];
rzz(0) q[2],q[10];
u3(pi,4.565362444196688,-4.565362444196688) q[2];
rzz(0) q[2],q[11];
u3(pi,0.34494687336415925,-0.34494687336415925) q[2];
rzz(0) q[2],q[12];
u3(pi/2,2.9468139090672256,-2.9468139090672256) q[2];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[13];
rzz(-pi/2) q[2],q[14];
rzz(pi/2) q[2],q[14];
rzz(pi/2) q[2],q[15];
rzz(pi/2) q[2],q[15];
rzz(-pi/2) q[2],q[16];
rzz(pi/2) q[2],q[16];
rzz(pi/2) q[2],q[17];
rzz(pi/2) q[2],q[17];
rzz(pi/2) q[2],q[18];
rzz(-pi/2) q[2],q[18];
rzz(pi/2) q[2],q[19];
rzz(pi/2) q[2],q[19];
rzz(pi/2) q[2],q[20];
rzz(-pi/2) q[2],q[20];
rzz(pi/2) q[2],q[21];
rzz(pi/2) q[2],q[21];
rzz(pi/2) q[2],q[22];
rzz(pi/2) q[2],q[22];
u3(pi/2,2.9468139090672256,-2.9468139090672256) q[2];
rzz(0) q[2],q[23];
u3(pi/2,0.19603538158400308,-0.19603538158400308) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,0.5233893360880595,-0.5233893360880595) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,1.7668317083788996,-1.7668317083788996) q[3];
rzz(0) q[3],q[11];
u3(pi,5.335680962856904,-5.335680962856904) q[3];
rzz(0) q[3],q[12];
u3(pi,3.047973192512817,-3.047973192512817) q[3];
rzz(0) q[3],q[13];
u3(pi,0.7602654221687299,-0.7602654221687299) q[3];
rzz(0) q[3],q[14];
u3(pi,4.755742959004229,-4.755742959004229) q[3];
rzz(0) q[3],q[15];
u3(pi,2.4680351886601413,-2.4680351886601413) q[3];
rzz(0) q[3],q[16];
u3(pi,0.18032741831605412,-0.18032741831605412) q[3];
rzz(0) q[3],q[17];
u3(pi,4.175804955151553,-4.175804955151553) q[3];
rzz(0) q[3],q[18];
u3(pi,1.8880971848074657,-1.8880971848074657) q[3];
rzz(0) q[3],q[19];
u3(pi,5.883574721642964,-5.883574721642964) q[3];
rzz(0) q[3],q[20];
u3(pi,3.5958669512988775,-3.5958669512988775) q[3];
rzz(0) q[3],q[21];
u3(pi,1.307530862424072,-1.307530862424072) q[3];
rzz(0) q[3],q[22];
u3(pi/2,4.876380116902077,-4.876380116902077) q[3];
rzz(-pi/2) q[3],q[23];
rzz(pi/2) q[3],q[23];
u3(pi/2,0.8834158541894498,-0.8834158541894498) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
u3(pi/2,4.883919939270692,-4.883919939270692) q[7];
rzz(pi/2) q[4],q[7];
u3(pi/2,3.313123612475796,-3.313123612475796) q[7];
u3(pi/2,6.258680884481586,-6.258680884481586) q[7];
rzz(pi/2) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,5.554335811546754,-5.554335811546754) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,2.4542121809843462,-2.4542121809843462) q[4];
rzz(0) q[4],q[12];
u3(pi,5.472654402553419,-5.472654402553419) q[4];
rzz(0) q[4],q[13];
u3(pi,2.0841325663914687,-2.0841325663914687) q[4];
rzz(0) q[4],q[14];
u3(pi/2,5.102574787960542,-5.102574787960542) q[4];
rzz(-pi/2) q[4],q[15];
rzz(pi/2) q[4],q[15];
rzz(pi/2) q[4],q[16];
rzz(pi/2) q[4],q[16];
rzz(-pi/2) q[4],q[17];
rzz(pi/2) q[4],q[17];
rzz(pi/2) q[4],q[18];
rzz(pi/2) q[4],q[18];
rzz(-pi/2) q[4],q[19];
rzz(pi/2) q[4],q[19];
rzz(pi/2) q[4],q[20];
rzz(pi/2) q[4],q[20];
rzz(pi/2) q[4],q[21];
rzz(pi/2) q[4],q[21];
rzz(pi/2) q[4],q[22];
rzz(pi/2) q[4],q[22];
u3(pi/2,1.9609821343707488,-1.9609821343707488) q[4];
rzz(0) q[4],q[23];
u3(pi/2,1.6198051721908973,-1.6198051721908973) q[5];
rzz(0.7853830994606742) q[5],q[6];
u3(pi/2,4.687884557686689,-4.687884557686689) q[7];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,3.190601498985794,-3.190601498985794) q[5];
rzz(0) q[5],q[13];
u3(pi,6.232291506191432,-6.232291506191432) q[5];
rzz(0) q[5],q[14];
u3(pi,2.890893559833328,-2.890893559833328) q[5];
rzz(0) q[5],q[15];
u3(pi,5.832052602124092,-5.832052602124092) q[5];
rzz(0) q[5],q[16];
u3(pi,2.490654655765988,-2.490654655765988) q[5];
rzz(0) q[5],q[17];
u3(pi,5.431813698056753,-5.431813698056753) q[5];
rzz(0) q[5],q[18];
u3(pi,2.0904157516986484,-2.0904157516986484) q[5];
rzz(0) q[5],q[19];
u3(pi,5.03220311252013,-5.03220311252013) q[5];
rzz(0) q[5],q[20];
u3(pi,1.6901768476313088,-1.6901768476313088) q[5];
rzz(0) q[5],q[21];
u3(pi,4.631964208452791,-4.631964208452791) q[5];
rzz(0) q[5],q[22];
u3(pi/2,1.3904689084788424,-1.3904689084788424) q[5];
rzz(pi/2) q[5],q[23];
rzz(-pi/2) q[5],q[23];
u3(pi/2,1.0065662862101699,-1.0065662862101699) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,2.565424560921425,-2.565424560921425) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,5.718955266594859,-5.718955266594859) q[6];
rzz(0) q[6],q[14];
u3(pi,2.6037519912952205,-2.6037519912952205) q[6];
rzz(0) q[6],q[15];
u3(pi,5.798751719996041,-5.798751719996041) q[6];
rzz(0) q[6],q[16];
u3(pi/2,2.684176763227119,-2.684176763227119) q[6];
rzz(pi/2) q[6],q[17];
rzz(pi/2) q[6],q[17];
rzz(pi/2) q[6],q[18];
rzz(-pi/2) q[6],q[18];
rzz(pi/2) q[6],q[19];
rzz(pi/2) q[6],q[19];
rzz(pi/2) q[6],q[20];
rzz(pi/2) q[6],q[20];
rzz(-pi/2) q[6],q[21];
rzz(pi/2) q[6],q[21];
rzz(pi/2) q[6],q[22];
rzz(pi/2) q[6],q[22];
u3(pi/2,2.684176763227119,-2.684176763227119) q[6];
rzz(0) q[6],q[23];
u3(pi/2,5.6818844732825,-5.6818844732825) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,3.660583759962827,-3.660583759962827) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,0.9694954928978101,-0.9694954928978101) q[7];
rzz(0) q[7],q[15];
u3(pi,4.537716428845097,-4.537716428845097) q[7];
rzz(0) q[7],q[16];
u3(pi,2.2500086585010095,-2.2500086585010095) q[7];
rzz(0) q[7],q[17];
u3(pi,6.245486195336508,-6.245486195336508) q[7];
rzz(0) q[7],q[18];
u3(pi,3.9577784249924215,-3.9577784249924215) q[7];
rzz(0) q[7],q[19];
u3(pi,1.670070654648334,-1.670070654648334) q[7];
rzz(0) q[7],q[20];
u3(pi,5.665548191483833,-5.665548191483833) q[7];
rzz(0) q[7],q[21];
u3(pi,3.3778404211397453,-3.3778404211397453) q[7];
rzz(0) q[7],q[22];
u3(pi/2,0.6635043684381643,-0.6635043684381643) q[7];
rzz(-pi/2) q[7],q[23];
rzz(pi/2) q[7],q[23];
u3(pi/2,3.387265199100515,-3.387265199100515) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,0.24692918257215776,-0.24692918257215776) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,1.8164688723056186,-1.8164688723056186) q[8];
rzz(0) q[8],q[16];
u3(pi,4.15507044363786,-4.15507044363786) q[8];
rzz(0) q[8],q[17];
u3(pi,5.691309251243269,-5.691309251243269) q[8];
rzz(0) q[8],q[18];
u3(pi/2,1.747353833926643,-1.747353833926643) q[8];
rzz(pi/2) q[8],q[19];
rzz(pi/2) q[8],q[19];
rzz(-pi/2) q[8],q[20];
rzz(pi/2) q[8],q[20];
rzz(pi/2) q[8],q[21];
rzz(pi/2) q[8],q[21];
rzz(-pi/2) q[8],q[22];
rzz(pi/2) q[8],q[22];
u3(pi/2,1.747353833926643,-1.747353833926643) q[8];
rzz(0) q[8],q[23];
u3(pi/2,0.8771326688822703,-0.8771326688822703) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,0.5189911063730339,-0.5189911063730339) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,2.4479289956771666,-2.4479289956771666) q[9];
rzz(0) q[9],q[17];
u3(pi,6.016778250155172,-6.016778250155172) q[9];
rzz(0) q[9],q[18];
u3(pi,3.7290704798110847,-3.7290704798110847) q[9];
rzz(0) q[9],q[19];
u3(pi,1.441362709466997,-1.441362709466997) q[9];
rzz(0) q[9],q[20];
u3(pi,5.436840246302496,-5.436840246302496) q[9];
rzz(0) q[9],q[21];
u3(pi,3.1491324759584085,-3.1491324759584085) q[9];
rzz(0) q[9],q[22];
u3(pi/2,0.4341681047261094,-0.4341681047261094) q[9];
rzz(pi/2) q[9],q[23];
rzz(pi/2) q[9],q[23];
u3(pi/2,0.8375486014470389,-0.8375486014470389) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,1.2384158240450964,-1.2384158240450964) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi/2,5.5499375818317285,-5.5499375818317285) q[10];
rzz(0) q[10],q[18];
u3(pi,2.834973210599429,-2.834973210599429) q[10];
rzz(0) q[10],q[19];
u3(pi,0.547265440255342,-0.547265440255342) q[10];
rzz(0) q[10],q[20];
u3(pi/2,4.116114694733347,-4.116114694733347) q[10];
rzz(pi/2) q[10],q[21];
rzz(-pi/2) q[10],q[21];
rzz(pi/2) q[10],q[22];
rzz(pi/2) q[10],q[22];
u3(pi/2,0.9745220411435538,-0.9745220411435538) q[10];
rzz(0) q[10],q[23];
u3(pi/2,1.7819113531161308,-1.7819113531161308) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,3.660583759962827,-3.660583759962827) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,3.352707679911027,-3.352707679911027) q[11];
rzz(0) q[11],q[19];
u3(pi,1.2088848531013523,-1.2088848531013523) q[11];
rzz(0) q[11],q[20];
u3(pi,0.06346017160251381,-0.06346017160251381) q[11];
rzz(0) q[11],q[21];
u3(pi,5.201220797283261,-5.201220797283261) q[11];
rzz(0) q[11],q[22];
u3(pi/2,3.0573979704735863,-3.0573979704735863) q[11];
rzz(pi/2) q[11],q[23];
rzz(-pi/2) q[11],q[23];
u3(pi/2,4.127424428286271,-4.127424428286271) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/128) q[12],q[18];
u3(pi/2,4.627565978737765,-4.627565978737765) q[19];
rzz(0.012275082143518329) q[12],q[19];
u3(pi/2,2.5566281014913734,-2.5566281014913734) q[12];
rzz(0) q[12],q[20];
u3(pi,6.1254773559693785,-6.1254773559693785) q[12];
rzz(0) q[12],q[21];
u3(pi,3.837769585625291,-3.837769585625291) q[12];
rzz(0) q[12],q[22];
u3(pi,1.550061815281204,-1.550061815281204) q[12];
rzz(0) q[12],q[23];
u3(pi/2,0.4391946529718531,-0.4391946529718531) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/128) q[13],q[19];
u3(pi/2,3.660583759962827,-3.660583759962827) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi/2,2.0099909797667497,-2.0099909797667497) q[13];
rzz(0) q[13],q[21];
u3(pi,5.578840234244755,-5.578840234244755) q[13];
rzz(0) q[13],q[22];
u3(pi,3.2911324639006674,-3.2911324639006674) q[13];
rzz(0) q[13],q[23];
u3(pi/2,4.951150022057514,-4.951150022057514) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/128) q[14],q[20];
u3(pi/2,0.7627786962916018,-0.7627786962916018) q[21];
rzz(0.012275082143518329) q[14],q[21];
u3(pi/2,3.3803536952626176,-3.3803536952626176) q[14];
rzz(0) q[14],q[22];
u3(pi,1.3345485592449442,-1.3345485592449442) q[14];
rzz(0) q[14],q[23];
u3(pi/2,0.10367255756846318,-0.10367255756846318) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/32) q[15],q[19];
rzz(0.049100328574073315) q[15],q[20];
rzz(pi/128) q[15],q[21];
u3(pi/2,3.660583759962827,-3.660583759962827) q[22];
rzz(0.012275082143518329) q[15],q[22];
u3(pi/2,1.6744688843633597,-1.6744688843633597) q[15];
rzz(0) q[15],q[23];
u3(pi/2,2.8010440099406595,-2.8010440099406595) q[16];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/8) q[16],q[18];
rzz(0.19640131429629326) q[16],q[19];
rzz(pi/32) q[16],q[20];
rzz(0.049100328574073315) q[16],q[21];
rzz(pi/128) q[16],q[22];
u3(pi/2,3.660583759962827,-3.660583759962827) q[23];
rzz(0.012275082143518329) q[16],q[23];
u3(pi/2,1.590902519777871,-1.590902519777871) q[17];
rzz(0.7853830994606742) q[17],q[18];
rzz(pi/8) q[17],q[19];
rzz(0.19640131429629326) q[17],q[20];
rzz(pi/32) q[17],q[21];
rzz(0.049100328574073315) q[17],q[22];
rzz(pi/128) q[17],q[23];
u3(pi/2,3.048601511043535,-3.048601511043535) q[18];
rzz(0.7853830994606742) q[18],q[19];
rzz(pi/8) q[18],q[20];
rzz(0.19640131429629326) q[18],q[21];
rzz(pi/32) q[18],q[22];
rzz(0.049100328574073315) q[18],q[23];
u3(pi/2,3.5330350982270815,-3.5330350982270815) q[19];
rzz(0.7853830994606742) q[19],q[20];
rzz(pi/8) q[19],q[21];
rzz(0.19640131429629326) q[19],q[22];
rzz(pi/32) q[19],q[23];
u3(pi/2,5.466999535776958,-5.466999535776958) q[20];
rzz(0.7853830994606742) q[20],q[21];
rzz(pi/8) q[20],q[22];
rzz(0.19640131429629326) q[20],q[23];
u3(pi/2,4.018725322472063,-4.018725322472063) q[21];
rzz(0.7853830994606742) q[21],q[22];
rzz(pi/8) q[21],q[23];
u3(pi/2,4.500645635532738,-4.500645635532738) q[22];
rzz(0.7853830994606742) q[22],q[23];
u3(pi/2,0.8896990394966294,-0.8896990394966294) q[1];
u3(pi/2,6.088406562657019,-6.088406562657019) q[2];
u3(pi/2,5.102574787960542,-5.102574787960542) q[4];
u3(pi/2,5.825769416816913,-5.825769416816913) q[6];
u3(pi/2,4.888946487516436,-4.888946487516436) q[8];
u3(pi/2,4.116114694733347,-4.116114694733347) q[10];
u3(pi/2,5.118282751228491,-5.118282751228491) q[12];
u3(pi/2,0.5761680926683681,-0.5761680926683681) q[13];
u3(pi/2,5.571928730406857,-5.571928730406857) q[14];
u3(pi/2,4.816061537953153,-4.816061537953153) q[15];
u3(pi/2,0.15079644737231007,-0.15079644737231007) q[23];
u3(pi,2.127486545011008,-2.127486545011008) q[24];
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
