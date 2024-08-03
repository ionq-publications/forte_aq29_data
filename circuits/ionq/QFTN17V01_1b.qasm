OPENQASM 2.0;
include "qelib1.inc";
qreg q[17];
creg c[17];
u(pi/2,2.364990949622396,-2.364990949622396) q[16];
u(pi/2,0.008796459430051584,-0.008796459430051584) q[16];
rzz(-pi/2) q[15],q[16];
u(pi/2,-1.5619998673648452,1.5619998673648452) q[16];
rzz(pi/8) q[14],q[16];
u(pi/2,0.4209734155810323,-0.4209734155810323) q[15];
rzz(0.7853830994606742) q[14],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/8) q[13],q[15];
u(pi/2,-1.5368671261361269,1.5368671261361269) q[14];
rzz(0.7853830994606742) q[13],q[14];
u(pi,3.2195041513988194,-3.2195041513988194) q[16];
rzz(pi/32) q[12],q[16];
u(pi,2.159530790077624,-2.159530790077624) q[15];
rzz(0.19640131429629326) q[12],q[15];
u(pi,3.5977519068910313,-3.5977519068910313) q[14];
rzz(pi/8) q[12],q[14];
u(pi/2,3.2540616705883076,-3.2540616705883076) q[13];
rzz(0.7853830994606742) q[12],q[13];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/32) q[11],q[15];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/8) q[11],q[13];
u(pi/2,-1.4344512056290994,1.4344512056290994) q[12];
rzz(0.7853830994606742) q[11],q[12];
u(pi,1.3131857292005336,-1.3131857292005336) q[16];
rzz(pi/128) q[10],q[16];
u(pi,-1.0423804424610934,1.0423804424610934) q[15];
rzz(0.049100328574073315) q[10],q[15];
u(pi,1.7988759534455157,-1.7988759534455157) q[14];
rzz(pi/32) q[10],q[14];
u(pi,1.9546989490635696,-1.9546989490635696) q[13];
rzz(0.19640131429629326) q[10],q[13];
u(pi,1.1586193706439158,-1.1586193706439158) q[12];
rzz(pi/8) q[10],q[12];
u(pi/2,-1.1196636217394023,1.1196636217394023) q[11];
rzz(0.7853830994606742) q[10],q[11];
rzz(0.012275082143518329) q[9],q[16];
rzz(pi/128) q[9],q[15];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/32) q[9],q[13];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/8) q[9],q[11];
u(pi/2,2.1168051299888027,-2.1168051299888027) q[10];
rzz(0.7853830994606742) q[9],q[10];
u(pi/2,0.6138672045114455,-0.6138672045114455) q[8];
u(pi/2,4.050141249007961,-4.050141249007961) q[16];
rzz(0) q[8],q[16];
u(pi/2,3.755459858101239,-3.755459858101239) q[8];
rzz(0.012275082143518329) q[8],q[15];
rzz(pi/128) q[8],q[14];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/32) q[8],q[12];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/8) q[8],q[10];
u(pi/2,3.7309554354032386,-3.7309554354032386) q[9];
rzz(0.7853830994606742) q[8],q[9];
u(pi/2,-0.638371627209446,0.638371627209446) q[7];
rzz(0) q[7],q[16];
u(pi,1.2905662620946874,-1.2905662620946874) q[7];
u(pi/2,3.441928911272977,-3.441928911272977) q[15];
rzz(0) q[7],q[15];
u(pi/2,3.218875832868102,-3.218875832868102) q[7];
rzz(0.012275082143518329) q[7],q[14];
rzz(pi/128) q[7],q[13];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/32) q[7],q[11];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/8) q[7],q[9];
u(pi/2,3.755459858101239,-3.755459858101239) q[8];
rzz(0.7853830994606742) q[7],q[8];
u(pi/2,0.8834158541894497,-0.8834158541894497) q[6];
rzz(0) q[6],q[16];
u(pi,4.194026192542374,-4.194026192542374) q[6];
rzz(0) q[6],q[15];
u(pi,1.3904689084788422,-1.3904689084788422) q[6];
u(pi/2,4.290158927742222,-4.290158927742222) q[14];
rzz(0) q[6],q[14];
u(pi/2,4.701079246831767,-4.701079246831767) q[6];
u(pi,1.3113007736083797,-1.3113007736083797) q[13];
rzz(0.012275082143518329) q[6],q[13];
u(pi,-0.7250795844485243,0.7250795844485243) q[12];
rzz(pi/128) q[6],q[12];
u(pi,2.1909467166135217,-2.1909467166135217) q[11];
rzz(0.049100328574073315) q[6],q[11];
u(pi,0.8758760318208343,-0.8758760318208343) q[10];
rzz(pi/32) q[6],q[10];
u(pi,2.4309643953477824,-2.4309643953477824) q[9];
rzz(0.19640131429629326) q[6],q[9];
u(pi,2.455468818045783,-2.455468818045783) q[8];
rzz(pi/8) q[6],q[8];
u(pi/2,0.07728317927830886,-0.07728317927830886) q[7];
rzz(0.7853830994606742) q[6],q[7];
u(pi/2,2.1601591086083416,-2.1601591086083416) q[5];
rzz(0) q[5],q[16];
u(pi,-0.8130441787490384,0.8130441787490384) q[5];
rzz(0) q[5],q[15];
u(pi,2.6665838443670165,-2.6665838443670165) q[5];
rzz(0) q[5],q[14];
u(pi,-0.136973439696515,0.136973439696515) q[5];
u(pi/2,-1.1749556524425826,1.1749556524425826) q[13];
rzz(0) q[5],q[13];
u(pi/2,3.173636898656409,-3.173636898656409) q[5];
u(pi,-0.1589645882716435,0.1589645882716435) q[12];
rzz(0.012275082143518329) q[5],q[12];
u(pi,-0.2469291825721578,0.2469291825721578) q[11];
rzz(pi/128) q[5],q[11];
u(pi,1.2296193646150448,-1.2296193646150448) q[10];
rzz(0.049100328574073315) q[5],q[10];
u(pi,2.838114803253019,-2.838114803253019) q[9];
rzz(pi/32) q[5],q[9];
u(pi,2.7432387051146074,-2.7432387051146074) q[8];
rzz(0.19640131429629326) q[5],q[8];
u(pi,2.4429024474314236,-2.4429024474314236) q[7];
rzz(pi/8) q[5],q[7];
u(pi/2,4.701079246831767,-4.701079246831767) q[6];
rzz(0.7853830994606742) q[5],q[6];
u(pi/2,5*pi/8,-5*pi/8) q[4];
rzz(0) q[4],q[16];
u(pi,4.7004509283010485,-4.7004509283010485) q[4];
rzz(0) q[4],q[15];
u(pi,0.7495840071465247,-0.7495840071465247) q[4];
rzz(0) q[4],q[14];
u(pi,3.0825307117023053,-3.0825307117023053) q[4];
rzz(0) q[4],q[13];
u(pi,-0.8683362094522188,0.8683362094522188) q[4];
u(pi/2,-1.3992653679088938,1.3992653679088938) q[12];
rzz(0) q[4],q[12];
u(pi/2,1.868619310355209,-1.868619310355209) q[4];
rzz(0.012275082143518329) q[4],q[11];
rzz(pi/128) q[4],q[10];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/32) q[4],q[8];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/8) q[4],q[6];
u(pi/2,3.173636898656409,-3.173636898656409) q[5];
rzz(0.7853830994606742) q[4],q[5];
u(pi/2,pi/4,-pi/4) q[3];
rzz(0) q[3],q[16];
u(pi,4.033176648678577,-4.033176648678577) q[3];
rzz(0) q[3],q[15];
u(pi,1.1045839770021715,-1.1045839770021715) q[3];
rzz(0) q[3],q[14];
u(pi,4.458548293974634,-4.458548293974634) q[3];
rzz(0) q[3],q[13];
u(pi,1.5299556222982291,-1.5299556222982291) q[3];
rzz(0) q[3],q[12];
u(pi,-1.3992653679088938,1.3992653679088938) q[3];
u(pi/2,3.4293625406586177,-3.4293625406586177) q[11];
rzz(0) q[3],q[11];
u(pi/2,1.8485131173722342,-1.8485131173722342) q[3];
u(pi,-0.6358583530865742,0.6358583530865742) q[10];
rzz(0.012275082143518329) q[3],q[10];
u(pi,0.18472564803107971,-0.18472564803107971) q[9];
rzz(pi/128) q[3],q[9];
u(pi,0.8702211650443727,-0.8702211650443727) q[8];
rzz(0.049100328574073315) q[3],q[8];
u(pi,0.136973439696515,-0.136973439696515) q[7];
rzz(pi/32) q[3],q[7];
u(pi,0.8136724972797564,-0.8136724972797564) q[6];
rzz(0.19640131429629326) q[3],q[6];
u(pi,1.8742741771316704,-1.8742741771316704) q[5];
rzz(pi/8) q[3],q[5];
u(pi/2,-1.2729733432345842,1.2729733432345842) q[4];
rzz(0.7853830994606742) q[3],q[4];
u(pi/2,0,0) q[2];
rzz(0) q[2],q[16];
u(pi,2.823035158515788,-2.823035158515788) q[2];
rzz(0) q[2],q[15];
u(pi,-0.9550441666912971,0.9550441666912971) q[2];
rzz(0) q[2],q[14];
u(pi,1.5500618152812038,-1.5500618152812038) q[2];
rzz(0) q[2],q[13];
u(pi,4.055167797253705,-4.055167797253705) q[2];
rzz(0) q[2],q[12];
u(pi,0.27708847204661957,-0.27708847204661957) q[2];
rzz(0) q[2],q[11];
u(pi,2.7815661354884025,-2.7815661354884025) q[2];
u(pi/2,2.1871768054292136,-2.1871768054292136) q[10];
rzz(0) q[2],q[10];
u(pi/2,-0.6779556946446773,0.6779556946446773) q[2];
rzz(0.012275082143518329) q[2],q[9];
rzz(pi/128) q[2],q[8];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/32) q[2],q[6];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/8) q[2],q[4];
u(pi/2,-1.293079536217559,1.293079536217559) q[3];
rzz(0.7853830994606742) q[2],q[3];
u(pi/2,pi/2,-pi/2) q[1];
rzz(0) q[1],q[16];
u(pi,-1.4017786420317657,1.4017786420317657) q[1];
rzz(0) q[1],q[15];
u(pi,2.077849381084289,-2.077849381084289) q[1];
rzz(0) q[1],q[14];
u(pi,-0.7257079029792421,0.7257079029792421) q[1];
rzz(0) q[1],q[13];
u(pi,2.7539201201368124,-2.7539201201368124) q[1];
rzz(0) q[1],q[12];
u(pi,-0.049637163926718575,0.049637163926718575) q[1];
rzz(0) q[1],q[11];
u(pi,3.4299908591893367,-3.4299908591893367) q[1];
rzz(0) q[1],q[10];
u(pi,0.6264335751258048,-0.6264335751258048) q[1];
u(pi/2,2.1080086705587515,-2.1080086705587515) q[9];
rzz(0) q[1],q[9];
u(pi/2,3.937043913478729,-3.937043913478729) q[1];
rzz(0.012275082143518329) q[1],q[8];
rzz(pi/128) q[1],q[7];
rzz(0.049100328574073315) q[1],q[6];
rzz(pi/32) q[1],q[5];
rzz(0.19640131429629326) q[1],q[4];
rzz(pi/8) q[1],q[3];
u(pi/2,2.4636369589451155,-2.4636369589451155) q[2];
rzz(0.7853830994606742) q[1],q[2];
rzz(-pi/2) q[0],q[16];
rzz(-pi/2) q[0],q[15];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[11];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[9];
u(pi/2,3.693884642090879,-3.693884642090879) q[8];
rzz(pi/2) q[0],q[8];
u(pi/2,-1.3929821826017144,1.3929821826017144) q[7];
rzz(pi/2) q[0],q[7];
u(pi/2,3.2088227363766153,-3.2088227363766153) q[6];
rzz(-pi/2) q[0],q[6];
u(pi/2,0.5749114556069324,-0.5749114556069324) q[5];
rzz(pi/2) q[0],q[5];
u(pi/2,1.868619310355209,-1.868619310355209) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,1.8485131173722342,-1.8485131173722342) q[3];
rzz(pi/2) q[0],q[3];
u(pi/2,-0.6779556946446773,0.6779556946446773) q[2];
rzz(pi/2) q[0],q[2];
rzz(-pi/2) q[0],q[1];
u(pi,pi,-pi) q[0];
u(pi/2,3.937043913478729,-3.937043913478729) q[1];
rzz(-pi/2) q[0],q[1];
u(pi/2,4.034433285740012,-4.034433285740012) q[2];
u(pi/2,1.6782387955476676,-1.6782387955476676) q[2];
rzz(pi/2) q[0],q[2];
u(pi/2,0.2777167905773379,-0.2777167905773379) q[3];
u(pi/2,3.8120085258658554,-3.8120085258658554) q[3];
rzz(-pi/2) q[0],q[3];
u(pi/2,0.29782298356031234,-0.29782298356031234) q[4];
u(pi/2,3.6354510187341083,-3.6354510187341083) q[4];
rzz(pi/2) q[0],q[4];
u(pi/2,-0.9958848711879644,0.9958848711879644) q[5];
u(pi/2,2.244353791724548,-2.244353791724548) q[5];
rzz(-pi/2) q[0],q[5];
u(pi/2,-1.503566244008075,1.503566244008075) q[6];
u(pi/2,1.6870352549777188,-1.6870352549777188) q[6];
rzz(pi/2) q[0],q[6];
u(pi/2,3.3194067977829755,-3.3194067977829755) q[7];
u(pi/2,0.20231856689118266,-0.20231856689118266) q[7];
rzz(pi/2) q[0],q[7];
u(pi/2,2.123088315295982,-2.123088315295982) q[8];
u(pi/2,-1.0059379676794518,1.0059379676794518) q[8];
rzz(pi/2) q[0],q[8];
rzz(pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[10];
rzz(pi/2) q[0],q[11];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[14];
rzz(-pi/2) q[0],q[15];
rzz(pi/2) q[0],q[16];
u(pi/2,-0.7753450669059611,0.7753450669059611) q[1];
u(pi/2,0.10744246875277086,-0.10744246875277086) q[2];
rzz(0.7853830994606742) q[1],q[2];
u(pi/2,-0.9003804545188346,0.9003804545188346) q[3];
rzz(pi/8) q[1],q[3];
u(pi/2,2.064654691939212,-2.064654691939212) q[4];
rzz(0.19640131429629326) q[1],q[4];
u(pi/2,3.815150118519444,-3.815150118519444) q[5];
rzz(pi/32) q[1],q[5];
u(pi/2,0.11623892818282244,-0.11623892818282244) q[6];
rzz(0.049100328574073315) q[1],q[6];
u(pi/2,-1.368477759903714,1.368477759903714) q[7];
rzz(pi/128) q[1],q[7];
u(pi/2,3.706451012705238,-3.706451012705238) q[8];
rzz(0.012275082143518329) q[1],q[8];
u(pi/2,-0.7753450669059611,0.7753450669059611) q[1];
u(pi/2,-0.2607521902479528,0.2607521902479528) q[9];
u(pi/2,2.889636922771891,-2.889636922771891) q[9];
rzz(0) q[1],q[9];
u(pi,2.5352652714469626,-2.5352652714469626) q[1];
rzz(0) q[1],q[10];
u(pi,-0.2682920126165682,0.2682920126165682) q[1];
rzz(0) q[1],q[11];
u(pi,3.211336010499487,-3.211336010499487) q[1];
rzz(0) q[1],q[12];
u(pi,0.40777872643595514,-0.40777872643595514) q[1];
rzz(0) q[1],q[13];
u(pi,3.88740674955201,-3.88740674955201) q[1];
rzz(0) q[1],q[14];
u(pi,1.0838494654884787,-1.0838494654884787) q[1];
rzz(0) q[1],q[15];
u(pi,4.5634774886045335,-4.5634774886045335) q[1];
rzz(0) q[1],q[16];
u(pi/2,2.4636369589451155,-2.4636369589451155) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u(pi/2,2.1042387593744434,-2.1042387593744434) q[9];
rzz(0.012275082143518329) q[2],q[9];
u(pi/2,4.034433285740012,-4.034433285740012) q[2];
rzz(0) q[2],q[10];
u(pi,1.548805178219768,-1.548805178219768) q[2];
rzz(0) q[2],q[11];
u(pi,-0.2814867017616456,0.2814867017616456) q[2];
rzz(0) q[2],q[12];
u(pi,4.172035043967245,-4.172035043967245) q[2];
rzz(0) q[2],q[13];
u(pi,2.341743163985832,-2.341743163985832) q[2];
rzz(0) q[2],q[14];
u(pi,0.5120796025351364,-0.5120796025351364) q[2];
rzz(0) q[2],q[15];
u(pi,-1.3182122774462772,1.3182122774462772) q[2];
rzz(0) q[2],q[16];
u(pi/2,1.8485131173722342,-1.8485131173722342) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u(pi/2,-0.9694954928978102,0.9694954928978102) q[10];
rzz(0.012275082143518329) q[3],q[10];
u(pi/2,0.2777167905773379,-0.2777167905773379) q[3];
rzz(0) q[3],q[11];
u(pi,3.101380267623844,-3.101380267623844) q[3];
rzz(0) q[3],q[12];
u(pi,-0.6773273761139594,0.6773273761139594) q[3];
rzz(0) q[3],q[13];
u(pi,1.827778605858542,-1.827778605858542) q[3];
rzz(0) q[3],q[14];
u(pi,4.3328845878310425,-4.3328845878310425) q[3];
rzz(0) q[3],q[15];
u(pi,0.5548052626239572,-0.5548052626239572) q[3];
rzz(0) q[3],q[16];
u(pi/2,1.8679909918244912,-1.8679909918244912) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u(pi/2,3.415539532982823,-3.415539532982823) q[11];
rzz(0.012275082143518329) q[4],q[11];
u(pi/2,3.4387873186193874,-3.4387873186193874) q[4];
rzz(0) q[4],q[12];
u(pi,0.40338049672092935,-0.40338049672092935) q[4];
rzz(0) q[4],q[13];
u(pi,3.7579731322241106,-3.7579731322241106) q[4];
rzz(0) q[4],q[14];
u(pi,0.8287521420169877,-0.8287521420169877) q[4];
rzz(0) q[4],q[15];
u(pi,4.183344777520168,-4.183344777520168) q[4];
rzz(0) q[4],q[16];
u(pi/2,3.717132427727443,-3.717132427727443) q[5];
rzz(0.7853830994606742) q[5],q[6];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u(pi/2,1.7241060482900785,-1.7241060482900785) q[12];
rzz(0.012275082143518329) q[5],q[12];
u(pi/2,2.1463361009325466,-2.1463361009325466) q[5];
rzz(0) q[5],q[13];
u(pi,0.0948760981384118,-0.0948760981384118) q[5];
rzz(0) q[5],q[14];
u(pi,-0.8664512538600649,0.8664512538600649) q[5];
rzz(0) q[5],q[15];
u(pi,4.456035019851763,-4.456035019851763) q[5];
rzz(0) q[5],q[16];
u(pi/2,3.2088227363766153,-3.2088227363766153) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u(pi/2,-1.1919202527719674,1.1919202527719674) q[13];
rzz(0.012275082143518329) q[6],q[13];
u(pi/2,-1.503566244008075,1.503566244008075) q[6];
rzz(0) q[6],q[14];
u(pi,1.2962211288711485,-1.2962211288711485) q[6];
rzz(0) q[6],q[15];
u(pi,3.7529465839783676,-3.7529465839783676) q[6];
rzz(0) q[6],q[16];
u(pi/2,-1.3929821826017144,1.3929821826017144) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u(pi/2,1.1303450367616077,-1.1303450367616077) q[14];
rzz(0.012275082143518329) q[7],q[14];
u(pi/2,3.3194067977829755,-3.3194067977829755) q[7];
rzz(0) q[7],q[15];
u(pi,0.34683182895631326,-0.34683182895631326) q[7];
rzz(0) q[7],q[16];
u(pi/2,3.694512960621596,-3.694512960621596) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u(pi/2,0.28211502029236346,-0.28211502029236346) q[15];
rzz(0.012275082143518329) q[8],q[15];
u(pi/2,-1.0178760197630932,1.0178760197630932) q[8];
rzz(0) q[8],q[16];
u(pi/2,2.092300707290802,-2.092300707290802) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u(pi/2,4.0319200116171405,-4.0319200116171405) q[16];
rzz(0.012275082143518329) q[9],q[16];
u(pi/2,-0.9789202708585795,0.9789202708585795) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u(pi/2,0.2657787384936965,-0.2657787384936965) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
u(pi/2,-1.4243981091376121,1.4243981091376121) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
u(pi/2,-1.198831756609865,1.198831756609865) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
u(pi/2,1.124061851454428,-1.124061851454428) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
u(pi/2,0.2758318349851838,-0.2758318349851838) q[15];
rzz(0.7853830994606742) q[15],q[16];
u(pi,pi/2,-pi/2) q[0];
u(pi/2,1.5909025197778712,-1.5909025197778712) q[1];
u(pi/2,2.4793449222130644,-2.4793449222130644) q[2];
u(pi/2,3.3778404211397453,-3.3778404211397453) q[3];
u(pi/2,1.1479379556217104,-1.1479379556217104) q[4];
u(pi/2,2.404575017057628,-2.404575017057628) q[5];
u(pi/2,0.2695486496780042,-0.2695486496780042) q[6];
u(pi/2,3.657442167309237,-3.657442167309237) q[7];
u(pi/2,2.1237166338267,-2.1237166338267) q[8];
u(pi/2,4.025636826309961,-4.025636826309961) q[16];
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