OPENQASM 2.0;
include "qelib1.inc";
qreg q[22];
creg c[22];
u3(pi/2,3.7831058734528287,-3.7831058734528287) q[21];
u3(pi/2,6.139300363645174,-6.139300363645174) q[21];
rzz(pi/2) q[20],q[21];
u3(pi/2,1.426911383260484,-1.426911383260484) q[21];
rzz(pi/8) q[19],q[21];
u3(pi/2,5.335680962856904,-5.335680962856904) q[20];
rzz(0.7853830994606742) q[19],q[20];
rzz(0.19640131429629326) q[18],q[21];
rzz(pi/8) q[18],q[20];
u3(pi/2,5.708273851572654,-5.708273851572654) q[19];
rzz(0.7853830994606742) q[18],q[19];
u3(pi,0.22556635252774715,-0.22556635252774715) q[21];
rzz(pi/32) q[17],q[21];
u3(pi,1.293079536217559,-1.293079536217559) q[20];
rzz(0.19640131429629326) q[17],q[20];
u3(pi,4.3002120242337085,-4.3002120242337085) q[19];
rzz(pi/8) q[17],q[19];
u3(pi/2,5.634760583478653,-5.634760583478653) q[18];
rzz(0.7853830994606742) q[17],q[18];
u3(pi,0.6358583530865741,-0.6358583530865741) q[21];
rzz(0.049100328574073315) q[16],q[21];
u3(pi,2.55160155324563,-2.55160155324563) q[20];
rzz(pi/32) q[16],q[20];
u3(pi,4.701707565362484,-4.701707565362484) q[19];
rzz(0.19640131429629326) q[16],q[19];
u3(pi,1.6072388015765384,-1.6072388015765384) q[18];
rzz(pi/8) q[16],q[18];
u3(pi/2,0.8432034682235006,-0.8432034682235006) q[17];
rzz(0.7853830994606742) q[16],q[17];
u3(pi,0.8727344391672445,-0.8727344391672445) q[21];
rzz(pi/128) q[15],q[21];
u3(pi,3.7504333098554947,-3.7504333098554947) q[20];
rzz(0.049100328574073315) q[15],q[20];
u3(pi,5.073672135547516,-5.073672135547516) q[19];
rzz(pi/32) q[15],q[19];
u3(pi,2.902203293386251,-2.902203293386251) q[18];
rzz(0.19640131429629326) q[15],q[18];
u3(pi,5.712043762756962,-5.712043762756962) q[17];
rzz(pi/8) q[15],q[17];
u3(pi/2,0.54789375878606,-0.54789375878606) q[16];
rzz(0.7853830994606742) q[15],q[16];
rzz(0.012275082143518329) q[14],q[21];
rzz(pi/128) q[14],q[20];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/32) q[14],q[18];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/8) q[14],q[16];
u3(pi/2,3.3721855543632837,-3.3721855543632837) q[15];
rzz(0.7853830994606742) q[14],q[15];
u3(pi/2,0.9217432845632453,-0.9217432845632453) q[13];
u3(pi/2,5.781158801135938,-5.781158801135938) q[21];
rzz(0) q[13],q[21];
u3(pi/2,4.063335938153039,-4.063335938153039) q[13];
rzz(0.012275082143518329) q[13],q[20];
rzz(pi/128) q[13],q[19];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/32) q[13],q[17];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/8) q[13],q[15];
u3(pi/2,5.334424325795468,-5.334424325795468) q[14];
rzz(0.7853830994606742) q[13],q[14];
u3(pi/2,2.4881413816431164,-2.4881413816431164) q[12];
rzz(0) q[12],q[21];
u3(pi,4.8267429529753585,-4.8267429529753585) q[12];
u3(pi/2,5.930698611446811,-5.930698611446811) q[20];
rzz(0) q[12],q[20];
u3(pi/2,0.8827875356587319,-0.8827875356587319) q[12];
rzz(0.012275082143518329) q[12],q[19];
rzz(pi/128) q[12],q[18];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/32) q[12],q[16];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/8) q[12],q[14];
u3(pi/2,0.9217432845632453,-0.9217432845632453) q[13];
rzz(0.7853830994606742) q[12],q[13];
u3(pi/2,6.019291524278043,-6.019291524278043) q[11];
rzz(0) q[11],q[21];
u3(pi,2.7545484386675305,-2.7545484386675305) q[11];
rzz(0) q[11],q[20];
u3(pi,5.649211909685166,-5.649211909685166) q[11];
u3(pi/2,3.6360793372648263,-3.6360793372648263) q[19];
rzz(0) q[11],q[19];
u3(pi/2,2.383840505543935,-2.383840505543935) q[11];
u3(pi,0.5786813667912399,-0.5786813667912399) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi,3.8685571936304712,-3.8685571936304712) q[17];
rzz(pi/128) q[11],q[17];
u3(pi,5.578211915714037,-5.578211915714037) q[16];
rzz(0.049100328574073315) q[11],q[16];
u3(pi,0.6572211831309848,-0.6572211831309848) q[15];
rzz(pi/32) q[11],q[15];
u3(pi,3.920079313149344,-3.920079313149344) q[14];
rzz(0.19640131429629326) q[11],q[14];
u3(pi,3.9401855061323183,-3.9401855061323183) q[13];
rzz(pi/8) q[11],q[13];
u3(pi/2,0.8827875356587319,-0.8827875356587319) q[12];
rzz(0.7853830994606742) q[11],q[12];
u3(pi/2,3.669380219392878,-3.669380219392878) q[10];
rzz(0) q[10],q[21];
u3(pi,0.5541769440932395,-0.5541769440932395) q[10];
rzz(0) q[10],q[20];
u3(pi,3.749176672794059,-3.749176672794059) q[10];
rzz(0) q[10],q[19];
u3(pi,0.6609910943152925,-0.6609910943152925) q[10];
u3(pi/2,5.499672099374291,-5.499672099374291) q[18];
rzz(0) q[10],q[18];
u3(pi/2,3.829601444725958,-3.829601444725958) q[10];
rzz(0.012275082143518329) q[10],q[17];
rzz(pi/128) q[10],q[16];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/32) q[10],q[14];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/8) q[10],q[12];
u3(pi/2,5.525433159133728,-5.525433159133728) q[11];
rzz(0.7853830994606742) q[10],q[11];
u3(pi/2,2.184663531306342,-2.184663531306342) q[9];
rzz(0) q[9],q[21];
u3(pi,0.21614157456697777,-0.21614157456697777) q[9];
rzz(0) q[9],q[20];
u3(pi,5.7038756218576285,-5.7038756218576285) q[9];
rzz(0) q[9],q[19];
u3(pi,4.909052680499411,-4.909052680499411) q[9];
rzz(0) q[9],q[18];
u3(pi,4.1136014206104745,-4.1136014206104745) q[9];
u3(pi/2,0.2978229835603124,-0.2978229835603124) q[17];
rzz(0) q[9],q[17];
u3(pi/2,2.145707782401829,-2.145707782401829) q[9];
u3(pi,4.200937696380271,-4.200937696380271) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi,0.12692034320502762,-0.12692034320502762) q[15];
rzz(pi/128) q[9],q[15];
u3(pi,2.5327519973240915,-2.5327519973240915) q[14];
rzz(0.049100328574073315) q[9],q[14];
u3(pi,2.6621856146519907,-2.6621856146519907) q[13];
rzz(pi/32) q[9],q[13];
u3(pi,3.901229757227805,-3.901229757227805) q[12];
rzz(0.19640131429629326) q[9],q[12];
u3(pi,1.8447432061879268,-1.8447432061879268) q[11];
rzz(pi/8) q[9],q[11];
u3(pi/2,0.6880087911361646,-0.6880087911361646) q[10];
rzz(0.7853830994606742) q[9],q[10];
u3(pi/2,2.110521944681623,-2.110521944681623) q[8];
rzz(0) q[8],q[21];
u3(pi,5.679371199159628,-5.679371199159628) q[8];
rzz(0) q[8],q[20];
u3(pi,3.3916634288155403,-3.3916634288155403) q[8];
rzz(0) q[8],q[19];
u3(pi,1.1039556584714532,-1.1039556584714532) q[8];
rzz(0) q[8],q[18];
u3(pi,5.099433195306952,-5.099433195306952) q[8];
rzz(0) q[8],q[17];
u3(pi,2.811725424962865,-2.811725424962865) q[8];
u3(pi/2,0.9361946107697583,-0.9361946107697583) q[16];
rzz(0) q[8],q[16];
u3(pi/2,0.09676105373056564,-0.09676105373056564) q[8];
rzz(0.012275082143518329) q[8],q[15];
rzz(pi/128) q[8],q[14];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/32) q[8],q[12];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/8) q[8],q[10];
u3(pi/2,5.287300435991622,-5.287300435991622) q[9];
rzz(0.7853830994606742) q[8],q[9];
u3(pi/2,2.4542121809843462,-2.4542121809843462) q[7];
rzz(0) q[7],q[21];
u3(pi,6.023061435462352,-6.023061435462352) q[7];
rzz(0) q[7],q[20];
u3(pi,3.7353536651182644,-3.7353536651182644) q[7];
rzz(0) q[7],q[19];
u3(pi,1.4476458947741766,-1.4476458947741766) q[7];
rzz(0) q[7],q[18];
u3(pi,5.443123431609675,-5.443123431609675) q[7];
rzz(0) q[7],q[17];
u3(pi,3.155415661265588,-3.155415661265588) q[7];
rzz(0) q[7],q[16];
u3(pi,0.8670795723907829,-0.8670795723907829) q[7];
u3(pi/2,5.453176528101163,-5.453176528101163) q[15];
rzz(0) q[7],q[15];
u3(pi/2,4.435928826868787,-4.435928826868787) q[7];
rzz(0.012275082143518329) q[7],q[14];
rzz(pi/128) q[7],q[13];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/32) q[7],q[11];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/8) q[7],q[9];
u3(pi/2,3.2383537073203588,-3.2383537073203588) q[8];
rzz(0.7853830994606742) q[7],q[8];
u3(pi/2,2.159530790077624,-2.159530790077624) q[6];
rzz(0) q[6],q[21];
u3(pi,5.1779730116466975,-5.1779730116466975) q[6];
rzz(0) q[6],q[20];
u3(pi,1.7894511754847462,-1.7894511754847462) q[6];
rzz(0) q[6],q[19];
u3(pi,4.6847429650331,-4.6847429650331) q[6];
rzz(0) q[6],q[18];
u3(pi,1.2962211288711487,-1.2962211288711487) q[6];
rzz(0) q[6],q[17];
u3(pi,4.1908845998887845,-4.1908845998887845) q[6];
rzz(0) q[6],q[16];
u3(pi,0.8023627637268332,-0.8023627637268332) q[6];
rzz(0) q[6],q[15];
u3(pi,3.697026234744469,-3.697026234744469) q[6];
u3(pi/2,5.701362347734756,-5.701362347734756) q[14];
rzz(0) q[6],q[14];
u3(pi/2,0.43228314913395555,-0.43228314913395555) q[6];
u3(pi,3.691371367968007,-3.691371367968007) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi,5.936981796753991,-5.936981796753991) q[12];
rzz(pi/128) q[6],q[12];
u3(pi,4.323459809870274,-4.323459809870274) q[11];
rzz(0.049100328574073315) q[6],q[11];
u3(pi,3.289875826839231,-3.289875826839231) q[10];
rzz(pi/32) q[6],q[10];
u3(pi,5.539884485340242,-5.539884485340242) q[9];
rzz(0.19640131429629326) q[6],q[9];
u3(pi,6.280043714525997,-6.280043714525997) q[8];
rzz(pi/8) q[6],q[8];
u3(pi/2,4.435928826868787,-4.435928826868787) q[7];
rzz(0.7853830994606742) q[6],q[7];
u3(pi/2,9*pi/8,-9*pi/8) q[5];
rzz(0) q[5],q[21];
u3(pi,0.26954864967800424,-0.26954864967800424) q[5];
rzz(0) q[5],q[20];
u3(pi,3.1642121206956397,-3.1642121206956397) q[5];
rzz(0) q[5],q[19];
u3(pi,6.058875591713275,-6.058875591713275) q[5];
rzz(0) q[5],q[18];
u3(pi,2.670353755551324,-2.670353755551324) q[5];
rzz(0) q[5],q[17];
u3(pi,5.56501722656896,-5.56501722656896) q[5];
rzz(0) q[5],q[16];
u3(pi,2.1771237089377267,-2.1771237089377267) q[5];
rzz(0) q[5],q[15];
u3(pi,5.071787179955362,-5.071787179955362) q[5];
rzz(0) q[5],q[14];
u3(pi,1.6832653437934113,-1.6832653437934113) q[5];
u3(pi/2,2.734442245684556,-2.734442245684556) q[13];
rzz(0) q[5],q[13];
u3(pi/2,4.701707565362484,-4.701707565362484) q[5];
u3(pi,0.9129468251331939,-0.9129468251331939) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi,1.3929821826017144,-1.3929821826017144) q[11];
rzz(pi/128) q[5],q[11];
u3(pi,0.03581415625092364,-0.03581415625092364) q[10];
rzz(0.049100328574073315) q[5],q[10];
u3(pi,1.6952033958770523,-1.6952033958770523) q[9];
rzz(pi/32) q[5],q[9];
u3(pi,3.064937792842202,-3.064937792842202) q[8];
rzz(0.19640131429629326) q[5],q[8];
u3(pi,3.0278669995298424,-3.0278669995298424) q[7];
rzz(pi/8) q[5],q[7];
u3(pi/2,3.5738758027237485,-3.5738758027237485) q[6];
rzz(0.7853830994606742) q[5],q[6];
u3(pi/2,3*pi/4,-3*pi/4) q[4];
rzz(0) q[4],q[21];
u3(pi,5.924415426139632,-5.924415426139632) q[4];
rzz(0) q[4],q[20];
u3(pi,3.6367076557955444,-3.6367076557955444) q[4];
rzz(0) q[4],q[19];
u3(pi,1.348999885451457,-1.348999885451457) q[4];
rzz(0) q[4],q[18];
u3(pi,5.344477422286956,-5.344477422286956) q[4];
rzz(0) q[4],q[17];
u3(pi,3.0567696519428686,-3.0567696519428686) q[4];
rzz(0) q[4],q[16];
u3(pi,0.7690618815987813,-0.7690618815987813) q[4];
rzz(0) q[4],q[15];
u3(pi,4.76453941843428,-4.76453941843428) q[4];
rzz(0) q[4],q[14];
u3(pi,2.476831648090193,-2.476831648090193) q[4];
rzz(0) q[4],q[13];
u3(pi,0.18912387774610553,-0.18912387774610553) q[4];
u3(pi/2,3.1535307056734343,-3.1535307056734343) q[12];
rzz(0) q[4],q[12];
u3(pi/2,3.75797313222411,-3.75797313222411) q[4];
u3(pi,5.296096895421673,-5.296096895421673) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi,4.031291693086422,-4.031291693086422) q[10];
rzz(pi/128) q[4],q[10];
u3(pi,2.626371458401067,-2.626371458401067) q[9];
rzz(0.049100328574073315) q[4],q[9];
u3(pi,5.61716766461855,-5.61716766461855) q[8];
rzz(pi/32) q[4],q[8];
u3(pi,3.5028758087526195,-3.5028758087526195) q[7];
rzz(0.19640131429629326) q[4],q[7];
u3(pi,1.0141061085787852,-1.0141061085787852) q[6];
rzz(pi/8) q[4],q[6];
u3(pi/2,4.701707565362484,-4.701707565362484) q[5];
rzz(0.7853830994606742) q[4],q[5];
u3(pi/2,3*pi/4,-3*pi/4) q[3];
rzz(0) q[3],q[21];
u3(pi,5.924415426139632,-5.924415426139632) q[3];
rzz(0) q[3],q[20];
u3(pi,3.6367076557955444,-3.6367076557955444) q[3];
rzz(0) q[3],q[19];
u3(pi,1.348999885451457,-1.348999885451457) q[3];
rzz(0) q[3],q[18];
u3(pi,5.344477422286956,-5.344477422286956) q[3];
rzz(0) q[3],q[17];
u3(pi,3.0567696519428686,-3.0567696519428686) q[3];
rzz(0) q[3],q[16];
u3(pi,0.7690618815987813,-0.7690618815987813) q[3];
rzz(0) q[3],q[15];
u3(pi,4.76453941843428,-4.76453941843428) q[3];
rzz(0) q[3],q[14];
u3(pi,2.476831648090193,-2.476831648090193) q[3];
rzz(0) q[3],q[13];
u3(pi,0.18912387774610553,-0.18912387774610553) q[3];
rzz(0) q[3],q[12];
u3(pi,4.184601414581604,-4.184601414581604) q[3];
u3(pi/2,2.581132524189374,-2.581132524189374) q[11];
rzz(0) q[3],q[11];
u3(pi/2,1.4696370433493051,-1.4696370433493051) q[3];
u3(pi,2.3021590965506005,-2.3021590965506005) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi,3.1522740686119985,-3.1522740686119985) q[9];
rzz(pi/128) q[3],q[9];
u3(pi,1.7360441003737197,-1.7360441003737197) q[8];
rzz(0.049100328574073315) q[3],q[8];
u3(pi,4.24429167499981,-4.24429167499981) q[7];
rzz(pi/32) q[3],q[7];
u3(pi,5.164778322501619,-5.164778322501619) q[6];
rzz(0.19640131429629326) q[3],q[6];
u3(pi,3.311238656883642,-3.311238656883642) q[5];
rzz(pi/8) q[3],q[5];
u3(pi/2,0.6163804786343174,-0.6163804786343174) q[4];
rzz(0.7853830994606742) q[3],q[4];
u3(pi/2,pi,-pi) q[2];
rzz(0) q[2],q[21];
u3(pi,5.74345968929286,-5.74345968929286) q[2];
rzz(0) q[2],q[20];
u3(pi,1.5230441184603318,-1.5230441184603318) q[2];
rzz(0) q[2],q[19];
u3(pi,3.5858138548073897,-3.5858138548073897) q[2];
rzz(0) q[2],q[18];
u3(pi,5.648583591154448,-5.648583591154448) q[2];
rzz(0) q[2],q[17];
u3(pi,1.4275397017912022,-1.4275397017912022) q[2];
rzz(0) q[2],q[16];
u3(pi,3.4903094381382602,-3.4903094381382602) q[2];
rzz(0) q[2],q[15];
u3(pi,5.553079174485318,-5.553079174485318) q[2];
rzz(0) q[2],q[14];
u3(pi,1.3326636036527904,-1.3326636036527904) q[2];
rzz(0) q[2],q[13];
u3(pi,3.3948050214691303,-3.3948050214691303) q[2];
rzz(0) q[2],q[12];
u3(pi,5.457574757816189,-5.457574757816189) q[2];
rzz(0) q[2],q[11];
u3(pi,1.2371591869836605,-1.2371591869836605) q[2];
u3(pi/2,0.14576989912656638,-0.14576989912656638) q[10];
rzz(0) q[2],q[10];
u3(pi/2,3.8390262226867273,-3.8390262226867273) q[2];
rzz(0.012275082143518329) q[2],q[9];
rzz(pi/128) q[2],q[8];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/32) q[2],q[6];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/8) q[2],q[4];
u3(pi/2,1.4696370433493051,-1.4696370433493051) q[3];
rzz(0.7853830994606742) q[2],q[3];
u3(pi/2,0,0) q[1];
rzz(0) q[1],q[21];
u3(pi,3.568220935947287,-3.568220935947287) q[1];
rzz(0) q[1],q[20];
u3(pi,1.2805131656031998,-1.2805131656031998) q[1];
rzz(0) q[1],q[19];
u3(pi,5.2759907024386985,-5.2759907024386985) q[1];
rzz(0) q[1],q[18];
u3(pi,2.988282932094611,-2.988282932094611) q[1];
rzz(0) q[1],q[17];
u3(pi,0.7005751617505239,-0.7005751617505239) q[1];
rzz(0) q[1],q[16];
u3(pi,4.696052698586023,-4.696052698586023) q[1];
rzz(0) q[1],q[15];
u3(pi,2.4083449282419354,-2.4083449282419354) q[1];
rzz(0) q[1],q[14];
u3(pi,0.12063715789784804,-0.12063715789784804) q[1];
rzz(0) q[1],q[13];
u3(pi,4.116114694733347,-4.116114694733347) q[1];
rzz(0) q[1],q[12];
u3(pi,1.8284069243892596,-1.8284069243892596) q[1];
rzz(0) q[1],q[11];
u3(pi,5.8238844612247584,-5.8238844612247584) q[1];
rzz(0) q[1],q[10];
u3(pi,3.536176690880671,-3.536176690880671) q[1];
u3(pi/2,1.790079494015464,-1.790079494015464) q[9];
rzz(0) q[1],q[9];
u3(pi/2,0.8212123196483719,-0.8212123196483719) q[1];
u3(pi,0.8519999276535519,-0.8519999276535519) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi,3.0027342583011243,-3.0027342583011243) q[7];
rzz(pi/128) q[1],q[7];
u3(pi,4.105433279711142,-4.105433279711142) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi,1.3056459068319182,-1.3056459068319182) q[5];
rzz(pi/32) q[1],q[5];
u3(pi,5.437468564833214,-5.437468564833214) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi,pi/20,-pi/20) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,3.8390262226867273,-3.8390262226867273) q[2];
rzz(0.7853830994606742) q[1],q[2];
rzz(pi/2) q[0],q[21];
rzz(-pi/2) q[0],q[20];
rzz(pi/2) q[0],q[19];
rzz(pi/2) q[0],q[18];
rzz(pi/2) q[0],q[17];
rzz(-pi/2) q[0],q[16];
rzz(pi/2) q[0],q[15];
rzz(pi/2) q[0],q[14];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[11];
rzz(pi/2) q[0],q[10];
rzz(pi/2) q[0],q[9];
u3(pi/2,0.09110618695410401,-0.09110618695410401) q[8];
rzz(pi/2) q[0],q[8];
u3(pi/2,6.044424265506762,-6.044424265506762) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,2.6188316360324517,-2.6188316360324517) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,3.830858081787394,-3.830858081787394) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,3.9753713438525247,-3.9753713438525247) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,5.127079210658542,-5.127079210658542) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,0.6974335690969341,-0.6974335690969341) q[2];
rzz(pi/2) q[0],q[2];
rzz(pi/2) q[0],q[1];
u3(pi,pi,-pi) q[0];
u3(pi/2,0.8212123196483719,-0.8212123196483719) q[1];
rzz(pi/2) q[0],q[1];
u3(pi/2,2.2682298958918308,-2.2682298958918308) q[2];
u3(pi/2,4.624424386084176,-4.624424386084176) q[2];
rzz(-pi/2) q[0],q[2];
u3(pi/2,3.5562828838636453,-3.5562828838636453) q[3];
u3(pi/2,0.02199114857512855,-0.02199114857512855) q[3];
rzz(pi/2) q[0],q[3];
u3(pi/2,5.546167670647421,-5.546167670647421) q[4];
u3(pi/2,2.2085396354736244,-2.2085396354736244) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,5.40165440858229,-5.40165440858229) q[5];
u3(pi/2,2.162044064200496,-2.162044064200496) q[5];
rzz(pi/2) q[0],q[5];
u3(pi/2,4.189627962827347,-4.189627962827347) q[6];
u3(pi/2,0.9990264638415542,-0.9990264638415542) q[6];
rzz(pi/2) q[0],q[6];
u3(pi/2,4.473627938711865,-4.473627938711865) q[7];
u3(pi/2,1.307530862424072,-1.307530862424072) q[7];
rzz(pi/2) q[0],q[7];
u3(pi/2,1.6619025137490007,-1.6619025137490007) q[8];
u3(pi/2,4.791557115255152,-4.791557115255152) q[8];
rzz(pi/2) q[0],q[8];
rzz(pi/2) q[0],q[9];
rzz(pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[11];
rzz(pi/2) q[0],q[12];
rzz(pi/2) q[0],q[13];
rzz(pi/2) q[0],q[14];
rzz(-pi/2) q[0],q[15];
rzz(pi/2) q[0],q[16];
rzz(pi/2) q[0],q[17];
rzz(pi/2) q[0],q[18];
rzz(pi/2) q[0],q[19];
rzz(-pi/2) q[0],q[20];
rzz(pi/2) q[0],q[21];
u3(pi/2,5.5336013000330615,-5.5336013000330615) q[1];
u3(pi/2,6.195220712879072,-6.195220712879072) q[2];
rzz(0.7853830994606742) q[1],q[2];
u3(pi/2,4.734380128959818,-4.734380128959818) q[3];
rzz(pi/8) q[1],q[3];
u3(pi/2,0.6377433086787281,-0.6377433086787281) q[4];
rzz(0.19640131429629326) q[1],q[4];
u3(pi/2,0.5912477374055991,-0.5912477374055991) q[5];
rzz(pi/32) q[1],q[5];
u3(pi/2,5.711415444226244,-5.711415444226244) q[6];
rzz(0.049100328574073315) q[1],q[6];
u3(pi/2,6.019919842808761,-6.019919842808761) q[7];
rzz(pi/128) q[1],q[7];
u3(pi/2,3.2207607884602556,-3.2207607884602556) q[8];
rzz(0.012275082143518329) q[1],q[8];
u3(pi/2,5.5336013000330615,-5.5336013000330615) q[1];
rzz(0) q[1],q[9];
u3(pi,1.85291134708726,-1.85291134708726) q[1];
rzz(0) q[1],q[10];
u3(pi,3.915681083434318,-3.915681083434318) q[1];
rzz(0) q[1],q[11];
u3(pi,5.9778225012506585,-5.9778225012506585) q[1];
rzz(0) q[1],q[12];
u3(pi,1.7574069304181303,-1.7574069304181303) q[1];
rzz(0) q[1],q[13];
u3(pi,3.8201766667651884,-3.8201766667651884) q[1];
rzz(0) q[1],q[14];
u3(pi,5.882318084581529,-5.882318084581529) q[1];
rzz(0) q[1],q[15];
u3(pi,1.6619025137490007,-1.6619025137490007) q[1];
rzz(0) q[1],q[16];
u3(pi,3.7246722500960585,-3.7246722500960585) q[1];
rzz(0) q[1],q[17];
u3(pi,5.787441986443117,-5.787441986443117) q[1];
rzz(0) q[1],q[18];
u3(pi,1.5663980970798708,-1.5663980970798708) q[1];
rzz(0) q[1],q[19];
u3(pi,3.629167833426929,-3.629167833426929) q[1];
rzz(0) q[1],q[20];
u3(pi,5.691937569773987,-5.691937569773987) q[1];
rzz(0) q[1],q[21];
u3(pi/2,2.2682298958918308,-2.2682298958918308) q[2];
rzz(0.7853830994606742) q[2],q[3];
rzz(pi/8) q[2],q[4];
rzz(0.19640131429629326) q[2],q[5];
rzz(pi/32) q[2],q[6];
rzz(0.049100328574073315) q[2],q[7];
rzz(pi/128) q[2],q[8];
u3(pi/2,1.7775131234011048,-1.7775131234011048) q[9];
rzz(0.012275082143518329) q[2],q[9];
u3(pi/2,0.6974335690969341,-0.6974335690969341) q[2];
rzz(0) q[2],q[10];
u3(pi,3.2999289233307185,-3.2999289233307185) q[2];
rzz(0) q[2],q[11];
u3(pi,5.362070341147059,-5.362070341147059) q[2];
rzz(0) q[2],q[12];
u3(pi,1.1416547703145308,-1.1416547703145308) q[2];
rzz(0) q[2],q[13];
u3(pi,3.204424506661589,-3.204424506661589) q[2];
rzz(0) q[2],q[14];
u3(pi,5.267194243008648,-5.267194243008648) q[2];
rzz(0) q[2],q[15];
u3(pi,1.046150353645401,-1.046150353645401) q[2];
rzz(0) q[2],q[16];
u3(pi,3.108920089992459,-3.108920089992459) q[2];
rzz(0) q[2],q[17];
u3(pi,5.1716898263395175,-5.1716898263395175) q[2];
rzz(0) q[2],q[18];
u3(pi,0.9512742555069894,-0.9512742555069894) q[2];
rzz(0) q[2],q[19];
u3(pi,3.0134156733233297,-3.0134156733233297) q[2];
rzz(0) q[2],q[20];
u3(pi,5.076185409670387,-5.076185409670387) q[2];
rzz(0) q[2],q[21];
u3(pi/2,1.200088393671301,-1.200088393671301) q[3];
rzz(0.7853830994606742) q[3],q[4];
rzz(pi/8) q[3],q[5];
rzz(0.19640131429629326) q[3],q[6];
rzz(pi/32) q[3],q[7];
rzz(0.049100328574073315) q[3],q[8];
rzz(pi/128) q[3],q[9];
u3(pi/2,0.136973439696515,-0.136973439696515) q[10];
rzz(0.012275082143518329) q[3],q[10];
u3(pi/2,2.7708847204661975,-2.7708847204661975) q[3];
rzz(0) q[3],q[11];
u3(pi,0.05592034923389832,-0.05592034923389832) q[3];
rzz(0) q[3],q[12];
u3(pi,4.051397886069397,-4.051397886069397) q[3];
rzz(0) q[3],q[13];
u3(pi,1.76369011572531,-1.76369011572531) q[3];
rzz(0) q[3],q[14];
u3(pi,5.759167652560809,-5.759167652560809) q[3];
rzz(0) q[3],q[15];
u3(pi,3.4714598822167213,-3.4714598822167213) q[3];
rzz(0) q[3],q[16];
u3(pi,1.1837521118726342,-1.1837521118726342) q[3];
rzz(0) q[3],q[17];
u3(pi,5.179229648708133,-5.179229648708133) q[3];
rzz(0) q[3],q[18];
u3(pi,2.8915218783640455,-2.8915218783640455) q[3];
rzz(0) q[3],q[19];
u3(pi,0.6038141080199583,-0.6038141080199583) q[3];
rzz(0) q[3],q[20];
u3(pi,4.599291644855457,-4.599291644855457) q[3];
rzz(0) q[3],q[21];
u3(pi/2,0.44107960856400696,-0.44107960856400696) q[4];
rzz(0.7853830994606742) q[4],q[5];
rzz(pi/8) q[4],q[6];
rzz(0.19640131429629326) q[4],q[7];
rzz(pi/32) q[4],q[8];
rzz(0.049100328574073315) q[4],q[9];
rzz(pi/128) q[4],q[10];
u3(pi/2,5.705760577449782,-5.705760577449782) q[11];
rzz(0.012275082143518329) q[4],q[11];
u3(pi/2,5.153468588948697,-5.153468588948697) q[4];
rzz(0) q[4],q[12];
u3(pi,2.4391325362471155,-2.4391325362471155) q[4];
rzz(0) q[4],q[13];
u3(pi,0.15079644737231007,-0.15079644737231007) q[4];
rzz(0) q[4],q[14];
u3(pi,4.146273984207809,-4.146273984207809) q[4];
rzz(0) q[4],q[15];
u3(pi,1.8585662138637216,-1.8585662138637216) q[4];
rzz(0) q[4],q[16];
u3(pi,5.85404375069922,-5.85404375069922) q[4];
rzz(0) q[4],q[17];
u3(pi,3.566335980355133,-3.566335980355133) q[4];
rzz(0) q[4],q[18];
u3(pi,1.2786282100110458,-1.2786282100110458) q[4];
rzz(0) q[4],q[19];
u3(pi,5.274105746846545,-5.274105746846545) q[4];
rzz(0) q[4],q[20];
u3(pi,2.9863979765024573,-2.9863979765024573) q[4];
rzz(0) q[4],q[21];
u3(pi/2,3.6348227002033906,-3.6348227002033906) q[5];
rzz(0.7853830994606742) q[5],q[6];
rzz(pi/8) q[5],q[7];
rzz(0.19640131429629326) q[5],q[8];
rzz(pi/32) q[5],q[9];
rzz(0.049100328574073315) q[5],q[10];
rzz(pi/128) q[5],q[11];
u3(pi/2,3.1371944238747673,-3.1371944238747673) q[12];
rzz(0.012275082143518329) q[5],q[12];
u3(pi/2,5.205619026998288,-5.205619026998288) q[5];
rzz(0) q[5],q[13];
u3(pi,2.490654655765988,-2.490654655765988) q[5];
rzz(0) q[5],q[14];
u3(pi,0.20294688542190065,-0.20294688542190065) q[5];
rzz(0) q[5],q[15];
u3(pi,4.1984244222574,-4.1984244222574) q[5];
rzz(0) q[5],q[16];
u3(pi,1.910716651913312,-1.910716651913312) q[5];
rzz(0) q[5],q[17];
u3(pi,5.9061941887488105,-5.9061941887488105) q[5];
rzz(0) q[5],q[18];
u3(pi,3.6184864184047236,-3.6184864184047236) q[5];
rzz(0) q[5],q[19];
u3(pi,1.3307786480606363,-1.3307786480606363) q[5];
rzz(0) q[5],q[20];
u3(pi,5.326256184896136,-5.326256184896136) q[5];
rzz(0) q[5],q[21];
u3(pi/2,5.662406598830243,-5.662406598830243) q[6];
rzz(0.7853830994606742) q[6],q[7];
rzz(pi/8) q[6],q[8];
rzz(0.19640131429629326) q[6],q[9];
rzz(pi/32) q[6],q[10];
rzz(0.049100328574073315) q[6],q[11];
rzz(pi/128) q[6],q[12];
u3(pi/2,5.857813661883529,-5.857813661883529) q[13];
rzz(0.012275082143518329) q[6],q[13];
u3(pi/2,4.091610272035346,-4.091610272035346) q[6];
rzz(0) q[6],q[14];
u3(pi,0.8268671864248335,-0.8268671864248335) q[6];
rzz(0) q[6],q[15];
u3(pi,3.7215306574424694,-3.7215306574424694) q[6];
rzz(0) q[6],q[16];
u3(pi,0.33300882128051806,-0.33300882128051806) q[6];
rzz(0) q[6],q[17];
u3(pi,3.2276722922981538,-3.2276722922981538) q[6];
rzz(0) q[6],q[18];
u3(pi,6.122964081846507,-6.122964081846507) q[6];
rzz(0) q[6],q[19];
u3(pi,2.734442245684556,-2.734442245684556) q[6];
rzz(0) q[6],q[20];
u3(pi,5.629105716702192,-5.629105716702192) q[6];
rzz(0) q[6],q[21];
u3(pi/2,2.853822766520968,-2.853822766520968) q[7];
rzz(0.7853830994606742) q[7],q[8];
rzz(pi/8) q[7],q[9];
rzz(0.19640131429629326) q[7],q[10];
rzz(pi/32) q[7],q[11];
rzz(0.049100328574073315) q[7],q[12];
rzz(pi/128) q[7],q[13];
u3(pi/2,5.683141110343936,-5.683141110343936) q[14];
rzz(0.012275082143518329) q[7],q[14];
u3(pi/2,4.424619093315865,-4.424619093315865) q[7];
rzz(0) q[7],q[15];
u3(pi,1.7102830406142833,-1.7102830406142833) q[7];
rzz(0) q[7],q[16];
u3(pi,5.705132258919065,-5.705132258919065) q[7];
rzz(0) q[7],q[17];
u3(pi,3.4174244885749774,-3.4174244885749774) q[7];
rzz(0) q[7],q[18];
u3(pi,1.1297167182308896,-1.1297167182308896) q[7];
rzz(0) q[7],q[19];
u3(pi,5.125194255066388,-5.125194255066388) q[7];
rzz(0) q[7],q[20];
u3(pi,2.837486484722301,-2.837486484722301) q[7];
rzz(0) q[7],q[21];
u3(pi/2,0.06660176425610362,-0.06660176425610362) q[8];
rzz(0.7853830994606742) q[8],q[9];
rzz(pi/8) q[8],q[10];
rzz(0.19640131429629326) q[8],q[11];
rzz(pi/32) q[8],q[12];
rzz(0.049100328574073315) q[8],q[13];
rzz(pi/128) q[8],q[14];
u3(pi/2,2.293362637120549,-2.293362637120549) q[15];
rzz(0.012275082143518329) q[8],q[15];
u3(pi/2,4.778990744640794,-4.778990744640794) q[8];
rzz(0) q[8],q[16];
u3(pi,2.623229865747477,-2.623229865747477) q[8];
rzz(0) q[8],q[17];
u3(pi,1.4520441244892024,-1.4520441244892024) q[8];
rzz(0) q[8],q[18];
u3(pi,0.2814867017616455,-0.2814867017616455) q[8];
rzz(0) q[8],q[19];
u3(pi,5.394114586213675,-5.394114586213675) q[8];
rzz(0) q[8],q[20];
u3(pi,4.2229288449554,-4.2229288449554) q[8];
rzz(0) q[8],q[21];
u3(pi/2,4.907167724907257,-4.907167724907257) q[9];
rzz(0.7853830994606742) q[9],q[10];
rzz(pi/8) q[9],q[11];
rzz(0.19640131429629326) q[9],q[12];
rzz(pi/32) q[9],q[13];
rzz(0.049100328574073315) q[9],q[14];
rzz(pi/128) q[9],q[15];
u3(pi/2,4.059566026968731,-4.059566026968731) q[16];
rzz(0.012275082143518329) q[9],q[16];
u3(pi/2,0.19477874452256716,-0.19477874452256716) q[9];
rzz(0) q[9],q[17];
u3(pi,3.934530639355857,-3.934530639355857) q[9];
rzz(0) q[9],q[18];
u3(pi,1.988628149722339,-1.988628149722339) q[9];
rzz(0) q[9],q[19];
u3(pi,0.042725660088821185,-0.042725660088821185) q[9];
rzz(0) q[9],q[20];
u3(pi,4.380636796165608,-4.380636796165608) q[9];
rzz(0) q[9],q[21];
u3(pi/2,0.12754866173574558,-0.12754866173574558) q[10];
rzz(0.7853830994606742) q[10],q[11];
rzz(pi/8) q[10],q[12];
rzz(0.19640131429629326) q[10],q[13];
rzz(pi/32) q[10],q[14];
rzz(0.049100328574073315) q[10],q[15];
rzz(pi/128) q[10],q[16];
u3(pi/2,3.4211943997592846,-3.4211943997592846) q[17];
rzz(0.012275082143518329) q[10],q[17];
u3(pi/2,4.839937642120435,-4.839937642120435) q[10];
rzz(0) q[10],q[18];
u3(pi,2.124973270888136,-2.124973270888136) q[10];
rzz(0) q[10],q[19];
u3(pi,6.120450807723635,-6.120450807723635) q[10];
rzz(0) q[10],q[20];
u3(pi,3.8327430373795477,-3.8327430373795477) q[10];
rzz(0) q[10],q[21];
u3(pi/2,5.698220755081167,-5.698220755081167) q[11];
rzz(0.7853830994606742) q[11],q[12];
rzz(pi/8) q[11],q[13];
rzz(0.19640131429629326) q[11],q[14];
rzz(pi/32) q[11],q[15];
rzz(0.049100328574073315) q[11],q[16];
rzz(pi/128) q[11],q[17];
u3(pi/2,2.339858208393678,-2.339858208393678) q[18];
rzz(0.012275082143518329) q[11],q[18];
u3(pi/2,0.987716730288631,-0.987716730288631) q[11];
rzz(0) q[11],q[19];
u3(pi,5.114512840044183,-5.114512840044183) q[11];
rzz(0) q[11],q[20];
u3(pi,3.9439554173166265,-3.9439554173166265) q[11];
rzz(0) q[11],q[21];
u3(pi/2,6.271875573626663,-6.271875573626663) q[12];
rzz(0.7853830994606742) q[12],q[13];
rzz(pi/8) q[12],q[14];
rzz(0.19640131429629326) q[12],q[15];
rzz(pi/32) q[12],q[16];
rzz(0.049100328574073315) q[12],q[17];
rzz(pi/128) q[12],q[18];
u3(pi/2,3.6178580998740055,-3.6178580998740055) q[19];
rzz(0.012275082143518329) q[12],q[19];
u3(pi/2,4.701707565362484,-4.701707565362484) q[12];
rzz(0) q[12],q[20];
u3(pi,2.158274153016188,-2.158274153016188) q[12];
rzz(0) q[12],q[21];
u3(pi/2,2.709937822986556,-2.709937822986556) q[13];
rzz(0.7853830994606742) q[13],q[14];
rzz(pi/8) q[13],q[15];
rzz(0.19640131429629326) q[13],q[16];
rzz(pi/32) q[13],q[17];
rzz(0.049100328574073315) q[13],q[18];
rzz(pi/128) q[13],q[19];
u3(pi/2,5.912477374055991,-5.912477374055991) q[20];
rzz(0.012275082143518329) q[13],q[20];
u3(pi/2,4.280734149781452,-4.280734149781452) q[13];
rzz(0) q[13],q[21];
u3(pi/2,2.535265271446963,-2.535265271446963) q[14];
rzz(0.7853830994606742) q[14],q[15];
rzz(pi/8) q[14],q[16];
rzz(0.19640131429629326) q[14],q[17];
rzz(pi/32) q[14],q[18];
rzz(0.049100328574073315) q[14],q[19];
rzz(pi/128) q[14],q[20];
u3(pi/2,5.762937563745116,-5.762937563745116) q[21];
rzz(0.012275082143518329) q[14],q[21];
u3(pi/2,5.428672105403162,-5.428672105403162) q[15];
rzz(0.7853830994606742) q[15],q[16];
rzz(pi/8) q[15],q[17];
rzz(0.19640131429629326) q[15],q[18];
rzz(pi/32) q[15],q[19];
rzz(0.049100328574073315) q[15],q[20];
rzz(pi/128) q[15],q[21];
u3(pi/2,4.053282841661551,-4.053282841661551) q[16];
rzz(0.7853830994606742) q[16],q[17];
rzz(pi/8) q[16],q[18];
rzz(0.19640131429629326) q[16],q[19];
rzz(pi/32) q[16],q[20];
rzz(0.049100328574073315) q[16],q[21];
u3(pi/2,0.273318560862312,-0.273318560862312) q[17];
rzz(0.7853830994606742) q[17],q[18];
rzz(pi/8) q[17],q[19];
rzz(0.19640131429629326) q[17],q[20];
rzz(pi/32) q[17],q[21];
u3(pi/2,2.3335750230864982,-2.3335750230864982) q[18];
rzz(0.7853830994606742) q[18],q[19];
rzz(pi/8) q[18],q[20];
rzz(0.19640131429629326) q[18],q[21];
u3(pi/2,3.611574914566826,-3.611574914566826) q[19];
rzz(0.7853830994606742) q[19],q[20];
rzz(pi/8) q[19],q[21];
u3(pi/2,2.764601535159018,-2.764601535159018) q[20];
rzz(0.7853830994606742) q[20],q[21];
u3(pi/2,2.0106192982974678,-2.0106192982974678) q[1];
u3(pi/2,1.3948671381938682,-1.3948671381938682) q[2];
u3(pi/2,1.884327273623158,-1.884327273623158) q[3];
u3(pi/2,0.27206192380087607,-0.27206192380087607) q[4];
u3(pi/2,2.6112918136638363,-2.6112918136638363) q[5];
u3(pi/2,2.3643626310916783,-2.3643626310916783) q[6];
u3(pi/2,0.12315043202071989,-0.12315043202071989) q[7];
u3(pi/2,2.067167966062084,-2.067167966062084) q[8];
u3(pi/2,1.836575065288593,-1.836575065288593) q[9];
u3(pi/2,1.1184069846779663,-1.1184069846779663) q[10];
u3(pi/2,1.7875662198925921,-1.7875662198925921) q[11];
u3(pi/2,5.898026047849477,-5.898026047849477) q[12];
u3(pi/2,1.1391414961916588,-1.1391414961916588) q[13];
u3(pi/2,5.756654378437937,-5.756654378437937) q[21];
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
