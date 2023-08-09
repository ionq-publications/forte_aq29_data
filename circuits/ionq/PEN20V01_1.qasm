OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
creg c[20];
u3(pi/2,3.666238626739289,-3.666238626739289) q[18];
u3(pi/2,3.9137961278421645,-3.9137961278421645) q[19];
rzz(-pi/2) q[18],q[19];
u3(pi/2,2.342999801047268,-2.342999801047268) q[19];
u3(pi/2,2.6427077401997336,-2.6427077401997336) q[19];
rzz(pi/2) q[18],q[19];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[17];
rzz(pi/2) q[17],q[19];
u3(pi/2,5.784300393789527,-5.784300393789527) q[19];
u3(pi/2,2.0439201804255194,-2.0439201804255194) q[19];
rzz(pi/2) q[17],q[19];
u3(pi/2,2.9197962122463537,-2.9197962122463537) q[16];
rzz(pi/2) q[16],q[19];
u3(pi/2,2.0439201804255194,-2.0439201804255194) q[19];
u3(pi/2,3.9879377144668835,-3.9879377144668835) q[19];
rzz(-pi/2) q[16],q[19];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[15];
rzz(pi/2) q[15],q[19];
u3(pi/2,0.8463450608770902,-0.8463450608770902) q[19];
u3(pi/2,1.5921591568393072,-1.5921591568393072) q[19];
rzz(pi/2) q[15],q[19];
u3(pi/2,3.822689940888061,-3.822689940888061) q[14];
rzz(pi/2) q[14],q[19];
u3(pi/2,4.7337518104291,-4.7337518104291) q[19];
u3(pi/2,0.09990264638415543,-0.09990264638415543) q[19];
rzz(-pi/2) q[14],q[19];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[13];
rzz(pi/2) q[13],q[19];
u3(pi/2,3.241495299973949,-3.241495299973949) q[19];
u3(pi/2,3.39794661412272,-3.39794661412272) q[19];
rzz(pi/2) q[13],q[19];
u3(pi/2,1.1498229112138643,-1.1498229112138643) q[12];
rzz(pi/2) q[12],q[19];
u3(pi/2,0.25635396053292714,-0.25635396053292714) q[19];
u3(pi/2,3.085043985825177,-3.085043985825177) q[19];
rzz(pi/2) q[12],q[19];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[11];
rzz(-pi/2) q[11],q[19];
u3(pi/2,6.22663663941497,-6.22663663941497) q[19];
u3(pi/2,2.45923872923009,-2.45923872923009) q[19];
rzz(pi/2) q[11],q[19];
u3(pi/2,6.1694596531196355,-6.1694596531196355) q[10];
rzz(pi/2) q[10],q[19];
u3(pi/2,2.45923872923009,-2.45923872923009) q[19];
u3(pi/2,4.349220869629709,-4.349220869629709) q[19];
rzz(pi/2) q[10],q[19];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
rzz(-pi/2) q[9],q[19];
u3(pi/2,1.2076282160399165,-1.2076282160399165) q[19];
u3(pi/2,1.8453715247186446,-1.8453715247186446) q[19];
rzz(pi/2) q[9],q[19];
u3(pi/2,1.1165220290858124,-1.1165220290858124) q[8];
rzz(pi/2) q[8],q[19];
u3(pi/2,4.986964178308438,-4.986964178308438) q[19];
u3(pi/2,0.5692565888304705,-0.5692565888304705) q[19];
rzz(pi/2) q[8],q[19];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
rzz(pi/2) q[7],q[19];
u3(pi/2,0.5692565888304705,-0.5692565888304705) q[19];
u3(pi/2,1.1579910521131977,-1.1579910521131977) q[19];
rzz(-pi/2) q[7],q[19];
u3(pi/2,6.037512761668864,-6.037512761668864) q[6];
rzz(pi/2) q[6],q[19];
u3(pi/2,1.1579910521131977,-1.1579910521131977) q[19];
u3(pi/2,3.1214864606068184,-3.1214864606068184) q[19];
rzz(pi/2) q[6],q[19];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(pi/2) q[5],q[19];
u3(pi/2,3.1214864606068184,-3.1214864606068184) q[19];
u3(pi/2,3.9068846240042667,-3.9068846240042667) q[19];
rzz(-pi/2) q[5],q[19];
u3(pi/2,0.5893627818134451,-0.5893627818134451) q[4];
u3(pi/2,5.477680950799163,-5.477680950799163) q[19];
rzz(pi/2) q[4],q[19];
u3(pi/2,pi,-pi) q[1];
u3(pi/2,5*pi/4,-5*pi/4) q[1];
rzz(-pi/2) q[0],q[1];
rzz(pi/2) q[0],q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
u3(pi/2,13*pi/8,-13*pi/8) q[2];
rzz(pi/2) q[0],q[2];
rzz(-pi/2) q[0],q[3];
u3(pi/2,7*pi/8,-7*pi/8) q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
rzz(-pi/2) q[0],q[3];
u3(pi/2,3.730955435403238,-3.730955435403238) q[4];
rzz(-pi/2) q[0],q[4];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[4];
u3(pi/2,5.203105752875415,-5.203105752875415) q[4];
rzz(pi/2) q[0],q[4];
u3(pi/2,0.6873804726054468,-0.6873804726054468) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,5.399769452990137,-5.399769452990137) q[5];
u3(pi/2,2.2091679540043425,-2.2091679540043425) q[5];
rzz(-pi/2) q[0],q[5];
u3(pi/2,2.8959201080790713,-2.8959201080790713) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,1.3251237812841747,-1.3251237812841747) q[6];
u3(pi/2,4.442212012175967,-4.442212012175967) q[6];
rzz(-pi/2) q[0],q[6];
u3(pi/2,6.037512761668864,-6.037512761668864) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,4.466716434873968,-4.466716434873968) q[7];
u3(pi/2,1.3131857292005336,-1.3131857292005336) q[7];
rzz(-pi/2) q[0],q[7];
u3(pi/2,4.258114682675606,-4.258114682675606) q[8];
rzz(-pi/2) q[0],q[8];
rzz(-pi/2) q[0],q[8];
u3(pi/2,0.5340707511102649,-0.5340707511102649) q[9];
rzz(pi/2) q[0],q[9];
rzz(-pi/2) q[0],q[9];
u3(pi/2,3.0278669995298424,-3.0278669995298424) q[10];
rzz(-pi/2) q[0],q[10];
rzz(-pi/2) q[0],q[10];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[11];
rzz(pi/2) q[0],q[11];
rzz(-pi/2) q[0],q[11];
u3(pi/2,4.291415564803658,-4.291415564803658) q[12];
rzz(-pi/2) q[0],q[12];
rzz(-pi/2) q[0],q[12];
u3(pi/2,3.669380219392878,-3.669380219392878) q[13];
rzz(-pi/2) q[0],q[13];
rzz(pi/2) q[0],q[13];
u3(pi/2,3.822689940888061,-3.822689940888061) q[14];
rzz(-pi/2) q[0],q[14];
rzz(-pi/2) q[0],q[14];
u3(pi/2,3.669380219392878,-3.669380219392878) q[15];
rzz(-pi/2) q[0],q[15];
rzz(pi/2) q[0],q[15];
u3(pi/2,2.9197962122463537,-2.9197962122463537) q[16];
rzz(-pi/2) q[0],q[16];
rzz(-pi/2) q[0],q[16];
u3(pi/2,3.669380219392878,-3.669380219392878) q[17];
rzz(-pi/2) q[0],q[17];
rzz(pi/2) q[0],q[17];
u3(pi/2,3.666238626739289,-3.666238626739289) q[18];
rzz(-pi/2) q[0],q[18];
rzz(-pi/2) q[0],q[18];
u3(pi/2,7*pi/4,-7*pi/4) q[1];
u3(pi/2,pi/2,-pi/2) q[1];
rzz(pi/2) q[1],q[2];
u3(pi/2,13*pi/8,-13*pi/8) q[2];
u3(pi/2,3*pi/8,-3*pi/8) q[2];
rzz(-pi/2) q[1],q[2];
rzz(pi/2) q[1],q[3];
u3(pi/2,5.693822525366141,-5.693822525366141) q[3];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[3];
rzz(pi/2) q[1],q[3];
rzz(pi/2) q[1],q[4];
u3(pi/2,2.061513099285622,-2.061513099285622) q[4];
u3(pi/2,5.007070371291412,-5.007070371291412) q[4];
rzz(pi/2) q[1],q[4];
rzz(-pi/2) q[1],q[5];
u3(pi/2,5.350760607594136,-5.350760607594136) q[5];
u3(pi/2,2.110521944681623,-2.110521944681623) q[5];
rzz(pi/2) q[1],q[5];
rzz(pi/2) q[1],q[6];
u3(pi/2,4.442212012175967,-4.442212012175967) q[6];
u3(pi/2,1.2516105131901736,-1.2516105131901736) q[6];
rzz(pi/2) q[1],q[6];
rzz(-pi/2) q[1],q[7];
u3(pi/2,4.454778382790327,-4.454778382790327) q[7];
u3(pi/2,1.288681306502533,-1.288681306502533) q[7];
rzz(pi/2) q[1],q[7];
rzz(pi/2) q[1],q[8];
u3(pi/2,5.82325614269404,-5.82325614269404) q[8];
u3(pi/2,2.6690971184898884,-2.6690971184898884) q[8];
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
rzz(pi/2) q[1],q[14];
rzz(pi/2) q[1],q[14];
rzz(pi/2) q[1],q[15];
rzz(-pi/2) q[1],q[15];
rzz(pi/2) q[1],q[16];
rzz(pi/2) q[1],q[16];
rzz(pi/2) q[1],q[17];
rzz(pi/2) q[1],q[17];
rzz(pi/2) q[1],q[18];
rzz(pi/2) q[1],q[18];
u3(pi/2,7*pi/8,-7*pi/8) q[2];
u3(pi/2,3*pi/2,-3*pi/2) q[2];
rzz(pi/2) q[2],q[3];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[3];
u3(pi/2,4.516353598800687,-4.516353598800687) q[3];
rzz(pi/2) q[2],q[3];
rzz(-pi/2) q[2],q[4];
u3(pi/2,1.865477717701619,-1.865477717701619) q[4];
u3(pi/2,4.614371289592689,-4.614371289592689) q[4];
rzz(pi/2) q[2],q[4];
rzz(pi/2) q[2],q[5];
u3(pi/2,2.110521944681623,-2.110521944681623) q[5];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[5];
rzz(pi/2) q[2],q[5];
rzz(-pi/2) q[2],q[6];
u3(pi/2,4.393203166779967,-4.393203166779967) q[6];
u3(pi/2,1.1535928223981722,-1.1535928223981722) q[6];
rzz(pi/2) q[2],q[6];
rzz(pi/2) q[2],q[7];
u3(pi/2,1.288681306502533,-1.288681306502533) q[7];
u3(pi/2,4.381265114696325,-4.381265114696325) q[7];
rzz(pi/2) q[2],q[7];
rzz(pi/2) q[2],q[8];
u3(pi/2,2.6690971184898884,-2.6690971184898884) q[8];
u3(pi/2,5.7861853493816815,-5.7861853493816815) q[8];
rzz(-pi/2) q[2],q[8];
rzz(pi/2) q[2],q[9];
u3(pi/2,5.237034953534185,-5.237034953534185) q[9];
u3(pi/2,2.082875929330033,-2.082875929330033) q[9];
rzz(pi/2) q[2],q[9];
rzz(pi/2) q[2],q[10];
rzz(-pi/2) q[2],q[10];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[11];
rzz(pi/2) q[2],q[12];
rzz(pi/2) q[2],q[12];
rzz(-pi/2) q[2],q[13];
rzz(pi/2) q[2],q[13];
rzz(pi/2) q[2],q[14];
rzz(-pi/2) q[2],q[14];
rzz(pi/2) q[2],q[15];
rzz(pi/2) q[2],q[15];
rzz(pi/2) q[2],q[16];
rzz(pi/2) q[2],q[16];
rzz(-pi/2) q[2],q[17];
rzz(pi/2) q[2],q[17];
rzz(pi/2) q[2],q[18];
rzz(pi/2) q[2],q[18];
u3(pi/2,6.087149925595583,-6.087149925595583) q[3];
u3(pi/2,3*pi/8,-3*pi/8) q[3];
rzz(-pi/2) q[3],q[4];
u3(pi/2,1.472778636002895,-1.472778636002895) q[4];
u3(pi/2,3.82897312619524,-3.82897312619524) q[4];
rzz(pi/2) q[3],q[4];
rzz(pi/2) q[3],q[5];
u3(pi/2,5.0560792166874124,-5.0560792166874124) q[5];
u3(pi/2,1.5217874813988959,-1.5217874813988959) q[5];
rzz(pi/2) q[3],q[5];
rzz(pi/2) q[3],q[6];
u3(pi/2,1.1535928223981722,-1.1535928223981722) q[6];
u3(pi/2,4.098521775873244,-4.098521775873244) q[6];
rzz(-pi/2) q[3],q[6];
rzz(pi/2) q[3],q[7];
u3(pi/2,4.381265114696325,-4.381265114696325) q[7];
u3(pi/2,1.141026451783813,-1.141026451783813) q[7];
rzz(pi/2) q[3],q[7];
rzz(pi/2) q[3],q[8];
u3(pi/2,2.644592695791888,-2.644592695791888) q[8];
u3(pi/2,5.737176503985681,-5.737176503985681) q[8];
rzz(-pi/2) q[3],q[8];
rzz(pi/2) q[3],q[9];
u3(pi/2,2.082875929330033,-2.082875929330033) q[9];
u3(pi/2,5.199964160221826,-5.199964160221826) q[9];
rzz(pi/2) q[3],q[9];
rzz(pi/2) q[3],q[10];
u3(pi/2,1.4495308503663304,-1.4495308503663304) q[10];
u3(pi/2,4.579185451872482,-4.579185451872482) q[10];
rzz(pi/2) q[3],q[10];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[11];
rzz(pi/2) q[3],q[12];
rzz(pi/2) q[3],q[12];
rzz(-pi/2) q[3],q[13];
rzz(pi/2) q[3],q[13];
rzz(pi/2) q[3],q[14];
rzz(pi/2) q[3],q[14];
rzz(pi/2) q[3],q[15];
rzz(pi/2) q[3],q[15];
rzz(pi/2) q[3],q[16];
rzz(pi/2) q[3],q[16];
rzz(-pi/2) q[3],q[17];
rzz(pi/2) q[3],q[17];
rzz(pi/2) q[3],q[18];
rzz(pi/2) q[3],q[18];
u3(pi/2,5.399769452990137,-5.399769452990137) q[4];
u3(pi/2,2.1601591086083416,-2.1601591086083416) q[4];
rzz(-pi/2) q[4],q[5];
u3(pi/2,4.663380134988689,-4.663380134988689) q[5];
u3(pi/2,0.7363893180014475,-0.7363893180014475) q[5];
rzz(pi/2) q[4],q[5];
rzz(pi/2) q[4],q[6];
u3(pi/2,0.956929122283451,-0.956929122283451) q[6];
u3(pi/2,3.70582269417452,-3.70582269417452) q[6];
rzz(pi/2) q[4],q[6];
rzz(pi/2) q[4],q[7];
u3(pi/2,1.141026451783813,-1.141026451783813) q[7];
u3(pi/2,4.0865837237896026,-4.0865837237896026) q[7];
rzz(-pi/2) q[4],q[7];
rzz(pi/2) q[4],q[8];
u3(pi/2,2.595583850395887,-2.595583850395887) q[8];
u3(pi/2,5.639158813193679,-5.639158813193679) q[8];
rzz(pi/2) q[4],q[8];
rzz(pi/2) q[4],q[9];
u3(pi/2,5.199964160221826,-5.199964160221826) q[9];
u3(pi/2,2.0093626612360316,-2.0093626612360316) q[9];
rzz(-pi/2) q[4],q[9];
rzz(pi/2) q[4],q[10];
u3(pi/2,4.579185451872482,-4.579185451872482) q[10];
u3(pi/2,1.4130883755846888,-1.4130883755846888) q[10];
rzz(pi/2) q[4],q[10];
rzz(pi/2) q[4],q[11];
u3(pi/2,5.231380086757723,-5.231380086757723) q[11];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[11];
rzz(pi/2) q[4],q[11];
rzz(-pi/2) q[4],q[12];
rzz(pi/2) q[4],q[12];
rzz(pi/2) q[4],q[13];
rzz(pi/2) q[4],q[13];
rzz(-pi/2) q[4],q[14];
rzz(pi/2) q[4],q[14];
rzz(pi/2) q[4],q[15];
rzz(pi/2) q[4],q[15];
rzz(-pi/2) q[4],q[16];
rzz(pi/2) q[4],q[16];
rzz(pi/2) q[4],q[17];
rzz(pi/2) q[4],q[17];
rzz(-pi/2) q[4],q[18];
rzz(pi/2) q[4],q[18];
u3(pi/2,2.3071856447963444,-2.3071856447963444) q[5];
u3(pi/2,3.0435749627977917,-3.0435749627977917) q[5];
rzz(pi/2) q[5],q[6];
u3(pi/2,3.70582269417452,-3.70582269417452) q[6];
u3(pi/2,6.062017184366865,-6.062017184366865) q[6];
rzz(pi/2) q[5],q[6];
rzz(pi/2) q[5],q[7];
u3(pi/2,0.9449910701998098,-0.9449910701998098) q[7];
u3(pi/2,3.6938846420908784,-3.6938846420908784) q[7];
rzz(-pi/2) q[5],q[7];
rzz(pi/2) q[5],q[8];
u3(pi/2,5.639158813193679,-5.639158813193679) q[8];
u3(pi/2,2.3009024594891647,-2.3009024594891647) q[8];
rzz(pi/2) q[5],q[8];
rzz(pi/2) q[5],q[9];
u3(pi/2,5.150955314825825,-5.150955314825825) q[9];
u3(pi/2,1.9113449704440304,-1.9113449704440304) q[9];
rzz(-pi/2) q[5],q[9];
rzz(pi/2) q[5],q[10];
u3(pi/2,1.4130883755846888,-1.4130883755846888) q[10];
u3(pi/2,4.505043865247763,-4.505043865247763) q[10];
rzz(pi/2) q[5],q[10];
rzz(pi/2) q[5],q[11];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[11];
u3(pi/2,5.194309293445364,-5.194309293445364) q[11];
rzz(pi/2) q[5],q[11];
rzz(pi/2) q[5],q[12];
u3(pi/2,2.7124510971094273,-2.7124510971094273) q[12];
u3(pi/2,5.841477380084861,-5.841477380084861) q[12];
rzz(pi/2) q[5],q[12];
rzz(pi/2) q[5],q[13];
rzz(pi/2) q[5],q[13];
rzz(-pi/2) q[5],q[14];
rzz(pi/2) q[5],q[14];
rzz(pi/2) q[5],q[15];
rzz(pi/2) q[5],q[15];
rzz(-pi/2) q[5],q[16];
rzz(pi/2) q[5],q[16];
rzz(pi/2) q[5],q[17];
rzz(pi/2) q[5],q[17];
rzz(-pi/2) q[5],q[18];
rzz(pi/2) q[5],q[18];
u3(pi/2,4.491220857571968,-4.491220857571968) q[6];
u3(pi/2,4.9084243619686925,-4.9084243619686925) q[6];
rzz(pi/2) q[6],q[7];
u3(pi/2,0.5522919885010856,-0.5522919885010856) q[7];
u3(pi/2,2.90848647869343,-2.90848647869343) q[7];
rzz(pi/2) q[6],q[7];
rzz(-pi/2) q[6],q[8];
u3(pi/2,5.442495113078958,-5.442495113078958) q[8];
u3(pi/2,1.9082033777904406,-1.9082033777904406) q[8];
rzz(pi/2) q[6],q[8];
rzz(pi/2) q[6],q[9];
u3(pi/2,5.052937624033824,-5.052937624033824) q[9];
u3(pi/2,1.7146812703293088,-1.7146812703293088) q[9];
rzz(pi/2) q[6],q[9];
rzz(pi/2) q[6],q[10];
u3(pi/2,4.505043865247763,-4.505043865247763) q[10];
u3(pi/2,1.2654335208659686,-1.2654335208659686) q[10];
rzz(-pi/2) q[6],q[10];
rzz(pi/2) q[6],q[11];
u3(pi/2,5.194309293445364,-5.194309293445364) q[11];
u3(pi/2,2.00370779445957,-2.00370779445957) q[11];
rzz(pi/2) q[6],q[11];
rzz(pi/2) q[6],q[12];
u3(pi/2,5.841477380084861,-5.841477380084861) q[12];
u3(pi/2,2.675380303797068,-2.675380303797068) q[12];
rzz(-pi/2) q[6],q[12];
rzz(pi/2) q[6],q[13];
u3(pi/2,5.228238494104134,-5.228238494104134) q[13];
u3(pi/2,2.0740794698999814,-2.0740794698999814) q[13];
rzz(pi/2) q[6],q[13];
rzz(pi/2) q[6],q[14];
rzz(pi/2) q[6],q[14];
rzz(-pi/2) q[6],q[15];
rzz(pi/2) q[6],q[15];
rzz(pi/2) q[6],q[16];
rzz(pi/2) q[6],q[16];
rzz(-pi/2) q[6],q[17];
rzz(pi/2) q[6],q[17];
rzz(pi/2) q[6],q[18];
rzz(pi/2) q[6],q[18];
u3(pi/2,4.479282805488327,-4.479282805488327) q[7];
u3(pi/2,0.3436902363027234,-0.3436902363027234) q[7];
rzz(pi/2) q[7],q[8];
u3(pi/2,1.9082033777904406,-1.9082033777904406) q[8];
u3(pi/2,4.264397867982785,-4.264397867982785) q[8];
rzz(-pi/2) q[7],q[8];
rzz(pi/2) q[7],q[9];
u3(pi/2,1.7146812703293088,-1.7146812703293088) q[9];
u3(pi/2,4.463574842220378,-4.463574842220378) q[9];
rzz(pi/2) q[7],q[9];
rzz(-pi/2) q[7],q[10];
u3(pi/2,1.2654335208659686,-1.2654335208659686) q[10];
u3(pi/2,4.210990792871759,-4.210990792871759) q[10];
rzz(pi/2) q[7],q[10];
rzz(pi/2) q[7],q[11];
u3(pi/2,2.00370779445957,-2.00370779445957) q[11];
u3(pi/2,5.047282757257362,-5.047282757257362) q[11];
rzz(pi/2) q[7],q[11];
rzz(pi/2) q[7],q[12];
u3(pi/2,5.81697295738686,-5.81697295738686) q[12];
u3(pi/2,2.626371458401067,-2.626371458401067) q[12];
rzz(-pi/2) q[7],q[12];
rzz(pi/2) q[7],q[13];
u3(pi/2,2.0740794698999814,-2.0740794698999814) q[13];
u3(pi/2,5.191167700791774,-5.191167700791774) q[13];
rzz(pi/2) q[7],q[13];
rzz(pi/2) q[7],q[14];
u3(pi/2,2.2405838805402403,-2.2405838805402403) q[14];
u3(pi/2,5.3696101635156746,-5.3696101635156746) q[14];
rzz(-pi/2) q[7],q[14];
rzz(pi/2) q[7],q[15];
rzz(pi/2) q[7],q[15];
rzz(pi/2) q[7],q[16];
rzz(pi/2) q[7],q[16];
rzz(pi/2) q[7],q[17];
rzz(pi/2) q[7],q[17];
rzz(pi/2) q[7],q[18];
rzz(pi/2) q[7],q[18];
u3(pi/2,5.8351941947776815,-5.8351941947776815) q[8];
u3(pi/2,2.4052033355883453,-2.4052033355883453) q[8];
rzz(-pi/2) q[8],q[9];
u3(pi/2,1.321982188630585,-1.321982188630585) q[9];
u3(pi/2,3.67817667882293,-3.67817667882293) q[9];
rzz(pi/2) q[8],q[9];
rzz(pi/2) q[8],q[10];
u3(pi/2,4.210990792871759,-4.210990792871759) q[10];
u3(pi/2,0.6766990575832414,-0.6766990575832414) q[10];
rzz(pi/2) q[8],q[10];
rzz(-pi/2) q[8],q[11];
u3(pi/2,1.9056901036675686,-1.9056901036675686) q[11];
u3(pi/2,4.850619057142641,-4.850619057142641) q[11];
rzz(pi/2) q[8],q[11];
rzz(pi/2) q[8],q[12];
u3(pi/2,5.76796411199086,-5.76796411199086) q[12];
u3(pi/2,2.5283537676090653,-2.5283537676090653) q[12];
rzz(pi/2) q[8],q[12];
rzz(pi/2) q[8],q[13];
u3(pi/2,5.191167700791774,-5.191167700791774) q[13];
u3(pi/2,2.0005662018059804,-2.0005662018059804) q[13];
rzz(-pi/2) q[8],q[13];
rzz(pi/2) q[8],q[14];
u3(pi/2,2.2280175099258814,-2.2280175099258814) q[14];
u3(pi/2,5.345105740817674,-5.345105740817674) q[14];
rzz(pi/2) q[8],q[14];
rzz(-pi/2) q[8],q[15];
u3(pi/2,2.0891591146372126,-2.0891591146372126) q[15];
u3(pi/2,5.218813716143364,-5.218813716143364) q[15];
rzz(pi/2) q[8],q[15];
rzz(pi/2) q[8],q[16];
rzz(pi/2) q[8],q[16];
rzz(pi/2) q[8],q[17];
rzz(-pi/2) q[8],q[17];
rzz(pi/2) q[8],q[18];
rzz(pi/2) q[8],q[18];
u3(pi/2,5.2489730056178265,-5.2489730056178265) q[9];
u3(pi/2,6.178884431080405,-6.178884431080405) q[9];
rzz(pi/2) q[9],q[10];
u3(pi/2,0.6766990575832414,-0.6766990575832414) q[10];
u3(pi/2,3.0328935477755863,-3.0328935477755863) q[10];
rzz(-pi/2) q[9],q[10];
rzz(pi/2) q[9],q[11];
u3(pi/2,4.850619057142641,-4.850619057142641) q[11];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[11];
rzz(pi/2) q[9],q[11];
rzz(pi/2) q[9],q[12];
u3(pi/2,2.5283537676090653,-2.5283537676090653) q[12];
u3(pi/2,5.4732827210841375,-5.4732827210841375) q[12];
rzz(pi/2) q[9],q[12];
rzz(-pi/2) q[9],q[13];
u3(pi/2,2.0005662018059804,-2.0005662018059804) q[13];
u3(pi/2,5.044141164603771,-5.044141164603771) q[13];
rzz(pi/2) q[9],q[13];
rzz(pi/2) q[9],q[14];
u3(pi/2,5.345105740817674,-5.345105740817674) q[14];
u3(pi/2,2.15450424183188,-2.15450424183188) q[14];
rzz(pi/2) q[9],q[14];
rzz(-pi/2) q[9],q[15];
u3(pi/2,2.0772210625535714,-2.0772210625535714) q[15];
u3(pi/2,5.194309293445364,-5.194309293445364) q[15];
rzz(pi/2) q[9],q[15];
rzz(pi/2) q[9],q[16];
u3(pi/2,4.482424398141917,-4.482424398141917) q[16];
u3(pi/2,1.3282653739377646,-1.3282653739377646) q[16];
rzz(pi/2) q[9],q[16];
rzz(pi/2) q[9],q[17];
rzz(-pi/2) q[9],q[17];
rzz(pi/2) q[9],q[18];
rzz(pi/2) q[9],q[18];
u3(pi/2,4.603689874570483,-4.603689874570483) q[10];
u3(pi/2,4.924132325236641,-4.924132325236641) q[10];
rzz(pi/2) q[10],q[11];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[11];
u3(pi/2,3.6725218120464684,-3.6725218120464684) q[11];
rzz(-pi/2) q[10],q[11];
rzz(pi/2) q[10],q[12];
u3(pi/2,5.4732827210841375,-5.4732827210841375) q[12];
u3(pi/2,1.9389909857956202,-1.9389909857956202) q[12];
rzz(pi/2) q[10],q[12];
rzz(pi/2) q[10],q[13];
u3(pi/2,5.044141164603771,-5.044141164603771) q[13];
u3(pi/2,1.7058848108992577,-1.7058848108992577) q[13];
rzz(-pi/2) q[10],q[13];
rzz(pi/2) q[10],q[14];
u3(pi/2,2.15450424183188,-2.15450424183188) q[14];
u3(pi/2,5.198079204629672,-5.198079204629672) q[14];
rzz(pi/2) q[10],q[14];
rzz(pi/2) q[10],q[15];
u3(pi/2,5.194309293445364,-5.194309293445364) q[15];
u3(pi/2,2.00370779445957,-2.00370779445957) q[15];
rzz(-pi/2) q[10],q[15];
rzz(pi/2) q[10],q[16];
u3(pi/2,1.3282653739377646,-1.3282653739377646) q[16];
u3(pi/2,4.445353604829557,-4.445353604829557) q[16];
rzz(pi/2) q[10],q[16];
rzz(pi/2) q[10],q[17];
u3(pi/2,2.0891591146372126,-2.0891591146372126) q[17];
u3(pi/2,5.218813716143364,-5.218813716143364) q[17];
rzz(pi/2) q[10],q[17];
rzz(pi/2) q[10],q[18];
rzz(pi/2) q[10],q[18];
u3(pi/2,5.2433181388413646,-5.2433181388413646) q[11];
u3(pi/2,6.18956584610261,-6.18956584610261) q[11];
rzz(pi/2) q[11],q[12];
u3(pi/2,1.9389909857956202,-1.9389909857956202) q[12];
u3(pi/2,4.295185475987965,-4.295185475987965) q[12];
rzz(pi/2) q[11],q[12];
rzz(-pi/2) q[11],q[13];
u3(pi/2,1.7058848108992577,-1.7058848108992577) q[13];
u3(pi/2,4.454778382790327,-4.454778382790327) q[13];
rzz(pi/2) q[11],q[13];
rzz(pi/2) q[11],q[14];
u3(pi/2,5.198079204629672,-5.198079204629672) q[14];
u3(pi/2,1.8598228509251575,-1.8598228509251575) q[14];
rzz(pi/2) q[11],q[14];
rzz(-pi/2) q[11],q[15];
u3(pi/2,2.00370779445957,-2.00370779445957) q[15];
u3(pi/2,5.046654438726644,-5.046654438726644) q[15];
rzz(pi/2) q[11],q[15];
rzz(pi/2) q[11],q[16];
u3(pi/2,4.445353604829557,-4.445353604829557) q[16];
u3(pi/2,1.2547521058437634,-1.2547521058437634) q[16];
rzz(pi/2) q[11],q[16];
rzz(pi/2) q[11],q[17];
u3(pi/2,5.218813716143364,-5.218813716143364) q[17];
u3(pi/2,2.052088321324853,-2.052088321324853) q[17];
rzz(-pi/2) q[11],q[17];
rzz(pi/2) q[11],q[18];
u3(pi/2,5.228238494104134,-5.228238494104134) q[18];
u3(pi/2,2.0740794698999814,-2.0740794698999814) q[18];
rzz(pi/2) q[11],q[18];
u3(pi/2,2.7243891491930685,-2.7243891491930685) q[12];
u3(pi/2,3.98291116622114,-3.98291116622114) q[12];
rzz(pi/2) q[12],q[13];
u3(pi/2,4.454778382790327,-4.454778382790327) q[13];
u3(pi/2,0.5277875658030853,-0.5277875658030853) q[13];
rzz(-pi/2) q[12],q[13];
rzz(pi/2) q[12],q[14];
u3(pi/2,1.8598228509251575,-1.8598228509251575) q[14];
u3(pi/2,4.608716422816227,-4.608716422816227) q[14];
rzz(pi/2) q[12],q[14];
rzz(pi/2) q[12],q[15];
u3(pi/2,5.046654438726644,-5.046654438726644) q[15];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[15];
rzz(-pi/2) q[12],q[15];
rzz(pi/2) q[12],q[16];
u3(pi/2,1.2547521058437634,-1.2547521058437634) q[16];
u3(pi/2,4.298327068641555,-4.298327068641555) q[16];
rzz(pi/2) q[12],q[16];
rzz(pi/2) q[12],q[17];
u3(pi/2,5.193680974914646,-5.193680974914646) q[17];
u3(pi/2,2.003079475928852,-2.003079475928852) q[17];
rzz(-pi/2) q[12],q[17];
rzz(pi/2) q[12],q[18];
u3(pi/2,2.0740794698999814,-2.0740794698999814) q[18];
u3(pi/2,5.191167700791774,-5.191167700791774) q[18];
rzz(pi/2) q[12],q[18];
u3(pi/2,5.240176546187775,-5.240176546187775) q[13];
u3(pi/2,0.684238879951857,-0.684238879951857) q[13];
rzz(pi/2) q[13],q[14];
u3(pi/2,4.608716422816227,-4.608716422816227) q[14];
u3(pi/2,0.6817256058289851,-0.6817256058289851) q[14];
rzz(pi/2) q[13],q[14];
rzz(-pi/2) q[13],q[15];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[15];
u3(pi/2,4.457919975443916,-4.457919975443916) q[15];
rzz(pi/2) q[13],q[15];
rzz(pi/2) q[13],q[16];
u3(pi/2,4.298327068641555,-4.298327068641555) q[16];
u3(pi/2,0.9600707149370408,-0.9600707149370408) q[16];
rzz(pi/2) q[13],q[16];
rzz(-pi/2) q[13],q[17];
u3(pi/2,2.003079475928852,-2.003079475928852) q[17];
u3(pi/2,5.046654438726644,-5.046654438726644) q[17];
rzz(pi/2) q[13],q[17];
rzz(pi/2) q[13],q[18];
u3(pi/2,5.191167700791774,-5.191167700791774) q[18];
u3(pi/2,2.0005662018059804,-2.0005662018059804) q[18];
rzz(pi/2) q[13],q[18];
u3(pi/2,5.394114586213675,-5.394114586213675) q[14];
u3(pi/2,2.1746104348148547,-2.1746104348148547) q[14];
rzz(pi/2) q[14],q[15];
u3(pi/2,4.457919975443916,-4.457919975443916) q[15];
u3(pi/2,0.5309291584566751,-0.5309291584566751) q[15];
rzz(-pi/2) q[14],q[15];
rzz(pi/2) q[14],q[16];
u3(pi/2,0.9600707149370408,-0.9600707149370408) q[16];
u3(pi/2,3.70896428682811,-3.70896428682811) q[16];
rzz(pi/2) q[14],q[16];
rzz(pi/2) q[14],q[17];
u3(pi/2,5.046654438726644,-5.046654438726644) q[17];
u3(pi/2,1.7090264035528475,-1.7090264035528475) q[17];
rzz(-pi/2) q[14],q[17];
rzz(pi/2) q[14],q[18];
u3(pi/2,2.0005662018059804,-2.0005662018059804) q[18];
u3(pi/2,5.044141164603771,-5.044141164603771) q[18];
rzz(pi/2) q[14],q[18];
u3(pi/2,5.2433181388413646,-5.2433181388413646) q[15];
u3(pi/2,6.067672051143327,-6.067672051143327) q[15];
rzz(pi/2) q[15],q[16];
u3(pi/2,3.70896428682811,-3.70896428682811) q[16];
u3(pi/2,6.0651587770204545,-6.0651587770204545) q[16];
rzz(pi/2) q[15],q[16];
rzz(pi/2) q[15],q[17];
u3(pi/2,4.850619057142641,-4.850619057142641) q[17];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[17];
rzz(pi/2) q[15],q[17];
rzz(pi/2) q[15],q[18];
u3(pi/2,5.044141164603771,-5.044141164603771) q[18];
u3(pi/2,1.7065131294299756,-1.7065131294299756) q[18];
rzz(-pi/2) q[15],q[18];
u3(pi/2,4.494362450225558,-4.494362450225558) q[16];
u3(pi/2,4.867583657472026,-4.867583657472026) q[16];
rzz(pi/2) q[16],q[17];
u3(pi/2,1.3163273218541232,-1.3163273218541232) q[17];
u3(pi/2,3.6725218120464684,-3.6725218120464684) q[17];
rzz(pi/2) q[16],q[17];
rzz(pi/2) q[16],q[18];
u3(pi/2,4.8481057830197685,-4.8481057830197685) q[18];
u3(pi/2,1.3138140477312514,-1.3138140477312514) q[18];
rzz(pi/2) q[16],q[18];
u3(pi/2,2.101725485251572,-2.101725485251572) q[17];
u3(pi/2,3.0737342522722537,-3.0737342522722537) q[17];
rzz(pi/2) q[17],q[18];
u3(pi/2,1.3138140477312514,-1.3138140477312514) q[18];
u3(pi/2,3.670008537923596,-3.670008537923596) q[18];
rzz(pi/2) q[17],q[18];
u3(pi,0,0) q[1];
u3(pi,pi,-pi) q[2];
u3(pi,5*pi/4,-5*pi/4) q[3];
u3(pi,15*pi/8,-15*pi/8) q[4];
u3(pi,2.552229871776348,-2.552229871776348) q[5];
u3(pi,0.7363893180014475,-0.7363893180014475) q[6];
u3(pi,2.184663531306342,-2.184663531306342) q[7];
u3(pi,3.3256899830901547,-3.3256899830901547) q[8];
u3(pi,2.7118227785787092,-2.7118227785787092) q[9];
u3(pi,0.8344070087934491,-0.8344070087934491) q[10];
u3(pi,4.941096925566026,-4.941096925566026) q[11];
u3(pi,1.7856812643004385,-1.7856812643004385) q[12];
u3(pi,5.880433128989375,-5.880433128989375) q[13];
u3(pi,2.413999795018397,-2.413999795018397) q[14];
u3(pi,2.271999807076138,-2.271999807076138) q[15];
u3(pi,5.323742910773264,-5.323742910773264) q[16];
u3(pi,0.1715309588860027,-0.1715309588860027) q[17];
u3(pi/2,2.0992122111287,-2.0992122111287) q[18];
u3(pi/2,3.37030059877113,-3.37030059877113) q[18];
u3(pi,2.7834510910805568,-2.7834510910805568) q[19];
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
