OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u3(0,-1.40658297059163,1.4065829705916304) q[0];
u3(pi/2,-pi/4,-pi) q[1];
cx q[0],q[1];
u3(0,0,pi/4) q[1];
cx q[0],q[1];
u3(pi/2,0,-pi) q[1];
u3(pi/2,-pi/8,-pi) q[2];
cx q[0],q[2];
u3(0,0,pi/8) q[2];
cx q[0],q[2];
u3(0,0,-pi/4) q[2];
cx q[1],q[2];
u3(0,0,pi/4) q[2];
cx q[1],q[2];
u3(pi/2,0,-pi) q[2];
u3(pi/2,-pi/16,-pi) q[3];
cx q[0],q[3];
u3(0,0,pi/16) q[3];
cx q[0],q[3];
u3(0,0,-pi/8) q[3];
cx q[1],q[3];
u3(0,0,pi/8) q[3];
cx q[1],q[3];
u3(0,0,-pi/4) q[3];
cx q[2],q[3];
u3(0,0,pi/4) q[3];
cx q[2],q[3];
u3(pi/2,0,-pi) q[3];
u3(pi/2,-pi,-pi) q[4];
u3(pi/2,-pi/2,-pi) q[5];
u3(pi/2,-pi/4,-pi) q[6];
u3(pi/2,7*pi/8,-pi) q[7];
u3(pi/2,-9*pi/16,-pi) q[8];
u3(pi/2,2.258019719767665,-pi) q[9];
u3(pi/2,-2.0125827937059606,-pi) q[10];
u3(pi/2,-1.0062913968529799,-pi) q[11];
u3(pi/2,-0.5031456984264899,-pi) q[12];
u3(pi/2,2.890019804376548,-pi) q[13];
u3(pi/2,1.445009902188274,-pi) q[14];
u3(-pi,-pi/2,pi/2) q[15];
cx q[14],q[15];
u3(0,0,-1.44500990218827) q[15];
cx q[14],q[15];
u3(0,0,-9.58737992428526e-5) q[14];
u3(0,0,1.44500990218827) q[15];
cx q[13],q[15];
u3(0,0,-2.89001980437655) q[15];
cx q[13],q[15];
u3(0,0,-0.000191747598485705) q[13];
u3(0,0,2.89001980437655) q[15];
cx q[12],q[15];
u3(0,0,-5.7800396087531) q[15];
cx q[12],q[15];
u3(0,0,-0.00038349519697141) q[12];
u3(0,0,5.7800396087531) q[15];
cx q[11],q[15];
u3(0,0,-11.5600792175062) q[15];
cx q[11],q[15];
u3(0,0,-0.000766990393942821) q[11];
u3(0,0,11.5600792175062) q[15];
cx q[10],q[15];
u3(0,0,-23.1201584350124) q[15];
cx q[10],q[15];
u3(0,0,-0.00153398078788564) q[10];
u3(0,0,23.1201584350124) q[15];
cx q[9],q[15];
u3(0,0,-46.2403168700248) q[15];
cx q[9],q[15];
u3(0,0,46.2403168700248) q[15];
cx q[8],q[15];
u3(0,0,-92.4806337400495) q[15];
cx q[8],q[15];
u3(0,0,92.4806337400495) q[15];
cx q[7],q[15];
u3(0,0,-184.961267480099) q[15];
cx q[7],q[15];
u3(0,0,184.961267480099) q[15];
cx q[6],q[15];
u3(0,0,-369.922534960198) q[15];
cx q[6],q[15];
u3(0,0,369.922534960198) q[15];
cx q[5],q[15];
u3(0,0,-739.845069920396) q[15];
cx q[5],q[15];
u3(0,0,739.845069920396) q[15];
cx q[4],q[15];
u3(0,0,-471*pi) q[15];
cx q[4],q[15];
u3(0,0,471*pi) q[15];
u3(0,0,-pi/32) q[4];
cx q[0],q[4];
u3(0,0,pi/32) q[4];
cx q[0],q[4];
u3(0,0,-pi/16) q[4];
cx q[1],q[4];
u3(0,0,pi/16) q[4];
cx q[1],q[4];
u3(0,0,-pi/8) q[4];
cx q[2],q[4];
u3(0,0,pi/8) q[4];
cx q[2],q[4];
u3(0,0,-pi/4) q[4];
cx q[3],q[4];
u3(0,0,pi/4) q[4];
cx q[3],q[4];
u3(pi/2,0,-pi) q[4];
u3(0,0,-pi/64) q[5];
cx q[0],q[5];
u3(0,0,pi/64) q[5];
cx q[0],q[5];
u3(0,0,-pi/32) q[5];
cx q[1],q[5];
u3(0,0,pi/32) q[5];
cx q[1],q[5];
u3(0,0,-pi/16) q[5];
cx q[2],q[5];
u3(0,0,pi/16) q[5];
cx q[2],q[5];
u3(0,0,-pi/8) q[5];
cx q[3],q[5];
u3(0,0,pi/8) q[5];
cx q[3],q[5];
u3(0,0,-pi/4) q[5];
cx q[4],q[5];
u3(0,0,pi/4) q[5];
cx q[4],q[5];
u3(pi/2,0,-pi) q[5];
u3(0,0,-pi/128) q[6];
cx q[0],q[6];
u3(0,0,pi/128) q[6];
cx q[0],q[6];
u3(0,0,-pi/64) q[6];
cx q[1],q[6];
u3(0,0,pi/64) q[6];
cx q[1],q[6];
u3(0,0,-pi/32) q[6];
cx q[2],q[6];
u3(0,0,pi/32) q[6];
cx q[2],q[6];
u3(0,0,-pi/16) q[6];
cx q[3],q[6];
u3(0,0,pi/16) q[6];
cx q[3],q[6];
u3(0,0,-pi/8) q[6];
cx q[4],q[6];
u3(0,0,pi/8) q[6];
cx q[4],q[6];
u3(0,0,-pi/4) q[6];
cx q[5],q[6];
u3(0,0,pi/4) q[6];
cx q[5],q[6];
u3(pi/2,0,-pi) q[6];
u3(0,0,-pi/256) q[7];
cx q[0],q[7];
u3(0,0,pi/256) q[7];
cx q[0],q[7];
u3(0,0,-pi/128) q[7];
cx q[1],q[7];
u3(0,0,pi/128) q[7];
cx q[1],q[7];
u3(0,0,-pi/64) q[7];
cx q[2],q[7];
u3(0,0,pi/64) q[7];
cx q[2],q[7];
u3(0,0,-pi/32) q[7];
cx q[3],q[7];
u3(0,0,pi/32) q[7];
cx q[3],q[7];
u3(0,0,-pi/16) q[7];
cx q[4],q[7];
u3(0,0,pi/16) q[7];
cx q[4],q[7];
u3(0,0,-pi/8) q[7];
cx q[5],q[7];
u3(0,0,pi/8) q[7];
cx q[5],q[7];
u3(0,0,-pi/4) q[7];
cx q[6],q[7];
u3(0,0,pi/4) q[7];
cx q[6],q[7];
u3(pi/2,0,-pi) q[7];
u3(0,0,-pi/512) q[8];
cx q[0],q[8];
u3(0,0,pi/512) q[8];
cx q[0],q[8];
u3(0,0,-pi/256) q[8];
cx q[1],q[8];
u3(0,0,pi/256) q[8];
cx q[1],q[8];
u3(0,0,-pi/128) q[8];
cx q[2],q[8];
u3(0,0,pi/128) q[8];
cx q[2],q[8];
u3(0,0,-pi/64) q[8];
cx q[3],q[8];
u3(0,0,pi/64) q[8];
cx q[3],q[8];
u3(0,0,-pi/32) q[8];
cx q[4],q[8];
u3(0,0,pi/32) q[8];
cx q[4],q[8];
u3(0,0,-pi/16) q[8];
cx q[5],q[8];
u3(0,0,pi/16) q[8];
cx q[5],q[8];
u3(0,0,-pi/8) q[8];
cx q[6],q[8];
u3(0,0,pi/8) q[8];
cx q[6],q[8];
u3(0,0,-pi/4) q[8];
cx q[7],q[8];
u3(0,0,pi/4) q[8];
cx q[7],q[8];
u3(pi/2,0,-pi) q[8];
u3(0,0,-pi/1024) q[9];
cx q[0],q[9];
u3(0,0,pi/1024) q[9];
cx q[0],q[9];
cx q[0],q[10];
u3(0,0,0.00153398078788564) q[10];
cx q[0],q[10];
cx q[0],q[11];
u3(0,0,-pi/1024) q[10];
u3(0,0,0.000766990393942821) q[11];
cx q[0],q[11];
cx q[0],q[12];
u3(0,0,-0.00153398078788564) q[11];
u3(0,0,0.00038349519697141) q[12];
cx q[0],q[12];
cx q[0],q[13];
u3(0,0,-0.000766990393942821) q[12];
u3(0,0,0.000191747598485705) q[13];
cx q[0],q[13];
cx q[0],q[14];
u3(0,0,-0.00038349519697141) q[13];
u3(0,0,9.58737992428526e-5) q[14];
cx q[0],q[14];
u3(0,0,-0.000191747598485705) q[14];
u3(0,0,-pi/512) q[9];
cx q[1],q[9];
u3(0,0,pi/512) q[9];
cx q[1],q[9];
cx q[1],q[10];
u3(0,0,pi/1024) q[10];
cx q[1],q[10];
cx q[1],q[11];
u3(0,0,-pi/512) q[10];
u3(0,0,0.00153398078788564) q[11];
cx q[1],q[11];
cx q[1],q[12];
u3(0,0,-pi/1024) q[11];
u3(0,0,0.000766990393942821) q[12];
cx q[1],q[12];
cx q[1],q[13];
u3(0,0,-0.00153398078788564) q[12];
u3(0,0,0.00038349519697141) q[13];
cx q[1],q[13];
cx q[1],q[14];
u3(0,0,-0.000766990393942821) q[13];
u3(0,0,0.000191747598485705) q[14];
cx q[1],q[14];
u3(0,0,-0.00038349519697141) q[14];
u3(0,0,-pi/256) q[9];
cx q[2],q[9];
u3(0,0,pi/256) q[9];
cx q[2],q[9];
cx q[2],q[10];
u3(0,0,pi/512) q[10];
cx q[2],q[10];
u3(0,0,-pi/256) q[10];
cx q[2],q[11];
u3(0,0,pi/1024) q[11];
cx q[2],q[11];
u3(0,0,-pi/512) q[11];
cx q[2],q[12];
u3(0,0,0.00153398078788564) q[12];
cx q[2],q[12];
u3(0,0,-pi/1024) q[12];
cx q[2],q[13];
u3(0,0,0.000766990393942821) q[13];
cx q[2],q[13];
u3(0,0,-0.00153398078788564) q[13];
cx q[2],q[14];
u3(0,0,0.00038349519697141) q[14];
cx q[2],q[14];
u3(0,0,-0.000766990393942821) q[14];
u3(0,0,-pi/128) q[9];
cx q[3],q[9];
u3(0,0,pi/128) q[9];
cx q[3],q[9];
cx q[3],q[10];
u3(0,0,pi/256) q[10];
cx q[3],q[10];
u3(0,0,-pi/128) q[10];
cx q[3],q[11];
u3(0,0,pi/512) q[11];
cx q[3],q[11];
u3(0,0,-pi/256) q[11];
cx q[3],q[12];
u3(0,0,pi/1024) q[12];
cx q[3],q[12];
u3(0,0,-pi/512) q[12];
cx q[3],q[13];
u3(0,0,0.00153398078788564) q[13];
cx q[3],q[13];
u3(0,0,-pi/1024) q[13];
cx q[3],q[14];
u3(0,0,0.000766990393942821) q[14];
cx q[3],q[14];
u3(0,0,-0.00153398078788564) q[14];
u3(0,0,-pi/64) q[9];
cx q[4],q[9];
u3(0,0,pi/64) q[9];
cx q[4],q[9];
cx q[4],q[10];
u3(0,0,pi/128) q[10];
cx q[4],q[10];
u3(0,0,-pi/64) q[10];
cx q[4],q[11];
u3(0,0,pi/256) q[11];
cx q[4],q[11];
u3(0,0,-pi/128) q[11];
cx q[4],q[12];
u3(0,0,pi/512) q[12];
cx q[4],q[12];
u3(0,0,-pi/256) q[12];
cx q[4],q[13];
u3(0,0,pi/1024) q[13];
cx q[4],q[13];
u3(0,0,-pi/512) q[13];
cx q[4],q[14];
u3(0,0,0.00153398078788564) q[14];
cx q[4],q[14];
u3(0,0,-pi/1024) q[14];
u3(0,0,-pi/32) q[9];
cx q[5],q[9];
u3(0,0,pi/32) q[9];
cx q[5],q[9];
cx q[5],q[10];
u3(0,0,pi/64) q[10];
cx q[5],q[10];
u3(0,0,-pi/32) q[10];
cx q[5],q[11];
u3(0,0,pi/128) q[11];
cx q[5],q[11];
u3(0,0,-pi/64) q[11];
cx q[5],q[12];
u3(0,0,pi/256) q[12];
cx q[5],q[12];
u3(0,0,-pi/128) q[12];
cx q[5],q[13];
u3(0,0,pi/512) q[13];
cx q[5],q[13];
u3(0,0,-pi/256) q[13];
cx q[5],q[14];
u3(0,0,pi/1024) q[14];
cx q[5],q[14];
u3(0,0,-pi/512) q[14];
u3(0,0,-pi/16) q[9];
cx q[6],q[9];
u3(0,0,pi/16) q[9];
cx q[6],q[9];
cx q[6],q[10];
u3(0,0,pi/32) q[10];
cx q[6],q[10];
u3(0,0,-pi/16) q[10];
cx q[6],q[11];
u3(0,0,pi/64) q[11];
cx q[6],q[11];
u3(0,0,-pi/32) q[11];
cx q[6],q[12];
u3(0,0,pi/128) q[12];
cx q[6],q[12];
u3(0,0,-pi/64) q[12];
cx q[6],q[13];
u3(0,0,pi/256) q[13];
cx q[6],q[13];
u3(0,0,-pi/128) q[13];
cx q[6],q[14];
u3(0,0,pi/512) q[14];
cx q[6],q[14];
u3(0,0,-pi/256) q[14];
u3(0,0,-pi/8) q[9];
cx q[7],q[9];
u3(0,0,pi/8) q[9];
cx q[7],q[9];
cx q[7],q[10];
u3(0,0,pi/16) q[10];
cx q[7],q[10];
u3(0,0,-pi/8) q[10];
cx q[7],q[11];
u3(0,0,pi/32) q[11];
cx q[7],q[11];
u3(0,0,-pi/16) q[11];
cx q[7],q[12];
u3(0,0,pi/64) q[12];
cx q[7],q[12];
u3(0,0,-pi/32) q[12];
cx q[7],q[13];
u3(0,0,pi/128) q[13];
cx q[7],q[13];
u3(0,0,-pi/64) q[13];
cx q[7],q[14];
u3(0,0,pi/256) q[14];
cx q[7],q[14];
u3(0,0,-pi/128) q[14];
u3(0,0,-pi/4) q[9];
cx q[8],q[9];
u3(0,0,pi/4) q[9];
cx q[8],q[9];
cx q[8],q[10];
u3(0,0,pi/8) q[10];
cx q[8],q[10];
u3(0,0,-pi/4) q[10];
cx q[8],q[11];
u3(0,0,pi/16) q[11];
cx q[8],q[11];
u3(0,0,-pi/8) q[11];
cx q[8],q[12];
u3(0,0,pi/32) q[12];
cx q[8],q[12];
u3(0,0,-pi/16) q[12];
cx q[8],q[13];
u3(0,0,pi/64) q[13];
cx q[8],q[13];
u3(0,0,-pi/32) q[13];
cx q[8],q[14];
u3(0,0,pi/128) q[14];
cx q[8],q[14];
u3(0,0,-pi/64) q[14];
u3(pi/2,0,-pi) q[9];
cx q[9],q[10];
u3(0,0,pi/4) q[10];
cx q[9],q[10];
u3(pi/2,0,-pi) q[10];
cx q[9],q[11];
u3(0,0,pi/8) q[11];
cx q[9],q[11];
u3(0,0,-pi/4) q[11];
cx q[10],q[11];
u3(0,0,pi/4) q[11];
cx q[10],q[11];
u3(pi/2,0,-pi) q[11];
cx q[9],q[12];
u3(0,0,pi/16) q[12];
cx q[9],q[12];
u3(0,0,-pi/8) q[12];
cx q[10],q[12];
u3(0,0,pi/8) q[12];
cx q[10],q[12];
u3(0,0,-pi/4) q[12];
cx q[11],q[12];
u3(0,0,pi/4) q[12];
cx q[11],q[12];
u3(pi/2,0,-pi) q[12];
cx q[9],q[13];
u3(0,0,pi/32) q[13];
cx q[9],q[13];
u3(0,0,-pi/16) q[13];
cx q[10],q[13];
u3(0,0,pi/16) q[13];
cx q[10],q[13];
u3(0,0,-pi/8) q[13];
cx q[11],q[13];
u3(0,0,pi/8) q[13];
cx q[11],q[13];
u3(0,0,-pi/4) q[13];
cx q[12],q[13];
u3(0,0,pi/4) q[13];
cx q[12],q[13];
u3(pi/2,0,-pi) q[13];
cx q[9],q[14];
u3(0,0,pi/64) q[14];
cx q[9],q[14];
u3(0,0,-pi/32) q[14];
cx q[10],q[14];
u3(0,0,pi/32) q[14];
cx q[10],q[14];
u3(0,0,-pi/16) q[14];
cx q[11],q[14];
u3(0,0,pi/16) q[14];
cx q[11],q[14];
u3(0,0,-pi/8) q[14];
cx q[12],q[14];
u3(0,0,pi/8) q[14];
cx q[12],q[14];
u3(0,0,-pi/4) q[14];
cx q[13],q[14];
u3(0,0,pi/4) q[14];
cx q[13],q[14];
u3(pi/2,0,-pi) q[14];
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
