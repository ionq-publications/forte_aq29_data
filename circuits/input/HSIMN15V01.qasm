OPENQASM 2.0;
include "qelib1.inc";
qreg q[15];
creg c[15];
rx(-pi) q[0];
rx(-pi) q[2];
rx(-pi) q[4];
rx(-pi) q[6];
rx(-pi) q[8];
rx(-pi) q[10];
rx(-pi) q[12];
rx(-pi) q[14];
rx(-0.025608672983337577) q[0];
rz(0.00919202934923197) q[0];
rx(0.0025855305740396517) q[1];
rz(-0.05787143632140813) q[1];
rx(-0.04270063376987876) q[10];
rz(0.019752160323741552) q[10];
rx(0.005176834476686487) q[11];
rz(0.02290549031720591) q[11];
rx(-0.028964586051869468) q[12];
rz(0.057543142376577025) q[12];
rx(0.034791982594946536) q[13];
rz(0.022665078055759496) q[13];
rx(-0.030707607780540158) q[14];
rz(-0.035716731816981415) q[14];
rx(0.03577302153839668) q[2];
rz(0.051748089320033504) q[2];
rx(0.038562764822781315) q[3];
rz(0.05728775695654175) q[3];
rx(0.049408274996776136) q[4];
rz(0.023958938850150217) q[4];
rx(-0.04161048111705812) q[5];
rz(0.005903442189977411) q[5];
rx(-0.0307326339177445) q[6];
rz(-0.041620087712266635) q[6];
rx(-0.0005077145349519085) q[7];
rz(0.015768461624668113) q[7];
rx(0.031882899615161) q[8];
rz(-0.01678296421528314) q[8];
rx(-0.04067306778111446) q[9];
rz(0.057810955769588945) q[9];
rz(1.5708) q[1];
cx q[1],q[0];
rz(-1.560328) q[0];
ry(1.560328) q[1];
cx q[0],q[1];
ry(-1.560328) q[1];
cx q[1],q[0];
rz(-1.5708) q[0];
rz(1.5708) q[11];
cx q[11],q[10];
rz(-1.560328) q[10];
ry(1.560328) q[11];
cx q[10],q[11];
ry(-1.560328) q[11];
cx q[11],q[10];
rz(1.5708) q[13];
cx q[13],q[12];
rz(-1.560328) q[12];
ry(1.560328) q[13];
cx q[12],q[13];
ry(-1.560328) q[13];
cx q[13],q[12];
cx q[12],q[11];
rz(-1.560328) q[11];
ry(1.560328) q[12];
cx q[11],q[12];
ry(-1.560328) q[12];
cx q[12],q[11];
rz(-1.5708) q[11];
rz(1.5708) q[14];
cx q[14],q[13];
rz(-1.560328) q[13];
ry(1.560328) q[14];
cx q[13],q[14];
ry(-1.560328) q[14];
cx q[14],q[13];
rz(-1.5708) q[13];
rz(1.5708) q[3];
cx q[3],q[2];
rz(-1.560328) q[2];
ry(1.560328) q[3];
cx q[2],q[3];
ry(-1.560328) q[3];
cx q[3],q[2];
cx q[2],q[1];
rz(-1.560328) q[1];
ry(1.560328) q[2];
cx q[1],q[2];
ry(-1.560328) q[2];
cx q[2],q[1];
rz(-1.5708) q[1];
rz(1.5708) q[5];
cx q[5],q[4];
rz(-1.560328) q[4];
ry(1.560328) q[5];
cx q[4],q[5];
ry(-1.560328) q[5];
cx q[5],q[4];
cx q[4],q[3];
rz(-1.560328) q[3];
ry(1.560328) q[4];
cx q[3],q[4];
ry(-1.560328) q[4];
cx q[4],q[3];
rz(-1.5708) q[3];
rz(1.5708) q[7];
cx q[7],q[6];
rz(-1.560328) q[6];
ry(1.560328) q[7];
cx q[6],q[7];
ry(-1.560328) q[7];
cx q[7],q[6];
cx q[6],q[5];
rz(-1.560328) q[5];
ry(1.560328) q[6];
cx q[5],q[6];
ry(-1.560328) q[6];
cx q[6],q[5];
rz(-1.5708) q[5];
rz(1.5708) q[9];
cx q[9],q[8];
rz(-1.560328) q[8];
ry(1.560328) q[9];
cx q[8],q[9];
ry(-1.560328) q[9];
cx q[9],q[8];
cx q[10],q[9];
ry(1.560328) q[10];
cx q[8],q[7];
rz(-1.560328) q[7];
ry(1.560328) q[8];
cx q[7],q[8];
ry(-1.560328) q[8];
cx q[8],q[7];
rz(-1.5708) q[7];
rz(-1.560328) q[9];
cx q[9],q[10];
ry(-1.560328) q[10];
cx q[10],q[9];
rz(-1.5708) q[9];
rx(-0.025608672983337577) q[0];
rz(0.00919202934923197) q[0];
rx(0.0025855305740396517) q[1];
rz(-0.05787143632140813) q[1];
rx(-0.04270063376987876) q[10];
rz(0.019752160323741552) q[10];
rx(0.005176834476686487) q[11];
rz(0.02290549031720591) q[11];
rx(-0.028964586051869468) q[12];
rz(0.057543142376577025) q[12];
rx(0.034791982594946536) q[13];
rz(0.022665078055759496) q[13];
rx(-0.030707607780540158) q[14];
rz(-0.035716731816981415) q[14];
rx(0.03577302153839668) q[2];
rz(0.051748089320033504) q[2];
rx(0.038562764822781315) q[3];
rz(0.05728775695654175) q[3];
rx(0.049408274996776136) q[4];
rz(0.023958938850150217) q[4];
rx(-0.04161048111705812) q[5];
rz(0.005903442189977411) q[5];
rx(-0.0307326339177445) q[6];
rz(-0.041620087712266635) q[6];
rx(-0.0005077145349519085) q[7];
rz(0.015768461624668113) q[7];
rx(0.031882899615161) q[8];
rz(-0.01678296421528314) q[8];
rx(-0.04067306778111446) q[9];
rz(0.057810955769588945) q[9];
rz(1.5708) q[1];
cx q[1],q[0];
rz(-1.560328) q[0];
ry(1.560328) q[1];
cx q[0],q[1];
ry(-1.560328) q[1];
cx q[1],q[0];
rz(-1.5708) q[0];
rz(1.5708) q[11];
cx q[11],q[10];
rz(-1.560328) q[10];
ry(1.560328) q[11];
cx q[10],q[11];
ry(-1.560328) q[11];
cx q[11],q[10];
rz(1.5708) q[13];
cx q[13],q[12];
rz(-1.560328) q[12];
ry(1.560328) q[13];
cx q[12],q[13];
ry(-1.560328) q[13];
cx q[13],q[12];
cx q[12],q[11];
rz(-1.560328) q[11];
ry(1.560328) q[12];
cx q[11],q[12];
ry(-1.560328) q[12];
cx q[12],q[11];
rz(-1.5708) q[11];
rz(1.5708) q[14];
cx q[14],q[13];
rz(-1.560328) q[13];
ry(1.560328) q[14];
cx q[13],q[14];
ry(-1.560328) q[14];
cx q[14],q[13];
rz(-1.5708) q[13];
rz(1.5708) q[3];
cx q[3],q[2];
rz(-1.560328) q[2];
ry(1.560328) q[3];
cx q[2],q[3];
ry(-1.560328) q[3];
cx q[3],q[2];
cx q[2],q[1];
rz(-1.560328) q[1];
ry(1.560328) q[2];
cx q[1],q[2];
ry(-1.560328) q[2];
cx q[2],q[1];
rz(-1.5708) q[1];
rz(1.5708) q[5];
cx q[5],q[4];
rz(-1.560328) q[4];
ry(1.560328) q[5];
cx q[4],q[5];
ry(-1.560328) q[5];
cx q[5],q[4];
cx q[4],q[3];
rz(-1.560328) q[3];
ry(1.560328) q[4];
cx q[3],q[4];
ry(-1.560328) q[4];
cx q[4],q[3];
rz(-1.5708) q[3];
rz(1.5708) q[7];
cx q[7],q[6];
rz(-1.560328) q[6];
ry(1.560328) q[7];
cx q[6],q[7];
ry(-1.560328) q[7];
cx q[7],q[6];
cx q[6],q[5];
rz(-1.560328) q[5];
ry(1.560328) q[6];
cx q[5],q[6];
ry(-1.560328) q[6];
cx q[6],q[5];
rz(-1.5708) q[5];
rz(1.5708) q[9];
cx q[9],q[8];
rz(-1.560328) q[8];
ry(1.560328) q[9];
cx q[8],q[9];
ry(-1.560328) q[9];
cx q[9],q[8];
cx q[10],q[9];
ry(1.560328) q[10];
cx q[8],q[7];
rz(-1.560328) q[7];
ry(1.560328) q[8];
cx q[7],q[8];
ry(-1.560328) q[8];
cx q[8],q[7];
rz(-1.5708) q[7];
rz(-1.560328) q[9];
cx q[9],q[10];
ry(-1.560328) q[10];
cx q[10],q[9];
rz(-1.5708) q[9];
rx(-0.025608672983337577) q[0];
rz(0.00919202934923197) q[0];
rx(0.0025855305740396517) q[1];
rz(-0.05787143632140813) q[1];
rx(-0.04270063376987876) q[10];
rz(0.019752160323741552) q[10];
rx(0.005176834476686487) q[11];
rz(0.02290549031720591) q[11];
rx(-0.028964586051869468) q[12];
rz(0.057543142376577025) q[12];
rx(0.034791982594946536) q[13];
rz(0.022665078055759496) q[13];
rx(-0.030707607780540158) q[14];
rz(-0.035716731816981415) q[14];
rx(0.03577302153839668) q[2];
rz(0.051748089320033504) q[2];
rx(0.038562764822781315) q[3];
rz(0.05728775695654175) q[3];
rx(0.049408274996776136) q[4];
rz(0.023958938850150217) q[4];
rx(-0.04161048111705812) q[5];
rz(0.005903442189977411) q[5];
rx(-0.0307326339177445) q[6];
rz(-0.041620087712266635) q[6];
rx(-0.0005077145349519085) q[7];
rz(0.015768461624668113) q[7];
rx(0.031882899615161) q[8];
rz(-0.01678296421528314) q[8];
rx(-0.04067306778111446) q[9];
rz(0.057810955769588945) q[9];
rz(1.5708) q[1];
cx q[1],q[0];
rz(-1.560328) q[0];
ry(1.560328) q[1];
cx q[0],q[1];
ry(-1.560328) q[1];
cx q[1],q[0];
rz(-1.5708) q[0];
rz(1.5708) q[11];
cx q[11],q[10];
rz(-1.560328) q[10];
ry(1.560328) q[11];
cx q[10],q[11];
ry(-1.560328) q[11];
cx q[11],q[10];
rz(1.5708) q[13];
cx q[13],q[12];
rz(-1.560328) q[12];
ry(1.560328) q[13];
cx q[12],q[13];
ry(-1.560328) q[13];
cx q[13],q[12];
cx q[12],q[11];
rz(-1.560328) q[11];
ry(1.560328) q[12];
cx q[11],q[12];
ry(-1.560328) q[12];
cx q[12],q[11];
rz(-1.5708) q[11];
rz(1.5708) q[14];
cx q[14],q[13];
rz(-1.560328) q[13];
ry(1.560328) q[14];
cx q[13],q[14];
ry(-1.560328) q[14];
cx q[14],q[13];
rz(-1.5708) q[13];
rz(1.5708) q[3];
cx q[3],q[2];
rz(-1.560328) q[2];
ry(1.560328) q[3];
cx q[2],q[3];
ry(-1.560328) q[3];
cx q[3],q[2];
cx q[2],q[1];
rz(-1.560328) q[1];
ry(1.560328) q[2];
cx q[1],q[2];
ry(-1.560328) q[2];
cx q[2],q[1];
rz(-1.5708) q[1];
rz(1.5708) q[5];
cx q[5],q[4];
rz(-1.560328) q[4];
ry(1.560328) q[5];
cx q[4],q[5];
ry(-1.560328) q[5];
cx q[5],q[4];
cx q[4],q[3];
rz(-1.560328) q[3];
ry(1.560328) q[4];
cx q[3],q[4];
ry(-1.560328) q[4];
cx q[4],q[3];
rz(-1.5708) q[3];
rz(1.5708) q[7];
cx q[7],q[6];
rz(-1.560328) q[6];
ry(1.560328) q[7];
cx q[6],q[7];
ry(-1.560328) q[7];
cx q[7],q[6];
cx q[6],q[5];
rz(-1.560328) q[5];
ry(1.560328) q[6];
cx q[5],q[6];
ry(-1.560328) q[6];
cx q[6],q[5];
rz(-1.5708) q[5];
rz(1.5708) q[9];
cx q[9],q[8];
rz(-1.560328) q[8];
ry(1.560328) q[9];
cx q[8],q[9];
ry(-1.560328) q[9];
cx q[9],q[8];
cx q[10],q[9];
ry(1.560328) q[10];
cx q[8],q[7];
rz(-1.560328) q[7];
ry(1.560328) q[8];
cx q[7],q[8];
ry(-1.560328) q[8];
cx q[8],q[7];
rz(-1.5708) q[7];
rz(-1.560328) q[9];
cx q[9],q[10];
ry(-1.560328) q[10];
cx q[10],q[9];
rz(-1.5708) q[9];
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
