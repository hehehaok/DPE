close all;
clear all;
clc;

%% ECEF <-> LLA
truePosInLLA = [60.161086788889, 24.545448080556, 54.1640000026673];
truePosInECEF = wgslla2xyz(truePosInLLA(1), truePosInLLA(2), truePosInLLA(3));
[lat, lon, alt] = wgsxyz2lla(truePosInECEF);

%% ECEF <-> ENU
testPosInECEF = [2894080.172; 1321678.105; 5509498.837];
refState = truePosInECEF;
curState = testPosInECEF;
coorInENU = ECEF2ENU(refState, curState);

%% corr_pos≤‚ ‘

