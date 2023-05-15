%% 载入芬兰接收机解算的星历参数
clear;
close all;
fclose('all');
clc;
eph = load('eph.mat');

%% 星历参数名称
attr = {'weekNumber', 'accuracy', 'health', 'T_GD', 't_oc', 'a_f2',...
        'a_f1', 'a_f0', 'C_rs', 'deltan', 'M_0', 'C_uc', 'e', 'C_us',...
        'sqrtA', 't_oe', 'C_ic', 'omega_0', 'C_is', 'i_0',...
        'IODC', 'C_rc', 'omega', 'omegaDot', 'iDot'}; % 'IODE'和 'IODE_sf2' 'IODE_sf3'有什么关系呢
for ii = [13 15 17 24 28]
    ephChannel = eph.eph(ii);
    load(['channel' num2str(ii) '.mat']);
    % 变量名统一
    weekNumber = weeknumber;
    deltan = delta_n;
    sqrtA = sqrt_A;
    omega_0 = OMEGA_0;
    omegaDot = OMEGADOT;
    iDot = IDOT;
    clearvars weeknumber delta_n sqrt_A OMEGA_0 OMEGADOT IDOT;
    for jj = 1 : length(attr)
        attrTemp1 = ephChannel.(attr{jj});
        attrTemp2 = eval(attr{jj});
        attrDiff = attrTemp1 - attrTemp2;
        fprintf('通道%d的星历参数%s的比较结果为%f,%f,两者差值为%f\n', ii, attr{jj}, attrTemp1, attrTemp2, attrDiff);
    end
    fprintf('-------------------------------------------------------------------------------------\n');
end