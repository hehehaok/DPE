%% ����������ջ��������������
clear;
close all;
fclose('all');
clc;
eph = load('eph.mat');

%% ������������
attr = {'weekNumber', 'accuracy', 'health', 'T_GD', 't_oc', 'a_f2',...
        'a_f1', 'a_f0', 'C_rs', 'deltan', 'M_0', 'C_uc', 'e', 'C_us',...
        'sqrtA', 't_oe', 'C_ic', 'omega_0', 'C_is', 'i_0',...
        'IODC', 'C_rc', 'omega', 'omegaDot', 'iDot'}; % 'IODE'�� 'IODE_sf2' 'IODE_sf3'��ʲô��ϵ��
for ii = [13 15 17 24 28]
    ephChannel = eph.eph(ii);
    load(['channel' num2str(ii) '.mat']);
    % ������ͳһ
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
        fprintf('ͨ��%d����������%s�ıȽϽ��Ϊ%f,%f,���߲�ֵΪ%f\n', ii, attr{jj}, attrTemp1, attrTemp2, attrDiff);
    end
    fprintf('-------------------------------------------------------------------------------------\n');
end