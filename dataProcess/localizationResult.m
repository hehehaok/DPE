clear;
close all;
fclose('all');
clc;

%% ����ѡ��
dataFlag = 2; % 0-���� 1-OAKBAT 2-TEXBAT 3-CPNTBAT
switch dataFlag
    case 0
        dataDir = 'finland_cleanStatic/';
        truePosInLLA = [60.161086788889, 24.545448080556, 54.1640000026673]; % ��ʵλ��(γ��(�����ϸ�)����(��������) �߶�) ����
        secDir = 'test0907_1_skip0s_proc40s/';
    case 1
        dataDir = 'oak_cleanStatic/';
        truePosInLLA = [35.930544444, -84.310652778, 248.6000]; % ��ʵλ��(γ��(�����ϸ�) ����(��������) �߶�) OAKBAT
        secDir = 'test10912_skip50s_proc40s/';
    case 2
        dataDir = 'texbat_ds4/';
        truePosInLLA = [30.287528140, -97.735720004, 166.1898]; % ��ʵλ��(γ��(�����ϸ�)����(��������) �߶�) TEXBAT
%         secDir = 'test20913_skip20s_proc40s/';
        secDir = 'test10921_skip120s_proc40s/'; 
    case 3
%         dataDir = 'cpntbat_cs08/';
        dataDir = 'CPNTBAT_cleanStatic/';
        truePosInLLA = [22.803448, 113.953065, 50.0000]; % ��ʵλ��(γ��(�����ϸ�) ����(��������) �߶�) CPNTBAT
        secDir = '0914test1_skip20s_proc40s/';
    otherwise
end

%% ��ȡDPE��λ����ļ�
predir   = ['./pre-simulator/' dataDir];
postdir  = ['./post-simulator/' dataDir];

% ����λ�ý�����
csvFile = 'usrp.csv';
csvDir = [postdir secDir csvFile];
csvData = importdata(csvDir);

% ����������ؽ��
receiveFile = 'end-of-dp/receiver.mat';
rxDir = [postdir secDir receiveFile];
load(rxDir);

% ��������ʱ����ز���
dpeRunTimeFile = 'dpe_runtime.mat';
load([postdir secDir dpeRunTimeFile]);

%% ��ʼ������
data = csvData.data; % ��λ�������
truePosInECEF = wgslla2xyz(truePosInLLA(1), truePosInLLA(2), truePosInLLA(3)); % ��ʵλ�� LLA -> ECEF
trueUtmZone = findUtmZone(truePosInLLA(1), truePosInLLA(2));
truePosInUTM = cart2utm(truePosInECEF(1), truePosInECEF(2), truePosInECEF(3), trueUtmZone); % ��ʵλ�� ECEF -> UTM ENU
DPEresultInUTM = zeros(size(data,1), 3); 
for ii = 1:size(data, 1)
    tempX = data(ii, 4); % ������ݵ�4��ΪECEF-X
    tempY = data(ii, 5); % ������ݵ�5��ΪECEF-Y
    tempZ = data(ii, 6); % ������ݵ�6��ΪECEF-Z
    tempLat = data(ii, 10); % ������ݵ�10��Ϊγ��
    tempLon = data(ii, 11); % ������ݵ�11��Ϊ����
    tempUtmZone = findUtmZone(tempLat, tempLon);
    DPEresultInUTM(ii,:) = cart2utm(tempX, tempY, tempZ, tempUtmZone); % DPE����λ�� ECEF -> UTM ENU
    % �����ENU��DPE��������UTM�����еĽ�����������Ϊ��������
end

DPEresultInENU = zeros(size(data,1), 3); % DPE����λ��(ENU)
for ii = 1:size(data, 1)
    if ii == 1
        refState = data(1, 4:6)';
    else
        refState = data(ii-1, 4:6)';
    end
    curState = data(ii, 4:6)';
    DPEresultInENU(ii,:) = ECEF2ENU(refState, curState)'; % DPE����λ�� ECEF -> ENU 
    % ����DPEÿ�ν��������ǰһ�ν�����Ϊ���Ľ�������
    % ��������ENU�Ǳ��ε�DPE�����������ǰһ�ε�DPE��������ENU����
end
DPEresultInECEF = data(:, [4 5 6]); % DPE����λ�� ECEF 
timeAxis = double(DPE_start_time) + (1:size(data,1))*DPE_interval; % DPE������ʱ��������

%%  �������ͼ
figure(100);
ax100 = subplot(1,1,1);
h1 = plot(ax100, timeAxis, fliplr(DPEresultInUTM-truePosInUTM), 'LineWidth', 1.5); % DPE�������������ʵλ�õ�ENU���
nameArr = {'color'};
valueArr = {'#5d70ea'; '#e34f26'; '#fbbc05'};
set(h1, nameArr, valueArr);
title('�������(ENU)');
xlabel('ʱ��/s');
ylabel('���(m)');
grid on;
axis('tight'); 
legend('U','N','E');
a = axis;
axis([a(1) a(2) a(3)-1 a(4)+1])

figure(101);
ax101 = subplot(1,1,1);
h2 = plot(ax101, timeAxis, fliplr(DPEresultInECEF-truePosInECEF'), 'LineWidth', 1.5); % DPE�������������ʵλ�õ�ECEF���
nameArr = {'color'};
valueArr = {'#5d70ea'; '#e34f26'; '#fbbc05'};
set(h2, nameArr, valueArr);
title('�������(ECEF)');
xlabel('ʱ��/s');
ylabel('���(m)');
grid on;
axis('tight'); 
legend('Z','Y','X');
a = axis;
axis([a(1) a(2) a(3)-1 a(4)+1])

%%  DPE�����������ͼ
figure(200);
ax200 = subplot(1,1,1);
color = zeros(size(DPEresultInUTM));
color(:, 2) = linspace(0.3, 1, size(DPEresultInUTM,1))';

% ��ʼ����Ϊԭ��
% % ������ʼ��
% scatter3(0, 0, 0, 50, 'r', 'filled'); % ��ʼλ�ã�����ʼλ����Ϊԭ�㣩
% text(0, 0, 0, '��ʼλ��');
% hold('on');
% 
% % �����յ�
% endE = DPEresultInUTM(end,1) - DPEresultInUTM(1,1);
% endN = DPEresultInUTM(end,2) - DPEresultInUTM(1,2);
% endU = DPEresultInUTM(end,3) - DPEresultInUTM(1,3); % ����λ���������ʼλ�õ�ENU�仯
% scatter3 (endE, endN, endU, 50, 'b', 'filled');
% text(endE, endN, endU, '����λ��');
% 
% % �������н����
% scatter3 (DPEresultInUTM(:,1) - DPEresultInUTM(1,1), ...
%           DPEresultInUTM(:,2) - DPEresultInUTM(1,2), ... 
%           DPEresultInUTM(:,3) - DPEresultInUTM(1,3), 20, color, '+'); % ������λ���������ʼλ�õ�ENU�仯
% 
% % ������ʵλ��
% trueE = truePosInUTM(1) - DPEresultInUTM(1,1);
% trueN = truePosInUTM(2) - DPEresultInUTM(1,2);
% trueU = truePosInUTM(3) - DPEresultInUTM(1,3); % ��ʵλ���������ʼλ�õ�ENU
% scatter3 (trueE, trueN, trueU, 'r+');
% text(trueE, trueN, trueU, '��ʵλ��');

% ��ʵλ����Ϊԭ��
% �������н����
scatter3 (ax200, DPEresultInUTM(2:end-1,1) - truePosInUTM(1), ...
                 DPEresultInUTM(2:end-1,2) - truePosInUTM(2), ... 
                 DPEresultInUTM(2:end-1,3) - truePosInUTM(3), 30, color(2:end-1,:), '+'); % ������λ���������ʵλ�õ�ENU�仯
hold('on');

% ������ʼ��
beginE = DPEresultInUTM(1,1) - truePosInUTM(1);
beginN = DPEresultInUTM(1,2) - truePosInUTM(2);
beginU = DPEresultInUTM(1,3) - truePosInUTM(3); % ��ʼλ��
scatter3 (ax200, beginE, beginN, beginU, 50, 'r', 'filled');
text(beginE, beginN, beginU, '��ʼλ��');

% �����յ�
endE = DPEresultInUTM(end,1) - truePosInUTM(1);
endN = DPEresultInUTM(end,2) - truePosInUTM(2);
endU = DPEresultInUTM(end,3) - truePosInUTM(3); % ����λ��
scatter3 (ax200, endE, endN, endU, 50, 'b', 'filled');
text(endE, endN, endU, '����λ��');

% ������ʵλ��
scatter3(ax200, 0, 0, 0, 'r+'); % ��ʵλ��??
text(0, 0, 0, '��ʵλ��');

xlabel('X');ylabel('Y');zlabel('Z');
title('DPE������ͼ')
hold('off');
view(0, 90); % ֻ������ͼ
axis('equal');
grid('minor');
a = axis;
axis([a(1)-5 a(2)+5 a(3)-5 a(4)+5])

%% ����ˮƽλ����ĸ�������
% -------------------------------------------------------------------------
% corr_pos�����ֲ����ɽ��ܣ�
% ����X/Y/Z/dT���������㣬�ֱ�Ϊ0��1������2^4=16��ֵ���ֲ�����
% X   |       0               1       |       0               1       |
% Y   |   0       1   |   0       1   |   0       1   |   0       1   |
% Z   | 0   1 | 0   1 | 0   1 | 0   1 | 0   1 | 0   1 | 0   1 | 0   1 | 
% dT  |0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|    
% Ŀǰ����X/Y��������ͼ�������Ҫ��Z��dTȡ����ȡֵ�е����ֵ
% ��corr_pos��Ӧ����25������25^4=390625�������
% -------------------------------------------------------------------------
corr_shape = size(corr_pos); % 500*390625
numOfEpoch = corr_shape(1);

% ����
dX = [-22, -19, -16, -13, -10, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 10, 13, 16, 19, 22]*3;
dY = [-22, -19, -16, -13, -10, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 10, 13, 16, 19, 22]*3;
org_indx = find(dX == 0);
org_indy = find(dY == 0);
[DPE_x, DPE_y] = meshgrid(dX, dY); % ��ˮƽ������
DPE_y = flipud(DPE_y);

f300 = figure(300); 
ax300 = subplot(1,1,1); 
f301 = figure(301);
ax301 = subplot(1,1,1);

for epoch = 1 : numOfEpoch    
% for epoch = 20
    cla(ax300, ax301);
    
    curr_corr_pos = corr_pos(epoch, :); % 1*390625
    corr_xy_zdt = reshape(curr_corr_pos, [625, 625]); % 625*625
    corr_xy = max(corr_xy_zdt); % 1*625
    curr_corr_xy_2d = flipud(reshape(corr_xy, [25, 25])); % 25*25
    [max_corr, max_idx] = max(curr_corr_xy_2d,[],'all','linear');
    curr_corr_xy_2d = curr_corr_xy_2d ./ max_corr; % ��һ��
    
    % ���Ƶ�ǰ��Ԫ�ĸ�������ͼ
    meshc(ax300, DPE_x, DPE_y, curr_corr_xy_2d);hold(ax300, 'on');
    
    scatter3(ax300, DPE_x(max_idx), DPE_y(max_idx), curr_corr_xy_2d(max_idx),  50, 'r', 'filled'); 
    text(ax300, DPE_x(max_idx), DPE_y(max_idx), curr_corr_xy_2d(max_idx), '�������ֵ');
    
    scatter3(ax300, dX(org_indx), dY(org_indy), curr_corr_xy_2d(org_indx, org_indy), 50, 'b', 'filled');
    text(ax300, dX(org_indx), dY(org_indy),curr_corr_xy_2d(org_indx, org_indy), '����ԭ��');
    
    xlabel(ax300, 'X'); ylabel(ax300, 'Y'); zlabel(ax300, 'Corr Score');
    grid(ax300, 'minor'); axis(ax300, 'tight');
    current_time = double(DPE_start_time) + epoch*DPE_corr_save_interval; % DPE������ʱ��������
    title(ax300, ['��' num2str(current_time) 's������ͼ']);
    view(ax300, -37.5, 30); 
    hold(ax300, 'off');
    
    % ����λ����֤ͼ
    
    scatter(ax301, DPE_x(:), DPE_y(:), 10, 'p', 'filled'); hold(ax301, 'on');
    
    scatter(ax301, DPEresultInENU(epoch,1), ...
            DPEresultInENU(epoch,2), 50, 'r', 'filled');
    text(ax301, DPEresultInENU(epoch,1), DPEresultInENU(epoch,2), 'DPE������');

    scatter(ax301, 0, 0, 50, 'b', 'filled');
    text(ax301, 0, 0, '����ԭ��');

    xlabel(ax301, 'X');ylabel(ax301, 'Y');
    title(ax301, ['��' num2str(current_time) 's������ͼ']);
    grid(ax301, 'minor'); axis(ax301, 'tight'); hold(ax301, 'off');
    
    pause(0.5);
end

%% save fig
figDir = ['saveFig/' dataDir secDir];
mkdir(figDir);
saveas(figure(100), [figDir 'ENU']);
saveas(figure(101), [figDir 'ECEF']);
saveas(figure(200), [figDir 'DPE_result']);
saveas(figure(300), [figDir 'manifold']);
saveas(figure(301), [figDir 'curr_epoch_grid']);




