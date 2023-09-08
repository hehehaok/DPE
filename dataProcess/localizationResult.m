clear;
close all;
fclose('all');
clc;

%% ����ѡ��
dataFlag = 0; % 0-�������� 1-OAKBAT���� 2-�������ݵڶ������� 3-OAKBAT���ݵڶ�������
switch dataFlag
    case 0
        dataDir = 'finland_cleanStatic/';
        truePosInLLA = [60.161086788889, 24.545448080556, 54.1640000026673]; % ��ʵλ��(γ��(�����ϸ�)����(��������) �߶�) ����
        secDir = 'test0907_1_skip0s_proc40s/';
    case 1
        dataDir = 'oak_cleanStatic/';
        truePosInLLA = [35.930544444, -84.310652778, 248.6000]; % ��ʵλ��(γ��(�����ϸ�) ����(��������) �߶�) oak��̬
        secDir = 'test0907_1_skip50s_proc40s/';
    case 2
        dataDir = 'finland_cleanStatic/';
        truePosInLLA = [60.161086788889, 24.545448080556, 54.1640000026673]; % ��ʵλ��(γ��(�����ϸ�)����(��������) �߶�) ����
        secDir = 'test1proc40s/';     
    case 3
        dataDir = 'oak_cleanStatic/';
        truePosInLLA = [35.930544444, -84.310652778, 248.6000]; % ��ʵλ��(γ��(�����ϸ�) ����(��������) �߶�) oak��̬
        secDir = 'test3proc60s/';
    otherwise
end

%% ��ȡDPE��λ����ļ�
predir   = ['./pre-simulator/' dataDir];
postdir  = ['./post-simulator/' dataDir];
% ����λ�ý�����
csvFile = 'usrp.csv';
csvDir = [postdir secDir csvFile];
csvData = importdata(csvDir);
% ������ؽ��
receiveFile = 'end-of-dp/receiver.mat';
rxDir = [postdir secDir receiveFile];
load(rxDir);

%% ��ʼ������
const.EARTH_SEMIMAJORAXIS = 6378137;
const.EARTH_FLATTENING = 1/298.257223563;
header = csvData.colheaders; % ��ͷ
data = csvData.data; % ��λ�������
truePosInXYZ = wgslla2xyz(const, truePosInLLA(1), truePosInLLA(2), truePosInLLA(3)); % ��ʵλ��(WGS84 XYZ)
trueUtmZone = findUtmZone(truePosInLLA(1), truePosInLLA(2));
truePosInENU = cart2utm(truePosInXYZ(1), truePosInXYZ(2), truePosInXYZ(3), trueUtmZone); % ��ʵλ��(UTM ENU)
DPEresultInENU = zeros(size(data,1), 3); % DPE����λ��(UTM ENU)
for ii = 1:size(data, 1)
    tempX = data(ii, 4);
    tempY = data(ii, 5);
    tempZ = data(ii, 6);
    tempLat = data(ii, 10);
    tempLon = data(ii, 11);
    tempUtmZone = findUtmZone(tempLat, tempLon);
    DPEresultInENU(ii,:) = cart2utm(tempX, tempY, tempZ, tempUtmZone);
end

%% ��ͼ

%%  �������ͼ
figure(100);
plot(DPEresultInENU-truePosInENU, 'LineWidth', 1.2);
title('�������(ENU)');
xlabel(['��������: ', num2str(round((data(2,3)-data(1,3))*1000)), 'ms']);
ylabel('��� (m)');
grid on;
axis('tight'); 
legend('E','N','U');
a = axis;
axis([a(1) a(2) a(3)-1 a(4)+1])

%%  �����ͼ
figure(200);
color = zeros(size(DPEresultInENU));
color(:, 2) = linspace(0.3, 1, size(DPEresultInENU,1))';

% �������н����
% scatter3 (DPEresultInENU(:,1) - truePosInENU(1), ...
%           DPEresultInENU(:,2) - truePosInENU(2), ... 
%           DPEresultInENU(:,3) - truePosInENU(3), 10, color, '+');
% hold('on');
scatter3 (DPEresultInENU(2,1) - truePosInENU(1), ...
          DPEresultInENU(2,2) - truePosInENU(2), ... 
          DPEresultInENU(2,3) - truePosInENU(3), 10, 'b+');
hold('on');

% ������ʼ��
beginE = DPEresultInENU(1,1) - truePosInENU(1);
beginN = DPEresultInENU(1,2) - truePosInENU(2);
beginU = DPEresultInENU(1,3) - truePosInENU(3); % ��ʼλ��
scatter3 (beginE, beginN, beginU, 50, 'r', 'filled');
text(beginE, beginN, beginU, '��ʼλ��');

% �����յ�
endE = DPEresultInENU(end,1) - truePosInENU(1);
endN = DPEresultInENU(end,2) - truePosInENU(2);
endU = DPEresultInENU(end,3) - truePosInENU(3); % ����λ��
scatter3 (endE, endN, endU, 50, 'b', 'filled');
text(endE, endN, endU, '����λ��');

% ������ʵλ��
scatter3(0, 0, 0, 'r+'); % ��ʵλ��
text(0, 0, 0, '��ʵλ��');

xlabel('X');ylabel('Y');zlabel('Z');
hold('off');
view(0, 90); % ֻ������ͼ
axis('equal');
grid('minor');

%% ����ˮƽλ����ĸ�������
% -------------------------------------------------------------------------
% ��Ҫ˼·��
% ����ˮƽλ����������Ӧ�����
% ����λ������������25������25^4=390625�������
% X    25^3��һ��ֵ
% Y    25^2��һ��ֵ
% Z    25^1��һ��ֵ
% dT   25^0��һ��ֵ
% Ŀǰ������X/Y��Z��dTȡ����ȡֵ�е����ֵ
% -------------------------------------------------------------------------
corr_shape = size(corr_pos);
numOfEpoch = corr_shape(1);
axis_del = [4, 5]; % ������X/Y����˺���ȡ���ֵ��ʱ����Z��dT��Ӧ����ȡ
corr_xyzdt = reshape(corr_pos, [numOfEpoch, 25, 25, 25, 25]);
corr_xy_2d = max(corr_xyzdt, [], axis_del);

% ����
dX = [-22, -19, -16, -13, -10, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 10, 13, 16, 19, 22]*3;
dY = [-22, -19, -16, -13, -10, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 10, 13, 16, 19, 22]*3;
org_indx = find(dX == 0);
org_indy = find(dY == 0);

figure(300);
% for epoch = 1 : numOfEpoch    
for epoch = 1
    curr_corr_xy_2d = corr_xy_2d(epoch, : , :);
    curr_corr_xy_2d = squeeze(curr_corr_xy_2d); % ����ά��(1,25,25) -> (25,25)
    [max_corr, max_idx] = max(curr_corr_xy_2d,[],'all','linear');
    curr_corr_xy_2d = curr_corr_xy_2d ./ max_corr; % ��һ��
    
    % ��ˮƽ������
%     axis_x = DPEresultInENU(epoch,1) + dX;
%     axis_y = DPEresultInENU(epoch,2) + dY;
    axis_x = dX;
    axis_y = dY;
    [DPE_x, DPE_y] = meshgrid(axis_x, axis_y);
    
    % ���Ƶ�ǰ��Ԫ�ĸ�������ͼ
    surf(DPE_x, DPE_y, curr_corr_xy_2d); hold on;
    scatter3(DPE_x(max_idx), DPE_y(max_idx), curr_corr_xy_2d(max_idx),  50, 'r', 'filled');
    scatter3(axis_x(org_indx), axis_y(org_indy), curr_corr_xy_2d(org_indx, org_indy), 50, 'b', 'filled');
    xlabel('X');ylabel('Y');zlabel('Corr Score');
    grid('minor');hold off; axis('tight');
end
