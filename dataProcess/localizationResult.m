%% 
clear;
close all;
fclose('all');
clc;

%% ��ȡDPE��λ����ļ�
predir   = './pre-simulator/';
postdir  = './post-simulator/';
secDir = 'Talbot_rooftop_4hr_T_big=1_1USRPproc40s_6605_7200_1/';
csvFile = 'usrp6.csv';
runTime = 40;
dir = [postdir secDir csvFile];
csvData = importdata(csvFile);

%% ��ʼ������
const.EARTH_SEMIMAJORAXIS = 6378137;
const.EARTH_FLATTENING = 1/298.257223563;
header = csvData.colheaders; % ��ͷ
data = csvData.data; % ��λ�������
truePosInLLA = [60.161086788889, 24.545448080556, 54.1640000026673]; % ��ʵλ��(γ��(�����ϸ�) ����(��������) �߶�)
truePosInXYZ = wgslla2xyz(const, truePosInLLA(1), truePosInLLA(2), truePosInLLA(3)); % ��ʵλ��(WGS84 XYZ)
trueUtmZone = findUtmZone(truePosInLLA(1), truePosInLLA(2));
truePosInENU = cart2utm(truePosInXYZ(1), truePosInXYZ(2), truePosInXYZ(3), trueUtmZone); % ��ʵλ��(UTM ENU)
DPEresultInENU = zeros(size(data,1), 3);
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
figure(1);
plot(DPEresultInENU-truePosInENU);
title('�������(UTM����ϵ)');
legend('E', 'N', 'U');
xlabel(['��������: ', num2str(round((data(2,3)-data(1,3))*1000)), 'ms']);
ylabel('��� (m)');
grid on;
axis('tight'); 

figure(2);
% plot3 (DPEresultInENU(:,1) - truePosInENU(1), ...
%        DPEresultInENU(:,2) - truePosInENU(2), ... 
%        DPEresultInENU(:,3) - truePosInENU(3), '+');
% color = linspace(1,2,size(DPEresultInENU, 1));
% color = repmat([0 0 1], size(DPEresultInENU, 1), 1);
color = zeros(size(DPEresultInENU));
color(:, 2) = linspace(0.3, 1, size(DPEresultInENU,1))';
scatter3 (DPEresultInENU(:,1) - truePosInENU(1), ...
       DPEresultInENU(:,2) - truePosInENU(2), ... 
       DPEresultInENU(:,3) - truePosInENU(3), 10, color, '+');
hold('on');
% plot3(0, 0, 0, 'r+', 'LineWidth', 1.5, 'MarkerSize', 10);
scatter3 (DPEresultInENU(1,1) - truePosInENU(1), ... % ��ʼλ��
          DPEresultInENU(1,2) - truePosInENU(2), ... 
          DPEresultInENU(1,3) - truePosInENU(3), 50, 'r', 'filled');
scatter3 (DPEresultInENU(end,1) - truePosInENU(1), ... % ����λ��
          DPEresultInENU(end,2) - truePosInENU(2), ... 
          DPEresultInENU(end,3) - truePosInENU(3), 50, 'b', 'filled');
scatter3(0, 0, 0, 'r+'); % ��ʵλ��
hold('off');
view(0, 90);
axis('equal');
grid('minor');






