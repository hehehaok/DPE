%% 
clear;
close all;
fclose('all');
clc;

%% ����ѡ��
dataFlag = 3; % 0-�������� 1-OAKBAT���� 2-�������ݵڶ������� 3-OAKBAT���ݵڶ�������
switch dataFlag
    case 0
        dataDir = 'finland_cleanStatic/';
        truePosInLLA = [60.161086788889, 24.545448080556, 54.1640000026673]; % ��ʵλ��(γ��(�����ϸ�)����(��������) �߶�) ����
        secDir = 'test1__1USRPproc40s_6605_7200/';
    case 1
        dataDir = 'oak_cleanStatic/';
        truePosInLLA = [35.930544444, -84.310652778, 248.6000]; % ��ʵλ��(γ��(�����ϸ�) ����(��������) �߶�) oak��̬
        secDir = 'test1__1USRPproc60s_6605_7200/';
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
csvFile = 'usrp6.csv';
dir = [postdir secDir csvFile];
csvData = importdata(dir);

%% ��ʼ������
const.EARTH_SEMIMAJORAXIS = 6378137;
const.EARTH_FLATTENING = 1/298.257223563;
header = csvData.colheaders; % ��ͷ
data = csvData.data; % ��λ�������
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
% �������ͼ
figure(100);
plot(DPEresultInENU(:,1)-truePosInENU(1));
title('�����������(UTM����ϵ)');
xlabel(['��������: ', num2str(round((data(2,3)-data(1,3))*1000)), 'ms']);
ylabel('��� (m)');
grid on;
axis('tight'); 
a = axis;
axis([a(1) a(2) a(3)-1 a(4)+1])

figure(101);
plot(DPEresultInENU(:,2)-truePosInENU(2));
title('�����������(UTM����ϵ)');
xlabel(['��������: ', num2str(round((data(2,3)-data(1,3))*1000)), 'ms']);
ylabel('��� (m)');
grid on;
axis('tight'); 
a = axis;
axis([a(1) a(2) a(3)-1 a(4)+1])

figure(102);
plot(DPEresultInENU(:,3)-truePosInENU(3));
title('�����������(UTM����ϵ)');
xlabel(['��������: ', num2str(round((data(2,3)-data(1,3))*1000)), 'ms']);
ylabel('��� (m)');
grid on;
axis('tight'); 
a = axis;
axis([a(1) a(2) a(3)-1 a(4)+1])

% �����ͼ
figure(200);
color = zeros(size(DPEresultInENU));
color(:, 2) = linspace(0.3, 1, size(DPEresultInENU,1))';
scatter3 (DPEresultInENU(:,1) - truePosInENU(1), ...
       DPEresultInENU(:,2) - truePosInENU(2), ... 
       DPEresultInENU(:,3) - truePosInENU(3), 10, color, '+');
hold('on');
beginE = DPEresultInENU(1,1) - truePosInENU(1);
beginN = DPEresultInENU(1,2) - truePosInENU(2);
beginU = DPEresultInENU(1,3) - truePosInENU(3); % ��ʼλ��
scatter3 (beginE, beginN, beginU, 50, 'r', 'filled');
text(beginE, beginN, beginU, '��ʼλ��');
endE = DPEresultInENU(end,1) - truePosInENU(1);
endN = DPEresultInENU(end,2) - truePosInENU(2);
endU = DPEresultInENU(end,3) - truePosInENU(3); % ����λ��
scatter3 (endE, endN, endU, 50, 'b', 'filled');
text(endE, endN, endU, '����λ��');
scatter3(0, 0, 0, 'r+'); % ��ʵλ��
text(0, 0, 0, '��ʵλ��');
hold('off');
view(0, 90);
axis('equal');
grid('minor');






