%% 
clear;
close all;
fclose('all');
clc;

%% 数据选择
dataFlag = 3; % 0-芬兰数据 1-OAKBAT数据 2-芬兰数据第二次运行 3-OAKBAT数据第二次运行
switch dataFlag
    case 0
        dataDir = 'finland_cleanStatic/';
        truePosInLLA = [60.161086788889, 24.545448080556, 54.1640000026673]; % 真实位置(纬度(北正南负)经度(东正西负) 高度) 芬兰
        secDir = 'test1__1USRPproc40s_6605_7200/';
    case 1
        dataDir = 'oak_cleanStatic/';
        truePosInLLA = [35.930544444, -84.310652778, 248.6000]; % 真实位置(纬度(北正南负) 经度(东正西负) 高度) oak静态
        secDir = 'test1__1USRPproc60s_6605_7200/';
    case 2
        dataDir = 'finland_cleanStatic/';
        truePosInLLA = [60.161086788889, 24.545448080556, 54.1640000026673]; % 真实位置(纬度(北正南负)经度(东正西负) 高度) 芬兰
        secDir = 'test1proc40s/';     
    case 3
        dataDir = 'oak_cleanStatic/';
        truePosInLLA = [35.930544444, -84.310652778, 248.6000]; % 真实位置(纬度(北正南负) 经度(东正西负) 高度) oak静态
        secDir = 'test3proc60s/';
    otherwise
end

%% 读取DPE定位结果文件
predir   = ['./pre-simulator/' dataDir];
postdir  = ['./post-simulator/' dataDir];
csvFile = 'usrp6.csv';
dir = [postdir secDir csvFile];
csvData = importdata(dir);

%% 初始化数据
const.EARTH_SEMIMAJORAXIS = 6378137;
const.EARTH_FLATTENING = 1/298.257223563;
header = csvData.colheaders; % 表头
data = csvData.data; % 定位结果数据
truePosInXYZ = wgslla2xyz(const, truePosInLLA(1), truePosInLLA(2), truePosInLLA(3)); % 真实位置(WGS84 XYZ)
trueUtmZone = findUtmZone(truePosInLLA(1), truePosInLLA(2));
truePosInENU = cart2utm(truePosInXYZ(1), truePosInXYZ(2), truePosInXYZ(3), trueUtmZone); % 真实位置(UTM ENU)
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

%% 画图
% 坐标误差图
figure(100);
plot(DPEresultInENU(:,1)-truePosInENU(1));
title('东向坐标误差(UTM坐标系)');
xlabel(['测算周期: ', num2str(round((data(2,3)-data(1,3))*1000)), 'ms']);
ylabel('误差 (m)');
grid on;
axis('tight'); 
a = axis;
axis([a(1) a(2) a(3)-1 a(4)+1])

figure(101);
plot(DPEresultInENU(:,2)-truePosInENU(2));
title('北向坐标误差(UTM坐标系)');
xlabel(['测算周期: ', num2str(round((data(2,3)-data(1,3))*1000)), 'ms']);
ylabel('误差 (m)');
grid on;
axis('tight'); 
a = axis;
axis([a(1) a(2) a(3)-1 a(4)+1])

figure(102);
plot(DPEresultInENU(:,3)-truePosInENU(3));
title('天向坐标误差(UTM坐标系)');
xlabel(['测算周期: ', num2str(round((data(2,3)-data(1,3))*1000)), 'ms']);
ylabel('误差 (m)');
grid on;
axis('tight'); 
a = axis;
axis([a(1) a(2) a(3)-1 a(4)+1])

% 网格点图
figure(200);
color = zeros(size(DPEresultInENU));
color(:, 2) = linspace(0.3, 1, size(DPEresultInENU,1))';
scatter3 (DPEresultInENU(:,1) - truePosInENU(1), ...
       DPEresultInENU(:,2) - truePosInENU(2), ... 
       DPEresultInENU(:,3) - truePosInENU(3), 10, color, '+');
hold('on');
beginE = DPEresultInENU(1,1) - truePosInENU(1);
beginN = DPEresultInENU(1,2) - truePosInENU(2);
beginU = DPEresultInENU(1,3) - truePosInENU(3); % 起始位置
scatter3 (beginE, beginN, beginU, 50, 'r', 'filled');
text(beginE, beginN, beginU, '起始位置');
endE = DPEresultInENU(end,1) - truePosInENU(1);
endN = DPEresultInENU(end,2) - truePosInENU(2);
endU = DPEresultInENU(end,3) - truePosInENU(3); % 结束位置
scatter3 (endE, endN, endU, 50, 'b', 'filled');
text(endE, endN, endU, '结束位置');
scatter3(0, 0, 0, 'r+'); % 真实位置
text(0, 0, 0, '真实位置');
hold('off');
view(0, 90);
axis('equal');
grid('minor');






