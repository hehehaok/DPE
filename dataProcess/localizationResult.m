clear;
close all;
fclose('all');
clc;

%% 数据选择
dataFlag = 0; % 0-芬兰数据 1-OAKBAT数据 2-芬兰数据第二次运行 3-OAKBAT数据第二次运行
switch dataFlag
    case 0
        dataDir = 'finland_cleanStatic/';
        truePosInLLA = [60.161086788889, 24.545448080556, 54.1640000026673]; % 真实位置(纬度(北正南负)经度(东正西负) 高度) 芬兰
        secDir = 'test0907_1_skip0s_proc40s/';
    case 1
        dataDir = 'oak_cleanStatic/';
        truePosInLLA = [35.930544444, -84.310652778, 248.6000]; % 真实位置(纬度(北正南负) 经度(东正西负) 高度) oak静态
        secDir = 'test0907_1_skip50s_proc40s/';
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
% 导入位置结果表格
csvFile = 'usrp.csv';
csvDir = [postdir secDir csvFile];
csvData = importdata(csvDir);
% 导入相关结果
receiveFile = 'end-of-dp/receiver.mat';
rxDir = [postdir secDir receiveFile];
load(rxDir);

%% 初始化数据
const.EARTH_SEMIMAJORAXIS = 6378137;
const.EARTH_FLATTENING = 1/298.257223563;
header = csvData.colheaders; % 表头
data = csvData.data; % 定位结果数据
truePosInXYZ = wgslla2xyz(const, truePosInLLA(1), truePosInLLA(2), truePosInLLA(3)); % 真实位置(WGS84 XYZ)
trueUtmZone = findUtmZone(truePosInLLA(1), truePosInLLA(2));
truePosInENU = cart2utm(truePosInXYZ(1), truePosInXYZ(2), truePosInXYZ(3), trueUtmZone); % 真实位置(UTM ENU)
DPEresultInENU = zeros(size(data,1), 3); % DPE解算位置(UTM ENU)
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

%%  坐标误差图
figure(100);
plot(DPEresultInENU-truePosInENU, 'LineWidth', 1.2);
title('坐标误差(ENU)');
xlabel(['测算周期: ', num2str(round((data(2,3)-data(1,3))*1000)), 'ms']);
ylabel('误差 (m)');
grid on;
axis('tight'); 
legend('E','N','U');
a = axis;
axis([a(1) a(2) a(3)-1 a(4)+1])

%%  网格点图
figure(200);
color = zeros(size(DPEresultInENU));
color(:, 2) = linspace(0.3, 1, size(DPEresultInENU,1))';

% 绘制所有解算点
% scatter3 (DPEresultInENU(:,1) - truePosInENU(1), ...
%           DPEresultInENU(:,2) - truePosInENU(2), ... 
%           DPEresultInENU(:,3) - truePosInENU(3), 10, color, '+');
% hold('on');
scatter3 (DPEresultInENU(2,1) - truePosInENU(1), ...
          DPEresultInENU(2,2) - truePosInENU(2), ... 
          DPEresultInENU(2,3) - truePosInENU(3), 10, 'b+');
hold('on');

% 绘制起始点
beginE = DPEresultInENU(1,1) - truePosInENU(1);
beginN = DPEresultInENU(1,2) - truePosInENU(2);
beginU = DPEresultInENU(1,3) - truePosInENU(3); % 起始位置
scatter3 (beginE, beginN, beginU, 50, 'r', 'filled');
text(beginE, beginN, beginU, '起始位置');

% 绘制终点
endE = DPEresultInENU(end,1) - truePosInENU(1);
endN = DPEresultInENU(end,2) - truePosInENU(2);
endU = DPEresultInENU(end,3) - truePosInENU(3); % 结束位置
scatter3 (endE, endN, endU, 50, 'b', 'filled');
text(endE, endN, endU, '结束位置');

% 绘制真实位置
scatter3(0, 0, 0, 'r+'); % 真实位置
text(0, 0, 0, '真实位置');

xlabel('X');ylabel('Y');zlabel('Z');
hold('off');
view(0, 90); % 只看俯视图
axis('equal');
grid('minor');

%% 生成水平位置域的概率流形
% -------------------------------------------------------------------------
% 主要思路：
% 生成水平位置域的网格对应的序号
% 对于位置域，网格点各向25个，共25^4=390625个网格点
% X    25^3变一次值
% Y    25^2变一次值
% Z    25^1变一次值
% dT   25^0变一次值
% 目前仅考虑X/Y，Z和dT取所有取值中的最大值
% -------------------------------------------------------------------------
corr_shape = size(corr_pos);
numOfEpoch = corr_shape(1);
axis_del = [4, 5]; % 仅考虑X/Y，因此后续取最大值的时候在Z和dT对应轴上取
corr_xyzdt = reshape(corr_pos, [numOfEpoch, 25, 25, 25, 25]);
corr_xy_2d = max(corr_xyzdt, [], axis_del);

% 网格
dX = [-22, -19, -16, -13, -10, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 10, 13, 16, 19, 22]*3;
dY = [-22, -19, -16, -13, -10, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 10, 13, 16, 19, 22]*3;
org_indx = find(dX == 0);
org_indy = find(dY == 0);

figure(300);
% for epoch = 1 : numOfEpoch    
for epoch = 1
    curr_corr_xy_2d = corr_xy_2d(epoch, : , :);
    curr_corr_xy_2d = squeeze(curr_corr_xy_2d); % 矩阵维度(1,25,25) -> (25,25)
    [max_corr, max_idx] = max(curr_corr_xy_2d,[],'all','linear');
    curr_corr_xy_2d = curr_corr_xy_2d ./ max_corr; % 归一化
    
    % 求水平坐标面
%     axis_x = DPEresultInENU(epoch,1) + dX;
%     axis_y = DPEresultInENU(epoch,2) + dY;
    axis_x = dX;
    axis_y = dY;
    [DPE_x, DPE_y] = meshgrid(axis_x, axis_y);
    
    % 绘制当前历元的概率流型图
    surf(DPE_x, DPE_y, curr_corr_xy_2d); hold on;
    scatter3(DPE_x(max_idx), DPE_y(max_idx), curr_corr_xy_2d(max_idx),  50, 'r', 'filled');
    scatter3(axis_x(org_indx), axis_y(org_indy), curr_corr_xy_2d(org_indx, org_indy), 50, 'b', 'filled');
    xlabel('X');ylabel('Y');zlabel('Corr Score');
    grid('minor');hold off; axis('tight');
end
