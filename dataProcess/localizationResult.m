function [] = localizationResult(dataDir, secDir, truePosInLLA)

%% 读取DPE定位结果文件
predir   = ['./pre-simulator/' dataDir];
postdir  = ['./post-simulator/' dataDir];

% 导入位置结果表格
csvFile = 'usrp.csv';
csvDir  = [postdir secDir csvFile];
csvData = importdata(csvDir);

% 导入网格相关结果
receiveFile = 'end-of-dp/costFunction.h5';
scoreDir = [postdir secDir receiveFile];
dpe_plan = h5read(scoreDir,'/dpe_plan');
dpe_plan = dpe_plan{1};
if strcmp(dpe_plan, 'GRID')
    corr_pos = h5read(scoreDir,'/GRID/corr_pos')';
elseif strcmp(dpe_plan, 'ARS')
    cost_score = h5read(scoreDir,'/ARS/costScore')';
    N_Iter = h5read(scoreDir,'/ARS/N_Iter')';
    deltaxt = permute(h5read(scoreDir,'/ARS/deltaxt'), [3 2 1]);
    delta_in_iteration = permute(h5read(scoreDir,'/ARS/delta_in_iteration'), [3 2 1]);
end

% 导入运行时间相关参数
dpeRunTimeFile = 'dpe_runtime.mat';
load([postdir secDir dpeRunTimeFile]);

%% 初始化数据
data = csvData.data; % 定位结果数据
truePosInECEF = wgslla2xyz(truePosInLLA(1), truePosInLLA(2), truePosInLLA(3)); % 真实位置 LLA -> ECEF
trueUtmZone = findUtmZone(truePosInLLA(1), truePosInLLA(2));
truePosInUTM = cart2utm(truePosInECEF(1), truePosInECEF(2), truePosInECEF(3), trueUtmZone); % 真实位置 ECEF -> UTM ENU
DPEresultInUTM = zeros(size(data,1), 3); 
for ii = 1:size(data, 1)
    tempX = data(ii, 4); % 表格数据第4列为ECEF-X
    tempY = data(ii, 5); % 表格数据第5列为ECEF-Y
    tempZ = data(ii, 6); % 表格数据第6列为ECEF-Z
    tempLat = data(ii, 10); % 表格数据第10列为纬度
    tempLon = data(ii, 11); % 表格数据第11列为经度
    tempUtmZone = findUtmZone(tempLat, tempLon);
    DPEresultInUTM(ii,:) = cart2utm(tempX, tempY, tempZ, tempUtmZone); % DPE解算位置 ECEF -> UTM ENU
    % 这里的ENU是DPE解算结果在UTM网格中的结果，可以理解为绝对坐标
end

DPEresultInENU = zeros(size(data,1), 3); % DPE解算位置(ENU)
for ii = 1:size(data, 1)
    if ii == 1
        refState = data(1, 4:6)';
    else
        refState = data(ii-1, 4:6)';
    end
    curState = data(ii, 4:6)';
    DPEresultInENU(ii,:) = ECEF2ENU(refState, curState)'; % DPE解算位置 ECEF -> ENU 
    % 由于DPE每次解算均是以前一次解算结果为中心建立网格
    % 因此这里的ENU是本次的DPE解算结果相对于前一次的DPE解算结果的ENU坐标
end
DPEresultInECEF = data(:, [4 5 6]); % DPE解算位置 ECEF 
timeAxis = double(DPE_start_time) + (1:size(data,1))*double(DPE_interval); % DPE解算结果时间坐标轴

%%  坐标误差图
figure(100);
ax100 = subplot(1,1,1);
h1 = plot(ax100, timeAxis, fliplr(DPEresultInUTM-truePosInUTM), 'LineWidth', 1.5); % DPE解算结果相对于真实位置的ENU误差
nameArr = {'color'};
valueArr = {'#5d70ea'; '#e34f26'; '#fbbc05'};
set(h1, nameArr, valueArr);
title('坐标误差(ENU)');
xlabel('时间/s');
ylabel('误差(m)');
grid on;
axis('tight'); 
legend('U','N','E');
a = axis;
axis([a(1) a(2) a(3)-1 a(4)+1])

figure(101);
ax101 = subplot(1,1,1);
h2 = plot(ax101, timeAxis, fliplr(DPEresultInECEF-truePosInECEF'), 'LineWidth', 1.5); % DPE解算结果相对于真实位置的ECEF误差
nameArr = {'color'};
valueArr = {'#5d70ea'; '#e34f26'; '#fbbc05'};
set(h2, nameArr, valueArr);
title('坐标误差(ECEF)');
xlabel('时间/s');
ylabel('误差(m)');
grid on;
axis('tight'); 
legend('Z','Y','X');
a = axis;
axis([a(1) a(2) a(3)-1 a(4)+1])

%%  DPE解算结果网格点图
figure(200);
ax200 = subplot(1,1,1);
color = zeros(size(DPEresultInUTM));
color(:, 2) = linspace(0.3, 1, size(DPEresultInUTM,1))';

% 起始点作为原点
% % 绘制起始点
% scatter3(0, 0, 0, 50, 'r', 'filled'); % 起始位置（以起始位置作为原点）
% text(0, 0, 0, '起始位置');
% hold('on');
% 
% % 绘制终点
% endE = DPEresultInUTM(end,1) - DPEresultInUTM(1,1);
% endN = DPEresultInUTM(end,2) - DPEresultInUTM(1,2);
% endU = DPEresultInUTM(end,3) - DPEresultInUTM(1,3); % 结束位置相对于起始位置的ENU变化
% scatter3 (endE, endN, endU, 50, 'b', 'filled');
% text(endE, endN, endU, '结束位置');
% 
% % 绘制所有解算点
% scatter3 (DPEresultInUTM(:,1) - DPEresultInUTM(1,1), ...
%           DPEresultInUTM(:,2) - DPEresultInUTM(1,2), ... 
%           DPEresultInUTM(:,3) - DPEresultInUTM(1,3), 20, color, '+'); % 各解算位置相对于起始位置的ENU变化
% 
% % 绘制真实位置
% trueE = truePosInUTM(1) - DPEresultInUTM(1,1);
% trueN = truePosInUTM(2) - DPEresultInUTM(1,2);
% trueU = truePosInUTM(3) - DPEresultInUTM(1,3); % 真实位置相对于起始位置的ENU
% scatter3 (trueE, trueN, trueU, 'r+');
% text(trueE, trueN, trueU, '真实位置');

% 真实位置作为原点

% 绘制真实位置
scatter3(ax200, 0, 0, 0, 'r+'); % 真实位置
text(0, 0, 0, '真实位置');
hold('on');

% 绘制起始点
beginE = DPEresultInUTM(1,1) - truePosInUTM(1);
beginN = DPEresultInUTM(1,2) - truePosInUTM(2);
beginU = DPEresultInUTM(1,3) - truePosInUTM(3); % 起始位置
scatter3 (ax200, beginE, beginN, beginU, 50, 'r', 'filled');
text(beginE, beginN, beginU, '起始位置');

% 绘制所有解算点
scatter3 (ax200, DPEresultInUTM(2:end-1,1) - truePosInUTM(1), ...
                 DPEresultInUTM(2:end-1,2) - truePosInUTM(2), ... 
                 DPEresultInUTM(2:end-1,3) - truePosInUTM(3), 30, color(2:end-1,:), '+'); % 各解算位置相对于真实位置的ENU变化


% 绘制终点
endE = DPEresultInUTM(end,1) - truePosInUTM(1);
endN = DPEresultInUTM(end,2) - truePosInUTM(2);
endU = DPEresultInUTM(end,3) - truePosInUTM(3); % 结束位置
scatter3 (ax200, endE, endN, endU, 50, 'b', 'filled');
text(endE, endN, endU, '结束位置');

xlabel('X');ylabel('Y');zlabel('Z');
title('DPE解算结果图')
hold('off');
view(0, 90); % 只看俯视图
axis('equal');
grid('minor');
a = axis;
axis([a(1)-5 a(2)+5 a(3)-5 a(4)+5])

%% 
if strcmp(dpe_plan, 'GRID')
    %% 生成水平位置域的概率流形
    % -------------------------------------------------------------------------
    % corr_pos网格点分布规律介绍：
    % 假设X/Y/Z/dT各向两个点，分别为0和1，则共有2^4=16个值，分布如下
    % X   |       0               1       |       0               1       |
    % Y   |   0       1   |   0       1   |   0       1   |   0       1   |
    % Z   | 0   1 | 0   1 | 0   1 | 0   1 | 0   1 | 0   1 | 0   1 | 0   1 | 
    % dT  |0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|0 1|    
    % 目前仅就X/Y两个方向画图，因此需要对Z和dT取所有取值中的最大值
    % 而corr_pos对应各向25个网格共25^4=390625个网格点
    % -------------------------------------------------------------------------
    corr_shape = size(corr_pos); % numOfEpoch*numOfGrid^4
    numOfEpoch = corr_shape(1);

    % 网格
%     a = -100 : 10 : 101;
%     pos_b = 125 : 25 : 501;
%     neg_b = -pos_b(end:-1:1);
%     pos_c = 550 : 50 : 1001;
%     neg_c = -pos_c(end:-1:1);
%     dtmp = [neg_c, neg_b, a, pos_b, pos_c];
    dtmp = [-22, -19, -16, -13, -10, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 10, 13, 16, 19, 22]*3;

    numOfGrid = length(dtmp);
    dX = dtmp;
    dY = dtmp;
    org_indx = find(dX == 0);
    org_indy = find(dY == 0);
    [DPE_x, DPE_y] = meshgrid(dX, dY); % 求水平坐标面
    DPE_y = flipud(DPE_y);

    f300 = figure(300); 
    ax300 = subplot(1,1,1); 
    f301 = figure(301);
    ax301 = subplot(1,1,1);

    for epoch = 1 : numOfEpoch    
%     for epoch = 2
        cla(ax300, ax301);

        curr_corr_pos = corr_pos(epoch, :); % 1*390625
        corr_xy_zdt = reshape(curr_corr_pos, [numOfGrid^2, numOfGrid^2]); % 625*625
        corr_xy = max(corr_xy_zdt); % 1*625
        curr_corr_xy_2d = flipud(reshape(corr_xy, [numOfGrid, numOfGrid])); % 25*25
        [max_corr, max_idx] = max(curr_corr_xy_2d,[],'all','linear'); 
        curr_corr_xy_2d = curr_corr_xy_2d ./ max_corr; % 归一化

        % 绘制当前历元的概率流型图
        meshc(ax300, DPE_x, DPE_y, curr_corr_xy_2d);hold(ax300, 'on');

        scatter3(ax300, DPE_x(max_idx), DPE_y(max_idx), curr_corr_xy_2d(max_idx),  50, 'r', 'filled'); 
        text(ax300, DPE_x(max_idx), DPE_y(max_idx), curr_corr_xy_2d(max_idx), '网格最大值');

        scatter3(ax300, dX(org_indx), dY(org_indy), curr_corr_xy_2d(org_indx, org_indy), 50, 'b', 'filled');
        text(ax300, dX(org_indx), dY(org_indy),curr_corr_xy_2d(org_indx, org_indy), '网格原点');

        xlabel(ax300, 'X'); ylabel(ax300, 'Y'); zlabel(ax300, 'Corr Score');
        grid(ax300, 'minor'); axis(ax300, 'tight');
        current_time = double(DPE_start_time) + epoch*DPE_corr_save_interval; % DPE解算结果时间坐标轴
        title(ax300, ['第' num2str(current_time) 's的流型图']);
        view(ax300, -37.5+90, 30-20); 
%         view(ax300, -37.5, 30); 
        hold(ax300, 'off');

        % 网格位置验证图

        scatter(ax301, DPE_x(:), DPE_y(:), 10, 'p', 'filled'); hold(ax301, 'on');

        scatter(ax301, DPEresultInENU(epoch,1), ...
                DPEresultInENU(epoch,2), 50, 'r', 'filled');
        text(ax301, DPEresultInENU(epoch,1), DPEresultInENU(epoch,2), 'DPE解算结果');

        scatter(ax301, 0, 0, 50, 'b', 'filled');
        text(ax301, 0, 0, '网格原点');

        xlabel(ax301, 'X');ylabel(ax301, 'Y');
        title(ax301, ['第' num2str(current_time) 's的网格图']);
        grid(ax301, 'minor'); axis(ax301, 'tight'); hold(ax301, 'off');

        pause(0.5);
    end
elseif strcmp(dpe_plan, 'ARS')
    fig400 = figure(400);
    ax400 = subplot(1,1,1);
%     fig401 = figure(401);
%     ax401 = subplot(1,1,1);
    corr_shape = size(cost_score);
    numOfEpoch = corr_shape(1);
%     for epoch = 1 : numOfEpoch    
    for epoch = 1
        cla(ax400);
        yyaxis(ax400, 'left');
        plot(ax400, 1:N_Iter, squeeze(delta_in_iteration(epoch, 1, :)), 'LineWidth', 1.0);
        ylabel(ax400, '$\delta x$ in ECEF(m)', 'Interpreter', 'latex');
        yyaxis(ax400, 'right');
        plot(ax400, 1:N_Iter, cost_score(epoch, :), 'LineWidth', 1.0);
        xlabel(ax400, 'ARS Iteration Step'); ylabel(ax400, 'cost function value');
        current_time = double(DPE_start_time) + epoch*DPE_corr_save_interval;
        title(ax400, ['第' num2str(current_time) 's的代价函数图']);
        grid(ax400, 'minor'); axis(ax400, 'tight'); hold(ax400, 'off');
        
        pause(1);
    end
end

end