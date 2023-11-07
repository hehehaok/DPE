function [] = plotManifold(dataDir, secDir, disp_axis)

%% ��ȡDPE��λ����ļ�
postdir  = ['./post-simulator/' dataDir];

% ����������ؽ��
receiveFile = 'end-of-dp/costFunction.h5';
scoreDir = [postdir secDir receiveFile];
corr_pos = h5read(scoreDir,'/GRID/corr_pos')';
axis_1d = h5read(scoreDir,'/GRID/axis_1d')';	
corr_pos_prn = permute(h5read(scoreDir,'/GRID/corr_pos_prn'), [3 2 1]);
dpe_prn = sort(h5read(scoreDir,'/GRID/dpe_prn')');

% ��������ʱ����ز���
dpeRunTimeFile = 'dpe_runtime.mat';
load([postdir secDir dpeRunTimeFile]);

% ����
prn_idx = 1;
curr_epoch = 1;

%% ����ˮƽλ����ĸ�������
f300 = figure(300); 
ax300 = subplot(1,1,1); 
f500 = figure(500); 
ax500 = subplot(1,1,1);

corr_shape = size(corr_pos); 
numOfEpoch = corr_shape(1);

% ����
dtmp = double(axis_1d);	
numOfGrid = length(dtmp);
dX = dtmp;
dY = dtmp;
[DPE_x, DPE_y] = meshgrid(dX, dY); % ��ˮƽ������
DPE_y = flipud(DPE_y);

% for epoch = 1 : numOfEpoch    
for epoch = curr_epoch
    cla(ax300);
    curr_corr_pos = corr_pos(epoch, :); 
    [curr_corr_2d, ~] = grid2dplane(curr_corr_pos, numOfGrid, 1, disp_axis);
    % ���Ƶ�ǰ��Ԫ�ĸ�������ͼ
    surf(ax300, DPE_x, DPE_y, curr_corr_2d);
    colormap(ax300, slanCM('prism2'));
    
    disp_axis = sort(disp_axis);
    xlabel(ax300, disp_axis(1)); ylabel(ax300, disp_axis(2)); zlabel(ax300, 'Corr Score');
    grid(ax300, 'minor'); axis(ax300, 'tight');
    current_time = double(DPE_start_time) + epoch*DPE_corr_save_interval; % DPE������ʱ��������
    title(ax300, ['��' num2str(current_time) 's������ͼ']);
    view(ax300, -37.5, 30); 
    pause(0.2);
    hold(ax300, 'off');
    
    % ���Ƶ������ǵĸ�������ͼ
    axis_idx = ax500;
    curr_corr_pos_prn = corr_pos_prn(prn_idx, epoch, :);
    [curr_corr_xy_2d_prn, ~] = grid2dPrnPlane(curr_corr_pos_prn, numOfGrid, 0, disp_axis);
%     curr_corr_xy_2d_prn = curr_corr_xy_2d_prn ./ max(max(curr_corr_pos));
    curr_corr_xy_2d_prn = curr_corr_xy_2d_prn ./ max(curr_corr_xy_2d_prn);
    surf(axis_idx, DPE_x, DPE_y, curr_corr_xy_2d_prn);
    axis(axis_idx, [dtmp(1), dtmp(end), dtmp(1), dtmp(end), 0, 1]);
    colormap(axis_idx, slanCM('parula'));
    title(axis_idx, ['PRN' num2str(dpe_prn(prn_idx)) '��' num2str(current_time) 's������ͼ']);
    view(axis_idx, 0, 90); 
    hold(axis_idx, 'on');

end

