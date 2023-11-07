function [] = plotCorrScore(dataDir, secDir)

%% ��ȡDPE��λ����ļ�
postdir  = ['./post-simulator/' dataDir];

% ����������ؽ��
receiveFile = 'end-of-dp/costFunction.h5';
scoreDir = [postdir secDir receiveFile];
dpe_prn = sort(h5read(scoreDir,'/GRID/dpe_prn')');
code_idc = h5read(scoreDir,'/PRN/code_idc')';	
corr_code = permute(h5read(scoreDir,'/PRN/code_corr_prn'), [3 2 1]); 
corr_code_rc = permute(h5read(scoreDir,'/PRN/corr_pos_prn_rc'), [3 2 1]);

% ��������ʱ����ز���
dpeRunTimeFile = 'dpe_runtime.mat';
load([postdir secDir dpeRunTimeFile]);

% ����
prn_idx = 1;
curr_epoch = 1;

%% ��֤��������������
f400 = figure(400);
ax400 = subplot(1,1,1);
ax401 = axes('Position', [0.6 0.6 0.28 0.28]);
ax402 = axes('Position', [0.2 0.6 0.28 0.28]);

corr_code_prn = reshape(corr_code(prn_idx, curr_epoch, :), [], 20); % �������ĵ�������źŵ����ֵ
corr_code_prn = fftshift(sum(corr_code_prn, 2));
% corr_code_prn = sum(corr_code_prn, 2);
corr_code_prn = corr_code_prn/max(corr_code_prn);

corr_code_rc_prn = reshape(corr_code_rc(prn_idx, curr_epoch, :), 1, []); % ��������ͬ������ֵ-����λ
corr_code_rc_prn = floor(mod(corr_code_rc_prn, length(code_idc))) + 1;
idx1 = find(corr_code_rc_prn > length(code_idc)/2);
idx2 = find(corr_code_rc_prn < length(code_idc)/2);
corr_code_rc_prn(idx1) = corr_code_rc_prn(idx1) - length(code_idc)/2;
corr_code_rc_prn(idx2) = corr_code_rc_prn(idx2) + length(code_idc)/2; 

plot(ax400, code_idc, corr_code_prn);hold(ax400, 'on');
scatter(ax400, code_idc(corr_code_rc_prn), corr_code_prn(corr_code_rc_prn))
axis(ax400, 'tight');

idx = find(code_idc > -2.1 & code_idc < 2.1);
plot(ax401, code_idc(idx), corr_code_prn(idx));hold(ax401, 'on');
scatter(ax401, code_idc(corr_code_rc_prn), corr_code_prn(corr_code_rc_prn));
axis(ax401, 'tight');

edges = [-10 -1.1:0.2:1.1 10];
% N = histcounts();
histogram(ax402, code_idc(corr_code_rc_prn), edges);

