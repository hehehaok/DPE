clear all;
close all;
fclose('all');
clc;

%% �����ļ�
trackFile = 'D:\study\��һ\DPESDR\DPEpython\DPE\pre-simulator\finland_cleanStatic\test1_skip0_start0\end-of-40_usrp6\channel_15.mat';
load(trackFile);

%% ������������
T = 1; % һ����Ԫʱ��/ms
timeAxis = T : T : T*double(channel__cpcount); % ʱ����
K = 20;
M = 5;
figure(100);
ax1 = subplot(3, 3, 1); % ɢ��ͼ
ax2 = subplot(3, 3, [2 3]); % ��������
ax3 = subplot(3, 3, [4 5 6]); % ����ɻ���ͼ
ax4 = subplot(3, 3, [7 8 9]); % �����
figure(200);
ax5 = gca; % ��խ�����ʱȷ��������

%% ��ͼ
% ɢ��ͼ
scatter(ax1, channel_array_iP(1:end-1), ...
             channel_array_qP(1:end-1), '.');
grid(ax1);
axis(ax1, 'equal');
title(ax1, 'Discrete-Time Scatter Plot');
xlabel(ax1, 'I prompt');
ylabel(ax1, 'Q prompt');

% ��������
plot(ax2, timeAxis/1000, channel_array_iP(1:end-1), ...
          timeAxis/1000, channel_array_qP(1:end-1));
grid(ax2);
title(ax2, 'In-phase (I_P) and Quad-phase (Q_P) component of the received signal');
xlabel(ax2, 'Time (s)');
axis(ax2, 'tight');
legendHandle = legend(ax2, '${I_P}$','${Q_P}$');
set(legendHandle, 'Interpreter', 'Latex'); 

% ����ɻ���ͼ
plot(ax3, timeAxis/1000, ...
        [sqrt(channel_array_iE(1:end-1).^2 + ...
              channel_array_qE(1:end-1).^2)', ...
         sqrt(channel_array_iP(1:end-1).^2 + ...
              channel_array_qP(1:end-1).^2)', ...
         sqrt(channel_array_iL(1:end-1).^2 + ...
              channel_array_qL(1:end-1).^2)'],'-*');

grid  (ax3);
title (ax3, 'Correlation results');
xlabel(ax3, 'Time (s)');
axis  (ax3, 'tight');
legendHandle = legend(ax3, '$\sqrt{I_{E}^2 + Q_{E}^2}$', ...
                           '$\sqrt{I_{P}^2 + Q_{P}^2}$', ...
                           '$\sqrt{I_{L}^2 + Q_{L}^2}$');
set(legendHandle, 'Interpreter', 'Latex');

% �����
plot(ax4, timeAxis/1000, channel_array_snr(1:end-1));
grid(ax4);
title(ax4, 'CNo');
xlabel(ax4, 'Time (s)');
axis(ax4, 'tight');

% ����һ�������ͼ��
dataCNo = CNoFromTrack(channel_array_iP(1:end-1), channel_array_qP(1:end-1), T/1000, K, M);
timeAxisInTight = (1:length(dataCNo))*K*M/1000;
plot(ax5, timeAxisInTight, dataCNo); 
grid(ax5);
title(ax5, 'CNo');
xlabel(ax5, 'Time (s)');
axis(ax5, 'tight');
% axis([0 450 37 52]);

%% ���������
function dataCNo = CNoFromTrack(ipResults, qpResults, T, K, M)
        CNoInterval = K * M;    
        dataCNo  = zeros(1, floor(length(ipResults)/CNoInterval));
        normalizedPower = 0;
        % ���������
        for loopCnt = 1:length(ipResults)
            if (rem(loopCnt, M)==0)
                wide = sum(ipResults(loopCnt-M+1:loopCnt).^2 + qpResults(loopCnt-M+1:loopCnt).^2);
                narrow = sum(abs(ipResults(loopCnt-M+1:loopCnt)))^2 + sum(abs(qpResults(loopCnt-M+1:loopCnt)))^2; 
                normalizedPower = narrow / wide + normalizedPower;
                if (rem(loopCnt, CNoInterval)==0)
                    CNoCnt = loopCnt/CNoInterval;
                    meanNormalizedPower = normalizedPower / K;
                    dataCNo(CNoCnt) = abs(10*log10((1/T) * (meanNormalizedPower-1) / (M-meanNormalizedPower)));
                    normalizedPower = 0;
                end
            end
        end 
end



