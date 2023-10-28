clear all;
close all;
fclose('all');
clc;

%% 保存位置
dir = 'D:\academic\研究生\研1\FGI-GSRx\test\CPNTBAT\processResult\GPS L1 trackResult\';   % 结果存放路径
refName = 'trackResult_GPSL1_cpntbat_cs01.mat';
channelIdx = 6;
load(refName);
[CNo_cs01, timeAxisInMs] = CNoFromTrack(trackData01, settings);

%% cs05
csName = 'trackResult_GPSL1_cpntbat_cs05.mat';
load(csName);
[CNo_cs05, ~] = CNoFromTrack(trackData05, settings);
figure(500);
plot(timeAxisInMs/1000, CNo_cs01(channelIdx, :), 'LineWidth', 2); 
hold on; grid on; 
plot(timeAxisInMs/1000, CNo_cs05(channelIdx, :), 'LineWidth', 2); 
axis([0 450 37 52]);
legend('cs01(无欺骗)', 'cs05(静态小功率位置欺骗)');
saveas(500, [dir 'cs05载噪比.png']);

% %% cs06
% csName = 'trackResult_GPSL1_cpntbat_cs06.mat';
% load(csName);
% [CNo_cs06, ~] = CNoFromTrack(trackData06, settings);
% figure(600);
% plot(timeAxisInMs/1000, CNo_cs01(channelIdx, :), 'LineWidth', 2); 
% hold on; grid on; 
% plot(timeAxisInMs/1000, CNo_cs06(channelIdx, :), 'LineWidth', 2); 
% axis([0 450 37 47]);
% legend('cs01(无欺骗)', 'cs06(静态小功率位置欺骗)');
% saveas(600, [dir 'cs06载噪比.png']);
% 
% %% cs07
% csName = 'trackResult_GPSL1_cpntbat_cs07.mat';
% load(csName);
% [CNo_cs07, ~] = CNoFromTrack(trackData07, settings);
% figure(700);
% plot(timeAxisInMs/1000, CNo_cs01(channelIdx, :), 'LineWidth', 2); 
% hold on; grid on; 
% plot(timeAxisInMs/1000, CNo_cs07(channelIdx, :), 'LineWidth', 2); 
% axis([0 450 37 47]);
% legend('cs01(无欺骗)', 'cs07(静态中等功率位置欺骗)');
% saveas(700, [dir 'cs07载噪比.png']);
% 
% %% cs08
% csName = 'trackResult_GPSL1_cpntbat_cs08.mat';
% load(csName);
% [CNo_cs08, ~] = CNoFromTrack(trackData08, settings);
% figure(800);
% plot(timeAxisInMs/1000, CNo_cs01(channelIdx, :), 'LineWidth', 2); 
% hold on; grid on; 
% plot(timeAxisInMs/1000, CNo_cs08(channelIdx, :), 'LineWidth', 2); 
% axis([0 450 37 47]);
% legend('cs01(无欺骗)', 'cs08(静态大功率位置欺骗)');
% saveas(800, [dir 'cs08载噪比.png']);
% 
% %% cs09
% csName = 'trackResult_GPSL1_cpntbat_cs09.mat';
% load(csName);
% [CNo_cs09, ~] = CNoFromTrack(trackData09, settings);
% figure(900);
% plot(timeAxisInMs/1000, CNo_cs01(channelIdx, :), 'LineWidth', 2); 
% hold on; grid on; 
% plot(timeAxisInMs/1000, CNo_cs09(channelIdx, 11:end), 'LineWidth', 2); % cs09的处理时间为0-420s，改为10-420s与其他保持一致
% axis([0 450 37 50]);
% legend('cs01(无欺骗)', 'cs09(动态大功率位置欺骗)');
% saveas(900, [dir 'cs09载噪比.png']);
% 
% %% cs10
% csName = 'trackResult_GPSL1_cpntbat_cs10.mat';
% load(csName);
% [CNo_cs10, ~] = CNoFromTrack(trackData10, settings);
% figure(1000);
% plot(timeAxisInMs/1000, CNo_cs01(channelIdx, :), 'LineWidth', 2); 
% hold on; grid on; 
% plot(timeAxisInMs/1000, CNo_cs10(channelIdx, 11:end), 'LineWidth', 2); % cs10的处理时间为0-420s，改为10-420s与其他保持一致
% axis([0 450 34 47]);
% legend('cs01(无欺骗)', 'cs10(动态小功率位置欺骗)');
% saveas(1000, [dir 'cs10载噪比.png']);

%% 计算载噪比
function [DataCNo, timeAxisInMs] = CNoFromTrack(tR, allSettings)
    for signalNr = 1:allSettings.sys.nrOfSignals
        signal = allSettings.sys.enabledSignals{signalNr};
        NumToProcess =  round(allSettings.sys.msToProcess/1000/allSettings.(signal).Nc);
        K = allSettings.(signal).K;
        M = allSettings.(signal).M;
        CNoInterval = K * M;    
        DataCNo  = zeros(tR.(signal).nrObs, floor(NumToProcess/CNoInterval));
        sampleSpacing = allSettings.(signal).Nc*1000;   
        % 计算载噪比
        for ii = 1 : tR.(signal).nrObs
            ipResults = tR.(signal).channel(ii).I_P(sampleSpacing:sampleSpacing:end);
            qpResults = tR.(signal).channel(ii).Q_P(sampleSpacing:sampleSpacing:end);
            normalizedPower = 0;
            for loopCnt = 1:NumToProcess
                if (rem(loopCnt, M)==0)
                    wide = sum(ipResults(loopCnt-M+1:loopCnt).^2 + qpResults(loopCnt-M+1:loopCnt).^2);
                    narrow = sum(abs(ipResults(loopCnt-M+1:loopCnt)))^2 + sum(abs(qpResults(loopCnt-M+1:loopCnt)))^2; 
                    normalizedPower = narrow / wide + normalizedPower;
                    if (rem(loopCnt, CNoInterval)==0)
                        CNoCnt = loopCnt/CNoInterval;
                        meanNormalizedPower = normalizedPower / K;
                        DataCNo(ii, CNoCnt) = abs(10*log10((1/(allSettings.(signal).Nc)) * (meanNormalizedPower-1) / (M-meanNormalizedPower)));
                        normalizedPower = 0;
                    end
                end
            end 
            timeAxisInMs = sampleSpacing*CNoInterval:sampleSpacing*CNoInterval:length(tR.(signal).channel(ii).I_P);
        end
    end
end