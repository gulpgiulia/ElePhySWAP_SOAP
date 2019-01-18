% analyzeResultsOfRecordingSet.m
%
%   Copyright 2017 Giulia De Bonis @ INFN, Rome - Italy
%   Version: 1.0 - April 4, 2017
%


%% Load analysis params and information about recording dataset.
%
setParamsAndOptions
nPlot = 1; % number of plots to be produced for each observable

%% Load AnalysisSummary.mat for each channel in the RecordingSet
%
for crsRecSet = ChannelSet % modified for running on a sub-set of recording channels (a cortical area)
   
   workingdir = [AnalysisDir RecordingSet(crsRecSet).label];
   results(crsRecSet)=load([workingdir '/AnalysisSummary.mat']);
   duration(crsRecSet)=load([workingdir '/UpDownStateDuration.mat']);
end

save([AnalysisDir 'Results.mat'], 'results');
save([AnalysisDir 'Duration.mat'], 'duration');


%% PyyBaseline
%

Frequency = 1.0e+03*[0.1984 0.3968 0.5952 0.7937 0.9921 1.1905 1.3889];

figure 
hold on
for crsRecSet = ChannelSet
    plot(Frequency,results(crsRecSet).PyyBaseline, '.-b', 'LineWidth', 1.);
end

text(1600, results(25).PyyBaseline(end), ['Ch.25 '], ...
   'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

set(gca, 'XLim', [Frequency(1)*0.90 Frequency(end)*1.10], 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'TickDir', 'out', 'Layer', 'top', 'Box', 'off');

xlabel('\omega/2\pi (Hz)');
ylabel('LFP p.s.d. (\muV^2/Hz)');

% % ... to find out which channel is 'out of scale'
% for crsRecSet = ChannelSet
%     if(results(crsRecSet).PyyBaseline(1)>1e4) disp(crsRecSet); end;
% end

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
print('-deps2c', [AnalysisDir 'Summary/' 'PyyBaselineAllChannels.eps']);

%% PyyBaseline - Mean, Std, StdErr
%

PyyMatrix = [];
for crsRecSet = ChannelSet([1:24 26:end]) % exclude Ch.25 (out of scale)
PyyMatrix=[PyyMatrix results(crsRecSet).PyyBaseline];
end

x = Frequency;
y = mean(PyyMatrix,2); % mean value over the Channel Set for each frequency
err1 = std(PyyMatrix,0,2)/sqrt(numel(x)); % standard error
err2 = std(PyyMatrix,0,2); % standard deviation

figure
hold on

% PYYBaseline of the different channels on the background...
GRAY_CLR = [0.9 0.9 0.9];
for crsRecSet = ChannelSet([1:24 26:end])
    plot(Frequency,results(crsRecSet).PyyBaseline, 'color', GRAY_CLR, 'LineWidth', 1.);
end

n=1;
hlgn(n) = errorbar(x,y,err1,'.-r', 'LineWidth', 1.);
slgn{n} = 'errorbar = standard error';
n=2;
hlgn(n) = errorbar(x,y,err2,'.-b', 'LineWidth', 1.);
slgn{n} = 'errorbar = standard deviation';
n=3;
hlgn(n) = plot(x,y, '.-k', 'LineWidth', 1.);
slgn{n} = 'mean';
% errorbar(x,y,err,'.-r', 'LineWidth', 1.,'CapSize',18);
% CapSize is not available in release R2016

legend(hlgn, slgn);

set(gca, 'XLim', [Frequency(1)*0.90 Frequency(end)*1.10], 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'TickDir', 'out', 'Layer', 'top', 'Box', 'off');

xlabel('\omega/2\pi (Hz)');
ylabel('LFP p.s.d. (\muV^2/Hz)');

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
print('-deps2c', [AnalysisDir 'Summary/' 'PyyBaselineMeanErrorbars.eps']);

%% PyyBaseline - Median, Percentile
%

MedianOverChannels = prctile(PyyMatrix, 50, 2); 
% prctile(x,50) --> the median of x
% PyyMatrix, each column is a channel. DIM = 2 --> median over channels

figure
hold on

% PYYBaseline of the different channels on the background...
for crsRecSet = ChannelSet([1:24 26:end])
    plot(Frequency,results(crsRecSet).PyyBaseline, 'color', GRAY_CLR, 'LineWidth', 1.);
end


% Patches representing different percetiles...

XP = [Frequency fliplr(Frequency)];
for prc = 10:10:40 % (percentile)
   clr = [1-prc/100 1-prc/100 1]; % color (RGB triplet)
   YP = prctile(PyyMatrix,[prc 100-prc],2); 
   YP = [YP(:,1) flipud(YP(:,2))];
   YP = reshape(YP,1,numel(YP));
   patch(XP, YP, clr, 'EdgeColor', 'none');
end

hlgn = []; slgn = []; 
n=1;
hlgn(n) = plot(Frequency, MedianOverChannels, '.-b', 'LineWidth', 1.);
slgn{n} = 'median';
% n=2;
% hlgn(n) = plot(Frequency, y, '.-r', 'LineWidth', 1.);
% slgn{n} = 'mean';
legend(hlgn, slgn);

set(gca, 'XLim', [Frequency(1) Frequency(end)], 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'TickDir', 'out', 'Layer', 'top', 'Box', 'off');

xlabel('\omega/2\pi (Hz)');
ylabel('LFP p.s.d. (\muV^2/Hz)');

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
print('-deps2c', [AnalysisDir 'Summary/' 'PyyBaselineMedianPercentile.eps']);


%% Mode Params: Mu1
%

Mu1 = NaN(1,numel(RecordingSet));
for crsRecSet = ChannelSet
    Mu1(crsRecSet)=(results(crsRecSet).ModeParams.Mu1);
end

% resultsMu1 = plotSummaryResults(Mu1,nPlot);
resultsMu1 = plotVarCorticalAreas(Mu1,nPlot);
varName = varname(Mu1);
for i=[1:nPlot]; set(resultsMu1.fig(i), 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]); end;
% figure(resultsMu1.fig(1));
% ylabel('Mu1');
% figure(resultsMu1.fig(2));
% ylabel('Mu1');
% figure(resultsMu1.fig(3));
% xlabel('Mu1');
print(resultsMu1.fig(1),'-deps2c', [AnalysisDir 'Summary/' 'Mu1.eps']);
print(resultsMu1.fig(2),'-deps2c', [AnalysisDir 'Summary/' 'Mu1Percentile.eps']);
print(resultsMu1.fig(3),'-deps2c', [AnalysisDir 'Summary/' 'Mu1Histo.eps']);


%% Mode Params: Sigma1
%

Sigma1 = NaN(1,numel(RecordingSet));
for crsRecSet = ChannelSet
    Sigma1(crsRecSet)=(results(crsRecSet).ModeParams.Sigma1);
end

resultsSigma1 = plotSummaryResults(Sigma1);
for i=[1:3]; set(resultsSigma1.fig(i), 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]); end;
figure(resultsSigma1.fig(1));
ylabel('Sigma1');
figure(resultsSigma1.fig(2));
ylabel('Sigma1');
figure(resultsSigma1.fig(3));
xlabel('Sigma1');
print(resultsSigma1.fig(1),'-deps2c', [AnalysisDir 'Summary/' 'Sigma1.eps']);
print(resultsSigma1.fig(2),'-deps2c', [AnalysisDir 'Summary/' 'Sigma1Percentile.eps']);
print(resultsSigma1.fig(3),'-deps2c', [AnalysisDir 'Summary/' 'Sigma1Histo.eps']);


%% Mode Params: Ampl1
%

Ampl1 = NaN(1,numel(RecordingSet));
for crsRecSet = ChannelSet
    Ampl1(crsRecSet)=(results(crsRecSet).ModeParams.Ampl1);
end

resultsAmpl1 = plotSummaryResults(Ampl1);
for i=[1:3]; set(resultsAmpl1.fig(i), 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]); end;
figure(resultsAmpl1.fig(1));
ylabel('Ampl1');
figure(resultsAmpl1.fig(2));
ylabel('Ampl1');
figure(resultsAmpl1.fig(3));
xlabel('Ampl1');
print(resultsAmpl1.fig(1),'-deps2c', [AnalysisDir 'Summary/' 'Ampl1.eps']);
print(resultsAmpl1.fig(2),'-deps2c', [AnalysisDir 'Summary/' 'Ampl1Percentile.eps']);
print(resultsAmpl1.fig(3),'-deps2c', [AnalysisDir 'Summary/' 'Ampl1Histo.eps']);


%% Mode Params: Mu2
%

Mu2 = NaN(1,numel(RecordingSet));
for crsRecSet = ChannelSet
    Mu2(crsRecSet)=(results(crsRecSet).ModeParams.Mu2);
end

resultsMu2 = plotSummaryResults(Mu2);
for i=[1:3]; set(resultsMu2.fig(i), 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]); end;
figure(resultsMu2.fig(1));
ylabel('Mu2');
figure(resultsMu2.fig(2));
ylabel('Mu2');
figure(resultsMu2.fig(3));
xlabel('Mu2');
print(resultsMu2.fig(1),'-deps2c', [AnalysisDir 'Summary/' 'Mu2.eps']);
print(resultsMu2.fig(2),'-deps2c', [AnalysisDir 'Summary/' 'Mu2Percentile.eps']);
print(resultsMu2.fig(3),'-deps2c', [AnalysisDir 'Summary/' 'Mu2Histo.eps']);


%% Mode Params: M2 -- cortical areas
%
CorticalArea = ['V' 'R' 'P' 'S' 'M'];
meanMu2 = [];
stderrMu2 = [];

for area = CorticalArea(1:end)
    if area == 'V'     % V = Visual Cortex
    subset = [18 19 20 22 23 24 25 27 28 29 30 31 32];
    Mu2V = Mu2(subset);
    resultsMu2V = plotSummaryResults(Mu2V, subset);   
    figure(resultsMu2V.fig(1));
    ylabel('Mu2V');    
    for i = 2:3; close(resultsMu2V.fig(i)); end;  
    meanMu2 = [meanMu2 resultsMu2V.mean];
    stderrMu2 = [stderrMu2 resultsMu2V.stderr];
elseif area == 'R' % R = Retrosplenial Cortex 
    subset = [13 17 21 26];
    Mu2R = Mu2(subset);
    resultsMu2R = plotSummaryResults(Mu2R, subset);   
    figure(resultsMu2R.fig(1));
    ylabel('Mu2R');
    for i = 2:3; close(resultsMu2R.fig(i)); end; 
    meanMu2 = [meanMu2 resultsMu2R.mean];
    stderrMu2 = [stderrMu2 resultsMu2R.stderr];
elseif area == 'P' % P = Parietal Association Area (PtA)
    subset = [14 15 16];
    Mu2P = Mu2(subset);
    resultsMu2P = plotSummaryResults(Mu2P, subset);    
    figure(resultsMu2P.fig(1));
    ylabel('Mu2P');
    for i = 2:3; close(resultsMu2P.fig(i)); end; 
    meanMu2 = [meanMu2 resultsMu2P.mean];
    stderrMu2 = [stderrMu2 resultsMu2P.stderr];
elseif area == 'S' % S = Somatosensory Cortex
    subset = [4 5 7 8 10 11 12];
    Mu2S = Mu2(subset);
    resultsMu2S = plotSummaryResults(Mu2S, subset);   
    figure(resultsMu2S.fig(1));
    ylabel('Mu2S');
    for i = 2:3; close(resultsMu2S.fig(i)); end;    
    meanMu2 = [meanMu2 resultsMu2S.mean];
    stderrMu2 = [stderrMu2 resultsMu2S.stderr];
elseif area == 'M' % M = Motor Cortex
    subset = [1 2 3 6 9];
    Mu2M = Mu2(subset);
    resultsMu2M = plotSummaryResults(Mu2M, subset);    
    figure(resultsMu2M.fig(1));
    ylabel('Mu2M');
    for i = 2:3; close(resultsMu2M.fig(i)); end;  
    meanMu2 = [meanMu2 resultsMu2M.mean];
    stderrMu2 = [stderrMu2 resultsMu2M.stderr];
end
end

%% Mu2 -- figure with subplots 
%

figure(resultsMu2.fig(1)); ax1 = gca; 
fig1 = get(ax1,'children'); xR1 = get(ax1, 'XLim'); 
xL1 = get(ax1, 'XTickLabel'); xT1 = get(ax1, 'XTick');
figure(resultsMu2V.fig(1)); ax2 = gca; 
fig2 = get(ax2,'children'); xR2 = get(ax2, 'XLim'); yR2 = get(ax2, 'YLim'); 
xL2 = get(ax2, 'XTickLabel'); xT2 = get(ax2, 'XTick');
figure(resultsMu2R.fig(1)); ax3 = gca; 
fig3 = get(ax3,'children'); xR3 = get(ax3, 'XLim'); yR3 = get(ax3, 'YLim'); 
xL3 = get(ax3, 'XTickLabel'); xT3 = get(ax3, 'XTick');
figure(resultsMu2P.fig(1)); ax4 = gca; 
fig4 = get(ax4,'children'); xR4 = get(ax4, 'XLim'); yR4 = get(ax4, 'YLim'); 
xL4 = get(ax4, 'XTickLabel'); xT4 = get(ax4, 'XTick');
figure(resultsMu2S.fig(1)); ax5 = gca; 
fig5 = get(ax5,'children'); xR5 = get(ax5, 'XLim'); yR5 = get(ax5, 'YLim'); 
xL5 = get(ax5, 'XTickLabel'); xT5 = get(ax5, 'XTick');
figure(resultsMu2M.fig(1)); ax6 = gca; 
fig6 = get(ax6,'children'); xR6 = get(ax6, 'XLim'); yR6 = get(ax6, 'YLim'); 
xL6 = get(ax6, 'XTickLabel'); xT6 = get(ax6, 'XTick');

figure; %create new figure
s1 = subplot(3,3,[1 2 4 5]); copyobj(fig1,s1); title('Mu2');
set(s1,'XLim',xR1,'XTickLabel',xL1,'XTick',xT1); 
s2 = subplot(3,3,3); copyobj(fig2,s2); title('Mu2V');
set(s2,'XLim',xR2,'YLim',yR2,'XTickLabel',xL2,'XTick',xT2);%,'FontSize', 4); 
s3 = subplot(3,3,6); copyobj(fig3,s3); title('Mu2R');
set(s3,'XLim',xR3,'YLim',yR3,'XTickLabel',xL3,'XTick',xT3); 
s4 = subplot(3,3,7); copyobj(fig4,s4); title('Mu2P');
set(s4,'XLim',xR4,'YLim',yR4,'XTickLabel',xL4,'XTick',xT4); 
s5 = subplot(3,3,8); copyobj(fig5,s5); title('Mu2S');
set(s5,'XLim',xR5,'YLim',yR5,'XTickLabel',xL5,'XTick',xT5); %,'FontSize', 8); 
s6 = subplot(3,3,9); copyobj(fig6,s6); title('Mu2M');
set(s6,'XLim',xR6,'YLim',yR6,'XTickLabel',xL6,'XTick',xT6); 

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
print('-deps2c', [AnalysisDir 'Summary/' 'Mu2_CorticalAreas.eps']);

%% Mu2Mean -- Cortical Areas
%

figure
hold on
err1Mu2 = [resultsMu2.mean + resultsMu2.stderr resultsMu2.mean + resultsMu2.stderr...
    resultsMu2.mean - resultsMu2.stderr resultsMu2.mean - resultsMu2.stderr];
err2Mu2 = [resultsMu2.mean + resultsMu2.std resultsMu2.mean + resultsMu2.std...
    resultsMu2.mean - resultsMu2.std resultsMu2.mean - resultsMu2.std];

n=3;
hlgn(n)=patch([0.5 5.5 5.5 0.5], err2Mu2, SHAD_RED_2, 'EdgeColor', 'none');
slgn{n} = 'standard deviation';
n=2;
hlgn(n)=patch([0.5 5.5 5.5 0.5], err1Mu2, SHAD_RED_1, 'EdgeColor', 'none');
slgn{n} = 'standard error';
errorbar([1:numel(CorticalArea)],meanMu2,stderrMu2,'ob', 'LineWidth', 1.);
set(gca,'xticklabels',[' ' cellstr(CorticalArea')' ' ']);
n=1;
hlgn(n)=plot([0.5 5.5], [resultsMu2.mean resultsMu2.mean], 'r-', 'LineWidth', 1.5);
slgn{n} = 'mean';

xlabel('Cortical Area');
ylabel('Mu2');

legend(hlgn, slgn, 'Location','southeast')

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
print('-deps2c', [AnalysisDir 'Summary/' 'Mu2Mean_CorticalAreas.eps']);

%var = [var sprintf('resultsMu2%c.mean ',area)]; 


%% Mode Params: Sigma2
%

Sigma2 = NaN(1,numel(RecordingSet));
for crsRecSet = ChannelSet
    Sigma2(crsRecSet)=(results(crsRecSet).ModeParams.Sigma2);
end

resultsSigma2 = plotSummaryResults(Sigma2);
for i=[1:3]; set(resultsSigma2.fig(i), 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]); end;
figure(resultsSigma2.fig(1));
ylabel('Sigma2');
figure(resultsSigma2.fig(2));
ylabel('Sigma2');
figure(resultsSigma2.fig(3));
xlabel('Sigma2');
print(resultsSigma2.fig(1),'-deps2c', [AnalysisDir 'Summary/' 'Sigma2.eps']);
print(resultsSigma2.fig(2),'-deps2c', [AnalysisDir 'Summary/' 'Sigma2Percentile.eps']);
print(resultsSigma2.fig(3),'-deps2c', [AnalysisDir 'Summary/' 'Sigma2Histo.eps']);


%% Mode Params: Ampl2
%

Ampl2 = NaN(1,numel(RecordingSet));
for crsRecSet = ChannelSet
    Ampl2(crsRecSet)=(results(crsRecSet).ModeParams.Ampl2);
end

resultsAmpl2 = plotSummaryResults(Ampl2);
for i=[1:3]; set(resultsAmpl2.fig(i), 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]); end;
figure(resultsAmpl2.fig(1));
ylabel('Ampl2');
figure(resultsAmpl2.fig(2));
ylabel('Ampl2');
figure(resultsAmpl2.fig(3));
xlabel('Ampl2');
print(resultsAmpl2.fig(1),'-deps2c', [AnalysisDir 'Summary/' 'Ampl2.eps']);
print(resultsAmpl2.fig(2),'-deps2c', [AnalysisDir 'Summary/' 'Ampl2Percentile.eps']);
print(resultsAmpl2.fig(3),'-deps2c', [AnalysisDir 'Summary/' 'Ampl2Histo.eps']);


%% Mode Params: Ampl1/Ampl2
%

AmplRatio = Ampl1./Ampl2;
resultsAmplRatio = plotSummaryResults(AmplRatio);
for i=[1:3]; set(resultsAmplRatio.fig(i), 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]); end;
figure(resultsAmplRatio.fig(1));
ylabel('Ampl1/Ampl2');
figure(resultsAmplRatio.fig(2));
ylabel('Ampl1/Ampl2');
figure(resultsAmplRatio.fig(3));
xlabel('Ampl1/Ampl2');
print(resultsAmplRatio.fig(1),'-deps2c', [AnalysisDir 'Summary/' 'AmplRatio.eps']);
print(resultsAmplRatio.fig(2),'-deps2c', [AnalysisDir 'Summary/' 'AmplRatioPercentile.eps']);
print(resultsAmplRatio.fig(3),'-deps2c', [AnalysisDir 'Summary/' 'AmplRatioHisto.eps']);


%% Mode Params: deltaMu = Mu2 - Mu1 ~ Mu2 (Mu1 ~ 0 for all the channels)
%
%% DownState
%

DownState = NaN(1,numel(RecordingSet));
for crsRecSet = ChannelSet
    DownState(crsRecSet)=(results(crsRecSet).DownState.mean_duration);
    DownStateStd(crsRecSet)=(results(crsRecSet).DownState.std_duration);
    DownStateErr(crsRecSet)=DownStateStd(crsRecSet)/sqrt(results(crsRecSet).DownState.numel_duration);
end

figure
hold on
hlgn = []; slgn = [];
err2 = std(DownState); % standard deviation
E2 = mean(DownState) + [err2 err2 -err2 -err2];
n=3;
hlgn(n) = patch([0.5 33 33 0.5], E2, SHAD_RED_2, 'EdgeColor', 'none');
slgn{n} = 'standard deviation';
err1 = std(DownState)/sqrt(numel(DownState)); % standard error
E1 = mean(DownState) + [err1 err1 -err1 -err1];
n=2;
hlgn(n) = patch([0.5 33 33 0.5], E1, SHAD_RED_1, 'EdgeColor', 'none');
slgn{n} = 'standard error';
errorbar([1:numel(ChannelSet)],DownState,DownStateErr,'ob', 'LineWidth', 1.);
hold on
n=1;
hlgn(n) = plot([0 numel(ChannelSet)], [mean(DownState) mean(DownState)], 'r-', 'LineWidth', 1.5);
slgn{n} = 'mean';
set(gca, 'XLim', [0 33]);
xlabel('Channels');
ylabel('Down State Duration [s]');
legend(hlgn, slgn, 'Location','northwest');

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
print('-deps2c', [AnalysisDir 'Summary/' 'DownStateErrorBar.eps']);

resultsDownState = plotSummaryResults(DownState,1);
for i=[1:3]; set(resultsDownState.fig(i), 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]); end;
figure(resultsDownState.fig(1));
ylabel('Down State Duration [s]');
figure(resultsDownState.fig(2));
ylabel('Down State Duration [s]');
figure(resultsDownState.fig(3));
xlabel('Down State Duration [s]');
print(resultsDownState.fig(1),'-deps2c', [AnalysisDir 'Summary/' 'DownState.eps']);
print(resultsDownState.fig(2),'-deps2c', [AnalysisDir 'Summary/' 'DownStatePercentile.eps']);
print(resultsDownState.fig(3),'-deps2c', [AnalysisDir 'Summary/' 'DownStateHisto.eps']);

%% Variable ordered by cortical areas

%% DownState -- Cortical Areas
%

resultsDownState = plotVarCorticalAreas(DownState);
resultsDownState.fig(1); yLab = get(get(gca,'YLabel'),'String'); ylabel([yLab ' [s]']);
print(resultsDownState.fig(1),'-deps2c', [AnalysisDir 'Summary/' 'DownState_CorticalAreas_Opt2.eps']);
