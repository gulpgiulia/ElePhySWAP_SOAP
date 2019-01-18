function plotLFPlogMUAandUD(LFP, LogMUA, UD, PeriodToPlot)
%
%  plotLFPlogMUAandUD(LFP, LogMUA, UD, PeriodToPlot)
%
%
%   Copyright 2008 Maurizio Mattia @ Ist. Super. Sanit�, Rome - Italy
%   Version: 1.0 - Oct. 19, 2008
%

Y_WINDOW = 8;

figure;
hold on;

n = 1;
ndx = find(LogMUA.time >= PeriodToPlot(1) & LogMUA.time <= PeriodToPlot(2));
hlgn(n) = plot(LogMUA.time(ndx), ...
   (LogMUA.value(ndx) - mean(LogMUA.value(ndx))) / std(LogMUA.value(ndx)), 'k'); % (normalised)
slgn{n} = 'log(MUA)';

if length(LFP) > 0
   n = n + 1;
   ndx = find(LFP.time >= PeriodToPlot(1) & LFP.time <= PeriodToPlot(2));
   hlgn(n) = plot(LFP.time(ndx), ...
      (LFP.value(ndx) - mean(LFP.value(ndx)))/std(LFP.value(ndx)) - Y_WINDOW, 'b');
   slgn{n} = 'LFP';
end

n = n + 1;
ndx = find(UD.time >= PeriodToPlot(1) & UD.time <= PeriodToPlot(2));
hlgn(n) = plot(UD.time(ndx), (UD.value(ndx)-0.5)*(Y_WINDOW+10)-Y_WINDOW/2, 'r');
slgn{n} = 'Up/Down';


set(gca, 'TickDir', 'out', 'Box', 'off', 'Layer', 'top', ...
   'XLim', PeriodToPlot, 'YLim', [-Y_WINDOW-6 +6]);
%legend(hlgn, slgn, -1); 
legend(hlgn, slgn); % modified for R2016a
ylabel('Y (a.u.)');
xlabel('Time (s)');

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0.0 3.5 8 4]);
print('-deps2c', ['LFPlogMUAandUD.eps']);
