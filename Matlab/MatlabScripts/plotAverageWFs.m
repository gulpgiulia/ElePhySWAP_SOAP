function [tWF, MeanWFs, StdWFs] = plotAverageWFs(WFs, SingleWFtoPlot, YRange)
%
%  [tWF, MeanWFs, StdWFs] = plotAverageWFs(WFs[, SingleWFtoPlot[, YRange]])
%
%
%   Copyright 2014 Maurizio Mattia @ Ist. Super. Sanit�, Rome - Italy
%   Version: 1.0 - Jan. 10, 2014
%

SAMPLED_WFs = 10;
GRAY_CLR = [0.9 0.9 0.9];

if exist('SingleWFtoPlot') == 0
   SingleWFtoPlot = SAMPLED_WFs;
end


MeanWFs = zeros(size(WFs(1).value));
StdWFs = zeros(size(WFs(1).value));
for n = 1:length(WFs)
   MeanWFs = MeanWFs + WFs(n).value;
   StdWFs = StdWFs + (WFs(n).value).^2;
end
MeanWFs = MeanWFs / length(WFs);
StdWFs = sqrt(StdWFs / length(WFs) - MeanWFs.^2);
tWF = WFs(1).time;
if size(tWF,1) > 1
   tWF = tWF';
end
if size(MeanWFs,1) > 1
   MeanWFs = MeanWFs';
   StdWFs = StdWFs';
end % (row vectors)

figure;
hold on;

XP = [tWF fliplr(tWF)];
YP = [MeanWFs+StdWFs fliplr(MeanWFs-StdWFs)];
patch(XP, YP, GRAY_CLR, 'EdgeColor', 'none');

clr = jet(min([length(WFs) SingleWFtoPlot]));
for n = 1:min([length(WFs) SingleWFtoPlot])
   hlgn(n) = plot(WFs(n).time, WFs(n).value, 'Color', clr(n,:));
   slgn{n} = ['WF ' num2str(n)];
end

hlgn(length(hlgn)+1) = plot(tWF, MeanWFs, 'k', 'LineWidth', 1);
slgn{length(hlgn)} = 'Mean WF';
%legend(hlgn, slgn, -1); 
legend(hlgn, slgn); % modified for R2016a

if exist('YRange') == 1
   YLim = YRange;
else
   YLim = get(gca, 'Ylim');
end
set(gca, 'YLim', YLim);

text(WFs(1).time(end), YLim(2), ['(n = ' num2str(length(WFs)) ')'], ...
   'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
