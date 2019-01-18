function [UpStateLen, DownStateLen] = plotStateDurationHist(Trans, CutInterval)
%function [UpStateLen, DownStateLen] = plotStateDurationHist(Trans, PeriodToRemove)
%
%  [UpStateLen, DownStateLen] = plotStateDurationHist(Trans[, PeriodToRemove])
%
%
%   Copyright 2009 Maurizio Mattia @ Ist. Super. Sanita', Rome - Italy
%   Version: 1.0 - Sep. 22, 2009
%

if exist('CutInterval')
   for idx = 1:1:size(CutInterval,2)
       PeriodToRemove(1)=CutInterval(1,idx);
       PeriodToRemove(2)=CutInterval(2,idx);
            
    
       ndx = find(Trans.time < PeriodToRemove(1));
       if numel(ndx) > 1
          StateLen = diff(Trans.time(ndx));
          StateVal = Trans.val(ndx(1:end-1));
       else
          StateLen = [];
          StateVal = [];
       end
       ndx = find(Trans.time > PeriodToRemove(2));
       if numel(ndx) > 1
          StateLen = [StateLen diff(Trans.time(ndx))];
          StateVal = [StateVal Trans.val(ndx(1:end-1))];
       end
   end
else
   StateLen = diff(Trans.time);
   StateVal = Trans.val(1:end-1);
end

% Down states...
DownStateLen = StateLen(find(StateVal == 0));

% Up states...
ndxU = find(Trans.val == 1);
if ndxU(end) > length(StateLen)
   ndxU = ndxU(1:end-1);
end
UpStateLen = StateLen(find(StateVal == 1));

figure;
hold on

Options.ValRange = [0 max(max([DownStateLen UpStateLen]))];
% Options.ValRange = [0 2.1];
Options.BinNum = 20;
% Options.BinNum = 22;
Options.EdgeColor = 'b';
Options.FaceColor = 'b';
Options.FaceAlpha = 0.5;

plotHistogramWithOptions(DownStateLen, Options);

Options.BinNum = 40;
Options.EdgeColor = 'r';
Options.FaceColor = 'r';
Options.FaceAlpha = 0.5;

plotHistogramWithOptions(UpStateLen, Options);

YLim = get(gca, 'YLim');
clear hlgn
hlgn(1) = plot([0 0]+mean(DownStateLen), YLim, 'b--', 'LineWidth', 1);
text(mean(DownStateLen), YLim(2)*0.75, ...
   [' ' num2str(mean(DownStateLen), 3) 's ' ...
   '(sd: ' num2str(std(DownStateLen), 3) 's, '...
   'n=' num2str(length(DownStateLen)) ')']);
hlgn(2) = plot([0 0]+mean(UpStateLen), YLim, 'r--', 'LineWidth', 1);
text(mean(UpStateLen), YLim(2)*0.5, ...
   [' ' num2str(mean(UpStateLen), 3) 's ' ...
   '(sd: ' num2str(std(UpStateLen), 3) 's, '...
   'n=' num2str(length(UpStateLen)) ')']);

set(gca, 'Layer', 'top', 'Box', 'off', 'TickDir', 'out');

legend(hlgn, {'Down', 'Up'});
xlabel('State duration (s)')
set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.5 3.5 5 4]);
% print('-deps2c', 'UpDownDuration.hist.eps');
print('-djpeg', '-r300', 'UpDownDuration.hist.jpg');
