function [UDcycleLen] = plotDurationHist(Trans, CutInterval)
%function [UDcycleLen] = plotDurationHist(Trans, PeriodToRemove)
%
%  [UpStateLen, DownStateLen] = plotDurationHist(Trans[, PeriodToRemove])
%  (from M.Mattia, "plotStateDurationHist") 
%   Superimpose a third histo. with the distribution of the UDcycle duration
%
%   Copyright 2017 Giulia De Bonis @ INFNs, Rome - Italy
%   Version: 1.0 - Mrch 30, 2017
%

% if exist('PeriodToRemove')
%    ndx = find(Trans.time < PeriodToRemove(1));
%    if numel(ndx) > 1
%       StateLen = diff(Trans.time(ndx));
%       StateVal = Trans.val(ndx(1:end-1));
%    else
%       StateLen = [];
%       StateVal = [];
%    end
%    ndx = find(Trans.time > PeriodToRemove(2));
%    if numel(ndx) > 1
%       StateLen = [StateLen diff(Trans.time(ndx))];
%       StateVal = [StateVal Trans.val(ndx(1:end-1))];
%    end
% else
%    StateLen = diff(Trans.time);
%    StateVal = Trans.val(1:end-1);
% end

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

% HERE the changes introduced to add the histogram of the UDcycle duration
nCycle=min(length(UpStateLen), length(DownStateLen));
% (takes into account the case in which the number of Up-states and Down-states differs)
UDcycleLen=DownStateLen(1:nCycle)+UpStateLen(1:nCycle); % duration of the UD cycle

Options.ValRange = [0 max(max([DownStateLen UpStateLen UDcycleLen]))];

Options.BinNum = 20; % 22;
Options.EdgeColor = 'b';
Options.FaceColor = 'b';
Options.FaceAlpha = 0.5;

plotHistogramWithOptions(DownStateLen, Options);

Options.BinNum = 40;
Options.EdgeColor = 'r';
Options.FaceColor = 'r';
Options.FaceAlpha = 0.5;

plotHistogramWithOptions(UpStateLen, Options);

Options.BinNum = 20;
Options.EdgeColor = 'g';
Options.FaceColor = 'g';
Options.FaceAlpha = 0.5;

plotHistogramWithOptions(UDcycleLen, Options);

YLim = get(gca, 'YLim');
clear hlgn
hlgn(1) = plot([0 0]+mean(DownStateLen), YLim, 'b--', 'LineWidth', 1);
text(mean(DownStateLen), YLim(2)*0.80, ...
   [' ' num2str(mean(DownStateLen), 3) 's ' ...
   '(sd: ' num2str(std(DownStateLen), 3) 's, '...
   'n=' num2str(length(DownStateLen)) ')']);
hlgn(2) = plot([0 0]+mean(UpStateLen), YLim, 'r--', 'LineWidth', 1);
text(mean(UpStateLen), YLim(2)*0.60, ...
   [' ' num2str(mean(UpStateLen), 3) 's ' ...
   '(sd: ' num2str(std(UpStateLen), 3) 's, '...
   'n=' num2str(length(UpStateLen)) ')']);
hlgn(3) = plot([0 0]+mean(UDcycleLen), YLim, 'g--', 'LineWidth', 1);
text(mean(UDcycleLen), YLim(2)*0.40, ...
   [' ' num2str(mean(UDcycleLen), 3) 's ' ...
   '(sd: ' num2str(std(UDcycleLen), 3) 's, '...
   'n=' num2str(length(UDcycleLen)) ')']);

set(gca, 'Layer', 'top', 'Box', 'off', 'TickDir', 'out');

legend(hlgn, {'Down', 'Up', 'UD'});
xlabel('Duration (s)')
set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.5 3.5 5 4]);
% print('-deps2c', 'UpDownDuration.hist.eps');
print('-djpeg', '-r300', 'UpDownCycleDuration.hist.jpg');
