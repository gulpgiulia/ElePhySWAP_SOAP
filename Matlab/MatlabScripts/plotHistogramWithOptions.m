function PeakCenters = plotHistogramWithOptions(Values, Options)
%
%  PeakCenters = plotHistogramWithOptions(Values[, ValRange[, Threshold]])
%
%
%   Copyright 2009 Maurizio Mattia @ Ist. Super. Sanità, Rome - Italy
%   Version: 1.0 - Apr. 25, 2009
%


%
% Adjust and control parameters...
%
if exist('Options','var') ~= 1
   Options.ValRange = [min(Values) max(Values)];
end
if isfield(Options, 'ValRange') == 0
   Options.ValRange = [min(Values) max(Values)];
end
if isfield(Options, 'BinNum') == 0
   Options.BinNum = 10;
end
if isfield(Options, 'EdgeColor') == 0
   Options.EdgeColor = 'k';
end
if isfield(Options, 'FaceColor') == 0
   Options.FaceColor = repmat(0.9, 1, 3);
end
if isfield(Options, 'FaceAlpha') == 0
   Options.FaceAlpha = 1;
end

%
%  Computes the histogram...
%
SampleSize = length(Values);
X = linspace(Options.ValRange(1), Options.ValRange(2), Options.BinNum);
Values = Values(Values >= Options.ValRange(1) & Values <= Options.ValRange(2));

%
%  Plots the histogram...
%
N = hist(Values, X) / SampleSize * 100;
dX = X(2) - X(1);
[XX, YY] = stairs(X - dX/2, N);
XX = [XX(1) XX' XX(end)+dX XX(end)+dX];
YY = [0 YY' YY(end) 0];

patch(XX, YY, Options.EdgeColor, 'FaceColor', Options.FaceColor, ...
   'EdgeColor', Options.EdgeColor, 'FaceAlpha', Options.FaceAlpha);
hold on;
hlgn = plot(XX, YY, 'Color', Options.EdgeColor, 'LineWidth', 1);

set(gca, 'Layer', 'top', 'Box', 'off', 'TickDir', 'out');
set(gca, 'XLim', Options.ValRange);

xlabel('Value (a.u)');
ylabel('Samples (%)');

%
% Shows peaks if a threshold is specified...
%
if isfield(Options, 'Threshold') == 1
   YLim = get(gca, 'YLim');
   ndx = find(N > Options.Threshold);
   
   ndxEnd = find(diff(ndx) > 2);
   if isempty(ndxEnd)
      ndxEnd = length(ndx);
      ndxStart = 1;
   else
      ndxStart = [1 ndxEnd+1];
      ndxEnd = [ndxEnd length(ndx)];
   end

   for k = 1:length(ndxStart)
      ndxPeak = ndx(ndxStart(k):ndxEnd(k));
      PeakCenters(k) = N(ndxPeak) * X(ndxPeak)'/sum(N(ndxPeak));
      plot(PeakCenters(k) + [0 0], [0 YLim(2)], 'r--');
      text(PeakCenters(k), Options.Threshold/2, [' ' num2str(PeakCenters(k), 3)]);
   end
   plot(Options.ValRange, Options.Threshold + [0 0], 'b:');
else
   PeakCenters = hlgn;
end
