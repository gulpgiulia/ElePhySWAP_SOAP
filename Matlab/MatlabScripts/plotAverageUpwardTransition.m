%function [t, MeanY, StdY] = plotAverageUpwardTransition(MUA, UD, Triggers, tRange, YRange)
function [t, MeanY, StdY] = plotAverageUpwardTransition(MUA, UD, Triggers, tRange, ...
    CutInterval, CutSamples, YRange)

%
%  [t, MeanY, StdY] = plotAverageUpwardTransition(MUA, UD, Triggers, tRange[, YRange])
%
%
%   Copyright 2014 Maurizio Mattia @ Ist. Super. Sanit�, Rome - Italy
%   Version: 1.0 - Jan. 10, 2014
%

TIME_STEP = 0.020;       % Time window within to compute the moving average
                         % and the selection of good samples.
MINIMUM_SAMPLES = 10;     % Minimum number of waveforms to compute a reliable
                         % average and st.dev. estimate.
DIST_FROM_TRANS = 0.100; % Time from the closest transitions acceptable to 
                         % say you are far enough from a transition.

GRAY_CLR = [0.9 0.9 0.9];

% if ~isempty(jList)
%     UD_WFs = extractEventCenteredWaveForms(UD.time, UD.value, Triggers, tRange, CutSamples, jList, tList);
%     MUA_WFs = extractEventCenteredWaveForms(MUA.time, MUA.value, Triggers, tRange, CutSamples, jList, tList);
% else
%     UD_WFs = extractEventCenteredWaveForms(UD.time, UD.value, Triggers, tRange);
%     MUA_WFs = extractEventCenteredWaveForms(MUA.time, MUA.value, Triggers, tRange);
% end

UD_WFs = extractEventCenteredWaveForms(UD.time, UD.value, Triggers, tRange, CutInterval, CutSamples);
MUA_WFs = extractEventCenteredWaveForms(MUA.time, MUA.value, Triggers, tRange, CutInterval, CutSamples);

ndx0 = find(MUA_WFs(1).time >= 0, 1, 'first');
dt = diff(MUA_WFs(1).time(1:2));

TimeSteps = round(TIME_STEP / dt); % (number of samples)
Samples = zeros(1, length(MUA_WFs(1).time));
t = zeros(1, length(MUA_WFs(1).time));
MeanY = zeros(1, length(MUA_WFs(1).time));
StdY = zeros(1, length(MUA_WFs(1).time));

%
%  MUAs average after upward transition...
%
UpStateEnd = zeros(1, length(UD_WFs));
for k = 1:length(UD_WFs) % (loop on the upward transitions)
   ndx = find(UD_WFs(k).time > 0 & UD_WFs(k).value' == 0, 1, 'first');
   
   if ndx<=0; disp([k ndx]); end;
   
   
   if length(ndx) == 1
      UpStateEnd(k) = UD_WFs(k).time(ndx) - DIST_FROM_TRANS;
   else
      UpStateEnd(k) = UD_WFs(k).time(end) - DIST_FROM_TRANS;
   end
end

for n = ndx0:TimeSteps:length(MUA_WFs(1).time)
   if n+TimeSteps-1 > length(MUA_WFs(1).time)
      ndx = n:length(MUA_WFs(1).time);
   else
      ndx = n:(n+TimeSteps-1); %  Time window within to compute the moving average
   end
   ndxSamples = find(UpStateEnd >= MUA_WFs(1).time(ndx(end))); % find 'long' transition, lasting more than time(ndx(end))
   if length(ndxSamples) > MINIMUM_SAMPLES
      t(ndx) = MUA_WFs(1).time(ndx);
      Samples(ndx) = length(ndxSamples);
      for k = ndxSamples
         MeanY(ndx) = MeanY(ndx) + MUA_WFs(k).value(ndx)';
         StdY(ndx) = StdY(ndx) + MUA_WFs(k).value(ndx)'.^2;
      end
   end
end

%
%  MUAs average before upward transition...
%
DownStateBegin = zeros(1, length(UD_WFs));
for k = 1:length(UD_WFs)
   ndx = find(fliplr(UD_WFs(k).time) < 0 & fliplr(UD_WFs(k).value') == 1, 1, 'first');
   if length(ndx) == 1
      ndx = length(UD_WFs(k).time) - ndx + 1;
      DownStateBegin(k) = UD_WFs(k).time(ndx) + DIST_FROM_TRANS;
   else
      DownStateBegin(k) = UD_WFs(k).time(1) + DIST_FROM_TRANS;
   end
end

for n = ndx0-1:-TimeSteps:1
   if n-TimeSteps+1 < 1
      ndx = 1:n;
   else
      ndx = (n-TimeSteps+1):n;
   end
   ndxSamples = find(DownStateBegin <= MUA_WFs(1).time(ndx(1)));
   if length(ndxSamples) > MINIMUM_SAMPLES
      t(ndx) = MUA_WFs(1).time(ndx);
      Samples(ndx) = length(ndxSamples);
      for k = ndxSamples
         MeanY(ndx) = MeanY(ndx) + MUA_WFs(k).value(ndx)';
         StdY(ndx) = StdY(ndx) + MUA_WFs(k).value(ndx)'.^2;
      end
   end
end

%
%  Computes the average MUA and its standard error...
%
ndx = find(Samples > 0);
t = t(ndx);
MeanY = MeanY(ndx) ./ Samples(ndx);
StdY = sqrt((StdY(ndx) ./ Samples(ndx) - MeanY.^2) ./ Samples(ndx));

%
% Plots the average downward transition dynamics computed
%

figure
hold on

XP = [t fliplr(t)];
YP = [MeanY+StdY fliplr(MeanY-StdY)];
patch(XP, YP, GRAY_CLR, 'EdgeColor', 'none');

plot(t, MeanY, 'b', 'LineWidth', 1);

if exist('YRange') == 1
   set(gca, 'YLim', YRange);
end

plot([0 0], get(gca, 'YLim'), 'k--');

set(gca, 'TickDir', 'out', 'Box', 'on', 'Layer', 'top', 'XLim', tRange);

xlabel('Time (s)');
ylabel('log(MUA)');

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1 3.5 6.7 4]);
print('-deps2c', 'LogMUAUpTrans.Average.eps')
