%function WFs = extractEventCenteredWaveForms(t, Y, Triggers, ...
%                                             WFTimeRange, BLTimeRange)
                                         
function WFs = extractEventCenteredWaveForms(t, Y, Triggers, ...
    WFTimeRange, CutInterval, CutSamples, BLTimeRange)                                        
%
%
%   Copyright 2008 Maurizio Mattia @ Ist. Super. Sanit�, Rome - Italy
%   Version: 1.0 - Jan. 21, 2008
%


SamplingFreq = 1 / (t(2) - t(1));
ndxRel = round(WFTimeRange(1)*SamplingFreq):round(WFTimeRange(2)*SamplingFreq);
%ndxSurge = round(-DELTA_T/2*SAMPLING_FREQ):round(DELTA_T/2*SAMPLING_FREQ);
if exist('BLTimeRange')
   ndxBaseline = round(BLTimeRange(1)*SamplingFreq):round(BLTimeRange(2)*SamplingFreq);
else
   ndxBaseline = [];
end
tWFs = ndxRel / SamplingFreq;

% ************************************************************************
% THR = NaN;
% if ~isempty(CutInterval)
%    if numel(tList) == 1 % the only element in jList is the threshold 
%             THR = tList;
%         elseif diff(jList)==1
%             THR = tList(end);
%         else display('WATCH OUT - this case has to be checked with care')    
%    end
% end
% ************************************************************************

for n = 1:length(Triggers)
    
    ndxOffset = round((Triggers(n) - t(1)) * SamplingFreq);
    
    if ~isempty(CutInterval)
        for idx = 1:1:size(CutInterval,2)
            if Triggers(n)>CutInterval(2,idx)
             ndxOffset = round((Triggers(n) - t(1)) * SamplingFreq)-CutSamples(idx);
            end
        end
    end
    
%     if ~isnan(THR)
%         if Triggers(n)>=THR
%             ndxOffset = round((Triggers(n) - t(1)) * SamplingFreq)-CutSamples;
%         end 
%     end

   ndx = ndxRel + ndxOffset;
   ndxMask = find(ndx>0 & ndx<=length(t));
    
% disp([t(ndxOffset) Triggers(n)])
% disp(t(ndx(end))-t(ndx(1)))
% disp([t(ndx(1)) t(ndx(end))])
   
   WFs(n).time = tWFs(ndxMask);
%   WFs(n).value = csY(ndx(ndxMask)) - mean(csY(ndxSurge + ndxOffset));
%   WFs(n).value = csY(ndx(ndxMask)) - mean(csY(ndxBaseline + ndxOffset));
   if isempty(ndxBaseline)
      WFs(n).value = Y(ndx(ndxMask))';
   else
      WFs(n).value = Y(ndx(ndxMask))' - mean(Y(ndxBaseline + ndxOffset));
   end
end
