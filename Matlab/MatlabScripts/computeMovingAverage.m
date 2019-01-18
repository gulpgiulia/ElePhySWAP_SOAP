function [t, MA] = computeMovingAverage(WF, WindowSize, StepSize)
%
%   [t, MA] = computeMovingAverage(WF, WindowSize, StepSize)
%
%
%   Copyright 2008 Maurizio Mattia @ Ist. Super. Sanit�, Rome - Italy
%   Version: 1.0 - Mar. 7, 2008
%

   SampleNum = floor((length(WF.value) - WindowSize)/StepSize) + 1;
   
   dt = WF.time(2) - WF.time(1);
   t = (0:SampleNum-1) * dt * StepSize + mean(WF.time(1:WindowSize)); % (the centre of the time interval)
   
   data = zeros(WindowSize, SampleNum);
   for k = 1:WindowSize
      data(k,:) = WF.value(k:StepSize:end-WindowSize+k);
   end
   MA = mean(data);