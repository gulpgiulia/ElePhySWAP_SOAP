function [t, Ys] = computeRasterOfWFs(WFs)
%
%  [t, Ys] = computeRasterOfWFs(WFs)
%
%
%   Copyright 2010 Maurizio Mattia @ Ist. Super. Sanità, Rome - Italy
%   Version: 1.0 - Sep. 22, 2010
%

%
% Evaluate the raster size...
%
MinSizeT = length(WFs(1).time);
for k = 2:length(WFs)
   if MinSizeT > length(WFs(k).time);
      MinSizeT = length(WFs(k).time);
   end
end

%
% Builds the raster data support...
%
t = zeros(1, MinSizeT);
Ys = zeros(length(WFs), MinSizeT);
for k = 1:length(WFs)
   t = t + WFs(k).time(1:MinSizeT);
   Ys(k,:) = WFs(k).value(1:MinSizeT);
end
t = t / length(WFs);
