function UD = computeBinaryUDState(MUA, Threshold, MinDuration)
%
%  UD = computeBinaryMUA(inMUA, Threshold, MinDuration)
%
%
%   Copyright 2008 Maurizio Mattia @ Ist. Super. Sanità, Rome - Italy
%   Version: 1.0 - Jul. 6, 2008
%

%
% Computes the binary raster...
UD.value = zeros(size(MUA.value));
UD.time = MUA.time;
UD.dt = diff(MUA.time(1:2));
UD.value(find(MUA.value >= Threshold)) = 1;

%
% Smooth raster removing small up and down states...
MinStateLen = floor(MinDuration / UD.dt);

updates = 1;
while updates > 0
   updates = 0;
   Trans.ndx = [1 find(diff(UD.value) ~= 0)+1 size(UD.value,2)+1];
   Trans.val = UD.value(Trans.ndx(1:end-1));
   Trans.val = [Trans.val abs(UD.value(end)-1)];
   StateLen = diff(Trans.ndx);
   ndx = find(StateLen < MinStateLen);
   for k = ndx
      change = 0;
      if k == 1
         change = StateLen(k+1) > StateLen(k);
      else
         if k == length(StateLen)
            change = StateLen(k-1) >= StateLen(k);
         else
            change = StateLen(k-1) >= StateLen(k) & ...
               StateLen(k+1) > StateLen(k);
         end
      end
      if change
         UD.value(Trans.ndx(k):Trans.ndx(k+1)-1) = Trans.val(k+1);
         updates = updates + 1;
      end
   end
   %       disp(['Updates: ' num2str(updates)]);
end
