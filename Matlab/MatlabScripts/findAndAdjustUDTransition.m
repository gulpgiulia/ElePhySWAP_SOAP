function [Trans, UD] = findAndAdjustUDTransition(MUA, UD, UDthreshold)
%
%  [Trans, UD] = findAndAdjustUDTransition(MUA, UD, UDthreshold)
%
%
%   Copyright 2008 Maurizio Mattia @ Ist. Super. Sanit�, Rome - Italy
%   Version: 1.0 - Sep. 23, 2008
%

FIT_WINDOW = [-50 50]/1000;

ndxFit = round(FIT_WINDOW/UD.dt);
ndxFit = ndxFit(1):ndxFit(2);

%
% Find upward and downward transitions...
Trans.ndx = find(diff(UD.value) ~= 0)+1;
Trans.val = UD.value(Trans.ndx);          % 1 for upward and 0 for downward transitions.
Trans.time = UD.time(Trans.ndx);

%
% Remove the transitions too close to the begin and the end of sampled
% period...
ndxTrans = find(Trans.time - UD.time(1) > -FIT_WINDOW(1) & ...
                UD.time(end) - Trans.time > FIT_WINDOW(2));
Trans.ndx = Trans.ndx(ndxTrans);
Trans.val = Trans.val(ndxTrans);
Trans.time = Trans.time(ndxTrans);

if numel(Trans.ndx) > 0

   % figure
   % hold on
   % DUshift = [];
   % UDshift = [];

   for k = 1:numel(Trans.ndx)
      p = polyfit(ndxFit, ...
         MUA.value(ndxFit+Trans.ndx(k)), 3); % degree=3 (cubic)
      p(end) = p(end) - UDthreshold; % p(end) is the constant term of the polynomial. 
      % "p(end) - UDthreshold" is a translation of the constant term centered at the Threshold
      %    r = roots(p)*MUA.dt;
      r = roots(p); %  Find polynomial roots.
      ndx = find(imag(r)==0); % select only real roots (1 or 3 for a 3-degree polynomial)
      
      % Upward transition...
      if Trans.val(k) == 1
         %       plot(ndxFit*MUA.dt, polyval(p,ndxFit),'r')
         [val, pos] = min(abs(r(ndx)));
         Shift = round(r(ndx(pos)));
         if Shift > 0
            UD.value(Trans.ndx(k):Trans.ndx(k)+Shift) = 0;
         else
            UD.value(Trans.ndx(k)+Shift:Trans.ndx(k)) = 1;
         end
         %       DUshift = [DUshift r(ndx(pos))];
         
      % Downward transition...
      else
         %       plot(ndxFit*MUA.dt, polyval(p,ndxFit),'b')
         [val, pos] = min(abs(r(ndx)));
         Shift = round(r(ndx(pos)));
         if Shift > 0
            UD.value(Trans.ndx(k):Trans.ndx(k)+Shift) = 1;
         else
            UD.value(Trans.ndx(k)+Shift:Trans.ndx(k)) = 0;
         end
         %       UDshift = [UDshift r(ndx(pos))];
      end
   end
   % figure; hist(DUshift,20)
   % figure; hist(UDshift,20)
   
   %
   % Find the updated upward and downward transitions...
   Trans.ndx = find(diff(UD.value) ~= 0)+1;
   Trans.val = UD.value(Trans.ndx);
   Trans.time = UD.time(Trans.ndx);
   
   %
   % Remove the transitions too close to the begin and the end of sampled
   % period...
   ndxTrans = find(Trans.time - UD.time(1) > -FIT_WINDOW(1) & ...
      UD.time(end) - Trans.time > FIT_WINDOW(2));
   Trans.ndx = Trans.ndx(ndxTrans);
   Trans.val = Trans.val(ndxTrans);
   Trans.time = Trans.time(ndxTrans);

end