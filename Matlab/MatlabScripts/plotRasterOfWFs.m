function ZRange = plotRasterOfWFs(WFs, ZRange, ZLabel)
%
%  ZRange = plotRasterOfWFs(WFs[, ZRange[, ZLabel]])
%
%
%   Copyright 2014 Maurizio Mattia @ Ist. Super. Sanit�, Rome - Italy
%   Version: 1.0 - Jan. 10, 2014
%

%
% Evaluate the raster size...
%
MinSizeT = length(WFs(1).time);
for k = 2:length(WFs)
   if MinSizeT > length(WFs(k).time);
      MinSizeT = length(WFs(k).time);
      %disp(k)
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

%
% Plots the raster...
%
figure;
imagesc(t, 1:size(Ys,1), Ys);

if exist('ZRange') == 0
   ZRange = [min(min(Ys)) max(max(Ys))];
else
   if isempty(ZRange)
      ZRange = [min(min(Ys)) max(max(Ys))];
   end
end
set(gca, 'CLim', ZRange);

h = colorbar;
set(gca, 'YDir', 'normal');
if exist('ZLabel') == 1
   set(get(h,'YLabel'),'String',ZLabel);
else
   set(get(h,'YLabel'),'String','Y (a.u.)');
end
set(h, 'YLim', ZRange);
xlabel('Time (s)');
ylabel('Trial');
set(gca,'Layer','top','TickDir','out','Box','on');

%
% Saves the plot...
%
set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1 3 6 5]);
print('-deps2c', 'RasterOfWFs.eps');
