function [ModeParams, MIN_2ND_PEAK_AREA,N,X,F1] = plotBimodalHistogram(Values, fid)
%function [ModeParams, MIN_2ND_PEAK_AREA,N,X,F1,F2] = plotBimodalHistogram(Values, fid)
%function [ModeParams, MIN_2ND_PEAK_AREA,N,X,F1,F2] = plotBimodalHistogram(Values, ValRange, BinNum, fid)

%  ModeParams = plotBimodalHistogram(Values[, ValRange[, BinNum]])
%
%
%   Copyright 2009 Maurizio Mattia @ Ist. Super. Sanit�, Rome - Italy
%   Version: 1.0 - Sep. 22, 2009
%   Version: 2.0 - May. 22, 2018 --> compute the "TAILS"
%   Version: 3.0 - June 2018 --> Skewness; fit of the "half" Gaussian

%*************************************************************************
GRAY_CLR = repmat(0.9, 1, 3);
BIN_NUM = 100;
PEAK_THRESHOLD = 3/5; % The level from the peak which is at about 1 st.dev from it.
TAIL_TO_NEGLECT = 0.005;
MIN_2ND_PEAK_AREA = 10.0;
%*************************************************************************
if ~exist('fid')
    fid=1;
end
%*************************************************************************

SampleSize = length(Values); %(for the normalization of the histograms)

if exist('BinNum')
   if BinNum>0
      BIN_NUM = BinNum;
   end
end

if exist('ValRange') ~= 1
   ValRange = [min(Values) max(Values)];
else
   if isempty(ValRange)
      ValRange = [min(Values) max(Values)];
   end
end
X = linspace(ValRange(1), ValRange(2), BIN_NUM);

Values = Values(find(Values >= ValRange(1) & Values <= ValRange(2)));

%*************************************************************************
% N  --> histrogram of the logMUA distribution
% >>> FIT the firstGaussian
% N2 --> N - firstGaussian; 
% sum(N2) = area of the firstTail [N.B. sum(N2>0)]
% >>> if firstTail > 10% >>> FIT the secondGaussian
% N3 --> N - firstGaussian - secondGaussian
% sum(N3) = area of the secondTail
% ...

% 'first' peak = largest peak (it can be left pr right...)

% FIT with a Gaussian
% >>> beta0 = initial values for the coeffients of the fit
% If the fit converges -->
% >>> beta1, if the peak is on the LEFT (MUA_DOWN)
% >>> beta2, if the peak is on the RIGHT (MUA_UP)
%*************************************************************************

% Computes histogram...
N = hist(Values, X) / SampleSize * 100; 
% N is a vector with length(X) elements, equal to the number of bins in the histogram.
% Each element of N stores the number of entries for each bin (normalised and in percentage)

% Smooth the histogram in order to find maxima...
knt = augknt(linspace(ValRange(1), ValRange(2), round(BIN_NUM/4)), 4);
spN = spap2(knt, 4, X, N); % Least Square Spline approx.
spN = spap2(newknt(spN), 4, X, N); % ... a better approximation
spN = spap2(newknt(spN), 4, X, N); % A second refinement of the fit...

% removes tails of the distribution for low and high values...
Y = cumsum(N);
Y = Y / Y(end); %[0 1]
a = X(find(Y > TAIL_TO_NEGLECT, 1, 'first'));
b = X(find(Y <= 1.0 - TAIL_TO_NEGLECT, 1, 'last'));

% finds extrema (minima and maxima) in the bulk of the distribution...
% (bulk = tails excluded)
XMaxMin = mean(fnzeros(fnder(spN), [a b]));
YMaxMin = fnval(spN, XMaxMin);

% Criterium 1 for peak selection (original):
[YMax, YMaxNdx] = max(YMaxMin);

% Criterium 2 for peak selection (18th Sep. 2009):
% ConcMaxMin = fnval(fnder(fnder(spN)),XMaxMin);
% [ConcMax, YMaxNdx] = min(ConcMaxMin);

% -------------------------------------------------------------------------
%    FIT the firstGaussian
% ---------------------------
% Gaussian fit on the first peak (N.B. 'first' = largest)
XFirstMax = XMaxMin(YMaxNdx);
YFirstMax = YMaxMin(YMaxNdx);
first=0;
% ... Is the first peak on left or on the right?
if XFirstMax-a < b-XFirstMax  % >>> firstPeak is on the LEFT (DOWN state)
   FirstLeft = 1;
   ndx = find(X >= XFirstMax, 1 , 'first');
   [val, pos] = max(N(ndx-1:ndx+1));
   ndx = ndx + pos - 2;
   ndxL = find(N(1:ndx)>val*PEAK_THRESHOLD, 1, 'first');
%    ndxL = find(N(1:ndx)>YFirstMax*PEAK_THRESHOLD, 1, 'first');
   
% --- beta0 = initial values for the coefficients of the fit 
   beta0(1) = X(ndx);
   beta0(2) = (ndx - ndxL)*(X(2) - X(1)); %(X is linspace) % first guess for SIGMA
   %    if beta0(2) == 0
   %       beta0(2) = X(2) - X(1);
   %    end
   beta0(3) = val;
   beta0(3) = YFirstMax; 
   %>>> Original
   %%%ndx = ndx + (ndx - ndxL - 1);
   %>>> 18th Sep. 2009
   %%% ndx = ndx + (ndx - ndxL - 1) + 1;
   %>>> [GDB] June 2018 --- HALF Gaussian
   ndx=ndx+1;

% --- GAUSSIAN FIT ---
   beta1 = nlinfit(X(1:ndx), N(1:ndx), @mygauss, beta0);
   ModeParams.Mu1 = beta1(1);
   ModeParams.Sigma1 = beta1(2);
   ModeParams.Ampl1 = beta1(3);
   first=1;
   
else  % >>> firstPeak is on the RIGHT (UP state)
   FirstLeft = 0;
   ndx = find(X <= XFirstMax, 1 , 'last');
   [val, pos] = max(N(ndx-1:ndx+1));
   ndx = ndx + pos - 2;
   ndxL = find(N(ndx:end)<val*PEAK_THRESHOLD, 1, 'first');
%    ndxL = find(N(ndx:end)<YFirstMax*PEAK_THRESHOLD, 1, 'first');

% --- beta0 = initial values for the coefficients of the fit 
   beta0(1) = X(ndx);
   beta0(2) = ndxL*(X(2) - X(1));
   %    if beta0(2) == 0
   %       beta0(2) = X(2) - X(1);
   %    end
   beta0(3) = val;
   beta0(3) = YFirstMax;
   %>>> Original
   %%%ndx = ndx - ndxL;
   %>>> 18th Sep. 2009
   %%%ndx = ndx - ndxL - 1;
   %>>> [GDB] June 2018 --- HALF Gaussian
   ndx=ndx-1;

% --- GAUSSIAN FIT ---
   beta2 = nlinfit(X(ndx:end), N(ndx:end), @mygauss, beta0);

   ModeParams.Mu2 = beta2(1);
   ModeParams.Sigma2 = beta2(2);
   ModeParams.Ampl2 = beta2(3);
   first=1;
end

if(first)
    fprintf(fid,'...first peak fitted\n');
    %[GDB] here, add a check if the fit procedure is not converging
end
if(FirstLeft == 0); fprintf(fid,'[WARNING. First peak is on the RIGHT]\n');end
ModeParams.FirstLeft = FirstLeft;

% -------------------------------------------------------------------------
%    FIT the secondGaussian
% ---------------------------
% Gaussian fit on the second peak (N.B. 'second' = second largest)
% N.B. if SecondPeakArea<MIN_2ND_PEAK_AREA --> the 2nd peak is NOT FITTED
second=0;
if FirstLeft == 0  % firstPeak is on the RIGHT
   N2 = N - mygauss(beta2, X);
   NN2=N2; %(back-up)
   
   N2(find(X>beta2(1)))=0;
   N2(find(N2<0))=0; % set at zero negative values, to interpret N2 as a distribution
%  N.B. estimate of the secondPeak: the sum is only on values N2>0
    
   SecondPeakArea = sum(N2(find(N2>0)));  
   fprintf(fid,['2nd peak area = ' num2str(SecondPeakArea) '%%\n']);
% %    if SecondPeakArea >= MIN_2ND_PEAK_AREA
% % % ------------------------------------------------------------------------
% %       % Criterium 0 for peak selection:
% %       %    [val, ndx] = max(N2);
% % 
% %       % Criterium 1 for peak selection (original):
% %       ndx = find(X>a);
% %       maxndx = find(diff(N2(ndx))<0, 1, 'first');
% % % [GDB] 'first' can cause "Warning: iteration limit", because a local
% % % (relative) maximum is selected as starting value for the fit
% % % --- to be corrected
% %       ndx = ndx(maxndx);
% %       val = N2(ndx);
% % 
% %       % Criterium 2 for peak selection (18th Sep. 2009):
% % %       ndx = find(X>=N2*X'/sum(N2), 1, 'first');
% % %       val = N2(ndx);
% % 
% %       ndxL = find(N2(1:ndx)>val*PEAK_THRESHOLD, 1, 'first');
% % 
% %       % Criterium 3 for peak selection (18th Sep. 2009):
% % %       spN2 = fnval(spN, X) - mygauss(beta2, X);
% % %       [val, ndx] = max(spN2);
% % %       
% % %       ndxL = find(spN2(1:ndx)>val*PEAK_THRESHOLD, 1, 'first');
% % % ------------------------------------------------------------------------
% % % --- beta0 = initial values for the coefficients of the fit 
% %       beta0(1) = X(ndx);
% %       beta0(2) = (ndx - ndxL)*(X(2) - X(1));
% %       beta0(3) = val;
% %       ndx = ndx + (ndx - ndxL - 1);
% % %       ndx = ndx + (ndx - ndxL - 1) + 1;
% %       ndx = min([ndx length(X)]);
% %       
% % % --- GAUSSIAN FIT ---
% %       beta1 = nlinfit(X(1:ndx), N2(1:ndx), @mygauss, beta0);
% %       ModeParams.Mu1 = beta1(1);
% %       ModeParams.Sigma1 = beta1(2);
% %       ModeParams.Ampl1 = beta1(3);
% %       second=1;
% %    end
else % firstPeak is on the LEFT >>> secondPeak is on the RIGHT (i.e. UP state)
   N2 = N - mygauss(beta1, X); % what is left after subtracting the first gaussian fit
   NN2=N2; %(back-up)

   N2(find(X<beta1(1)))=0;
   N2(find(N2<0))=0; % set at zero negative values, to interpret N2 as a distribution
%  N.B. estimate of the secondPeak: the sum is only on values N2>0
   
   SecondPeakArea = sum(N2(find(N2>0)));
   fprintf(fid,['2nd peak area = ' num2str(SecondPeakArea) '%%\n']);
% %    if SecondPeakArea >= MIN_2ND_PEAK_AREA
% % % ------------------------------------------------------------------------
% %       % Criterium 0 for peak selection:
% %       %    [val, ndx] = max(N2);
% % 
% %       % Criterium 1 for peak selection (original):
% %       ndx = fliplr(find(X<b));      
% % % [GDB] 'first' can cause "Warning: iteration limit", because a local
% % % (relative) maximum is selected as starting value for the fit
% %       %%% maxndx = find(diff(N2(ndx))<0, 1, 'first');
% %       maxndxALL = find(diff(N2(ndx))<0);
% %       [peak, idx]= max(N2(ndx(maxndxALL)));
% %       maxndx=maxndxALL(idx);
% %       
% %       ndx = ndx(maxndx);
% %       val = N2(ndx);
% % 
% %       % Criterium 2 for peak selection (18th Sep. 2009):
% % %       ndx = find(X>=N2*X'/sum(N2), 1, 'first');
% % %       val = N2(ndx);
% % 
% %       ndxL = find(N2(ndx:end)<val*PEAK_THRESHOLD, 1, 'first');
% % 
% %       % Criterium 3 for peak selection (18th Sep. 2009):
% % %       spN2 = fnval(spN, X) - mygauss(beta1, X);
% % %       [val, ndx] = max(spN2);
% % %       ndxL = find(spN2(ndx:end)<val*PEAK_THRESHOLD, 1, 'first');
% % % ------------------------------------------------------------------------
% % % --- beta0 = initial values for the coefficients of the fit 
% %       beta0(1) = X(ndx);
% %       beta0(2) = ndxL*(X(2) - X(1)); % guess value for sigma
% %       beta0(3) = val;
% %       ndx = ndx - ndxL - 1; % ...maybe this is to define an interval symmetrically around the peak
% % 
% %       % but ndx MUST be > 0!
% %       if ndx < 1; ndx = 1; end
% % 
% % % --- GAUSSIAN FIT ---
% %       beta2 = nlinfit(X(ndx:end), N2(ndx:end), @mygauss, beta0);
% %       ModeParams.Mu2 = beta2(1);
% %       ModeParams.Sigma2 = beta2(2);
% %       ModeParams.Ampl2 = beta2(3);
% %       second=1;
% %    end
end

if(second);fprintf(fid,'...second peak fitted\n');end

% Add SecodnPeakArea in the list of ModeParams
ModeParams.SecondPeakArea = SecondPeakArea;

% -------------------------------------------------------------------------
%    check the TAIL
% ---------------------------
% % if exist('beta1') && exist('beta2')
% %     N3 = N - mygauss(beta1, X) - mygauss(beta2, X);
% %     tail = sum(N3(find(N3>0)));
% %     ModeParams.Tail = tail;
% %     if tail >= MIN_2ND_PEAK_AREA
% %         disp(['WARNING. The tail of the distribution is ' num2str(SecondPeakArea) '%']);
% %         disp('...check the hypothesis of bimodality...');
% %     end
% % end


% -------------------------------------------------------------------------
% INDEX of ASYMMETRY (SKEWNESS)
% -----------------------------

%%%[ fi=N2/sum(N2), sum(N2/sum(N2))=1 ]

%%% ----
%%% TOT
%%% ----
MU=sum(X.*(N/sum(N)));
SIG=sqrt(sum((X-MU).^2.*N/sum(N)));
%...or sqrt(sum((X.*X).*(N2/sum(N2)))-MU*MU)
MED=X(find(cumsum(N/sum(N))>=0.5,1));
Q1=X(find(cumsum(N/sum(N))>=0.25,1));
Q3=X(find(cumsum(N/sum(N))>=0.75,1));
[a,b]=max(N/sum(N));
MOD=X(b);

alpha(1)=MU-MED;
alpha(2)=MU-MOD;
alpha(3)=Q3+Q1-2*MED;
gamma=sum((X-MU).^3.*N/sum(N))/SIG^3; % Fisher

ModeParams.SkewnessTot.MU=MU;
ModeParams.SkewnessTot.SIG=SIG;
ModeParams.SkewnessTot.MED=MED;
ModeParams.SkewnessTot.Q1=Q1;
ModeParams.SkewnessTot.Q3=Q3;
ModeParams.SkewnessTot.MOD=MOD;
ModeParams.SkewnessTot.alpha=alpha;
ModeParams.SkewnessTot.gamma=gamma;

%%% ----
%%% TAIL
%%% ----
MU=sum(X.*(N2/sum(N2)));
SIG=sqrt(sum((X-MU).^2.*N2/sum(N2)));
%...or sqrt(sum((X.*X).*(N2/sum(N2)))-MU*MU)
MED=X(find(cumsum(N2/sum(N2))>=0.5,1));
Q1=X(find(cumsum(N2/sum(N2))>=0.25,1));
Q3=X(find(cumsum(N2/sum(N2))>=0.75,1));
[a,b]=max(N2/sum(N2));
MOD=X(b);

alpha(1)=MU-MED;
alpha(2)=MU-MOD;
alpha(3)=Q3+Q1-2*MED;
gamma=sum((X-MU).^3.*N2/sum(N2))/SIG^3; % Fisher

ModeParams.SkewnessTail.MU=MU;
ModeParams.SkewnessTail.SIG=SIG;
ModeParams.SkewnessTail.MED=MED;
ModeParams.SkewnessTail.Q1=Q1;
ModeParams.SkewnessTail.Q3=Q3;
ModeParams.SkewnessTail.MOD=MOD;
ModeParams.SkewnessTail.alpha=alpha;
ModeParams.SkewnessTail.gamma=gamma;


% -------------------------------------------------------------------------
%    plot the HISTOGRAM
% ---------------------------
F1=figure; 
hold on;

[XX, YY] = stairs(X - (X(2) - X(1))/2, N);
XX = [XX(1) XX' XX(end)];
YY = [0 YY' 0];

% --- Histogram filled in grey (distribution of logMUA)
patch(XX, YY, GRAY_CLR, 'FaceColor', GRAY_CLR, 'EdgeColor', GRAY_CLR);
plot(XX, YY, 'k', 'LineWidth', 1);

FitHist = zeros(size(N)); % FitHist = total fit
% --- Fit of the LEFT peak (blue)
if exist('beta1')
   plot(X, mygauss(beta1, X), 'b', 'LineWidth', 1);
   FitHist = mygauss(beta1, X);
   plot(X, mygauss(beta1, X), 'b', 'LineWidth', 1);
else
% --- Fit of the tail (blue)
    plot(X,N2,'b:')
end
% --- Fit of the RIGHT peak (red)
if exist('beta2')
   plot(X, mygauss(beta2, X), 'r', 'LineWidth', 1);
   FitHist = FitHist + mygauss(beta2, X);
   plot(X, mygauss(beta2, X), 'r', 'LineWidth', 1);
else
% --- Fit of the tail (red)
    plot(X,N2,'r:')
end

% % % --- TOTAL FIT (green)
% % plot(X, FitHist, 'g', 'LineWidth', 1); 
% % if exist('N3')
% % % --- Fit of the tail (black)
% %     plot(X,N3,'k:')
% % end

YLim = get(gca, 'YLim');

% --- Vertical dashed lines to mark the position of Mu1, Mu2
if exist('beta1')
   % --- Fit of the LEFT peak (blue)
   plot(ModeParams.Mu1 + [0 0], [0 YLim(2)], 'b--');
   text(ModeParams.Mu1, ModeParams.Ampl1+0.5, ...
        [' ' num2str(ModeParams.Mu1, 3) ...
         ' (s.d.: ' num2str(ModeParams.Sigma1, 3) ')']);
end
if exist('beta2')
   % --- Fit of the RIGHT peak (red)
   plot(ModeParams.Mu2 + [0 0], [0 YLim(2)], 'r--');
   text(ModeParams.Mu2, ModeParams.Ampl2+0.5, ...
        [' ' num2str(ModeParams.Mu2, 3) ...
         ' (s.d.: ' num2str(ModeParams.Sigma2, 3) ')']);
end

set(gca, 'Layer', 'top', 'Box', 'off', 'TickDir', 'out');
set(gca, 'XLim', ValRange);
xlabel('log(MUA)');

xlabel('Value');
ylabel('Samples (%)');

%%%%%%%%%%%%%%%%%%%
% plot of the Tail
%%%%%%%%%%%%%%%%%%%
axes('Position',[.55 .55 .35 .35])
box on
hold on;
if exist('beta1')
    plot(X,N2,'r') % RIGHT tail is in RED
elseif exist('beta2')
    plot(X,N2,'b') % LEFT tail is in BLUE
end
    
plot(ModeParams.SkewnessTail.MU + [0 0], [0 max(N2)], 'k--');
plot(ModeParams.SkewnessTail.MOD + [0 0], [0 max(N2)], 'm--');
plot(ModeParams.SkewnessTail.MED + [0 0], [0 max(N2)], 'c--');
set(gca, 'YLim', [0 max(N2)]);
set(gca, 'XLim', ValRange);
xlabel('Value');
ylabel('Samples (%)');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scale-free dynamics at large firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % F2=figure;
% % % if(FirstLeft ~= 0)
% % %     N3=N2;
% % %     N3(find(X<ModeParams.SkewnessTail.MU))=0;
% % %     ndx = find(X>=ModeParams.SkewnessTail.MU,1);
% % %     %ndx3 = find(X>=ModeParams.SkewnessTail.Q3,1);
% % %     beta0=[-3,N2(ndx)/(ModeParams.SkewnessTail.MU)^(-3)];
% % %     %beta0=[-2,N3(ndx),ModeParams.SkewnessTail.MU];
% % %     beta3 = nlinfit(X(ndx:end), N2(ndx:end), @powerlaw, beta0);
% % %     FIT = powerlaw(beta3,X(ndx:end));
% % % 
% % %     subplot(2,1,1) 
% % %     plot(X,N2)
% % %     hold on
% % %     plot(X,N3)
% % %     plot(X(ndx:end),FIT)
% % %     %text(ModeParams.Mu1-2*ModeParams.Sigma1,0.8*max(N2),['\gamma = ' num2str(beta3(1))],'FontSize',14);
% % % 
% % %     subplot(2,1,2) 
% % %     plot(X,log(N2))
% % %     hold on
% % %     plot(X,log(N3))
% % %     plot(X(ndx:end),log(FIT))
% % %     text(ModeParams.SkewnessTail.MU*2,log(N3(ndx)),['\gamma = ' num2str(beta3(1))],'FontSize',14);
% % % 
% % %     ModeParams.BetaTail=beta3;
% % % end

% -------------------------------------------------------------------------
%
% Work out the error of the fit...
%
HistError = N/100.; % number of entries for each bin (normalised)
HistError = sqrt((HistError .* (1. - HistError)) / length(Values)); %takes into account of the number of elements for each bin
ndx = find(abs(N-FitHist)/100. > HistError); %difference between the distribution (N) and its fit
ModeParams.FitError = sum(abs(N(ndx)-FitHist(ndx))/100. - HistError(ndx));
