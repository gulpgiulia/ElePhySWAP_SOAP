function [t, MUA] = computeSpectralEstimateOfMUA(t, LFP, MUAFreqBand, PlotPowerSpectrum, PyyBaseline)
%
%  [t, MUA] = computeSpectralEstimateOfMUA(t, LFP, MUAFreqBand[, PlotPowerSpectrum[, PyyBaseline]])
%
%  e.g. MUAFreqBand = [200 1500];
%
%   Copyright 2014 Maurizio Mattia @ Ist. Super. Sanità, Rome - Italy
%   Version: 1.0 - Aug. 20, 2014
%

DETRENDING_ORDER = 2; % 0 - no detrending; 
                      % 1 - mean detrending; 
                      % 2 - mean and slope detrending (linear): this should
                      % be the default.
                      % 3 - quadratic detrending.
                      % 4 - cubic detrending.

dt = t(2) - t(1);
fs = 1 / dt; % Sampling frequency of the input signal.

% FFTWindowSize = floor(fs / MUAFreqBand(1));
FFTWindowSize = round(fs / MUAFreqBand(1));
sn = floor(length(LFP)/FFTWindowSize); % Samples number of the output. 

X = reshape(LFP(1:FFTWindowSize*sn), FFTWindowSize, sn);

% Polinomial detrending...
if DETRENDING_ORDER > 0 % 0-th order detrending...
   X = X - repmat(mean(X), FFTWindowSize, 1); 
end
if DETRENDING_ORDER > 1 % 1-st order detrending...
   Xx = repmat(linspace(-FFTWindowSize/2,FFTWindowSize/2,FFTWindowSize)', 1, sn);
   X = X -  Xx .* repmat(mean(diff(X,1,1)), FFTWindowSize, 1);
end
if DETRENDING_ORDER > 2 % 2-nd order detrending...
   X = X - 1/2 * (Xx.^2) .* repmat(mean(diff(X,2,1)), FFTWindowSize, 1);
end
if DETRENDING_ORDER > 3 % 3-rd order detrending...
   X = X - 1/6 * (Xx.^3) .* repmat(mean(diff(X,3,1)), FFTWindowSize, 1);
end

% Power spectrum...
Y = fft(X, FFTWindowSize);
Pyy =  Y.* conj(Y) / FFTWindowSize;
F = fs * (0:FFTWindowSize/2)' / FFTWindowSize;

% ndxF = find(F>=MUAFreqBand(1) & F<=MUAFreqBand(2));
ndxF = find(F<=MUAFreqBand(2));
ndxF = ndxF(2:end);

if exist('PyyBaseline') ~= 1
   PyyBaseline = mean(Pyy(ndxF,:), 2);
   PyyBaseline = PyyBaseline / PyyBaseline(ndxF(1));
end

if exist('PlotPowerSpectrum') ~= 1
   PlotPowerSpectrum = 0;
end
if PlotPowerSpectrum
   figure;
   plot(F(ndxF), PyyBaseline,'.-');
   set(gca, 'XLim', [F(2) F(end)], 'XScale', 'log', 'YScale', 'log');
   xlabel('\omega/2\pi (Hz)');
   ylabel('baseline Power spectrum (a.u.)');
end

NormPyy = Pyy(ndxF,:) ./ repmat(PyyBaseline, 1, sn);
MUA = mean(NormPyy);
% t = (FFTWindowSize/2:FFTWindowSize:length(LFP))*dt + t(1);
t = mean(reshape(t(1:FFTWindowSize*sn), FFTWindowSize, sn));


% plotHistogram(MUA.value);
% plotHistogram(sqrt(MUA.value))
% plotHistogram(log(MUA.value))
% 
% sqrtMUA.value = sqrt(MUA.value);
% sqrtMUA.value = sqrtMUA.value / std(sqrtMUA.value);
% sqrtMUA.time = MUA.time;
% 
% logMUA.value = log(MUA.value);
% logMUA.value = logMUA.value/std(logMUA.value);
% logMUA.time = MUA.time;
