% analyzeUpDownStateVsRecordingSet.m
%
%   Copyright 2015 Maurizio Mattia @ Ist. Super. Sanità, Rome - Italy
%   Version: 1.0 - Feb. 18, 2015
%


%% INITIALIZATION (path, parameters and option, debug level, counters...)
%
close('all');
%clear('all');
clearvars;

debug=0; % debug level
step='step1'; % (first step of the SWAP)

setPath;
setParamsAndOptions;

disp(['[' mouse ' mouse]']);
disp(['FileName: ' FileName]);

%% MUAdir for Cristiano
%
MUAdir = [BaseDir 'RESULTS/Nov2016/' emission '/' mouse '/' mouse '_MUAdir/' FileName '/'];
if exist(MUAdir,'dir') == 0
    mkdir(MUAdir);
end


%%
Disc=[];
WeakBi=[];
SkewP=[];
SkewN=[];
RightPeak=[];
FitWarning=[];
LargeUDTH=[];
FewTrans=[];

%% Open data file...
%
hSMR = openSONandLoadInfo(DataFile);
if isempty(hSMR)
   return
end

%% Create ResultsDir
%
Results = AnalysisDir;
if ~exist(Results,'dir') % if NOT exist --> CREATE
   mkdir(Results);
end

%% Analysis loop on recording set...

%% Get the information on the stimulation from EventArray...
%
if strcmpi(emission,'Evoked') 
   Stim = getSONEventData(hSMR, EventArray.port, Options.PeriodToAnalyze);
end

%%
try
    
%% START LOOP on the RecordingSet    
fprintf('\n');
%for crsRecSet = 1:length(RecordingSet)
for crsRecSet = ChannelSet % modified for running on a sub-set of recording channels (a cortical area)
   %% Makes the working directory and change the current dir to it...
   %
   tempdir = [AnalysisDir RecordingSet(crsRecSet).label];
   if exist(tempdir,'dir') == 0
      mkdir(tempdir);
   end
   cd(tempdir);
   disp('----------------------------------------');
   fprintf(['Recording set: ' RecordingSet(crsRecSet).label '\n']);

   %% Loads dataset...
   %
   RawSig = getSONAnalogData(hSMR, RecordingSet(crsRecSet).port, Options.PeriodToAnalyze);
   RawSig.dt = diff(RawSig.time(1:2));
   fprintf('Analyzed period: [%g,%g] s\n', RawSig.time(1), RawSig.time(end));        
   
   clear ndx
   close all
     
   
   %%  Compute baseline for spectral MUA estimate...
   %
   MovingWindowSize = 1/Options.LogMUA.FreqBand(1); % [s]
   DetrendingOrder = 0;
   [PyyFreq, PyyBaseline] = plotMedianPSDofLFP(RawSig, MovingWindowSize, DetrendingOrder);
   ndx = find(PyyFreq>=Options.LogMUA.FreqBand(1)*0.99 & PyyFreq<=Options.LogMUA.FreqBand(2)*1.01);
   PyyFreq = PyyFreq(ndx);
   PyyBaseline = PyyBaseline(ndx);
   
   
   %% Compute spectral estimate of MUA...
   %
   [MUA.time, MUA.value] = computeSpectralEstimateOfMUA(RawSig.time, RawSig.value, ...
         Options.LogMUA.FreqBand, 0, PyyBaseline);
   MUA.dt = diff(MUA.time(1:2));

   
   %% Compute smoothed log(MUA): rsMUA...
   %
   logMUA.value = log(MUA.value);
   logMUA.time = MUA.time;
   logMUA.dt = MUA.dt;
   
   [rsMUA.time, rsMUA.value] = computeMovingAverage(logMUA, ...
      round(Options.LogMUA.SmoothingWindow / logMUA.dt), 1);
   rsMUA.dt = diff(rsMUA.time(1:2));
   
   
   %% Check for discontinuities in the LogMUA (and correct/remove them)
   %
   REF=median(rsMUA.value);
   % PARAMETERS
   lowTHR=REF-2;
   gap = 20; % number of samples around the discontinuity, to isolate and cut the discontinuity
   
   % LOOP for looking for discontinuities in the LogMuA
   newloop = 1;
   firstloop = 1; i=0; A=[];
   while (newloop)
       j=i+1;
       if firstloop; from = 1; firstloop = 0; else; from=A(i); end
       if ~isempty(find(rsMUA.value(from:end)<lowTHR,1)) % the START of the discontinuity is found. Look for the END
           A(j)=find(rsMUA.value(from:end)<lowTHR,1);% find below THR
           if j~=1; A(j)=A(j)+A(j-1); end
           i=j;           
           j=i+1; from=A(i);
           if ~isempty(find(rsMUA.value(from:end)>lowTHR,1)) % the END of the discontinuity is found. Look for if there is ANOTHER one
               A(j)=find(rsMUA.value(from:end)>lowTHR,1); % find above THR
%              if j~=1; A(j)=A(j)+A(j-1);end; %(here the check is not needed, if here, j cannot be = 1)
               A(j)=A(j)+A(j-1);   
               newloop = 1; i=j;
           else % no END of the discontinuity found... some problems?
               disp('WARNING: logMUA below THR till the end of the recording time'); %break;
           end
       else; newloop = 0; 
       end % the START of the discontinuity is NOT found. There are NO/NO MORE discontinuities. End the loop      
   end % (wjile loop)

   CutInterval=[]; CutSamples = [];
   % Introduce a "safety GAP" to isolate the discontinuity   
   if ~mod(numel(A),2) % numel(A) is even --> short and limited discontinuities identified
        % DOWNWARD
        B=A(1:2:end)-gap;
        C=A(2:2:end)+gap;
        CutInterval = [rsMUA.time(B);rsMUA.time(C)];
   else; disp('WARNING: odd number of discontinuity borders');
   end
   
   if ~isempty(CutInterval)
       disp('DISCONTINUITY FOUND'); Disc=[Disc crsRecSet]; 
       disp(CutInterval)
   end
   
%    test=rsMUA;
%    % CUT away the discontinuity 
%    for idx = length(B):-1:1
%        test.value(B(idx):C(idx))=[];
%        test.time(B(idx):C(idx))=[];
%    end
%    
%    figure
%    plot(rsMUA.time,rsMUA.value)
%    hold on
%    plot(test.time,test.value)
%    
%    rsMUA=test;
   
   %% Subsample LFP...
   %
   WindowSize = round(1/RawSig.dt/Options.LFP.SubsamplingRate);
   [ssLFP.time,ssLFP.value] = computeMovingAverage(RawSig, WindowSize, WindowSize);
   ssLFP.dt = diff(ssLFP.time(1:2));
   

   %% High-pass filtering of raw signal (LFP), if required...
   %
   if isfield(RecordingSet(crsRecSet),'HPFCutOffFreq')
      if ~isempty(RecordingSet(crsRecSet).HPFCutOffFreq)
         disp('LFP high-pass filtering...')
         [lpLFP.time,lpLFP.value] = computeMovingAverage(ssLFP, ...
            round(1/RecordingSet(crsRecSet).HPFCutOffFreq/ssLFP.dt), ...
            round(1/RecordingSet(crsRecSet).HPFCutOffFreq/ssLFP.dt/10));
         ssLFP.value = ssLFP.value - interp1(lpLFP.time,lpLFP.value,ssLFP.time,'PCHIP'); % modified for R2016a
         clear lpLFP
      end
   end
   
   %% Low-pass filtering of LFP...
   %
   [ssLFP.time,ssLFP.value] = computeMovingAverage(ssLFP, ...
         round(Options.LFP.SmoothingWindow/ssLFP.dt), 1);

   
   %% Remove bad periods (Cut Intervals with Discontinuities)...

   if ~isempty(CutInterval)
       
%  *** Check if any overlapping intervals...
       %%% CutIntervalBackup = CutInterval;
       for idx = 1:1:size(CutInterval,2)-1
           if (CutInterval(2,idx)>=CutInterval(1,idx+1) && CutInterval(2,idx)<=CutInterval(2,idx+1))
               %disp([CutInterval(2,idx) CutInterval(1,idx+1)])
               CutInterval(2,idx)=CutInterval(2,idx+1);
               CutInterval(:,idx+1)=[];
           end
       end

% *** Evaluate CutSamples
       ndx=[];

       for idx = 1:1:size(CutInterval,2)
         ndx = [ndx find(rsMUA.time>=CutInterval(1,idx) & rsMUA.time<=CutInterval(2,idx))];
         CutSamples = [CutSamples numel(ndx)]; %CutSamples(end) = total number of CutSamples (cumulative)
       end
       %ndx=unique(ndx);
       mdx=1:1:length(rsMUA.value);
       jdx=setdiff(mdx,ndx);

       %%% rsMUAbackup = rsMUA;
       rsMUA.time = rsMUA.time(jdx);
       rsMUA.value = rsMUA.value(jdx);

       % CUT the discontinuity also in the LFP
       %%% ssLFPbackup = ssLFP;
       ssLFP.time = ssLFP.time(jdx);
       ssLFP.value = ssLFP.value(jdx);

       %CutSamples=numel(ndx);
       
   end
   
   figure
   plot(rsMUA.time,rsMUA.value)
   
   XPos=RecordingSet(crsRecSet).XPos; YPos=RecordingSet(crsRecSet).YPos;
   % Save rsMUA (for Cristiano)
   save([MUAdir RecordingSet(crsRecSet).label(1:end-4) 'rsMUA.mat'], 'rsMUA',...
       'XPos', 'YPos');
   
   
   %% Compute the bimodal distribution of resampled MUA (I)
   %
   clear('ModeParams');
   %[ModeParams, THR_2ndPeak,~,~,F1,F2] = plotBimodalHistogram(rsMUA.value);
   [ModeParams, THR_2ndPeak,~,~,F1] = plotBimodalHistogram(rsMUA.value);
   %close(gcf);

   figure(F1);
   xlabel('log(MUA)');
   set(F1.Children(2).Title,'String','logMUA')
   set(F1.Children(2).XLabel,'String','log(MUA)')
   %title('logMUA');
   set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.5 3.5 5 4]);
   print('-deps2c', ['rsLogMUA.bihist.SW_' ...
   num2str(Options.LogMUA.SmoothingWindow*1000) '.eps']);
   
% %    if (ModeParams.FirstLeft ~=0)
% %        figure(F2); 
% %        xlabel('log(MUA)');
% %        %title('Distribution of the log(MUA) - TAIL');
% %        set(F2.Children(2).Title,'String','Distribution of the log(MUA) - TAIL')
% %        set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.5 3.5 5 4]);
% %        print('-deps2c', 'rsLogMUA.scalefree.eps');
% %    else
% %        close(F2)
% %    end
    
% ModeParams = 
%            Mu1: 
%         Sigma1: 
%          Ampl1: 
%            Mu2: 
%         Sigma2: 
%          Ampl2: 
%       FitError:
% SecondPeakArea:

%[GDB] there should be added a check if the fit procedure is not converging...

   % --> check the RightPeak
   if ModeParams.FirstLeft==0; % the fisrt (i.e. largest) peak is on the RIGHT)
       RightPeak=[RightPeak crsRecSet]; 
   end
   
   %% STANDARDIZE logMUA to have the first (i.e. largest) peak in 0...
   %  [GDB] [was: SHIFT, is: STANDARDIZE]
   %  NOTE: this is not the final reference value for log(MUA)
   %  corresponding to no firing activity... 
   %
   
% %    if isfield(ModeParams,'Mu1')
   if ModeParams.FirstLeft
      LogMUAReference = ModeParams.Mu1;
      Delta = ModeParams.Sigma1;
   else
      LogMUAReference = ModeParams.Mu2;
      Delta = ModeParams.Sigma2;
   end
   % keep (and print) the absolute scale
   RecordingSet(crsRecSet).LogMUAReference = LogMUAReference;
   fprintf('LogMUAReference = %f\n', LogMUAReference);
   
% %    if isfield(ModeParams,'Sigma1')
% %       Delta = ModeParams.Sigma1;
% %    else
% %       Delta = ModeParams.Sigma2;
% %    end
   
   rsMUA.value = rsMUA.value - RecordingSet(crsRecSet).LogMUAReference;
   %[GDB] June 2018 NOT only shift, but also STANDARDIZE the distribution
   rsMUA.value = (rsMUA.value - RecordingSet(crsRecSet).LogMUAReference)/Delta;
   

   %% Compute the bimodal distribution of MUA (II)
   %
   %disp('--- STANDARDIZE logMUA and re-fit...');
   disp('--- SHIFT logMUA and re-fit...');
   %[ModeParams,~,hMUA,hX,F1,F2] = plotBimodalHistogram(rsMUA.value);
   [ModeParams,~,hMUA,hX,F1] = plotBimodalHistogram(rsMUA.value);

   % hMUA is the histogram content of the distribution of logMUA i.e.
   % frequency count of the nBins values (nBins = 100, centered in hX)
   
   figure(F1);
   xlabel('log(MUA)');
   %set(F1.Children(2).Title,'String','logMUA (standardized)')
   set(F1.Children(2).Title,'String','logMUA (shifted)')
   set(F1.Children(2).XLabel,'String','log(MUA)')
   set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.5 3.5 5 4]);
   %print('-deps2c', ['rsLogMUA.STANDARD.bihist.SW_' ...
   %   num2str(Options.LogMUA.SmoothingWindow*1000) '.eps']);
   print('-deps2c', ['rsLogMUA.SHIFT.bihist.SW_' ...
      num2str(Options.LogMUA.SmoothingWindow*1000) '.eps']);
  %(SW = Smoothing Window)
  
% %    if ModeParams.FirstLeft ~=0
% %        figure(F2); 
% %        xlabel('log(MUA)');
% %        %title('Distribution of the log(MUA) - TAIL (standardized)');
% %        set(F2.Children(2).Title,'String','Distribution of the log(MUA) - TAIL')
% %        set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.5 3.5 5 4]);
% %        %print('-deps2c', 'rsLogMUA.STANDARD.scalefree.eps');
% %        print('-deps2c', 'rsLogMUA.SHIFT.scalefree.eps');
% %     else
% %        close(F2)
% %     end
   
   
   %% Check ModeParams: fit parameters, skewness, comparison between the fit and the distribution indexes
   %
   disp(ModeParams);
   disp('Total Distribution:')
   disp(ModeParams.SkewnessTot);
   disp('Tail Distribution:')
   disp(ModeParams.SkewnessTail);
   
   if ModeParams.SecondPeakArea<THR_2ndPeak
       %fprintf('--- WeakBimodality (2nd peak is NOT FITTED)\n')
       fprintf('--- WeakBimodality\n')
       WeakBi=[WeakBi crsRecSet];
   end
   
   if (ModeParams.SkewnessTail.gamma)>1
       disp('--- [tail] StrongAsymmetry (positive)');
       SkewP=[SkewP crsRecSet];
   elseif (ModeParams.SkewnessTail.gamma)<-1
       disp('--- [tail] StrongAsymmetry (negative)');
       SkewN=[SkewN crsRecSet];
   end
   
   % ... some checks (considering the values of Mu2 and the indexes of N2 distribution)
%    if isfield(ModeParams,'Mu2') && isfield(ModeParams,'Sigma2')
%        if(abs(ModeParams.SkewnessTail.MU-ModeParams.Mu2))>ModeParams.Sigma2
%            disp('GaussianFit far from distribution...')
%        end
% %        list=[ModeParams.Mu2 ModeParams.Mu2+ModeParams.Sigma2...
% %            ModeParams.Mu2-ModeParams.Sigma2 ModeParams.Skewness.MU...
% %            ModeParams.Skewness.MU+ModeParams.Skewness.SIG...
% %            ModeParams.Skewness.MU-ModeParams.Skewness.SIG];      
%    end
       

   %%  Set UD_THRESHOLD
   %
% % %    % option-1: if the SD of the first peak is stable across the channels, 
% % %    %           use a fixed threshold (AbsoluteThreshold) 
% % %    % option-2: modulate the threshold for each channel considering the 
% % %    %           value of SD1 --> FalsePositive = 2.5% of the Gaussian (FirstPeak)
% % %    % option-3: modulate the threshold for each channel considering the 
% % %    %           value of Mu1 and Mu2
% % %    
% % %    clear('opt');
% % %    UD.ThresholdModulation = 0.6; 
% % %    opt(1)=RecordingSet(crsRecSet).AbsoluteThreshold;
% % %    %%%opt(1)=0.35;
% % %    if isfield(ModeParams,'Sigma1') && isfield(ModeParams,'Mu1')
% % %        opt(2)= ModeParams.Mu1+2.*ModeParams.Sigma1; % 2sigma --> above threshold is NOT DOWN
% % %    end    
% % %    if isfield(ModeParams,'Mu1') && isfield(ModeParams,'Mu2')
% % %        opt(3)= UD.ThresholdModulation*(ModeParams.Mu2-ModeParams.Mu1)+ ModeParams.Mu1;
% % %    end
% % %    % ... but the Gaussian Fit of the second peak may be meaningless if the
% % %    % distribution is strongly skewed...
% % %    
% % %    fprintf('options for UD_THRESHOLD:')
% % %    for i=1:numel(opt); fprintf(' %f ', opt(i));end
% % %    fprintf('\n');
% % %       
% % %    
% % %    % ------------
% % %    % UD_THRESHOLD
% % %    % ------------
% % %    %UD_THRESHOLD = mean(opt);
% % %    if (numel(opt)~=1) UD_THRESHOLD = opt(2); else UD_THRESHOLD = opt(1); end
% % %    fprintf('UD_THRESHOLD = %f\n', UD_THRESHOLD);
% % %    
% % %    
% % %    % ... some checks (considering the values of Mu2 and Sigma2, and the skewness)
% % %    if(abs(ModeParams.SkewnessTail.gamma)<1) % if not a strong asymmetry...
% % %        if isfield(ModeParams,'Mu2') && isfield(ModeParams,'Sigma2')
% % %            if(UD_THRESHOLD>ModeParams.Mu2) 
% % %                disp('WARNING(1). UD_THRESHOLD is larger than the SecondPeak (Mu2)!');
% % %                LargeUDTH=[LargeUDTH crsRecSet];
% % %            elseif (ModeParams.Mu2-ModeParams.Sigma2)<UD_THRESHOLD
% % %                disp('WARNING(2). UD_THRESHOLD is too close to the SecondPeak (Mu2)!'); 
% % %            elseif (ModeParams.Mu2-2*ModeParams.Sigma2)<UD_THRESHOLD
% % %                disp('WARNING(3). UD_THRESHOLD is not far from the SecondPeak (Mu2)!'); 
% % %            end
% % %        end
% % %    else % if strong asymmetry, comparison made with MU and not Mu2
% % %        if(UD_THRESHOLD>ModeParams.SkewnessTail.MU) 
% % %            disp('WARNING. UD_THRESHOLD is larger than the SecondPeak (MU)!');
% % %            LargeUDTH=[LargeUDTH crsRecSet];
% % %        end
% % %    end
   
   % ...for STANDARDIZED variable
   UD_THRESHOLD = 2.; % threshold set at 2*sigma from the first Gaussian
   
   if(UD_THRESHOLD>ModeParams.SkewnessTail.MU) 
       disp('WARNING. UD_THRESHOLD is larger than Tail.MU!');
       LargeUDTH=[LargeUDTH crsRecSet];
   end
   
   if ModeParams.FirstLeft
       if(abs(2*ModeParams.Sigma1-UD_THRESHOLD)>0.05) OR (abs(ModeParams.Mu1)>0.05)
           % a 5% difference with respect to the standard gaussian is admitted 
           disp('WARNING - check the results of the fit (Modeparams)');
           FitWarning=[FitWarning crsRecSet];
       end
   else
       if(abs(2*ModeParams.Sigma2-UD_THRESHOLD)>0.05) OR (abs(ModeParams.Mu2)>0.05)
           % a 5% difference with respect to the standard gaussian is admitted 
           disp('WARNING - check the results of the fit (ModeParams)');
           FitWarning=[FitWarning crsRecSet];
       end
   end

   % -------------------------------------------------------------------------
   % evaluate the area of the histogram (N=hMUA) beyond the first peak (UD_THRESHOLD)
   % --> estimate the percentage in UP state
   % (rough) linear interpolation of the histogram
   % ---------------------------
   
%    X0=find(hX<UD_THRESHOLD,1,'last');
%    X1=find(hX>UD_THRESHOLD,1);
%    N0=hMUA(X0); N1=hMUA(X1);
%    NTH=N0+(UD_THRESHOLD-hX(X0))*(N1-N0)/(hX(X1)-hX(X0));
      
   UPpercent=sum(hMUA(find(hX>=UD_THRESHOLD)));
   fprintf('histo_TH = %f\n', hX(find(hX>=UD_THRESHOLD,1)));
   fprintf('~~~ estimate of the percentage in the UP state: %f\n',UPpercent); 
   ModeParams.UPpercent = UPpercent;
   
   
%% Set MIN_STATE_LEN
%
   MIN_STATE_LEN = RecordingSet(crsRecSet).MinStateDuration;

% As an alternative, a modulation of the "UDcycle" (somehow determined...) 
% can be considered...
%%% MIN_STATE_LEN = someUDcycle * RecordingSet(crsRecSet).UDcycleFraction;

   
   %% The following part of is perfomed only if a threshold for MUA Up and
   %  Down state detection has been set..
   %
   LogMUAbaseline = NaN;
   UpState = [];
   DownState = [];
   UDcycle = [];
   
   if ~isnan(UD_THRESHOLD)
      disp('Up and Down states detectable...')

      %%  Comute Binary MUA
      %      
      UD = computeBinaryUDState(rsMUA, UD_THRESHOLD, MIN_STATE_LEN);
      
      
      %% Single out and adjust transitions times...
      %
      NoTrans=0;
      [Trans, UD] = findAndAdjustUDTransition(rsMUA, UD, UD_THRESHOLD);
      nTrans=length(Trans.val);

      if length(Trans.val)<3 % at least 3 Trans for a valid evaluation of the states duration
          NoTrans=1;
          disp(['WARNING: ' num2str(length(Trans.val)) ' transition(s) found (nTrans < 3)']);
          disp('... it is not possibile to evaluate the duration of Up and Down states.');
          FewTrans=[FewTrans crsRecSet];
      end
    
      
      %% Save if needed
      %
      if isfield(Options, 'SaveMUA')
         if Options.SaveMUA
            save('UD.mat', 'UD', 'Trans');
         end
      end
      
      
      %% Plots together a sample of logMUA and RawSig...
      %
      PERIOD_TO_PLOT = [-20 0] + RawSig.time(end);
      plotLFPlogMUAandUD(RawSig, rsMUA, UD, PERIOD_TO_PLOT);
      
      
      %%  Histograms of Up and Down state duration...
      %   ... and of UDcycle duration
      %
      if NoTrans~=1
  
          if ~isempty(CutInterval)
            for idx = 1:1:size(CutInterval,2)
                F=find(Trans.time>CutInterval(1,idx) & Trans.time<CutInterval(2,idx),1);
                % Find if any Trans are in the removed intervals
                if ~isempty(F); disp('WARNING: Transitions found in the removed intervals'); end
            end
          end
                             
          if ~isempty(CutInterval)
              [UpStateLen, DownStateLen] = plotStateDurationHist(Trans, CutInterval);
              [UDcycleLen] = plotDurationHist(Trans, CutInterval);
          else
              [UpStateLen, DownStateLen] = plotStateDurationHist(Trans);
              [UDcycleLen] = plotDurationHist(Trans);              
          end
          
          DownState.mean_duration = mean(DownStateLen);
          DownState.std_duration = std(DownStateLen);
          DownState.numel_duration = length(DownStateLen);
          UpState.mean_duration = mean(UpStateLen);
          UpState.std_duration = std(UpStateLen);
          UpState.numel_duration = length(UpStateLen);
          UDcycle.mean_duration = mean(UDcycleLen);
          UDcycle.std_duration = std(UDcycleLen);
          UDcycle.numel_duration = length(UDcycleLen);
          UDcyle.stderr_duration = std(UpStateLen)/sqrt(numel(UpStateLen) ...
              + std(DownStateLen)/sqrt(numel(DownStateLen)));
          
          fprintf('... nTrans = %d\n',UDcycle.numel_duration);
          save('UpDownStateDuration', 'UpStateLen', 'DownStateLen', 'UDcycleLen');
          
          disp(['<UD cycle> = ' num2str(mean(UpStateLen)+mean(DownStateLen),3) ...
              's (st.err.: ' num2str(std(UpStateLen)/sqrt(numel(UpStateLen)) ...
              + std(DownStateLen)/sqrt(numel(DownStateLen)),2) 's)']);
      
      %%  Rasters and PSTH of log(MUA) and UD around UPWARD transitions...
      %
          WF_TIME_RANGE = [-0.5 1.5]; % Time range around the trigger event.

          % Extract useful UPWARD transitions times, the trigger for raster and PSTH...
          ndx = Trans.val == 1;
          Triggers = Trans.time(ndx); % Triggers = times of the upward transitions
          SortStateLen = UpStateLen;
          Triggers = Triggers(1:numel(SortStateLen));
    % ************************************************************************     
          % Check if any Trigger is too much close to the file border...
          ndx = find(Triggers - UD.time(1) > -WF_TIME_RANGE(1) & ...
             UD.time(end) - Triggers > WF_TIME_RANGE(2));
          Triggers = Triggers(ndx);
          SortStateLen = SortStateLen(ndx);
    % ************************************************************************      
    % Check if any TriggerWindow OVERLAPS with any CUT INTERVALS...
          jList=[]; %tList=[];
          if ~isempty(CutInterval)
              for j = 1:length(Triggers)
                TriggerInterval = (Triggers(j)+WF_TIME_RANGE);
                for idx = 1:1:size(CutInterval,2)  

                    [timeline,P] = sort([TriggerInterval CutInterval(:,idx)']);
                    Q=diff(P);
    % [1 2 3 4], [3 4 1 2] --> NOT OVERLAPPING intervals
    % [1 3 4 2], [1 3 2 4],[3 1 2 4],[3 1 4 2] --> OVERLAPPING intervals 
                    if Q(1)~=1 % OVERLAPPING INTERVALS
    %                    disp([j Triggers(j) timeline P]);
                        jList = [jList j];
                    end                
                end            
              end
          end
          jList = unique(jList);
    % ************************************************************************      
          % Remove bad waveforms... 

          if ~isempty(jList)
              RecordingSet(crsRecSet).WFsToRemove = jList;
          end             

          if isfield(RecordingSet(crsRecSet), 'WFsToRemove')
             if numel(RecordingSet(crsRecSet).WFsToRemove) > 0
                ndx = 1:numel(Triggers);
                ndx = setdiff(ndx, RecordingSet(crsRecSet).WFsToRemove);
                Triggers = Triggers(ndx);
                SortStateLen = SortStateLen(ndx);
             end
          end
    % ************************************************************************       
          save('UpwardTransitions.mat', 'Triggers');

          WFs = extractEventCenteredWaveForms(UD.time, UD.value, ...
                  Triggers, WF_TIME_RANGE, CutInterval, CutSamples);

          plotRasterOfWFs(WFs, [0.0 1.0], 'Up/Down');
          delete('RasterOfWFs.eps')

          colormap((1-gray/2));
          YLim = get(gca, 'YLim');
          hold on
          plot([0 0], YLim, 'r-');
          xlabel('Time (s)');
          ylabel('Down-Up Transition');

          set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
          print('-deps2c', 'UD.Raster.eps');


          %%  PSTH of log(MUA)...
          %

          WFs = extractEventCenteredWaveForms(rsMUA.time, rsMUA.value, ...
             Triggers, WF_TIME_RANGE, CutInterval, CutSamples);

          [tWF, MeanWFs, StdWFs] = plotAverageWFs(WFs, 5);

          set(gca, 'TickDir', 'out', 'Box', 'on', 'Layer', 'top', ...
             'XLim', WF_TIME_RANGE);

          xlabel('Time (s)');
          ylabel('log(MUA)');

          set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0.25 3.5 7.5 4]);
          print('-deps2c', 'LogMUA.Average.eps')

          SampleSize = length(WFs);
          save('LogMUA.Average.mat', 'tWF', 'MeanWFs', 'StdWFs', 'SampleSize');


          %% Average Log(MUA) around upward transitions of pure sequencies DOWN-UP...
          %

          [t, MeanY, StdY] = plotAverageUpwardTransition(rsMUA, UD, Triggers, ...
              WF_TIME_RANGE, CutInterval, CutSamples);

          save('LogMUAUpTrans.Average.mat', 't', 'MeanY', 'StdY');


          %%  Rasters of log(MUA)...
          %
          plotRasterOfWFs(WFs, prctile(rsMUA.value,[2.5 100-2.5]), 'log(MUA)');
          [RasterUpward.t, RasterUpward.Ys] = computeRasterOfWFs(WFs);
          delete('RasterOfWFs.eps')

          xlabel('Time (s)');
          ylabel('Down-Up Transition');
          set(gca,'TickDir','out');

          set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
          print('-deps2c', 'LogMUA.Raster.eps');


          %%  Rasters of log(MUA) sorted by Up state length...
          %
          [val,ndx] = sort(SortStateLen);


          WFs = extractEventCenteredWaveForms(rsMUA.time, rsMUA.value, ...
             Triggers(ndx), WF_TIME_RANGE, CutInterval, CutSamples);

          plotRasterOfWFs(WFs, prctile(rsMUA.value,[2.5 100-2.5]), 'log(MUA)');
          delete('RasterOfWFs.eps')


          xlabel('Time (s)');
          ylabel('Down-Up Transition');

          set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
          print('-deps2c', 'LogMUA_Sorted.Raster.eps');


          %%  Rasters and PSTH of log(MUA) around DOWNWARD transitions...
          %
          WF_TIME_RANGE = [-0.5 1.5]; % Time range around the trigger event.

          % Extract useful downward transitions times, the trigger for raster and PSTH...
          Triggers = Trans.time(Trans.val == 0);
          SortStateLen = DownStateLen;
          Triggers = Triggers(1:numel(SortStateLen));
    % ************************************************************************     
    % Check if any Trigger is too much close to the file border...
          ndx = find(Triggers - UD.time(1) > -WF_TIME_RANGE(1) & ...
             UD.time(end) - Triggers > WF_TIME_RANGE(2));
          Triggers = Triggers(ndx);
          SortStateLen = SortStateLen(ndx);
    % ************************************************************************      
    % Check if any TriggerWindow OVERLAPS with any CUT INTERVALS...
          jList=[]; %tList=[];
          if ~isempty(CutInterval)
              for j = 1:length(Triggers)
                TriggerInterval = (Triggers(j)+WF_TIME_RANGE);
                for idx = 1:1:size(CutInterval,2)  

                    [timeline,P] = sort([TriggerInterval CutInterval(:,idx)']);
                    Q=diff(P);
    % [1 2 3 4], [3 4 1 2] --> NOT OVERLAPPING intervals
    % [1 3 4 2], [1 3 2 4],[3 1 2 4],[3 1 4 2] --> OVERLAPPING intervals 
                    if Q(1)~=1 % OVERLAPPING INTERVALS
    %                    disp([j Triggers(j) timeline P]);
                        jList = [jList j];
                    end                
                end            
              end
          end
          jList = unique(jList);
    % ************************************************************************ 
      % Remove bad waveforms...
          if ~isempty(jList)
              RecordingSet(crsRecSet).WFsToRemove = jList;
          end

          if isfield(RecordingSet(crsRecSet), 'WFsToRemove')
             if numel(RecordingSet(crsRecSet).WFsToRemove) > 0
                ndx = 1:numel(Triggers);
                ndx = setdiff(ndx, RecordingSet(crsRecSet).WFsToRemove);
                Triggers = Triggers(ndx);
                SortStateLen = SortStateLen(ndx);
             end
          end
    % ************************************************************************ 
          save('DownwardTransitions.mat', 'Triggers');

          WFs = extractEventCenteredWaveForms(rsMUA.time, rsMUA.value, ...
             Triggers, WF_TIME_RANGE, CutInterval, CutSamples);

          [tWF, MeanWFs, StdWFs] = plotAverageWFs(WFs, 5);

          set(gca, 'TickDir', 'out', 'Box', 'off', 'Layer', 'top', ...
             'XLim', WF_TIME_RANGE);

          xlabel('Time (s)');
          ylabel('log(MUA)');

          set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0.25 3.5 7.5 4]);
          print('-deps2c', 'LogMUA.Downward.Average.eps')

          SampleSize = length(WFs);
          save('LogMUA.Downward.Average.mat', 'tWF', 'MeanWFs', 'StdWFs', 'SampleSize');

          LogMUAbaseline = min(MeanWFs(tWF > 0));

          %  Rasters of log(MUA)...
          %
          plotRasterOfWFs(WFs, prctile(rsMUA.value,[2.5 100-2.5]), 'log(MUA)');
          %       plotRasterOfWFs(WFs, [-0.75 3.75], 'log(MUA)');
          [RasterDownward.t, RasterDownward.Ys] = computeRasterOfWFs(WFs);
          delete('RasterOfWFs.eps')

          xlabel('Time (s)');
          ylabel('Up-Down Transition');

          set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
          print('-deps2c', 'LogMUA.Downward.Raster.eps');


          %%  Rasters of log(MUA) sorted by Down state duration...
          %
          [val,ndx] = sort(SortStateLen);

          WFs = extractEventCenteredWaveForms(rsMUA.time, rsMUA.value, ...
             Triggers(ndx), WF_TIME_RANGE, CutInterval, CutSamples);

          plotRasterOfWFs(WFs, prctile(rsMUA.value,[2.5 100-2.5]), 'log(MUA)');
          delete('RasterOfWFs.eps')

          xlabel('Time (s)');
          ylabel('Up-Down Transition');

          set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
          print('-deps2c', 'LogMUA_Sorted.Downward.Raster.eps');


          %% Average Log(MUA) around downward transitions...
          %

          [t, MeanY, StdY] = plotAverageDownwardTransition(rsMUA, UD, Triggers, ...
                  WF_TIME_RANGE, CutInterval, CutSamples);

          save('LogMUADownTrans.Average.mat', 't', 'MeanY', 'StdY');

       end % if NoTrans~=1
   end % if ~isnan(UD_THRESHOLD) 
 
   %% Saves some analysis results...
   %

  save('AnalysisSummary.mat','ModeParams','LogMUAbaseline','PyyFreq',...
       'PyyBaseline','LogMUAReference','Delta','UD_THRESHOLD','nTrans',...
       'UpState','DownState','UDcycle');
   %...there was opt in the list of variables

   %% Returns to the root directory...
   %
   cd('..');
end % Analysis loop on recording set (channels)

%% Check if any channel has to be DISCARDED on the basis of 'quality parameters'
%
% QUALITY PARAMETERS:
% - SD = Standard Deviation of the Gaussian Peak that fits the logMUA distribution of the DownState
% - logMuaRef = Mu of the Gaussian Peak (before the standardization)
% - Frequency = 1/<UDcycle>
disp('>> >> >> >> >> >> >> >> >> >> >> >> >> >> >> >> >> >> >> >> >> >>');
OUT = []; % Contains the list of Outlier Channels (for any of the quality parameter: SD, logMuaRef, Freq)
%% QUALITY PARAMETER: SD and logMuaRef
%
SDs = NaN(1,numel(RecordingSet));
LogMuaRef = NaN(1,numel(RecordingSet));

% *** Analysis loop on recording set...
for crsRecSet = ChannelSet % modified for running on a sub-set of recording channels (a cortical area)
   workingdir = [AnalysisDir RecordingSet(crsRecSet).label];
   load([workingdir '/AnalysisSummary.mat'], 'Delta');
   SDs(crsRecSet)=Delta;
   load([workingdir '/AnalysisSummary.mat'], 'LogMUAReference');
   LogMuaRef(crsRecSet)=LogMUAReference;
end

% *** Check the stability of a MUA threshold across channels.
%  In other words, if the same threshold can be used for all the channels.
%  If yes, this should provide more reliable profiles of traveling
%  wavefronts.
% [GDB] ...but now we are using standardized variables, and UD_THRESHOLD is
% set at 2 sigma of the first (i.e. dominant) peak, fitted with a Gaussian

%%
%%%%%%%%%%%
%   SD 
%%%%%%%%%%%
figure
stem(SDs);
set(gca, 'XLim', [0 33])%, 'YLim',[0 0.5]) % limit set for the SD
xlabel('Channels')
ylabel('SD')

%[GDB] there should be added a check if the fit procedure is not converging...
% --> NaN values for logMuaRef and SD

SDsOK = ~isnan(SDs); % median and mean are computed only on non-NaN values
medSD = median(SDs(SDsOK));
meanSD=mean(SDs(SDsOK));
title(['[logMUA] SD of the LargestPeak. Median SD = ' num2str(medSD,3) ]);

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0 0 5 4]);
print('-deps2c', 'DownMUASDsVsChannels.eps');

% *** Find OUTLIERS in the SDs
Q1=prctile(SDs, 25); Q3=prctile(SDs, 75); IQR=Q3-Q1; ww=1.5*IQR;
% OUT=find(SDs>Q3+ww);
% if (numel(OUT)~=0); 
%     disp('WARNING: outlier channels');disp(OUT);

OUT_SD.low=find(SDs<Q1-ww);
if (numel(OUT_SD.low)~=0) 
    disp(['(OUT_SD.low: ' num2str(OUT_SD.low) ')']);end

%%%OUT_SD.high = find(SDs>0.35);
%%% ...not a proper outlier, but a value too much close to the UD_THRESHOLD = 0.4)
OUT_SD.high = find(SDs>Q3+ww);
if (numel(OUT_SD.high)~=0)
    disp('WARNING: SD outlier channels');
    disp(['(OUT_SD.high: ' num2str(OUT_SD.high) ')']);
end % These are the only ones that count for outliers

OUT = [OUT OUT_SD.high];


%%
%%%%%%%%%%%%%%%%
%   LogMuaRef 
%%%%%%%%%%%%%%%%
figure
stem(LogMuaRef);
set(gca, 'XLim', [0 33]) % limit set for the SD
xlabel('Channels')
ylabel('LogMuaRef')

%[GDB] there should be added a check if the fit procedure is not converging...
% --> NaN values for logMuaRef and SD

LogMuaRefsOK = ~isnan(LogMuaRef); % median and mean are computed only on non-NaN values
medRef = median(LogMuaRef(LogMuaRefsOK));
meanRef=mean(LogMuaRef(LogMuaRefsOK));
title(['[logMUA] Reference Value. Median = ' num2str(medRef,3) ]);

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0 0 5 4]);
print('-deps2c', 'LogMUARefVsChannels.eps');

% *** Find OUTLIERS in the LogmuaRef
Q1=prctile(LogMuaRef, 25); Q3=prctile(LogMuaRef, 75); IQR=Q3-Q1; ww=1.5*IQR;

OUT_Ref.low=find(LogMuaRef<Q1-ww);
if (numel(OUT_Ref.low)~=0) 
    disp(['(OUT_LogMuaRef.low: ' num2str(OUT_Ref.low) ')']);end

OUT_Ref.high = find(LogMuaRef>Q3+ww);
if (numel(OUT_Ref.high)~=0)
    disp(['(OUT_LogMuaRef.high: ' num2str(OUT_Ref.high) ')']);end

% [[LogMuaRef Outliers do not count as outliers]]


%% QUALITY PARAMETER: FREQUENCY
%
Freq = NaN(1,numel(RecordingSet));

for crsRecSet = ChannelSet % modified for running on a sub-set of recording channels (a cortical area)    
   workingdir = [AnalysisDir RecordingSet(crsRecSet).label];
   if exist([workingdir '/UpDownStateDuration.mat'], 'file') 
       load([workingdir '/UpDownStateDuration.mat']);
       Freq(crsRecSet) = 1/mean(UDcycleLen);
   else
       Freq(crsRecSet) = 0;
   end 
end

figure

stem(Freq)
set(gca, 'XLim', [0 33]); %set(gca, 'XLim', [0 33], 'YLim',[0 1.1*max(Freq)]) 
xlabel('Channels')
ylabel('Frequency of Transitions')

FreqOK = ~isnan(Freq);
medFreq = median(Freq(FreqOK));
meanFreq=mean(Freq(FreqOK));
title(['Median Frequency = ' num2str(medFreq,3) ]);

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0 0 5 4]);
print('-deps2c', 'FrequencyVsChannels.eps');

% *** Find OUTLIERS in the Frequency
Q1=prctile(Freq, 25); Q3=prctile(Freq, 75); IQR=Q3-Q1; ww=1.5*IQR;
OUT_F.high=find(Freq>Q3+ww);
if (numel(OUT_F.high)~=0)
    disp(['(OUT_F.high: ' num2str(OUT_F.high) ')']); end
OUT_F.low=find(Freq<Q1-ww);
if (numel(OUT_F.low)~=0)
    disp('WARNING: Frequency outlier channels');
    disp(['(OUT_F.low: ' num2str(OUT_F.low) ')']); end
OUT = [OUT OUT_F.low]; % These are the only ones that count for outliers


%%
OUT = [OUT RightPeak FewTrans]; % These are the only ones that count for outliers
OUT = unique(OUT);
% if isempty(OUT); OUT = [OUT 0]; end; % previously added for compatibility with Python io


%% [old version] QUALITY PARAMETER SD
% % % %
% % % SDs = NaN(1,numel(RecordingSet));
% % % UDTHR = NaN(1,numel(RecordingSet));
% % % 
% % % % *** Analysis loop on recording set...
% % % for crsRecSet = ChannelSet % modified for running on a sub-set of recording channels (a cortical area)
% % %    workingdir = [AnalysisDir RecordingSet(crsRecSet).label];
% % %    load([workingdir '/AnalysisSummary.mat'], 'ModeParams');
% % %    load([workingdir '/AnalysisSummary.mat'], 'UD_THRESHOLD');
% % %    PeakMean = Inf;
% % %    if isfield(ModeParams,'Mu1')
% % %        PeakMean = ModeParams.Mu1;
% % %        PeakStd = ModeParams.Sigma1;
% % %    elseif isfield(ModeParams,'Mu2')
% % %        PeakStd = ModeParams.Sigma2;
% % %    else
% % %        PeakStd = NaN;
% % %    end
% % % % THIS CAUSES 'fake' problems in the SDs stem when the second gaussian peak is not correctly recognised 
% % % %    if isfield(ModeParams,'Mu2')
% % % %       if PeakMean > ModeParams.Mu2
% % % %          PeakMean = ModeParams.Mu2;
% % % %          PeakStd = ModeParams.Sigma2;
% % % %       end
% % % %    end
% % %    SDs(crsRecSet) = PeakStd;
% % %    UDTHR(crsRecSet) = UD_THRESHOLD;
% % % end
% % % 
% % % % *** Check the stability of a MUA threshold across channels.
% % % %  In other words if the same threshold can be used for all the channels.
% % % %  If yes, this should provide more reliable profiles of traveling
% % % %  wavefronts.
% % % %
% % % figure
% % % stem(SDs);
% % % set(gca, 'XLim', [0 33], 'YLim',[0 0.35]) % limit set for the SD
% % % xlabel('Channels')
% % % ylabel('Down MUA SD')
% % % 
% % % SDsOK = ~isnan(SDs); % median and mean are computed only on non-NaN values
% % % MED = median(SDs(SDsOK));
% % % meanSD=mean(SDs(SDsOK));
% % % title(['Median SD = ' num2str(MED,3) ]);
% % % % [modified for running on a sub-set of recording channels (a cortical
% % % % area)]
% % % 
% % % % *** Find OUTLIERS in the SDs
% % % Q1=prctile(SDs, 25); Q3=prctile(SDs, 75); IQR=Q3-Q1; ww=1.5*IQR;
% % % % OUT=find(SDs>Q3+ww);
% % % % if (numel(OUT)~=0); 
% % % %     disp('WARNING: outlier channels');disp(OUT);
% % % 
% % % OUT_SD.low=find(SDs<Q1-ww);
% % % if (numel(OUT_SD.low)~=0) 
% % %     disp(['(OUT_SD.low: ' num2str(OUT_SD.low) ')']);end
% % % 
% % % %%%OUT_SD.high = find(SDs>0.35);
% % % %%% ...not a proper outlier, but a value too much close to the UD_THRESHOLD = 0.4)
% % % OUT_SD.high = find(SDs>Q3+ww);
% % % if (numel(OUT_SD.high)~=0)
% % %     disp('WARNING: SD outlier channels');
% % %     disp(['(OUT_SD.high: ' OUT_SD.high]);
% % % end % These are the only ones that count for outliers
% % % 
% % % OUT = [OUT OUT_SD.high];
% % % 
% % % set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0 0 5 4]);
% % % print('-deps2c', 'DownMUASDsVsChannels.eps');
% % % 
% % % %% Add the UD_THRESHOLD in the plot
% % % figure
% % % stem(SDs);
% % % hold on;
% % % stem(UDTHR);
% % % set(gca, 'XLim', [0 33], 'YLim',[0 max(UDTHR)+max(UDTHR)*0.1]) % limit set for the SD
% % % xlabel('Channels')
% % % ylabel('Down MUA SD and UD_THRESHOLD')
% % % set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0 0 5 4]);
% % % print('-deps2c', 'DownMUASD-and-UDTHR_vsChannels.eps');
% % % 
% % % %% QUALITY PARAMETER: FREQUENCY
% % % %
% % % Freq = NaN(1,numel(RecordingSet));
% % % 
% % % for crsRecSet = ChannelSet % modified for running on a sub-set of recording channels (a cortical area)    
% % %    workingdir = [AnalysisDir RecordingSet(crsRecSet).label];
% % %    if exist([workingdir '/UpDownStateDuration.mat'], 'file') 
% % %        load([workingdir '/UpDownStateDuration.mat']);
% % %        Freq(crsRecSet) = 1/mean(UDcycleLen);
% % %    else
% % %        Freq(crsRecSet) = 0;
% % %    end 
% % % end
% % % 
% % % figure
% % % 
% % % stem(Freq)
% % % set(gca, 'XLim', [0 33]); %set(gca, 'XLim', [0 33], 'YLim',[0 1.1*max(Freq)]) 
% % % xlabel('Channels')
% % % ylabel('Frequency of Transitions')
% % % 
% % % FreqOK = ~isnan(Freq);
% % % MED = median(Freq(FreqOK));
% % % meanFreq=mean(Freq(FreqOK));
% % % title(['Median Frequency = ' num2str(MED,3) ]);
% % % 
% % % % Find OUTLIERS in the Frequency
% % % Q1=prctile(Freq, 25); Q3=prctile(Freq, 75); IQR=Q3-Q1; ww=1.5*IQR;
% % % OUT_F.high=find(Freq>Q3+ww);
% % % if (numel(OUT_F.high)~=0)
% % %     disp(['(OUT_F.high: ' num2str(OUT_F.high) ')']); end
% % % OUT_F.low=find(Freq<Q1-ww);
% % % if (numel(OUT_F.low)~=0)
% % %     disp('WARNING: Frequency outlier channels');
% % %     disp(['(OUT_SD.high: ' OUT_F.low]);
% % % end
% % % OUT = [OUT OUT_F.low]; % These are the only ones that count for outliers
% % % 
% % % 
% % % % modified for running on a sub-set of recording channels (a cortical area)
% % % %title(['Median SD = ' num2str(median(SDs),3)])
% % % set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0 0 5 4]);
% % % print('-deps2c', 'FrequencyVsChannels.eps');
% % % 
% % % OUT = unique(OUT);
% % % 
% % % % if isempty(OUT); OUT = [OUT 0]; end; % previously added for compatibility with Python io


%% Save Results for each channel in the RecordingSet
%
for crsRecSet = ChannelSet % modified for running on a sub-set of recording channels (a cortical area)        
    workingdir = [AnalysisDir RecordingSet(crsRecSet).label];
%    disp(exist([workingdir '/AnalysisSummary.mat'], 'file'))
   if exist([workingdir '/AnalysisSummary.mat'], 'file')      
       results(crsRecSet)=load([workingdir '/AnalysisSummary.mat']);
   end
%   disp(exist([workingdir '/UpDownStateDuration.mat'], 'file'))
   if exist([workingdir '/UpDownStateDuration.mat'], 'file')
       duration(crsRecSet)=load([workingdir '/UpDownStateDuration.mat']);
   end   
   if exist([workingdir '/LogMUAUpTrans.Average.mat'], 'file')
       uptrans(crsRecSet)=load([workingdir '/LogMUAUpTrans.Average.mat']);
   end   
   if exist([workingdir '/LogMUADownTrans.Average.mat'], 'file')
       downtrans(crsRecSet)=load([workingdir '/LogMUADownTrans.Average.mat']);
   end   
end

save([AnalysisDir 'Results.mat'], 'results', 'OUT');
save([AnalysisDir 'Duration.mat'], 'duration', 'OUT');
save([AnalysisDir 'UpwardTrans.mat'], 'uptrans', 'OUT');
save([AnalysisDir 'DownwardTrans.mat'], 'downtrans', 'OUT');
save([MUAdir 'outliers.mat'], 'OUT');

disp('>> >> >> >> >> >> >> >> >> >> >> >> >> >> >> >> >> >> >> >> >> >>');
disp(['...a total of ' num2str(numel(Disc)) ' channels with DISCONTINUITY']); 
if (numel(Disc))~=0; disp(Disc); end

%%% logMUA distribution
disp(['...a total of ' num2str(numel(WeakBi)) ' channels with WeakBimodality']); 
if (numel(WeakBi))~=0; disp(WeakBi); end
disp(['...a total of ' num2str(numel(RightPeak)) ' channels with the largest peak on the RIGHT (negative asymmetry)']); 
if (numel(RightPeak))~=0; disp(RightPeak); end
disp(['...a total of ' num2str(numel(SkewP)) ' channels with positive StrongAsymmetry (second peak)']); 
if (numel(SkewP))~=0; disp(SkewP); end
disp(['...a total of ' num2str(numel(SkewN)) ' channels with negative StrongAsymmetry (second peak)']); 
if (numel(SkewN))~=0; disp(SkewN); end
disp(['...a total of ' num2str(numel(LargeUDTH)) ' channels with UD_THRESHOLD beyond the Tail.MU']); 
if (numel(LargeUDTH))~=0; disp(LargeUDTH); end
disp(['...a total of ' num2str(numel(FitWarning)) ' channels with FitWarning (standardized Gaussian)']); 
if (numel(FitWarning))~=0; disp(FitWarning); end

disp(['...a total of ' num2str(numel(FewTrans)) ' channels with FewTrans']); 
if (numel(FewTrans))~=0; disp(FewTrans); end

disp(['...a total of ' num2str(numel(OUT)) ' outlier channels found']); 
if (numel(OUT))~=0; disp(OUT); end


catch
   msg = lasterror;
   fprintf('%s', msg.message)
   cd('..');
   if ~isempty(hSMR)
      closeSON(hSMR);
   end
   rethrow(lasterror);
end


%% Close file

if ~isempty(hSMR)
   closeSON(hSMR);
end

exit
