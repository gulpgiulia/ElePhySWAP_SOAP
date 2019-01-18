% plotMapOfUpOnsetTimeLags_2D.m
%
%   Copyright 2017 Maurizio Mattia @ Ist. Super. Sanita', Rome - Italy
%   Version: 1.0 - Jan. 18, 2017
%
% [GDB]

% (v1) First version of the new release with improvements introduced to fit
% the experimental setup of the 32-electrode array. In partciular:
% - output path for results (files and plots)
% - treatment of outlier channels
% - WaveCollection and cumulative function of the number of channels
% (GLOBALITY check)
% - inspection of suppressed channels/transitions in non-unique waves and
% identificaton of CRITICAL WAVES
% - optimization of MIN_CH_NUM (~ half of the array) and of NUM_OF_CLUSTERS
% - borders of CorticalAreas in ContourPlots
%
% (v2) Solutions for fixing the CRITICAL WAVES.
% - optimize the IWI with the median of the most active channels
% In addition:
% - MIN_CH_NUM set at 8 (replace the "globality" requirement with the
% observation of the distribution of WaveCollection(2)) --> it should
% preserve "local" waves, i.e. waves propagating locally in cortical areas,
% for instance anterior areas (motor + somatosensory), or visual area
% - ...maybe NUM_OF_CLUSTERS set at 4 for all the experiments
% 


%% Load analysis params and information about recording dataset.
%
setParamsAndOptions
version = 'v2';


%% [GDB] Set the Path for the results
%
OutputDir = [AnalysisDir 'WAVES_' version '/'];
if exist(OutputDir,'dir') == 0
    mkdir(OutputDir);
end

FilesToDelete = [OutputDir 'ClusteredWave_*'];
if ~isempty(dir([OutputDir 'ClusteredWave_*']))
     delete(FilesToDelete);
    % (to be sure that the new execution replace the previous results)
end

fid = fopen([OutputDir FileName '.out'], 'wt');
fprintf(fid,'\n*** FileName: %s\n', FileName);

%% [GDB] EXCLUDE the OUTLIER CHANNELS
%
nCh=32;
%ChannelSet defined in SetParametersAndOptions
load([AnalysisDir 'Results.mat']);
if ~isempty(OUT)
    nCh=nCh-numel(OUT);
    ChannelSet= setdiff(ChannelSet, OUT);
end
disp(['Number of Channels: ' num2str(nCh)])
fprintf(fid,'Number of Channels: %d\n', nCh);


%% Setting additional constants used throughout the script (colormap)
%
MAX_ABS_TIMELAG = Options.UpTimeLags.MaxAbsTimeLag;

TL_RANGE = [-300 300];

RED_BLUE = 1;
CM = jet();
if RED_BLUE
   CM = gray(32);
   CM(:,3) = 1;
   Red = flipud(gray(32));
   Red(2:end,1) = 1;
   CM = [CM; Red(2:end,:)];
end
CM(1,:) = [1 1 1]*0.25;


%% LOAD the detected upward transitions for the channels...
%
UpTrans = [];
ChLabel = [];
%TransPerCh = zeros(1,numel(RecordingSet));
TransPerCh = zeros(1,numel(ChannelSet)); % to take into account the case of a subset of electrodes
for crsRecSet = ChannelSet
   disp(['Recording set: ' RecordingSet(crsRecSet).label]);
   fprintf(fid,'Recording set: %s\n',RecordingSet(crsRecSet).label);
   ldir = [AnalysisDir RecordingSet(crsRecSet).label '/'];
   load([ldir 'UpwardTransitions.mat']);
   TickLabels{crsRecSet} = RecordingSet(crsRecSet).label(1:2);
   
   UpTrans = [UpTrans Triggers];
   ChLabel = [ChLabel repmat(crsRecSet, size(Triggers))];
   TransPerCh(crsRecSet) = numel(Triggers);
end

[UpTrans,ndx] = sort(UpTrans);
ChLabel = ChLabel(ndx);

figure
plot(UpTrans,ChLabel,'k.');
set(gca, 'YLim',[0 33]); % (for a better-looking representation in case of a subset of electrodes)
xlabel('Time (s)');
ylabel('Channel');
print('-deps2c', [OutputDir 'SequenceOfUpTrans.eps']);  % path for the results


%% ExpectedTrans i.e. estimate of the Number of Waves
%ExpectedTrans used to estimate/optimize IWI

%ExpectedTrans = median(TransPerCh);
%"TransPerCh(find(TransPerCh))" to account for the case of a subset of electrodes

eExpectedTrans(1) = median(TransPerCh(find(TransPerCh)));

[Active,nA] = sort(TransPerCh(find(TransPerCh)),'descend');
eExpectedTrans(2) = median(Active(1:8)); % consider THE 8 MOST ACTIVE CHANNELS only

eExpectedTrans(3) = mean(TransPerCh(find(TransPerCh)))+std(TransPerCh(find(TransPerCh)));
% mix median and std? maybe mean...?

sel1=1;
ExpectedTrans=eExpectedTrans(sel1);

fprintf('\n--- Estimate of Expected Transitions ---');
for i =1:1:numel(eExpectedTrans)
    if i==sel1; comment='******'; else comment=''; end;
    fprintf('\nExpectedTrans(%d): %d\t%s', i,round(eExpectedTrans(i)),comment);
end
fprintf('\n\n')


%% Compute the time lags matrix looking for optimal MAX_ABS_TIMELAG...
%  (depends on the distance between electrodes in the array)
%
DeltaTL = diff(UpTrans);
OneMoreLoop = 1;
while OneMoreLoop
   WnW = DeltaTL<=MAX_ABS_TIMELAG; % parameter initialized in SetParametersAndOptions
   ndxBegin = find(WnW==1,1,'first');
   nw = 0;
   clear('Wave','WaveUnique','WaveSize','WaveTime');
   while ndxBegin < numel(DeltaTL)
      ndxEnd = find(WnW(ndxBegin:end)==0,1,'first')+ndxBegin-1; % isolated transitions are EXCLUDED
      if isempty(ndxEnd)
         ndxEnd = numel(UpTrans);
      end
      nw = nw + 1;
      Wave(nw).ndx = ndxBegin:ndxEnd;
      WaveUnique(nw) = numel(Wave(nw).ndx) == numel(unique(ChLabel(Wave(nw).ndx)));
      WaveSize(nw) = numel(Wave(nw).ndx);
      WaveTime(nw) = mean(UpTrans(Wave(nw).ndx));
      if ndxEnd == numel(UpTrans)
         ndxBegin = numel(DeltaTL);
      else
         ndxBegin = find(WnW(ndxEnd:end)==1,1,'first')+ndxEnd-1;
      end
   end
   OneMoreLoop = 0;
   if min(WaveUnique) == 0 % the first scan requests that all the waves are UNIQUE
      fprintf('MaxAbsTimelag too large: %f -> %f\n',MAX_ABS_TIMELAG,MAX_ABS_TIMELAG*0.75);
      fprintf(fid,'MaxAbsTimelag too large: %f -> %f\n',MAX_ABS_TIMELAG,MAX_ABS_TIMELAG*0.75);
      MAX_ABS_TIMELAG = MAX_ABS_TIMELAG*0.75;
      OneMoreLoop = 1;
   else
%       if max(WaveSize) < numel(RecordingSet)
%          fprintf('MaxAbsTimelag too small: %f -> %f\n',MAX_ABS_TIMELAG,MAX_ABS_TIMELAG*1.25);
%          MAX_ABS_TIMELAG = MAX_ABS_TIMELAG*1.125;
%          OneMoreLoop = 1;
%       end
   end
end % while OneMoreLoop


%% ...a closer look at isolated transitions (i.e. non-waves)
totTrans=numel(UpTrans);
totWaveTrans=sum(WaveSize); 
%(totWaveTrans is always less than totTrans: isolated transitions are excluded)

% %OR
% totWaveTrans=0;
% for i=1:1:numel(Wave)
%     totWaveTrans=totWaveTrans+numel(Wave(i).ndx)
% end


%% Wave Hunt -- step 1: Compute the Time Lags, Find Unique Waves --> WaveCollection1
% [GDB]
%
WaveColl(1).nWaves=length(WaveSize);

X=[1:1:nCh];
dX = X(2) - X(1);
for i=X
    WaveColl(1).NumberOf(i)=length(find(WaveSize==i));
end

figure

subplot(1,2,1)
%plot([1:nCh],WaveColl(1).NumberOf);
%[XX, YY] = stairs(X - dX/2, WaveColl(1).NumberOf/WaveColl(1).nWaves*100); %percentage
[XX, YY] = stairs(X - dX/2, WaveColl(1).NumberOf); %absolute NumberOf
XX = [XX(1) XX' XX(end)+dX XX(end)+dX];
YY = [0 YY' YY(end) 0];
patch(XX, YY, 'b', 'FaceColor', 'b', 'EdgeColor', 'b', 'FaceAlpha', 0.5);
set(gca,'XLim',[0 33]);
xlabel('Number of Channels involved');
ylabel('Number of Waves in the Collection');
hold on
yyaxis right
plot(cumsum(WaveColl(1).NumberOf/WaveColl(1).nWaves))
dim0 = [0.2175 0.825 0.1 0.1];
str0 = {['nWaves = ' num2str(WaveColl(1).nWaves)]};
annotation('textbox',dim0,'String',str0,'FitBoxToText','on','BackgroundColor','none');

subplot(1,2,2)
patch(XX, YY, 'b', 'FaceColor', 'b', 'EdgeColor', 'b', 'FaceAlpha', 0.5);
value=sum(WaveColl(1).NumberOf(ceil(nCh/2):nCh))/WaveColl(1).nWaves;
set(gca,'XLim',[ceil(nCh/2) 33]);
xlabel('Number of Channels involved');
ylabel('Number of Waves in the Collection');

dim1 = [0.65 0.825 0.1 0.1];
% Size and location, specified as a four-element vector of the form [x y w h]. 
% The first two elements specify the coordinates of the lower left corner of the text box.
% The second two elements specify the width and height of the annotation, respectively.
% By default, the units are normalized to the figure. 
% The lower left corner of the figure maps to (0,0) and the upper right corner maps to (1,1). 
str1 = {'Percentage of waves involving more', ['than half of the electrode array: ' num2str(value)]};
% Each element of the cell array displays on a separate line
annotation('textbox',dim1,'String',str1,'FitBoxToText','on','BackgroundColor','white');

dim2 = [0.1 0.9 0.1 0.1];
str2 = ('WavesCollection 1 - Sequence of UpTrans and Optimization of MaxAbsTimeLag');
annotation('textbox',dim2,'String',str2,'FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none');
print('-deps2c', [OutputDir 'WavesCollection1.eps']);  % path for the results

%title({'WaveCollection 1','Sequence of UpTrans and Optimization of MaxAbsTimeLag'});


%% Estimate of IWI (NOT USED)
% Using the distribution of DeltaTL i.e. time lag between consecutive transitions
[hTL,Xh] = hist(DeltaTL,100);

percentIW=ceil(ExpectedTrans*100/numel(DeltaTL)); % estimate (percentage) how many transitions are Inter-Waves

indexIW=find((cumsum(hTL)/numel(DeltaTL))>(100-percentIW)/100,1); 
% indexIW estimates at which bin of the hist(DeltaTL) enters the contribution of inter-wave-intervals
% (rather than inter-transition-intervals)

eIWI(1)=Xh(indexIW); % IWI= the first value at which the DeltaTL distribution starts to represent inter-waves

wmDelta=sum(Xh(1:indexIW-1).*hTL(1:indexIW-1))/sum(hTL(1:indexIW-1));
% weighted mean of the DeltaTL that should represent inter-transition-intervals 
% (rather than inter-wave-intervals)
eIWI(2)=wmDelta*8; % weighted average (assuming MIN_CH_NUM= 8)
eIWI(3)=Xh(1)*8; % minimum estimate 

%eIWI = 3 different estimates of "minimum" IWI

sel2=1;
EIWI=eIWI(sel2);

fprintf('\n--- Estimate of InterWaveInterval ---');
for i =1:1:numel(eIWI)
    if i==sel2; comment='******'; else comment=''; end;
    fprintf('\nIWI(%d): %f\t%s', i,eIWI(i),comment);
end
fprintf('\n\n')


%% IWI distribution and histogram
%
IWI = diff(WaveTime);          % IWI = Inter-Wave-Interval...

%%% histogram(X,edges) sorts X into bins with the bin edges specified by 
%%% the vector, edges. Each bin includes the left edge, but does not include 
%%% the right edge, except for the last bin which includes both edges.
figure
hIWI=histogram(IWI,[0:0.250:round(max(IWI))]); 
% 250ms is the MINIMUM IWI to be compatible with a frequency within 4Hz
eWC2=sum(hIWI.Values(2:end)); % estimates the number of slow waves (frequency <4Hz)
% This number should be comparable with the number of elements in WaveCollection(2)

%% Recollect small-waves in full waves (involving a wider area).

MAX_IWI = max(IWI)*0.5;
%MAX_IWI = median(IWI); % ... another estimate (NOT TAKING INTO ACCOUNT the 4Hz LIMIT)
medianIWI = median(IWI(find(IWI>=0.250))); % estimate TAKING INTO ACCOUNT the 4Hz LIMIT

%%ExpectedTrans=eWC2;
ACCEPTABLE_REJECTION_RATE = 0.1;

OneMoreLoop = 1;
while OneMoreLoop
   WnW = IWI<=MAX_IWI;
   ndxBegin = find(WnW==1,1,'first');
   nw = 0;
   clear('FullWave','FullWaveUnique','FullWaveSize','FullWaveTime','FullWaveUniqueSize');
   while ndxBegin < numel(IWI)
      ndxEnd = find(WnW(ndxBegin:end)==0,1,'first')+ndxBegin-1;
      if isempty(ndxEnd)
         ndxEnd = numel(WaveTime);
      end
      nw = nw + 1;
      FullWave(nw).ndx = [Wave(ndxBegin:ndxEnd).ndx];
      FullWaveUnique(nw) = numel(FullWave(nw).ndx) == numel(unique(ChLabel(FullWave(nw).ndx)));
      FullWaveSize(nw) = numel(FullWave(nw).ndx);
      FullWaveUniqueSize(nw) = numel(unique(ChLabel(FullWave(nw).ndx)));
      FullWaveTime(nw) = mean(UpTrans(FullWave(nw).ndx));
      if ndxEnd == numel(WaveTime)
         ndxBegin = numel(IWI);
      else
         ndxBegin = find(WnW(ndxEnd:end)==1,1,'first')+ndxEnd-1;
         for j = ndxEnd+1:ndxBegin-1
            nw = nw + 1;
            FullWave(nw).ndx = Wave(j).ndx;
            FullWaveUnique(nw) = numel(FullWave(nw).ndx) == numel(unique(ChLabel(FullWave(nw).ndx)));
            FullWaveSize(nw) = numel(FullWave(nw).ndx);
            FullWaveUniqueSize(nw) = numel(unique(ChLabel(FullWave(nw).ndx)));
            FullWaveTime(nw) = mean(UpTrans(FullWave(nw).ndx));
         end
      end
   end
   BadWavesNum = numel(find(FullWaveUnique==0));
   fprintf('Max wave size: %d, Bad waves: %d (%.3g%%), Good waves: %d\n',...
       max(FullWaveSize),BadWavesNum,BadWavesNum/numel(FullWaveUnique)*100,numel(FullWaveUnique)-BadWavesNum);
   fprintf('Num. waves: %d, max. non-unique ch.: %d\n',...
       numel(FullWaveUnique),max(FullWaveSize-FullWaveUniqueSize));
   fprintf(fid,'Max wave size: %d, Bad waves: %d (%.3g%%), Good waves: %d\n',...
       max(FullWaveSize),BadWavesNum,BadWavesNum/numel(FullWaveUnique)*100,numel(FullWaveUnique)-BadWavesNum);
   fprintf(fid,'Num. waves: %d, max. non-unique ch.: %d\n',...
       numel(FullWaveUnique),max(FullWaveSize-FullWaveUniqueSize));
   
   OneMoreLoop = 0;
   if numel(FullWaveUnique) <= ExpectedTrans % If not we have an artifactual amplification of small waves...
      if BadWavesNum/numel(FullWaveUnique) > ACCEPTABLE_REJECTION_RATE;
         if min(FullWaveUnique) == 0 % at lest a Wave non-unique
            fprintf('MaxIWI too large: %f -> %f\n',MAX_IWI,MAX_IWI*0.75);
            fprintf(fid,'MaxIWI too large: %f -> %f\n',MAX_IWI,MAX_IWI*0.75);   
            MAX_IWI = MAX_IWI*0.75;
            OneMoreLoop = 1;
         else % only unique waves
            %if max(WaveSize) < numel(RecordingSet) % at least one wave MUST BE GLOBAL (i.e. involving the full set of electrodes)
            if max(WaveSize) < numel(ChannelSet) % at least one wave MUST BE GLOBAL (i.e. involving the full set of electrodes)
               fprintf('MaxIWI too small: %f -> %f\n',MAX_IWI,MAX_IWI*1.25);
               fprintf(fid,'MaxIWI too small: %f -> %f\n',MAX_IWI,MAX_IWI*1.25); 
               MAX_IWI = MAX_IWI*1.25;
               OneMoreLoop = 1;
            end
         end
      end %else...no more loop
   end %else...no more loop
end % while OneMoreLoop

fprintf('\nExpected waves     : %d\n', round(ExpectedTrans));
fprintf('Reconstructed waves: %d\n', numel(FullWaveUnique)); 
fprintf(fid,'\nExpected waves     : %d\n', round(ExpectedTrans));
fprintf(fid,'Reconstructed waves: %d\n', numel(FullWaveUnique)); 
%... but still several non-unique waves, i.e. waves with repeated channels

totFullWaveTrans=sum(FullWaveSize); 


%% Remove from waves nonunique channels...
%

FullWaveBackup = FullWave;
TOTNonUnique=numel(find(FullWaveUnique==0));

fprintf('\tWaves to clean     : %d\n', numel(find(FullWaveUnique==0)));
fprintf('\t(Max. non-unique ch.: %d)\n', max(FullWaveSize-FullWaveUniqueSize));
fprintf(fid,'\tWaves to clean     : %d\n', numel(find(FullWaveUnique==0)));
fprintf(fid,'\t(Max. non-unique ch.: %d)\n', max(FullWaveSize-FullWaveUniqueSize));

listOfnw=[]; % List of "critical" waves: non-unique but possibly hiding multiple unique waves
nAdded=0; % total number of new waves added to the collection 
          % (useful also to point back to the original wave indexing in FullWaveBackup)
nNew=[];  
nSegments=[]; % (how many segments non-unique waves are made of)
nSize=[]; % (size of the segments)

alert(1:6)=0;

%THR = 3*wmDelta;

% *** TWO THRESHOLDS ***
% THR1 = maxRisingTime, estimated as 1/2 of mean(UpStateDuration)
% THR2 = 250ms = min IWI (frequency limit) 

meanUP=0;
for i = ChannelSet
    meanUP=meanUP+results(i).UpState.mean_duration;
end
meanUP=meanUP/numel(ChannelSet);
THR1 = 0.5*meanUP; %RisingTime
THR2 = 0.250; %4Hz LIMIT

% Save plot of "critical" non-unique waves, before & after the treatment
% Create the folder for the storage of plots if non-existing, and empty if existing
FilesToDelete = [OutputDir 'nW/nW*'];
if exist([OutputDir 'nW/'],'dir') == 0
    mkdir([OutputDir 'nW/']);
elseif ~isempty(dir([OutputDir 'nW/']))
     delete(FilesToDelete);    
end

nw=1;
NonUnique=0; % COUNTER for NonUnique Waves
while nw<=numel(FullWaveUnique)
%for nw = 1:numel(FullWaveUnique) 
    % numel(FullWaveUnique) = number of waves in the collection
    % FullWaveUnique==0 if the waves contains repeated channels
    
    NewFullWave=[];NewFullWaveUnique=[];NewFullWaveSize=[];NewFullWaveTime=[];seg=[]; % empty arrays
    nCHS=[];
    
    if FullWaveUnique(nw) == 0
        NonUnique=NonUnique+1;
              
        chs = ChLabel(FullWave(nw).ndx);
        nCHS(1)=numel(chs);
        wndx = 1:1:numel(chs);
% FullWave(nw).ndx = index of the transitions that make up the wave(nw)
% ChLabel() = the corresponding channel for each transition 
% --> in this list there are repeated channels
        nch = hist(chs,1:numel(RecordingSet)); % WARNING! not numel(ChannelSet) even in the case of a subset
 
% --- (1) First CHANNEL CLEANING: check the TIME DISTANCE BETWEEN REPETITIONS
        rep=find(nch>1);
        
        k=1;
        while k<=numel(rep)
        %for i=rep
            i=rep(k);
        
            timeStamp=UpTrans(FullWave(nw).ndx(chs==i));
            idx=FullWave(nw).ndx(chs==i);
            if diff(timeStamp)<0.125
                % open a time-window around each occurrence of i and check
                % how many members of the clan are in the time-window 
                fprintf('Channel %d has close repetitions\n',i);
                nClan=[];
                for j=1:1:nch(i)
                    window=[timeStamp(j)-0.125 timeStamp(j)+0.125];
                    %window=[timeStamp(j)-0.125 timeStamp(j)+0.125];
                    inWindow=FullWave(nw).ndx(UpTrans(FullWave(nw).ndx)>window(1) & UpTrans(FullWave(nw).ndx)<window(2)); 
                    inWindow=setdiff(inWindow,idx(j));
                    nClan(j)=numel(intersect(ChLabel(inWindow),arrayMask{i}{2}));
                end
                
                if ~isempty(find(nClan==0)) % to be more generic, consider the case min(nClan)
                    [FullWave(nw).ndx,a]=setdiff(FullWave(nw).ndx,idx(find(nClan==0)));
                    b=setdiff(wndx,a);
                    chs(b)=[];
                    fprintf('%d repetition(s) deleted\n',numel(b));
                    nch=hist(chs,1:numel(RecordingSet));
                    if(nch(i)>1) 
                        fprintf('...Now, re-do channel %d\n',i);
                    else
                        k=k+1;
                    end
                else
                    fprintf('CloseRepProblem for channel i = %d NOT SOLVED\n',i);
                    k=k+1;
                end
            else
                k=k+1;
            end
            
        end
        nCHS(2)=numel(chs); % update the number of transitions (to keep count of removed transitions)

        
        
        
        
%[GDB] Check the time difference between transitions 
        delta=diff(UpTrans(FullWave(nw).ndx));
        mD=mean(delta);
        stdD=std(delta);
       
% --- PLOT non-unique wave, in order to CHECK how it is treated ---
        f=figure('visible','off');
        
        % --- deltaHistogram
        subplot(1,4,1) 
        nBin=round(numel(delta)/3);
        hDelta=histogram(delta,nBin);
        xlabel('delta (s)');
        hold on
        line([mD mD],[0 max(hDelta.Values)],'Color','r')
        line([mD+stdD mD+stdD],[0 max(hDelta.Values)],'Color',[0.75 0 0],'LineStyle','--')
        line([mD+2*stdD mD+2*stdD],[0 max(hDelta.Values)],'Color',[0.5 0 0],'LineStyle','--')
        line([mD+3*stdD mD+3*stdD],[0 max(hDelta.Values)],'Color',[0.25 0 0],'LineStyle','--')
        hold off
        set(gca,'Position',[0.05 0.11 0.2 0.8150]); %[left bottom width height]
               
        % --- UpTrans
        subplot(1,4,2:4)        
        plot(UpTrans(FullWave(nw).ndx),ChLabel(FullWave(nw).ndx),'k.')
        set(gca, 'YLim',[0 33]);
        hold on
        xlabel('Time (s)');
        ylabel('Channel');
        set(gca,'YLim',[0 33]); 
        set(gca,'Position', [0.35 0.11 0.6 0.8150]); %[left bottom width height]
        
        dim0 = [0.735 0.9 0.1 0.1];
        str0 = {['nW = ' num2str(nw) '   [ ' num2str(NonUnique) '/' num2str(TOTNonUnique) ' ]']};
        annotation('textbox',dim0,'String',str0,'FitBoxToText','on','BackgroundColor','none');

%             XRange = [min(UpTrans(FullWaveBackup(i).ndx)) min(UpTrans(FullWaveBackup(i).ndx))+2*MAX_IWI];
%             XRange = XRange + [-1 +1]*diff(XRange)*0.1;
%             set(gca,'XLim',XRange)            
%             print('-dpdf', [OutputDir 'nW/nW' num2str(i)])
% -----------------------------------------------------------------               


% Look for CANDIDATE SEGMENTS  
        THR=mD+3*stdD;
        B=find(delta>THR); % larger sepration between transitions 
                           % --> possible BREAK of the non-unique wave
        
        if ~isempty(B) % --- IDENTIFY SEGMENTS
            % SEGMENTS i.e. candidate NEW waves
            segments=numel(B)+1;
            istart=1;
            i=1;
            while i<segments
                seg(i).ndx=(istart:B(i));
                seg(i).chs=chs(istart:B(i));
                istart=B(i)+1;
                i=i+1;
            end
            seg(i).ndx=(istart:numel(chs));
            seg(i).chs=chs(istart:end);

            % --- CLEAN SEGMENTS---
            delSeg=[];
            for i=1:1:segments % CHECK if SEGMENTS have to be cleaned
                
                % 1) non-unique segments
                if(numel(unique(seg(i).chs))~=numel(seg(i).chs)) % the subset of transition is NOT unique 
                    % --> 'CLEAN', i.e. scan the transition sequence and
                    % find repetition of channels
                    fprintf('ALERT(1) nw=%d, seg=%d\n', nw,i);
                    fprintf(fid,'ALERT(1) nw=%d, seg=%d\n', nw,i);
                    alert(1)=alert(1)+1;
                    
                    %j=1;
                    %while j<=numel(seg(i).chs) % while-loop beacause the numel(chs) is updated in the loop
                    delCh=[];
                    for j=1:numel(seg(i).chs)
                        listNdx=(find(seg(i).chs==seg(i).chs(j)));
                        if(numel(listNdx)~=1)
%                         % Keep the first occurrance, remove the others (interpreted as noise)
%                             seg(i).chs(listNdx(2:end))=[];
%                             seg(i).ndx(listNdx(2:end))=[];

%                         % Keep the occurrance which is the closest to the other occurances in the "clan"
                            t0Clan=0; nClan=0;
                            for k = arrayMask{seg(i).chs(j)}{2}
                                %disp(k)
                                %disp(UpTrans(FullWave(nw).ndx(seg(i).ndx(seg(i).chs==k))))
                                tClan=UpTrans(FullWave(nw).ndx(seg(i).ndx(seg(i).chs==k)));
                                if ~isempty(tClan)
                                    nClan=nClan+1;
                                    t0Clan=t0Clan+tClan;   
                                end                               
                            end
                            t0Clan=t0Clan/nClan; 
                            tCh=UpTrans(FullWave(nw).ndx(seg(i).ndx(seg(i).chs==seg(i).chs(j))));
                            [t0Ch,index]=min(abs(tCh-t0Clan));
                            delCh=[delCh setdiff(listNdx,listNdx(index))];                            
                        end 
                        %j=j+1;
                    end 
                    seg(i).chs(unique(delCh))=[];
                    seg(i).ndx(unique(delCh))=[];
                end 
                
                % 2) channels non-LocallyConnected (see arrayMask)
                if numel(seg(i).chs)<=5 % 5 is min(numel(arrayMask{:}{2})
                    %for j=seg(i).chs
                    delList=[];
                    for j=1:1:numel(seg(i).chs)
                        k=seg(i).chs(j);
                        %disp(['k = ' num2str(k)])
                        %disp(['intersect: ' num2str(intersect(arrayMask{k}{2},setdiff(seg(i).chs,k)))])
                        if isempty(intersect(arrayMask{k}{2},setdiff(seg(i).chs,k)))
                          delList=[delList j];  
                        end 
                    end
                    if ~isempty(delList)
                        fprintf('ALERT(2) nw=%d, seg=%d\n', nw,i);
                        alert(2)=alert(2)+1;
                        seg(i).chs(delList)=[];
                        seg(i).ndx(delList)=[];
                    end
                end
                
                % PREPARE TO REMOVE EMPTY SEGMENTS
                if isempty(seg(i).ndx)
                   fprintf('ALERT(3) nw=%d, seg=%d\n', nw,i);
                   fprintf(fid,'ALERT(3) nw=%d, seg=%d\n', nw,i);
                   alert(3)=alert(3)+1;
                   delSeg=[delSeg i];
                end
            
            end
                
            % REMOVE EMPTY SEGMENTS
            if ~isempty(delSeg)
                seg(delSeg)=[];
            end
            segments=numel(seg); % update the value in 'segments' = number of segments

            % coalescence of segments if no repetitions with adjacent(s) one(s)
            % N.B. a "small" intersection is admitted
            %for i=1:1:(segments-1)   
            i=1;
            while i<=(segments-1)
                %if isempty(intersect(seg(i).chs,seg(i+1).chs))
                if numel(intersect(seg(i).chs,seg(i+1).chs))<=floor(1/4*min(numel(seg(i).chs),numel(seg(i+1).chs)))
                    % CANDIDATE SEGMENTS for COALESCENCE
                    % check also if distance between segments'border is smaller than 250ms = 1/4Hz
                    distance=UpTrans(FullWave(nw).ndx(seg(i+1).ndx(1))) - UpTrans(FullWave(nw).ndx(seg(i).ndx(end)));
                    
                    if distance>=0.250
                        % FREQUENCY ALERT: distance compatible with SWA, the two segments should be kept separated
                        fprintf('>>> ALERT(4) (Frequency Alert)\n nw=%d, seg=%d AND seg=%d\n', nw,i,i+1 );
                        fprintf(fid,'>>> ALERT(4) (Frequency Alert)\n nw=%d, seg=%d AND seg=%d\n', nw,i,i+1 );
                        alert(4)=alert(4)+1;
                        disp(['--- ' num2str(seg(i).chs)]);
                        disp(['--- ' num2str(seg(i+1).chs)]);
                        fprintf('TOT Transitions: %d\n',numel(seg(i).chs)+numel(seg(i+1).chs));
                        disp(['intersect: ' num2str(intersect(seg(i).chs,seg(i+1).chs))]);
                        fprintf('numel(intesect): %d\n',numel(intersect(seg(i).chs,seg(i+1).chs)));
                        fprintf('threshold: %d\n',floor(1/4*min(numel(seg(i).chs),numel(seg(i+1).chs))));
                        i=i+1; %increment the pointer only if no coalescence is made
                    else
                        % COALESCENCE
                        % The two segments are close enough that can be merged into a single wave
                        % (consider them separated waves would mean the SWA frequency is larger than 4Hz)
                        
                        disp('*** Empty/small intersection for consecutive segments ==> coalescence of segments');
                        fprintf('ALERT(5) nw=%d, seg=%d AND seg=%d\n', nw,i,i+1);
                        fprintf(fid,'*** Empty/small intersection for consecutive segments ==> coalescence of segments\n');
                        fprintf(fid,'ALERT(5) nw=%d, seg=%d AND seg=%d\n', nw,i,i+1);
                        alert(5)=alert(5)+1;
                        disp(['--- ' num2str(seg(i).chs)]);
                        disp(['--- ' num2str(seg(i+1).chs)]);
                        fprintf('TOT Transitions: %d\n',numel(seg(i).chs)+numel(seg(i+1).chs));
                        disp(['intersect: ' num2str(intersect(seg(i).chs,seg(i+1).chs))]);
                        fprintf('numel(intesect): %d\n',numel(intersect(seg(i).chs,seg(i+1).chs)));
                        fprintf('threshold: %d\n',floor(1/4*min(numel(seg(i).chs),numel(seg(i+1).chs))));

                        % COALESCENCE of consecutive SEGMENTS
                        mergedCHS=[seg(i).chs seg(i+1).chs];
                        mergedNDX=[seg(i).ndx seg(i+1).ndx];
                        
                        % CHECK for REPETITIONS (and treat them as usual...
                        % looking at the meanTime in the Clan)
                        delCh=[];
                        for j=1:numel(mergedCHS)
                            listNdx=(find(mergedCHS==mergedCHS(j)));
                            if(numel(listNdx)~=1)
    %                         % Keep the first occurrance, remove the others (interpreted as noise)
    %                             seg(i).chs(listNdx(2:end))=[];
    %                             seg(i).ndx(listNdx(2:end))=[];

    %                         % Keep the occurrance which is the closest to the other occurances in the "clan"
                                t0Clan=0; nClan=0;
                                for k = arrayMask{mergedCHS(j)}{2}
                                    %disp(k)
                                    %disp(UpTrans(FullWave(nw).ndx(seg(i).ndx(seg(i).chs==k))))
                                    tClan=UpTrans(FullWave(nw).ndx(mergedNDX(mergedCHS==k)));
                                    if ~isempty(tClan)
                                        nClan=nClan+1;
                                        t0Clan=t0Clan+tClan;
                                    end
                                end
                                t0Clan=t0Clan/nClan;
                                tCh=UpTrans(FullWave(nw).ndx(mergedNDX(mergedCHS==mergedCHS(j))));
                                [t0Ch,index]=min(abs(tCh-t0Clan));
                                delCh=[delCh setdiff(listNdx,listNdx(index))];                            
                            end 
                        end 
                        
                        mergedCHS(unique(delCh))=[];
                        mergedNDX(unique(delCh))=[];
                        % N.B. unique(delCh) is empty if delCh is empty -->
                        % no cancellation is made in mergedCHS and mergedNDX
                        
                        seg(i).chs=mergedCHS; % COALESCENCE
                        seg(i).ndx=mergedNDX; % COALESCENCE  
                        seg(i+1)=[]; % coalesced segments are at index i, segment at index i+1 is REMOVED 
                        segments=segments-1;
                        
                        % if segments are merged, do not increment the pointer but
                        % check if the new segment can be merged with the following one
                                                                     
                    end % END IF (FREQUENCY ALERT)
                    
                else % consecutive segments intersect too much...
                    i=i+1; %increment the pointer only if no coalescence is made
                end
            end
                    
            if(segments~=numel(seg)) disp('ERROR - Number of Segments');end; % (for doubt's sake)
                
            % $$$$$ N.B. the number of segments has to be updated    
            for i=1:1:segments      
                NewFullWave(i).ndx = FullWave(nw).ndx(seg(i).ndx);
                NewFullWaveUnique(i) = 1; % update the logical value (...we are "cleaning" the waves)
                NewFullWaveSize(i) = numel(NewFullWave(i).ndx);
                NewFullWaveTime(i) = mean(UpTrans(NewFullWave(i).ndx));
                nSize=[nSize NewFullWaveSize(i)];
            end
            nSegments=[nSegments segments];
            
        else % NO SEGMENTS identified -->
             % repeated chiannels are due to noise (and not to the presence of more than one wave)
             % CLEAN the wave, i.e. keep only the first channel occurrance
            fprintf('ALERT(6) nw=%d\n', nw); 
            fprintf(fid,'ALERT(6) nw=%d\n', nw);
            alert(6)=alert(6)+1;
            
            % ..no more so sure that the while loop is correct when elemets are deleted (and not inserted) in the array
%             j=1;
%             while j<=numel(chs) % while-loop beacause the numel(chs) is updated in the loop 
%                 listNdx=(find(chs==chs(j)));
%                 if(numel(listNdx)~=1)
%                 % Keep the first occurrance, remove the others (interpreted as noise)
%                     chs(listNdx(2:end))=[];
%                     FullWave(nw).ndx(listNdx(2:end))=[];
%                 end
%                 j=j+1;
%             end 
            
            delCh=[];
            for j=1:numel(chs)
                listNdx=(find(chs==chs(j)));
                if(numel(listNdx)~=1)
% Keep the occurrance which is the closest to the other occurances in the "clan"
                    t0Clan=0; nClan=0;
                    for k = arrayMask{chs(j)}{2}
                        %disp(k)
                        tClan=UpTrans(FullWave(nw).ndx(chs==k));
                        if ~isempty(tClan)
                            %nClan=nClan+1;
                            %t0Clan=t0Clan+tClan;
                            nClan=nClan+numel(tClan); % take into account the case the CLAN has repetions
                            t0Clan=t0Clan+sum(tClan); % take into account the case the CLAN has repetions  
                        end                               
                    end
                    t0Clan=t0Clan/nClan; 
                    tCh=UpTrans(FullWave(nw).ndx(chs==chs(j)));
                    [t0Ch,index]=min(abs(tCh-t0Clan));
                    delCh=[delCh setdiff(listNdx,listNdx(index))];                            
                end 
            end 
            chs(unique(delCh))=[];
            FullWave(nw).ndx(unique(delCh))=[];
       
            % wave is "cleaned"; store and plot the updated wave
            NewFullWave.ndx = FullWave(nw).ndx;
            NewFullWaveUnique = 1; % update the logical value (...we are "cleaning" the waves)
            NewFullWaveSize = numel(NewFullWave.ndx);
            NewFullWaveTime = mean(UpTrans(NewFullWave.ndx));
            
            nSize=[nSize NewFullWaveSize];
            nSegments=[nSegments 1];
        end

        
        % --- REPLACE CurrentWave with NewWave(s) 
        % [its segments or its 'cleaned' version]
        if nw~=1 
            FullWave=[FullWave(1:nw-1) NewFullWave FullWave(nw+1:end)];
            FullWaveUnique = [FullWaveUnique(1:nw-1) NewFullWaveUnique FullWaveUnique(nw+1:end)];
            FullWaveSize = [FullWaveSize(1:nw-1) NewFullWaveSize FullWaveSize(nw+1:end)];
            FullWaveTime = [FullWaveTime(1:nw-1) NewFullWaveTime FullWaveTime(nw+1:end)];
        else
            FullWave=[NewFullWave FullWave(nw+1:end)];
            FullWaveUnique = [NewFullWaveUnique FullWaveUnique(nw+1:end)];
            FullWaveSize = [NewFullWaveSize FullWaveSize(nw+1:end)];
            FullWaveTime = [NewFullWaveTime FullWaveTime(nw+1:end)];
        end
        
% --- PLOT the wave after the treatment --- 
        mask={[1 0 0] [0 1 0] [0 0 1]}; % red, green, blue [alternance of 3 colors]
        % C={'r','b','g','m','c','y','k'};
        for j=1:1:numel(NewFullWave)
            c=mask{mod(j,3)+1}; % selection from a set of 3 colors)
            plot(UpTrans(FullWave(nw-1+j).ndx),ChLabel(FullWave(nw-1+j).ndx),...
                'color',c,'marker','o','linestyle','none');
        end
        print('-dpdf', [OutputDir 'nW/nW' num2str(nw)])
        close(f);
% -----------------------------------------        

        % --- INCREMENT the pointer
        if numel(NewFullWave)>1 % SEGMENTS ARE NEW WAVES          
            nAdded=nAdded+numel(NewFullWave)-1;
            nw=nw+numel(NewFullWave); % increment (point at the next wave)
        else % no segments identified, the current wave is a New Wave, because it has been cleaned
            nw=nw+1; % increment (point at the next wave)
        end
                  
        
    else nw=nw+1; % increment (point at the next wave) [current wave is already unique]      
    end
    
    if ~isempty(NewFullWave)
        nNew=[nNew numel(NewFullWave)-1];
    end
    
end % END WHILE LOOP  
fprintf('alert = [%d,%d,%d,%d,%d,%d]\n', alert); 
fprintf(fid,'alert = [%d,%d,%d,%d,%d,%d]\n', alert); 


%% Plot the distribution of nSize (i.e. number of channels involved) and nSegments for AddedWaves
%
    % X=[1:1:nCh];
    % dX = X(2) - X(1);
    % for i=X
    %     NumberOf(i)=length(find(nSize==i));
    % end

figure

subplot(1,2,1)
hist(nSize,[1:1:31]);
    % [XX, YY] = stairs(X - dX/2, NumberOf); %absolute NumberOf
    % XX = [XX(1) XX' XX(end)+dX XX(end)+dX];
    % YY = [0 YY' YY(end) 0];
    % patch(XX, YY, 'b', 'FaceColor', 'g', 'EdgeColor', 'g', 'FaceAlpha', 0.5);
set(gca,'XLim',[0 33]);
xlabel('Number of Channels involved');
ylabel('Number of New Waves');

subplot(1,2,2)
hist(nSegments,1:1:max(nSegments));
xlabel('Number of Segments in New Waves');
ylabel('Number of NonUnique Waves');

dim2 = [0.1 0.9 0.1 0.1];
str2 = ('NonUnique Waves... after the "cleaning"');
annotation('textbox',dim2,'String',str2,'FitBoxToText','on','FontSize',12,'FontWeight','bold',...
    'LineStyle','none');%,'HorizontalAlignment','center' );

dim0 = [0.75 0.875 0.1 0.1];
str0 = {['NonUnique Waves: ' num2str(NonUnique)],['Added Waves: ' num2str(nAdded)],...
    ['New Waves: ' num2str(numel(nSize))]};
annotation('textbox',dim0,'String',str0,'FitBoxToText','on','BackgroundColor','white');

print('-deps2c', [OutputDir 'nSize_nSegments.eps']);
%title({'WaveCollection 1','Sequence of UpTrans and Optimization of MaxAbsTimeLag'});


%% Remove from waves nonunique channels...
%
if exist('OLDv','var')
    FullWaveBackup = FullWave;

    fprintf('\tWaves to clean     : %d\n', numel(find(FullWaveUnique==0)));
    fprintf('\t(Max. non-unique ch.: %d)\n', max(FullWaveSize-FullWaveUniqueSize));
    fprintf(fid,'\tWaves to clean     : %d\n', numel(find(FullWaveUnique==0)));
    fprintf(fid,'\t(Max. non-unique ch.: %d)\n', max(FullWaveSize-FullWaveUniqueSize));

    listOfnw=[]; % List of "critical" waves: non-unique but possibly hiding multiple unique waves
    for nw = 1:numel(FullWaveUnique) 
        % numel(FullWaveUnique) = number of waves in the collection
        % FullWaveUnique==0 if the waves contains repeated channels
       if FullWaveUnique(nw) == 0
          chs = ChLabel(FullWave(nw).ndx);
    % FullWave(nw).ndx = index of the transitions that make up the wave(nw)
    % ChLabel() = the corresponding channel for each transition 
    % --> in this list there are repeated channels

          nch = hist(chs,1:numel(RecordingSet)); % WARNING! not numel(ChannelSet) even in the case of a subset

          ndx = find(nch>1); % non-unique channels

    %[GDB] CHECK the amount of channels and transitions suppressed from the non-unique wave      

          C=length(ndx); % amount of suppressed channels from non-unique wave
          L1=length(chs);
          [chs,ndx] = setdiff(chs,ndx); % unique channels
    % "ndx" now contains the index of transition in the wave corresponding to the ordered list of unique channels
          L2=length(chs);
          T=L1-L2; % amount of suppressed transitions from non-unique wave

          if C>ceil(nCh/2)
              fprintf('\n*** nW = %d\tsuppressed channels = %d \tsuppressed transitions = %d',nw,C,T);
              fprintf(fid,'\n*** nW = %d\tsuppressed channels = %d \tsuppressed transitions = %d',nw,C,T);    
              listOfnw=[listOfnw, nw];
          end

          if T>ceil(nCh/2)
              fprintf('\n### nW = %d\tsuppressed channels = %d \tsuppressed transitions = %d',nw,C,T);
              fprintf(fid,'\n### nW = %d\tsuppressed channels = %d \tsuppressed transitions = %d',nw,C,T);
              listOfnw=[listOfnw, nw];
          end

          FullWave(nw).ndx = FullWave(nw).ndx(ndx);
    % FullWave(nw) is now a "reduced" wave with only transitions occurred in unique (ie. non repeated) channels
    % N.B. REPEATED CHANNELS ARE SUPPRESSED ! ! !  (and their transitions removed)

          FullWaveUnique(nw) = 1; % update the logical value (...we are "cleaning" the waves)
          FullWaveSize(nw) = numel(FullWave(nw).ndx);
          FullWaveTime(nw) = mean(UpTrans(FullWave(nw).ndx));
       end
    end

    % Save plot of "critical" non-unique waves = hidden unique waves
    % Create the folder for the storage of plots if non-existing, and empty if existing
    FilesToDelete = [OutputDir 'nW/nW*'];
    if exist([OutputDir 'nW/'],'dir') == 0
        mkdir([OutputDir 'nW/']);
    elseif ~isempty(dir([OutputDir 'nW/']))
         delete(FilesToDelete);    
    end

    listOfnw=unique(listOfnw); % remove duplicate indexes
    fprintf('\n\nCritical waves: %d\n',numel(listOfnw));
    fprintf(fid,'\n\nCritical waves: %d\n',numel(listOfnw));
    if(~isempty(listOfnw))
        for i = listOfnw
          figure
          plot(UpTrans(FullWaveBackup(i).ndx),ChLabel(FullWaveBackup(i).ndx),'k.')
          hold on
          plot(UpTrans(FullWave(i).ndx),ChLabel(FullWave(i).ndx),'ro')
          xlabel('Time (s)');
          ylabel('Channel');
          dim0 = [0.2175 0.825 0.1 0.1];
          str0 = {['nW = ' num2str(i)]};
          annotation('textbox',dim0,'String',str0,'FitBoxToText','on','BackgroundColor','none');
    %       XRange = [min(UpTrans(FullWaveBackup(i).ndx)) min(UpTrans(FullWaveBackup(i).ndx))+2*MAX_IWI];
    %       XRange = XRange + [-1 +1]*diff(XRange)*0.1;
    %       set(gca,'XLim',XRange)
          print('-dpdf', [OutputDir 'nW/nW' num2str(i)])
        end
    else
        for i = [1:1:numel(FullWave)]
            if mod(i,20)==0
                figure
                plot(UpTrans(FullWaveBackup(i).ndx),ChLabel(FullWaveBackup(i).ndx),'k.')
                hold on
                plot(UpTrans(FullWave(i).ndx),ChLabel(FullWave(i).ndx),'ro')
                xlabel('Time (s)');
                ylabel('Channel');
                dim0 = [0.2175 0.825 0.1 0.1];
                str0 = {['nW = ' num2str(i)]};
                annotation('textbox',dim0,'String',str0,'FitBoxToText','on','BackgroundColor','none');
    %             XRange = [min(UpTrans(FullWaveBackup(i).ndx)) min(UpTrans(FullWaveBackup(i).ndx))+2*MAX_IWI];
    %             XRange = XRange + [-1 +1]*diff(XRange)*0.1;
    %             set(gca,'XLim',XRange)            
                print('-dpdf', [OutputDir 'nW/nW' num2str(i)])
            end
        end
    end
    % HOW TO SOLVE THE PROBLEM of CRITICAL NON-UNIQUE WAVES? Possible solutions:
    % - reduce IWI (for critical waves only, of for the full waves collection?)
    % - keep the firts occurrance of each channel (to preserve at least one wave)
    % - check the time distance between transitions at the same channel. 
    %   If larger than a threshold, try to "cut" the sequence in more than one wave

end

%% Wave Hunt -- step 2: Coalescence of 'Short' Waves, Rejection of 'Small' Waves --> WaveCollection2
% [GDB]
%
WaveColl(2).nWaves=length(FullWaveSize);

X=[1:1:nCh];
dX = X(2) - X(1);
for i=X
    WaveColl(2).NumberOf(i)=length(find(FullWaveSize==i));
end

figure

subplot(1,2,1)
%plot([1:nCh],WaveColl(1).NumberOf);
%[XX, YY] = stairs(X - dX/2, WaveColl(1).NumberOf/WaveColl(1).nWaves*100); %percentage
[XX, YY] = stairs(X - dX/2, WaveColl(2).NumberOf); %absolute NumberOf
XX = [XX(1) XX' XX(end)+dX XX(end)+dX];
YY = [0 YY' YY(end) 0];
patch(XX, YY, 'b', 'FaceColor', 'b', 'EdgeColor', 'b', 'FaceAlpha', 0.5);
set(gca,'XLim',[0 33]);
xlabel('Number of Channels involved');
ylabel('Number of Waves in the Collection');
hold on
yyaxis right
plot(cumsum(WaveColl(2).NumberOf/WaveColl(2).nWaves))
dim0 = [0.2175 0.825 0.1 0.1];
str0 = {['nWaves = ' num2str(WaveColl(2).nWaves)]};
annotation('textbox',dim0,'String',str0,'FitBoxToText','on','BackgroundColor','none');

subplot(1,2,2)
patch(XX, YY, 'b', 'FaceColor', 'b', 'EdgeColor', 'b', 'FaceAlpha', 0.5);
value=sum(WaveColl(2).NumberOf(ceil(nCh/2):nCh))/WaveColl(2).nWaves;
set(gca,'XLim',[ceil(nCh/2) 33]);
xlabel('Number of Channels involved');
ylabel('Number of Waves in the Collection');

dim1 = [0.65 0.825 0.1 0.1];
% Size and location, specified as a four-element vector of the form [x y w h]. 
% The first two elements specify the coordinates of the lower left corner of the text box.
% The second two elements specify the width and height of the annotation, respectively.
% By default, the units are normalized to the figure. 
% The lower left corner of the figure maps to (0,0) and the upper right corner maps to (1,1). 
str1 = {'Percentage of waves involving more', ['than half of the electrode array: ' num2str(value)]};
% Each element of the cell array displays on a separate line
annotation('textbox',dim1,'String',str1,'FitBoxToText','on','BackgroundColor','white');

dim2 = [0.1 0.9 0.1 0.1];
str2 = ('WavesCollection 2 - "FullWaves" and Optimization of IWI');
annotation('textbox',dim2,'String',str2,'FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none');
print('-deps2c', [OutputDir 'WavesCollection2.eps']);
%title({'WaveCollection 1','Sequence of UpTrans and Optimization of MaxAbsTimeLag'});

%% Find Maxima and Minima in the histogram of FullWaveSize

[hh,xx]=hist(FullWaveSize,(1:1:nCh));

for i=1:1:numel(hh)-1
    if abs(hh(i)-hh(i+1))==1
        hh(i+1)=hh(i);
    end   
end
RI=diff(hh);

bMax=[];
bMin=[];
for i=1:1:numel(RI)-1
    %disp([RI(i), RI(i+1)])
    if sign(RI(i))*sign(RI(i+1))~=1
        if sign(RI(i))~=0
            if sign(RI(i))>sign(RI(i+1))
                fprintf('BIN %d is a MAXIMUM\n',i+1);
                bMax=[bMax i+1];
            else
                fprintf('BIN %d is a MINIMUM\n',i+1);
                bMin=[bMin i+1];
            end
        else
            fprintf('BIN %d has the same value than %d\n',i+1,i);
            if ~isempty(intersect(i,bMin))
                bMin=[bMin i+1];
            end
            if ~isempty(intersect(i,bMax))
                bMin=[bMin i+1];
            end               
        end  
    end
end
if hh(end)>hh(end-1)
    fprintf('BIN %d is a MAXIMUM\n',numel(hh));
    bMax=[bMax i+1];
end


%% Remove small waves and those rejected...
%
%MIN_CH_NUM = 4; 
% If a sub-area is selected, this number could be too large...

if SelectedArea == 'R'
    MIN_CH_NUM = 3;
elseif SelectedArea == 'P'
    MIN_CH_NUM = 2;
% maybe a parameter ad hoc for each cortical area? ...to be OPTIMISED
else
    %MIN_CH_NUM = 4; % too few for 32 electrodes?
    %MIN_CH_NUM = 28; 
%--- [GDB] introduce a requirement of GLOBALITY
%-- at least half of the array:
    if value >= 0.5
        MIN_CH_NUM = ceil(nCh/2);
    % add a requirement concerning an ACCEPTABLE_REJECTION_RATE = 0.5 
    % i.e. at least 50% of the waves in the collection must involve the 50%
    % or more of the electrode array
    else
        R=[ceil(nCh/2):-1:1];
        for i=R
            value=sum(WaveColl(2).NumberOf(i:nCh))/WaveColl(2).nWaves;
            if value >= 0.5
                MIN_CH_NUM = i;
                break
            end
        end
    end    
end
    
if version=='v2'
    MIN_CH_NUM = 16;
end
    
ndx = find(FullWaveUnique==1 & FullWaveSize >= MIN_CH_NUM);
RejectedWaves = numel(FullWaveUnique) - numel(ndx); % rejected beacuse too small 
% (i.e. involving too few channels)
Wave = FullWave(ndx);
WaveUnique = FullWaveUnique(ndx);
WaveSize = FullWaveSize(ndx);
WaveTime = FullWaveTime(ndx);

fprintf('\nMinimum number of channels: %d\n',MIN_CH_NUM);
fprintf('Rejected waves: %d (too small, i.e. involve less than minimum number of channels)\n',RejectedWaves);
fprintf('Accepted waves: %d\n',numel(Wave));
fprintf(fid,'\nMinimum number of channels: %d\n',MIN_CH_NUM);
fprintf(fid,'Rejected waves: %d (too small, i.e. involve less than minimum number of channels)\n',RejectedWaves);
fprintf(fid,'Accepted waves: %d\n',numel(Wave));
WaveColl(3).nWaves=numel(Wave);

TimeLagMatrix = NaN(numel(Wave),numel(RecordingSet)); % TLM initialized with NaN
UpTransNdxMatrix = NaN(numel(Wave),numel(RecordingSet)); % UpTransMatrix initialized with NaN
for k = 1:numel(Wave)
   TLMs = UpTrans(Wave(k).ndx);
   CLs = ChLabel(Wave(k).ndx);
   TimeLagMatrix(k,CLs) = TLMs - mean(TLMs); % each wave is centered at the mean time
   UpTransNdxMatrix(k,CLs) = Wave(k).ndx;
end
TimeLagRange = [min(TimeLagMatrix(isnan(TimeLagMatrix)==0)) ...
                max(TimeLagMatrix(isnan(TimeLagMatrix)==0))]; % max duration of the waves
TLMtoPlot = TimeLagMatrix;
TLMtoPlot(isnan(TLMtoPlot)==1) = -1; % '-1' replace 'NaN' in the plot


% %% Recognize different local waves...
% %
% XPos = [RecordingSet(:).XPos];
% Ch1stCol = find(XPos<0.5);
% Ch2ndCol = find(XPos>=0.5 & XPos<=2.0);
% Ch3rdCol = find(XPos>2.0);
% 
% if isempty(Ch1stCol)
%    WaveType = ones(size(TimeLagMatrix,1),1);
% else
%    WaveType = min(isnan(TimeLagMatrix(:,Ch1stCol)),[],2);
% end
% if isempty(Ch2ndCol)
%    WaveType = [WaveType ones(size(TimeLagMatrix,1),1)];
% else
%    WaveType = [WaveType min(isnan(TimeLagMatrix(:,Ch2ndCol)),[],2)];
% end
% if isempty(Ch3rdCol)
%    WaveType = [WaveType ones(size(TimeLagMatrix,1),1)];
% else
%    WaveType = [WaveType min(isnan(TimeLagMatrix(:,Ch3rdCol)),[],2)];
% end
% % WaveType = [min(isnan(TimeLagMatrix(:,Ch1stCol)),[],2) ...
% %             min(isnan(TimeLagMatrix(:,Ch2ndCol)),[],2) ...
% %             min(isnan(TimeLagMatrix(:,Ch3rdCol)),[],2)];
% SortingIndex=WaveType*(2.^[2 1 0]');

%% Duration of Waves
%

WaveDuration=[];
for i = 1:1:numel(FullWave)
    WaveDuration = [WaveDuration UpTrans(FullWave(i).ndx(end))-UpTrans(FullWave(i).ndx(1))];
end

%% Pick a wave every 20 and plot it to check how it goes...

deltaX=round(max(WaveDuration)/2,2);
% 
for i = 1:1:numel(FullWave)
        if mod(i,20)==0
            figure
%            plot(UpTrans(FullWaveBackup(i).ndx),ChLabel(FullWaveBackup(i).ndx),'k.')
%            hold on
            plot(UpTrans(FullWave(i).ndx),ChLabel(FullWave(i).ndx),'k.')
            set(gca, 'YLim',[0 33]);
            set(gca, 'XLim',[FullWaveTime(i)-deltaX FullWaveTime(i)+deltaX]);
            xlabel('Time (s)');
            ylabel('Channel');
            dim0 = [0.2175 0.825 0.1 0.1];
            str0 = {['nW = ' num2str(i)]};
            annotation('textbox',dim0,'String',str0,'FitBoxToText','on','BackgroundColor','none');
%             XRange = [min(UpTrans(FullWaveBackup(i).ndx)) min(UpTrans(FullWaveBackup(i).ndx))+2*MAX_IWI];
%             XRange = XRange + [-1 +1]*diff(XRange)*0.1;
%             set(gca,'XLim',XRange)            
        end
end


%% Plot timelags matrix sorting waves by type...
%
WHOLE_SLICE_WAVES = 0;
% [SortedSI,ndx]=sort(SortingIndex);
ndx = 1:size(TimeLagMatrix,1); % index along the waves in the collection (number of rows in the TimeLagMatrix)
SortingIndex = zeros(size(ndx))+WHOLE_SLICE_WAVES;

figure
subplot(1,3,2:3) % RIGHT PLOT
imagesc(TimeLagMatrix(ndx,:)*1000,TL_RANGE); % *1000 --> ms
hcb = colorbar();
colormap(CM);
caxis(TL_RANGE);
set(get(hcb,'XLabel'),'String','\delta (ms)') 
set(gca,'Box','on','Layer','top','TickDir','out')
set(gca,'YDir','norm','YTickLabel',[])
title(['Rejected waves: ' num2str(RejectedWaves) ' out of ' num2str(numel(FullWaveUnique))])
% ... rejected beacuse too small i.e. involving  too few channels
% (plot range to be compared with TimeLagRange)

subplot(1,3,1) % LEFT PLOT
hold on
% SIRange = [0 max(SortingIndex)];
for k = 1:numel(ndx) %number of waves
%    ndxCM = floor((SortedSI(k)-SIRange(1))/(diff(SIRange)+1)*size(CM,1))+1;
%    plot(WaveSize(ndx(k)),k,'ko','MarkerFaceColor',CM(ndxCM,:),'MarkerSize',4);
   plot(WaveSize(ndx(k)),k,'ko','MarkerFaceColor',[1 1 1]*0.8,'MarkerSize',4);
end
% set(gca,'XLim',[0.5 numel(RecordingSet)+0.5],'YLim',[1 numel(ndx)])
set(gca,'XLim',[0.5 numel(ChannelSet)+0.5],'YLim',[1 numel(ndx)]) % to account for a SelectedArea of sub-set of channels
set(gca,'Box','on','Layer','top','TickDir','out')
xlabel('Size (Ch.)') % (number of channels involved by the wave)
ylabel('Down-Up Transition') % (number of waves in the collection)

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0 0 4 6]);
print('-deps2c', [OutputDir 'MapOfUpOnsetTimeLags_AllWaves.eps']); % path for the results


%% Estimate missing timelags to perform PCA... *** replace the NaN with a reasonable value ***
%
NEAR_NEIGH = 5;

D = zeros(1,numel(WaveTime)*(numel(WaveTime)-1)/2);
k = 0;
for r = 1:numel(WaveTime)
   for c = r+1:numel(WaveTime)
      k = k + 1;
      TwoRows = TimeLagMatrix([r c],:);
      D(k) = mean(abs(diff(TwoRows(:,sum(isnan(TwoRows)) == 0)))); % ! ! ! WARNING ! ! !
%     A bug can occur when too few electrodes are active in an array with a
%     large number of electrodes: a TwoRows with no channels in common
%     would result in a D(k) == 0 (see also the optimisation of MIN_CH_NUM)<-- DONE Jan2018
   end
end

FilledTLM = TimeLagMatrix; % 'Filled' i.e. NaN replaced with estimate
DistMatrix = squareform(D);
NeighborRadius = zeros(1,numel(WaveTime));
for k = 1:numel(WaveTime)
   [~,ndx] = sort(DistMatrix(k,:));
   
   NeighborRadius(k) = mean(DistMatrix(k,ndx((1:NEAR_NEIGH)+1))); 
   % mean of NEAR_NEIGH values (the first NEAR_NEIGH values in ascending order)
   
   ndxF = find(isnan(FilledTLM(k,:)));
   for j = ndxF
      n = 0;
      r = 1;
      RefVal = 0;
      while n < NEAR_NEIGH && r<size(TimeLagMatrix,1) % correction needed when a cortical area (subset of reconrding channel) is enabled
%      while n < NEAR_NEIGH
         r = r + 1;
         if ~isnan(TimeLagMatrix(ndx(r),j))             
            n = n + 1;
            RefVal = RefVal + TimeLagMatrix(ndx(r),j);
         end
      end
%       if RefVal~=0 
%           disp([n RefVal/n])
%       end
      FilledTLM(k,j) = RefVal/NEAR_NEIGH;
   end
end


%% Remove outliers wave...
%
QUARTILE_THRESHOLD = 2.0;

q = prctile(NeighborRadius,[25 50 75]);
ThresholdDistance = q(3)+QUARTILE_THRESHOLD*(q(3)-q(1)); % inter-quartile range

ndx = find(NeighborRadius <= ThresholdDistance);
RejectedWaves = numel(WaveSize) - numel(ndx);

WaveSize = WaveSize(ndx);
WaveTime = WaveTime(ndx);
TimeLagMatrix = TimeLagMatrix(ndx,:);
FilledTLM = FilledTLM(ndx,:);
FilledTLMbackup=FilledTLM;

WaveColl(4).nWaves=numel(WaveSize);

fprintf('\nOutlier waves to exclude: %d\n', RejectedWaves);
fprintf('\n==> TOTAL NUMBER of WAVES: %d\n',numel(WaveSize));
fprintf('(Reduction of the number of elements in the WaveCollection: %d --> %d --> %d --> %d)\n',...
    WaveColl(1).nWaves,WaveColl(2).nWaves,WaveColl(3).nWaves,WaveColl(4).nWaves);
fprintf(fid,'\nOutlier waves to exclude: %d\n', RejectedWaves);
fprintf(fid,'\n==> TOTAL NUMBER of WAVES: %d\n',numel(WaveSize));
fprintf(fid,'(Reduction of the number of elements in the WaveCollection: %d --> %d --> %d --> %d)\n',...
    WaveColl(1).nWaves,WaveColl(2).nWaves,WaveColl(3).nWaves,WaveColl(4).nWaves);


%% Principal Component Analysis (PCA)...
%
FilledTLM=FilledTLMbackup(:,ChannelSet); % remove empty columns if a ChannelSet (subsample of electrodes) is used
CovMat = cov(FilledTLM);
[pc,variances,explained] = pcacov(CovMat);

X1 = FilledTLM * pc(:,1);
X2 = FilledTLM * pc(:,2);

RED_BLUE = 1;
CM = jet();
if RED_BLUE
   CM = gray(32);
   CM(:,3) = 1;
   Red = flipud(gray(32));
   Red(2:end,1) = 1;
   CM = [CM; Red(2:end,:)];
end

figure
subplot(1,2,1)
plot(1:length(pc), explained,'bo-', 'MarkerFaceColor', 'w', 'LineWidth', 1.)
set(gca,'Box','on','Layer','top','TickDir','out')
xlabel('PC number');
ylabel('Variance explained (%)');

subplot(1,2,2)
hold on
for k = 1:numel(WaveTime)
   ndxClr = floor((mean(FilledTLM(k,end-2:end))*1000. - TL_RANGE(1))/diff(TL_RANGE)*size(CM,1)) + 1;
   ndxClr = max([1 ndxClr]);
   ndxClr = min([size(CM,1) ndxClr]);
   plot(X1(k),X2(k),'ko','MarkerFaceColor',CM(ndxClr,:));
end
XRange = [min(X1) max(X1)];
XRange = XRange + [-1 +1]*diff(XRange)*0.05;
YRange = [min(X2) max(X2)];
YRange = YRange + [-1 +1]*diff(YRange)*0.05;
set(gca,'XLim',XRange,'YLim',YRange)
set(gca,'Box','on','Layer','top','TickDir','out')
plotOpt=0;
if (explained(1)>=4*explained(2) && plotOpt==1) % if the PC1 is much more larger than PC2
    set(gca,'DataAspectRatio',[1 1 1]) % equal data unit lengths in all directions
    set(gca,'XLim',XRange,'YLim',XRange)
end
xlabel('PC_1');
ylabel('PC_2');

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0 0 7 3]);
print('-deps2c', [OutputDir 'PCAOfUpOnsetTimeLags.eps']);  % path for the results


%% Plot the cleaned and sorted time lag matrix...
%
X1_SORTING = 1;
if isfield(Options.UpTimeLags,'PC2')
   if Options.UpTimeLags.PC2
      X1_SORTING = 0;
   end
end
if X1_SORTING
   [SortedPC1,ndx] = sort(X1); % X1 = projection of the Time Lags in the PC1 direction (PC reference frame)
else
   [SortedPC2,ndx] = sort(X2);
end

figure
subplot(1,3,2:3)
imagesc(FilledTLM(ndx,:)*1000,TL_RANGE);
hcb = colorbar();
set(get(hcb,'XLabel'),'String','\delta (ms)')
set(gca,'Box','on','Layer','top','TickDir','out')
set(gca,'YDir','norm','YTickLabel',[])
title(['Outliers = ' num2str(RejectedWaves)])
xlabel('Channel')
%xticks=double(1:numel(ChannelSet));
if SelectedArea~='A'
    set(gca,'XTick',[1:numel(ChannelSet)],'xticklabels',num2str(ChannelSet(:))) 
    % when a reduced number of electrodes is considered, the xticklabels reflects the correct electrode numbering
end


colormap(CM);

subplot(1,3,1)
hold on
plot([0 0],[1 numel(ndx)],'k--') % vertical line at x=0
plot(mean(FilledTLM(ndx,end-2:end),2)*1000,1:numel(ndx),'b.'); % mean(X,2) = mean over rows (PC1-ordered)                                                        
% set(gca,'XLim',Options.PeriodToAnalyze,'YLim',[1 numel(ndx)])
set(gca,'YLim',[1 numel(ndx)])
set(gca,'Box','on','Layer','top','TickDir','out')
xlabel(['\langle\delta_{' num2str(size(FilledTLM,2)-2) '-' num2str(size(FilledTLM,2)) '}\rangle (ms)'])
ylabel('Down-Up Transition')
title(['n = ' num2str(numel(WaveSize))])

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0 0 4 6]);
print('-deps2c', [OutputDir 'MapOfUpOnsetTimeLags_CleanAndSortedWaves.eps']);  % path for the results


%% Save sorted time lag matrix and wave information...
%
WaveSize = WaveSize(ndx);
WaveTime = WaveTime(ndx);

FilledTLM = FilledTLM(ndx,:);
TimeLagMatrix = TimeLagMatrix(ndx,:);

% save('MapOfUpOnsetTimeLags.mat','Channels','WaveSize','WaveTime',... 
%                                 'FilledTLM','TimeLagMatrix');
% Error using save. Variable 'Channels' not found.
save([OutputDir 'MapOfUpOnsetTimeLags.mat'],'ChannelSet','WaveSize','WaveTime',...
                                'FilledTLM','TimeLagMatrix'); % path for the results
                             

%% (I) First estimate of the NUM_OF_CLUSTERS
%
NUM_OF_CLUSTERS = 4; %6; % 8; % 16 (INITIALIZATION VALUE)
% NUM_OF_CLUSTERS = max(MIN_NUM_OF_CLUSTERS,ceil(numel(ChannelSet)/2)); 
    % a (naive) algorithm to optimise the number of clusters taking into
    % account the number of electrodes in the array. 
    
    % BUT The number of waves in the collection (size(TimeLagMatrix,1)) 
    % should be taken into account... 
    % --> MAX_NUM_OF_CLUSTERS

MIN_NUM_OF_CLUSTERS = 2;
MIN_NUM_OF_WAVES = 10; % Minimum number of waves per cluster
MAX_NUM_OF_CLUSTERS = floor(size(TimeLagMatrix,1)/MIN_NUM_OF_WAVES);

% To estimate NUM_OF_CLUSTERS, consider the cumulative sum of PCA

cumsumPCA=cumsum(explained);
thPCA=90; % threshold
for  i = 1:numel(explained)
    if cumsumPCA(i)>=thPCA
        NUM_OF_CLUSTERS = i;
        break
    end
end

fprintf('\nMIN_NUM_OF_CLUSTERS: %d', MIN_NUM_OF_CLUSTERS);
fprintf('\nMAX_NUM_OF_CLUSTERS: %d',MAX_NUM_OF_CLUSTERS);
fprintf('\nNumber Of Clusters: %d\n',NUM_OF_CLUSTERS);
fprintf(fid,'\nMIN_NUM_OF_CLUSTERS: %d', MIN_NUM_OF_CLUSTERS);
fprintf(fid,'\nMAX_NUM_OF_CLUSTERS: %d',MAX_NUM_OF_CLUSTERS);
fprintf(fid,'\nNumber Of Clusters: %d\n',NUM_OF_CLUSTERS);

if NUM_OF_CLUSTERS < MIN_NUM_OF_CLUSTERS
    fprintf('\tWARNING: NUM_OF_CLUSTERS < MIN_NUM_OF_CLUSTERS');
    fprintf(fid,'\tWARNING: NUM_OF_CLUSTERS < MIN_NUM_OF_CLUSTERS');
    NUM_OF_CLUSTERS = MIN_NUM_OF_CLUSTERS;
    fprintf(' --> Number Of Clusters: %d',NUM_OF_CLUSTERS);
    fprintf(fid,' --> Number Of Clusters: %d',NUM_OF_CLUSTERS);
end

if NUM_OF_CLUSTERS > MAX_NUM_OF_CLUSTERS
    fprintf('\tWARNING: NUM_OF_CLUSTERS > MAX_NUM_OF_CLUSTERS');
    fprintf(fid,'\tWARNING: NUM_OF_CLUSTERS > MAX_NUM_OF_CLUSTERS');    
    NUM_OF_CLUSTERS = MAX_NUM_OF_CLUSTERS;
    fprintf(' --> Number Of Clusters: %d\n',NUM_OF_CLUSTERS);
    fprintf(fid,' --> Number Of Clusters: %d\n',NUM_OF_CLUSTERS);
end

if isfield(Options.UpTimeLags,'NumOfClusteredWaves')
   NUM_OF_CLUSTERS = Options.UpTimeLags.NumOfClusteredWaves;
end


%% CLUSTERING by k-means 
%  [GDB] while loop and optimization
% 

Set_NumOfClusters=0;
if Set_NumOfClusters
    NUM_OF_CLUSTERS = 4;
end

CLUSTER=1;
while CLUSTER==1
    [T,Centroids,sumd] = kmeans(FilledTLM,NUM_OF_CLUSTERS);
    CovMat = cov(Centroids);
    [pc,~,~] = pcacov(CovMat);
    % X1 = Centroids*pc(:,1);
    X1 = -Centroids*pc(:,1);
    X2 = Centroids*pc(:,2);
    [Angles,Rho] = cart2pol(X1,X2);
    [~,ndx] = sort(Angles);
    Centroids = Centroids(ndx,:);
    [~,ndx] = sort(ndx);
    T = ndx(T);
    [SortedT,ndxC] = sort(T);

% (II) Second estimate of NUM_OF_CLUSTERS 
% -- Identify the number of "small clusters", i.e. clusters with too few waves
%    (less than MIN_NUM_OF_WAVES/2)
% -- Reduce the NUM_OF_CLUSTERS


    if Set_NumOfClusters==0;
        nSmallClusters=0;
        for wc = 1:NUM_OF_CLUSTERS
            ndxCl= find(T == wc)';
            NumOfWaves = numel(ndxCl);
            if NumOfWaves < MIN_NUM_OF_WAVES/2
                nSmallClusters=nSmallClusters+1;
            end
        end

        if nSmallClusters~=0
            disp(['nSmallClusters = ' num2str(nSmallClusters) ' --> reduce the Number of Clusters and re-try'])
            fprintf(fid,'nSmallClusters = %d --> reduce the Number of Clusters and re-try\n',nSmallClusters); 
            CLUSTER=1;
        else
            CLUSTER=0; % no need to try a further clustering with a different number of clusters
        end

        NUM_OF_CLUSTERS = NUM_OF_CLUSTERS - nSmallClusters;
        if NUM_OF_CLUSTERS < MIN_NUM_OF_CLUSTERS
            disp('WARNING: NUM_OF_CLUSTERS < MIN_NUM_OF_CLUSTERS ! ! !')
            fprintf(fid,'WARNING: NUM_OF_CLUSTERS < MIN_NUM_OF_CLUSTERS ! ! !'); 
            NUM_OF_CLUSTERS = MIN_NUM_OF_CLUSTERS;
        end  
    else CLUSTER=0; end;
end

percentage = cumsumPCA(NUM_OF_CLUSTERS);
disp(['Number of Clusters: ' num2str(NUM_OF_CLUSTERS) ' (' num2str(percentage) ' percentage of PCA)'])
fprintf(fid,'Number of Clusters: %d (%f percentage of PCA)\n',NUM_OF_CLUSTERS,percentage); 

nWaves = size(FilledTLM,1);
save([OutputDir 'ClusterInfo.mat'],'nWaves','NUM_OF_CLUSTERS','percentage');


%% Plot FilledTimeLag Matrix sorted by cluster.
%
TLRange = [-150 150];

RED_BLUE = 1;
CM = jet();
if RED_BLUE
   CM = gray(32);
   CM(:,3) = 1;
   Red = flipud(gray(32));
   Red(2:end,1) = 1;
   CM = [CM; Red(2:end,:)];
end

figure

subplot(1,4,1)
plot(SortedT,1:numel(WaveTime),'k.')
set(gca,'YDir','norm','XTickLabel',[],'YLim',[0.5 numel(WaveTime)+0.5])
ylabel('Down-Up Transition');
set(gca,'Box','on','Layer','top','TickDir','out')

subplot(1,4,2:4)
imagesc(FilledTLM(ndxC,:)*1000.,TLRange);
colormap(CM)
hcb = colorbar();
set(get(hcb,'XLabel'),'String','\delta (s)')
set(gca,'YDir','norm','YTickLabel',[])
xlabel('Channel')
title(['k-means clusters: ' num2str(NUM_OF_CLUSTERS)])
set(gca,'Box','on','Layer','top','TickDir','out')
if SelectedArea~='A'
    set(gca,'XTick',[1:numel(ChannelSet)],'xticklabels',num2str(ChannelSet(:))) 
    % when a reduced number of electrodes is considered, the xticklabels reflects the correct electrode numbering
end

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0 0 4 6]);
print('-deps2c', [OutputDir 'MapOfUpOnsetTimeLags_KmeansClustered.eps']);  % path for the results


%% Plot propagation graphs and wave features for each cluster...
% Prepare the GRID with elctrode positions
%
% TLRange = [-150 150]; % in ms...
TLRange = [-100 100]; % in ms... % "short" front
TLStep = 10;          % in ms...
CLs = [fliplr(-TLStep:-TLStep:TLRange(1)) 0:TLStep:TLRange(2)];

RED_BLUE = 1;
NumOfCLs = (numel(CLs)-1)/2;
CM = jet(2*NumOfCLs+1);
if RED_BLUE
   CM = gray(NumOfCLs+1);
   CM(:,3) = 1;
   CM(end,:) = 0;
   Red = flipud(gray(NumOfCLs+1));
   Red(2:end,1) = 1;
   CM = [CM; Red(2:end,:)];
end

% ELECTRODE POSITIONS
% XPos = zeros(size(ChannelSet));
% YPos = zeros(size(ChannelSet));
XPos = zeros(size(RecordingSet));
YPos = zeros(size(RecordingSet));
% k = 0;
% for Ch = ChannelSet
%    k = k + 1;
% %    XPos(k) = RecordingSet(Ch).XPos;
% %    YPos(k) = RecordingSet(Ch).YPos;
%    XPos(Ch) = RecordingSet(Ch).XPos;
%    YPos(Ch) = RecordingSet(Ch).YPos;
% end

for Ch = 1:1:numel(RecordingSet)
   XPos(Ch) = RecordingSet(Ch).XPos;
   YPos(Ch) = RecordingSet(Ch).YPos;
end


XRange = [min(XPos) max(XPos)];
XRange = XRange + [-1 1]*diff(XRange)*0.05;
YRange = [min(YPos) max(YPos)];
YRange = YRange + [-1 1]*diff(YRange)*0.05;

% CREATION of the GRID: linearly spaced meshed grid
YSamples = 20;
YGrid = linspace(YRange(1), YRange(2), YSamples+1);
XGrid = XRange(1):diff(YGrid(1:2)):XRange(2);


%% Order the clusters according to the number of waves per cluster
% [GDB]
%
for wc = 1:NUM_OF_CLUSTERS
    ndxCl = find(T == wc)';
    WavesInCluster(wc) = numel(ndxCl);
end

[C,sortedwc]=sort(WavesInCluster, 'descend');


%%
for wc = 1:NUM_OF_CLUSTERS
   
   %% 0. Computes the features of the wavefront in the cluster...
   %
   TLGrid = zeros(numel(YGrid),numel(XGrid));
   [YG,XG] = ndgrid(YGrid,XGrid);
   ndxC = find(T == sortedwc(wc))'; % take the waves of the selected cluster
   NumOfWaves = numel(ndxC);
   for sw = ndxC
      ndx = find(~isnan(TimeLagMatrix(sw,:)));
      
      X = XPos(ndx);
      Y = YPos(ndx);
      TLs = TimeLagMatrix(sw,ndx);
      
      st = tpaps([X;Y],TLs,1);
      TLG = reshape(fnval(st,[XG(:)';YG(:)']), size(XG));

      TLGrid = TLGrid + TLG; % sum, to compute the mean
   end % for sw = ...
   TLGrid = TLGrid/numel(ndxC); % mean value
   

   %% 1. Plot contour lines at different time lags
   %
   figure
   hold('on')

   % compute velocity field...
   ssXG = XG(1:2:end,1:2:end);
   ssYG = YG(1:2:end,1:2:end);
   ssTLGrid = TLGrid(1:2:end,1:2:end); %% ss = 'sub-sampled' (--> estimate of the velocity field)
   st = tpaps([ssXG(:)';ssYG(:)'],ssTLGrid(:)',1); % Thin-plate smoothing spline (INTERPLOTATION)
   fx = fnder(st,[1,0]); % partial derivative (X)
   fy = fnder(st,[0,1]); % partial derivative (Y)
   TX = reshape(fnval(fx,[ssXG(:)';ssYG(:)']), size(ssXG));
   TY = reshape(fnval(fy,[ssXG(:)';ssYG(:)']), size(ssXG));
    
   AbsoluteSpeed = 1./sqrt(TX.^2+TY.^2);
   XSpeed = TX.*(AbsoluteSpeed.^2);
   YSpeed = TY.*(AbsoluteSpeed.^2);

   % plot velocity field as background...
   SpeedThreshold = prctile(AbsoluteSpeed(:),90); % ... not used?
   ndx = find(AbsoluteSpeed<=SpeedThreshold);
%    quiver(ssXG(ndx),ssYG(ndx),XSpeed(ndx),YSpeed(ndx),'k')

%    plot contour lines of Up state wave front at equally spaced times...
%    [C,h] = contour(XGrid,YGrid,TLGrid*1000.,CLs,'LineWidth',1.);
   contour(XGrid,YGrid,TLGrid*1000.,CLs,'LineWidth',1.);
   caxis(TLRange);
   hcb = colorbar();
   set(get(hcb,'YLabel'),'String','\delta (ms)')
   colormap(CM)
   
   % SUPERIMPOSE ELECTRODE POSITIONS
   plot(XPos,YPos,'ko','MarkerFaceColor',[0 0 0]+0.75,'MarkerSize',3.)
   % (includes the position of missing channels)
   
   eps = 0.1;
   step = 0.55;
   % SUPERIMPOSE borders of Cortical Areas
   %%% ('LineStyle',':')
   % -- Retrosplenial
   rectangle('Position',[XPos(13)-eps YPos(13)-eps 0.75*step 3.5*step],'LineWidth',1,'EdgeColor',[0.5 0.5 0.5])
   % -- Parietal
   rectangle('Position',[XPos(14)-eps YPos(14)-eps 2.75*step 0.5*step],'LineWidth',1,'EdgeColor',[0.5 0.5 0.5])
   % -- Visual
   v = [XPos(18)-eps YPos(18)-eps; XPos(30)-eps YPos(30)+eps; XPos(32)+0.5*step YPos(30)+eps;...
       XPos(32)+0.5*step YPos(25)+0.5*step; XPos(25)+eps YPos(25)+0.5*step; XPos(25)+eps YPos(25)-0.5*step;...
       XPos(32)+0.5*step YPos(25)-0.5*step; XPos(32)+0.5*step YPos(18)-eps; XPos(18)-eps YPos(18)-eps];
   f = [1 2 3 4 5 6 7 8 9];
   patch('Faces',f,'Vertices',v,'LineWidth',1,'FaceColor','none','EdgeColor',[0.5 0.5 0.5]);
   % -- Somatosensory
   v = [XPos(5)+eps YPos(5)-eps; XPos(5)+eps YPos(12)-0.5*step; XPos(12)+0.5*step YPos(12)-0.5*step;...
       XPos(12)+0.5*step YPos(12)+0.5*step; XPos(10)-eps YPos(12)+0.5*step; XPos(10)-eps YPos(5)-eps];
   f = [1 2 3 4 5 6];
   patch('Faces',f,'Vertices',v,'LineWidth',1,'FaceColor','none','EdgeColor',[0.5 0.5 0.5]);
   % -- Motor
   v = [XPos(9)-eps YPos(9)+eps; XPos(9)-eps YPos(1)-eps; XPos(2)+eps YPos(1)-eps;...
       XPos(2)+eps YPos(1)+0.5*step; XPos(9)+0.5*step YPos(1)+0.5*step; XPos(9)+0.5*step YPos(9)+eps];
   f = [1 2 3 4 5 6];
   patch('Faces',f,'Vertices',v,'LineWidth',1,'FaceColor','none','EdgeColor',[0.5 0.5 0.5]);
   
   xlabel('X (mm)');
   ylabel('Y (mm)');
   title(sprintf('Cluster %d (n = %d), Median |v| = %.3f mm/s', wc, numel(ndxC), median(AbsoluteSpeed(:))));
   
   NumOfWavePerCluster(wc) = numel(ndxC);
   SpeedOfWavePerCluster(wc) = median(AbsoluteSpeed(:));
   
   set(gca,'XLim',XRange,'YLim',YRange,'YDir','norm','DataAspectRatio',[1 1 1])
   set(gca,'Box','on','Layer','top','TickDir','out')

   
   %% 5. Save the figure as EPS and SaveClusterInfo
   %
   set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0 0 5 6.5]);
   print('-deps2c',[OutputDir 'ClusteredWave_' num2str(wc) '.eps']);

   Cluster(wc).NumOfWaves=NumOfWaves;
   Cluster(wc).MedianV=median(AbsoluteSpeed(:));
   Cluster(wc).Grid.X=XGrid;
   Cluster(wc).Grid.Y=YGrid;
   Cluster(wc).Grid.Z=TLGrid;
   
   
end % for wc = ...

save([OutputDir 'ClusterInfo.mat'], 'Cluster','-append');
fclose(fid);
source=[OutputDir FileName '.out'];
destination=[pwd '/Step2-Output/' version '/' FileName '.out'];
copyfile(source,destination);
close('all')


%% END

