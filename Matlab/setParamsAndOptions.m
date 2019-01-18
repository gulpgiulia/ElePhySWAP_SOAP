% setParamsAndOptions.m
%
%   Load analysis params and information about recording dataset.
%
%   Copyright 2015 Maurizio Mattia @ Ist. Super. Sanit�, Rome - Italy
%   Version: 1.0 - Feb. 18, 2015
%

%% Load analysis params and information about recording dataset.
Options.PeriodToAnalyze = [0 500];

Options.SaveMUA = 0;

Options.LogMUA.FreqBand = [200 1500];
Options.LogMUA.SmoothingWindow = 0.040; % It depends on the sampling rate.
Options.UD.ThresholdModulation = 0.6;   % Between 0 and 1. 0.5 is the mean of the Up and Down levels.

Options.LFP.SmoothingWindow = Options.LogMUA.SmoothingWindow; % s...
Options.LFP.SubsamplingRate = Options.LogMUA.FreqBand(1);   % Hz...

Options.LFPMUAplot.MUArange = [-0.75 3.25];
Options.LFPMUAplot.LFPrange = [-2000 2000];

% N.B. struct UpTimeLags depends on the electrode array
Options.UpTimeLags.MaxAbsTimeLag = 0.800;       % Maximum reasonable time lag between electrodes...
Options.UpTimeLags.PCNumOfSTD = 3;              % Num of st.dev. in the PC space to avoid outliers in the PCA.
% % Options.UpTimeLags.InvertElectrodePos = 1;      % If exist and is 1 invert the order of the electrodes position.
Options.UpTimeLags.NumOfClusters = 10;          % Number of clusters of wavefronts.
Options.UpTimeLags.WavefrontTimeGap = 20;       % Time gap between wavefront in ms.


%% MEA-32Ch
% -------------------------------------------------------------------------
% ELECTRODE ARRAY - 32 channels
% -------------------------------------------------------------------------
% --> distinguish between RIGHT (RH) and LEFT (LH) Hemisphere
% --> define a reference system
% --> enable the opportunity to select a specific Cortical Area:
%   V = Visual Cortex
%   R = Retrosplenial Cortex 
%   P = Parietal Association Area (PtA)
%   S = Somatosensory Cortex
%   M = Motor Cortex

hemisphere = FileName(26:27); % 'RH/LH'
SelectedArea = 'A'; % A = all; or 'V', or 'R', or 'P' ... (see above)


%% Numbering of the Electrodes
%
if strcmpi(hemisphere,'RH')
    
    n = 1;
    RecordingSet(n).label ='01.Ch27';
    RecordingSet(n).port = '27';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;    %High-Pass Filter cut-off frequency
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 0*0.550;

    n = n + 1;
    RecordingSet(n).label ='02.Ch18';
    RecordingSet(n).port = '18';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 0*0.550;

    n = n + 1;
    RecordingSet(n).label ='03.Ch28';
    RecordingSet(n).port = '28';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 1*0.550;

    n = n + 1;
    RecordingSet(n).label ='04.Ch19';
    RecordingSet(n).port = '19';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;    
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 1*0.550;

    n = n + 1; 
    RecordingSet(n).label ='05.Ch10';
    RecordingSet(n).port = '10';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 1*0.550;

    n = n + 1; 
    RecordingSet(n).label ='06.Ch29';
    RecordingSet(n).port = '29';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 2*0.550;

    n = n + 1; 
    RecordingSet(n).label ='07.Ch20';
    RecordingSet(n).port = '20';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 2*0.550;

    n = n + 1; 
    RecordingSet(n).label ='08.Ch11';
    RecordingSet(n).port = '11';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 2*0.550;

    n = n + 1;
    RecordingSet(n).label ='09.Ch30';
    RecordingSet(n).port = '30';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 3*0.550;
 
    n = n + 1;
    RecordingSet(n).label ='10.Ch21';
    RecordingSet(n).port = '21';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 3*0.550;

    n = n + 1;
    RecordingSet(n).label ='11.Ch12';
    RecordingSet(n).port = '12';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 3*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='12.Ch4';
    RecordingSet(n).port = '4';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 3*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='13.Ch31';
    RecordingSet(n).port = '31';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 4*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='14.Ch22';
    RecordingSet(n).port = '22';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 4*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='15.Ch13';
    RecordingSet(n).port = '13';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 4*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='16.Ch5';
    RecordingSet(n).port = '5';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 4*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='17.Ch32';
    RecordingSet(n).port = '32';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 5*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='18.Ch23';
    RecordingSet(n).port = '23';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 5*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='19.Ch14';
    RecordingSet(n).port = '14';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 5*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='20.Ch6';
    RecordingSet(n).port = '6';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 5*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='21.Ch33';
    RecordingSet(n).port = '33';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 6*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='22.Ch24';
    RecordingSet(n).port = '24';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 6*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='23.Ch15';
    RecordingSet(n).port = '15';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 6*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='24.Ch7';
    RecordingSet(n).port = '7';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 6*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='25.Ch3';
    RecordingSet(n).port = '3';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 4*0.550;
    RecordingSet(n).YPos = 6*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='26.Ch34';
    RecordingSet(n).port = '34';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 7*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='27.Ch25';
    RecordingSet(n).port = '25';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 7*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='28.Ch16';
    RecordingSet(n).port = '16';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 7*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='29.Ch8';
    RecordingSet(n).port = '8';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 7*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='30.Ch26';
    RecordingSet(n).port = '26';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 8*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='31.Ch17';
    RecordingSet(n).port = '17';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 8*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='32.Ch9';
    RecordingSet(n).port = '9';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 8*0.550;
            
elseif strcmpi(hemisphere,'LH')
    
    n = 1;
    RecordingSet(n).label ='01.Ch10';
    RecordingSet(n).port = '10';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 0*0.550;

    n = n + 1;
    RecordingSet(n).label ='02.Ch19';
    RecordingSet(n).port = '19';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 0*0.550;

    n = n + 1;
    RecordingSet(n).label ='03.Ch9';
    RecordingSet(n).port = '9';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 1*0.550;

    n = n + 1;
    RecordingSet(n).label ='04.Ch18';
    RecordingSet(n).port = '18';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;    
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 1*0.550;

    n = n + 1; 
    RecordingSet(n).label ='05.Ch27';
    RecordingSet(n).port = '27';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 1*0.550;

    n = n + 1; 
    RecordingSet(n).label ='06.Ch8';
    RecordingSet(n).port = '8';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 2*0.550;

    n = n + 1; 
    RecordingSet(n).label ='07.Ch17';
    RecordingSet(n).port = '17';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 2*0.550;

    n = n + 1; 
    RecordingSet(n).label ='08.Ch26';
    RecordingSet(n).port = '26';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 2*0.550;

    n = n + 1;
    RecordingSet(n).label ='09.Ch7';
    RecordingSet(n).port = '7';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 3*0.550;
 
    n = n + 1;
    RecordingSet(n).label ='10.Ch16';
    RecordingSet(n).port = '16';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 3*0.550;

    n = n + 1;
    RecordingSet(n).label ='11.Ch25';
    RecordingSet(n).port = '25';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 3*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='12.Ch33';
    RecordingSet(n).port = '33';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 3*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='13.Ch6';
    RecordingSet(n).port = '6';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 4*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='14.Ch15';
    RecordingSet(n).port = '15';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 4*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='15.Ch24';
    RecordingSet(n).port = '24';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 4*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='16.Ch32';
    RecordingSet(n).port = '32';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 4*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='17.Ch5';
    RecordingSet(n).port = '5';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 5*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='18.Ch14';
    RecordingSet(n).port = '14';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 5*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='19.Ch23';
    RecordingSet(n).port = '23';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 5*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='20.Ch31';
    RecordingSet(n).port = '31';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 5*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='21.Ch4';
    RecordingSet(n).port = '4';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 6*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='22.Ch13';
    RecordingSet(n).port = '13';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 6*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='23.Ch22';
    RecordingSet(n).port = '22';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 6*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='24.Ch30';
    RecordingSet(n).port = '30';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 6*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='25.Ch34';
    RecordingSet(n).port = '34';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 4*0.550;
    RecordingSet(n).YPos = 6*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='26.Ch3';
    RecordingSet(n).port = '3';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 7*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='27.Ch12';
    RecordingSet(n).port = '12';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 7*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='28.Ch21';
    RecordingSet(n).port = '21';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 7*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='29.Ch29';
    RecordingSet(n).port = '29';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 7*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='30.Ch11';
    RecordingSet(n).port = '11';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 8*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='31.Ch20';
    RecordingSet(n).port = '20';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 8*0.550;
    
    n = n + 1;
    RecordingSet(n).label ='32.Ch28';
    RecordingSet(n).port = '28';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 8*0.550;
       
else 
    error('Hemisphere is not correclty identified. Check your input FileName. Exit.')
end

% EVOKED -- Stimulation Port
if strcmpi(emission,'Evoked') 
    EventArray.port = 'Stim';
end

%% SelectedArea
% 
if SelectedArea == 'V'     % V = Visual Cortex
    ChannelSet = [18 19 20 22 23 24 25 27 28 29 30 31 32];
elseif SelectedArea == 'R' % R = Retrosplenial Cortex 
    ChannelSet = [13 17 21 26];
elseif SelectedArea == 'P' % P = Parietal Association Area (PtA)
    ChannelSet = [14 15 16];
elseif SelectedArea == 'S' % S = Somatosensory Cortex
    ChannelSet = [4 5 7 8 10 11 12];
elseif SelectedArea == 'M' % M = Motor Cortex
    ChannelSet = [1 2 3 6 9];
elseif SelectedArea == 'A' % A = All the electrodes
    ChannelSet = [1:32];
else
    disp('Please, specify the Selected Area of the Cortex')
end


%% LOCALITY criteria
%

% %%Geometry only
arrayMask={{1,[2 3 4]} {2,[1 3 4]} {3,[1 2 4 7 6]} {4,[2 1 3 6 7 8 5]}...
    {5,[2 4 7 8]} {6,[3 4 7 10 9]} {7,[3 4 5 8 11 10 9 6]} {8,[5 4 7 10 11 12]}...
    {9,[6 7 10 14 13]} {10,[9 6 7 8 11 15 14 13]} {11,[12 16 15 14 10 7 8]}...
    {12,[8 11 15 16]} {13,[9 10 14 18 17]} {14,[13 9 10 11 15 19 18 17]}...
    {15,[14 10 11 12 16 20 19 18]} {16,[20 19 15 11 12]} {17,[13 14 18 22 21]}...
    {18,[17 13 14 15 19 23 22 21]} {19,[18 14 15 16 20 24 23 22]} ...
    {20,[25 24 23 19 15 16]} {21,[17 18 22 27 26]} {22,[21 17 18 19 23 28 27 26]}...
    {23,[22 18 19 20 24 29 28 27]} {24,[29 28 23 19 20 25]} {25,[29 24 20]}...
    {26,[21 22 27 30]} {27,[26 21 22 23 28 31 30]} {28,[27 22 23 24 29 32 31 30]}...
    {29,[32 31 28 23 24 25]} {30,[26 27 28 31]} {31,[30 27 28 29 32]}...
    {32,[31 28 29]}};

% %%Geometry + Areas
% arrayMask={{1,[2 3 4 6 9]} {2,[1 3 4 5 6 9]} {3,[1 2 4 7 6 9]} {4,[2 1 3 6 7 8 5 10 11 12]}...
%     {5,[2 4 7 8 10 11 12]} {6,[3 4 7 10 9 1 2]} {7,[3 4 5 8 11 10 9 6 12]}...
%     {8,[5 4 7 10 11 12]} {9,[6 7 10 14 13 1 2 3 6]} {10,[9 6 7 8 11 15 14 13 4 5 12]}...
%     {11,[12 16 15 14 10 7 8 5 4]} {12,[8 11 15 16 5 4 7 10]} {13,[9 10 14 18 17 21 26]}...
%     {14,[13 9 10 11 15 19 18 17 16]} {15,[14 10 11 12 16 20 19 18]} {16,[20 19 15 11 12 14]}...
%     {17,[13 14 18 22 21 26]} {18,[17 13 14 15 19 23 22 21 20 24 25 27 28 29 30 31 32]}...
%     {19,[18 14 15 16 20 24 23 22 25 27 28 29 30 31 32]}...
%     {20,[25 24 23 19 15 16 18 22 27 28 29 30 31 32]} {21,[17 18 22 27 26 13]}...
%     {22,[21 17 18 19 23 28 27 26 20 24 25 29 30 31 32]}...
%     {23,[22 18 19 20 24 29 28 27 25 30 31 32]} {24,[29 28 23 19 20 25 18 22 27 30 31 32]}...
%     {25,[29 24 20 18 19 22 23 27 28 30 31 32]} {26,[21 22 27 30 17 13]} ...
%     {27,[26 21 22 23 28 31 30 18 19 20 24 25 29 32]} {28,[27 22 23 24 29 32 31 30 18 19 20 25]}...
%     {29,[32 31 28 23 24 25 30 27 22 18 19 20]} {30,[26 27 28 31 18 19 20 22 23 24 25 29 32]}...
%     {31,[30 27 28 29 32 18 19 20 22 23 24 25]} {32,[31 28 29 18 19 20 22 23 24 25 27 30]}};
