% setPath.m
%
%   v0: Setting of the workspace directories for multiple users
%   v1: setPath for SWAP 
%
%   Copyright 2017 Giulia De Bonis @ INFN, Rome - Italy
%   Version: 1.0 - Jan. 19, 2017
%   Version: 2.0 - May. 21, 2018

mouse = 'WT'; % 'KO' or 'WT'
%%% disp(['[' mouse ' mouse]']);
emission = 'Spontaneous'; % or: 'Evoked'

BaseDir = erase(userpath,'Matlab');
DataDir = [BaseDir 'DATA/WBS_' mouse '_Spont_Evoked/']; 

FileList = dir([DataDir '*' emission '*']);
nExp=size(FileList,1);

if(debug)
    fprintf('FileList:\n')
    for i=1:nExp
        fprintf([FileList(i).name '\n'])
    end
end

i=11; % <-- change this number to scroll the FileList
FileName = erase(FileList(i).name,'.smr');
%%% disp(['FileName: ' FileName]);

DataFile = [DataDir FileName '.smr']; 

AnalysisDir = [BaseDir 'RESULTS/Nov2016/' emission '/' mouse '/' FileName '/'];

%%% OutputFile 
% fid = fopen([userpath '/output/' mouse '/' step '/' FileName '.out'], 'wt');