%function [mxpos, mypos] = DensityPeakFinder(Density) 
%   DensityPeakFinder.m
%
%   Find peaks in the 2D Density map (logMUA-LFP plane). 
%   The peaks are related to the expression of bi-modality in the logMUA distribution. 
%   N.B. Very naive algorithm. It fails if more than one local maximum is on the
%   same column or row.
%
%   Copyright 2017 Giulia De Bonis @ INFN, Rome - Italy
%   Version: 1.0 - Apr. 3, 2017 --> val = max(Density)
%   Version: 2.0 - Apr. 4, 2017 --> val = sum(Density)
%

% D = 10000; 
% delta = 1/D;
% (parameters D and delta are set in analyzeUpDownStatevsRecordingSet)

% -- small bug --
% % COLUMNS (Xdim, i.e. logMUA)
%     val = max(Density)';
%     val = round(val*D)/D;
%     A = (diff(val) == 0); % A = 1 if consecutive values are the same
%     idx = find(~A); % keep only indexes for which there is a variation between consecutive values
%     B = val(idx); % clean values (repeated consecutive values are removed)
%     gridX = xx(idx);
       

% *** COLUMNS (Xdim, i.e. logMUA)
    if (algo==1); val = max(Density);
    else (algo==2); val = sum(Density,1); end;    
    val = round(val*D)/D;
    A = diff([0 val])~=0;
    idx = find(A); % keep only indexes for which there is a variation between consecutive values
    B = val(idx); % clean values (repeated consecutive values are removed)
%    gridX = xx(idx);
    
    mx = []; mxpos = [];
    
    for i=2:length(B)-1        
%        if ((B(i)-B(i-1))>=delta) && ((B(i)-B(i+1))>=delta) % local maximum
        if ((B(i)> B(i-1)) && ((B(i)> B(i+1)))) % local maximum
            mx = [mx B(i)]; mxpos = [mxpos idx(i)];
        end
    end
    nx=numel(mx);
    
% *** ROWS (Ydim, i.e. LFP)
    if (algo==1); val = max(Density,[],2)';
    elseif (algo==2); val = sum(Density,2)'; end;    
    val = round(val*D)/D;
    A = diff([0 val])~=0;
    idx = find(A);
    B = val(idx); % clean values
%    gridY = yy(idx);

    my = []; mypos = [];
   
    for i=2:length(B)-1        
%        if ((B(i)-B(i-1))>=delta) && ((B(i)-B(i+1))>=delta) % local maximum
        if ((B(i)> B(i-1)) && ((B(i)> B(i+1)))) % local maximum
            my = [my B(i)]; mypos = [mypos idx(i)];
        end
    end    
    ny=numel(my);
    
    % ... more checks on type/number of peaks found in xx and yy...
    % if numel(mx)!=numel(my)... work around to get the same numel
    if nx>ny; % 1) mxpos>mypos
        disp('(numel(mx)...correction in progress)')
        nDel = length(mxpos)-length(mypos);
        [S,ndx] = sort(mx,'ascend');
        mx(ndx(1:nDel))=[];
        mxpos(ndx(1:nDel))=[];
    elseif ny>nx; % 2) mypos>mxpos
        disp('(numel(my)...correction in progress)')
        nDel = length(mypos)-length(mxpos);
        [S,ndx] = sort(my,'ascend');
        my(ndx(1:nDel))=[];
        mypos(ndx(1:nDel))=[];        
    end
    
    if length(mxpos)==length(mypos)
        % plot(xx(mxpos),yy(mypos),'g*');
        if (algo==1); plot(xx(mxpos),yy(mypos),'g*');
        elseif (algo==2); plot(xx(mxpos),yy(mypos),'y*'); end;
        disp(['Number of Local Maxima: ' num2str(length(mxpos))])
    
        if (length(mxpos)==1); disp('*** weak bi-modality ***');
        elseif (length(mxpos)==2); disp('*** strong bi-modality ***');
        elseif (length(mxpos)>2); disp('*** noisy bi-modality ***'); 
        else disp('WARNING...no maxima detected'); end;
    else disp('WARNING...something unexpected is happening')
    end
    
%     xPeak=xx(mxpos);
%     yPeak=yy(mypos);
%end