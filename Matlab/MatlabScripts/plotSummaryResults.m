function [resultsV] = plotSummaryResults(Variable,DrawPlot,subset)
%
%  [resultsV] = plotSummaryResults(Variable [, DrawPlot][, subset])
%
%   Copyright 2017 Giulia De Bonis @ INFN, Rome - Italy
%   Version: 1.0 - Apr. 10, 2017
%

if DrawPlot ~0;

    xmax=33;
    if exist('subset') == 1
        xmax=numel(subset)+1;
    end


    %% Variable: mean, median, stderror, stddeviation
    % 

    if DrawPlot >=1;   
        hlgn = []; slgn = [];
        figure;
        hold on

        SHAD_RED_2 = [1 0.8 0.8];
        err2 = std(Variable); % standard deviation
        E2 = mean(Variable) + [err2 err2 -err2 -err2];
        n=3;
        hlgn(n) = patch([0.5 xmax-0.5 xmax-0.5 0.5], E2, SHAD_RED_2, 'EdgeColor', 'none', 'FaceAlpha', 0.75);
        slgn{n} = 'standard deviation';

        SHAD_RED_1 = [1 0.6 0.6];
        err1 = std(Variable)/sqrt(numel(Variable)); % standard error
        E1 = mean(Variable) + [err1 err1 -err1 -err1];
        n=2;
        hlgn(n) = patch([0.5 xmax-0.5 xmax-0.5 0.5], E1, SHAD_RED_1, 'EdgeColor', 'none', 'FaceAlpha', 0.75);
        slgn{n} = 'standard error';

        stem(Variable,'k');
        set(gca, 'XLim', [0 xmax], 'YLim',[min(Variable)-0.1*abs(min(Variable)) 1.1*max(Variable)]);
        xlabel('Channels');
        ylabel('Variable');

        if exist('subset') == 1
            set(gca,'xTick', [1:numel(subset)],'xticklabels',num2str(subset(:))); 
        end

        n=1;
        hlgn(n)= plot([0 xmax-0.5], [mean(Variable) mean(Variable)], 'r-', 'LineWidth', 1.5);
        slgn{n} = 'mean';

        % n=4;
        % hlgn(n) = plot([0 xmax], [median(Variable) median(Variable)], 'b-', 'LineWidth', 1.5);
        % slgn{n} = 'median';

        legend(hlgn, slgn, 'Location','southeast')

        fig(1) = gcf;
        %(gcf = current figure handle)
    end
    

    %% Variable: mean, median, percentile
    %

    if DrawPlot >=2;
        hlgn = []; slgn = [];
        figure;
        hold on

        for prc = 10:10:40 % (percentile)
           clr = [1-prc/100 1-prc/100 1]; % color (RGB triplet)
           YP = prctile(Variable,[prc 100-prc]);
           patch([0.5 xmax xmax 0.5],[YP(1) YP(1) YP(2) YP(2)], clr, 'EdgeColor', 'none');
        end

        stem(Variable,'k');
        set(gca, 'XLim', [0 xmax], 'YLim',[min(Variable)-0.1*abs(min(Variable)) 1.1*max(Variable)]);
        xlabel('Channels');
        ylabel('Variable');

        if exist('subset') == 1
            set(gca,'xTick', [1:numel(subset)],'xticklabels',num2str(subset(:))); 
        end

        n=1;
        % hlgn(n)= plot([0 xmax], [mean(Variable) mean(Variable)], 'r-', 'LineWidth', 1.5);
        % slgn{n} = 'mean';
        %  
        % n=2;
        hlgn(n) = plot([0 xmax], [median(Variable) median(Variable)], 'b-', 'LineWidth', 1.5);
        slgn{n} = 'median';

        legend(hlgn, slgn, 'Location','northeast')

        fig(2) = gcf;
    end

    %% Variable: histogram
    if DrawPlot ==3;
        figure;
        histogram(Variable,'Normalization', 'probability');

        set(gca, 'Layer', 'top', 'Box', 'off', 'TickDir', 'out');

        xlabel('Variable (a.u)');
        ylabel('Samples (%)');

        fig(3) = gcf;
    end
end
%% Results
%

resultsV.mean = mean(Variable);
resultsV.std = std(Variable);
resultsV.stderr = std(Variable)/sqrt(numel(Variable));
resultsV.median = median(Variable);
if DrawPlot ~0; resultsV.fig = fig; end;


