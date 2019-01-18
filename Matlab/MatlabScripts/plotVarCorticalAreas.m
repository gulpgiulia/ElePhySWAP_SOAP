function [resultsV] = plotVarCorticalAreas(Variable, nPlot)
%
%  [resultsV] = plotVarCorticalAreas(Variable, nPlot)
%
%   Copyright 2017 Giulia De Bonis @ INFN, Rome - Italy
%   Version: 1.0 - Apr. 13, 2017
%

%% COLORS

SHAD_RED_1 = [1 0.6 0.6]; SHAD_RED_2 = [1 0.8 0.8]; 
SHAD_BLUE_1 = [0.6 0.6 1]; SHAD_BLUE_2 = [0.4 0.4 1];
grey1 = [0.4,0.4,0.4]; grey2 = [0.6,0.6,0.6];

%%
Var(5).name = 'V'; Var(5).subset = [18 19 20 22 23 24 25 27 28 29 30 31 32];
Var(4).name = 'R'; Var(4).subset = [13 17 21 26];
Var(3).name = 'P'; Var(3).subset = [14 15 16];
Var(2).name = 'S'; Var(2).subset = [4 5 7 8 10 11 12];
Var(1).name = 'M'; Var(1).subset = [1 2 3 6 9];
% (cortical areas ordered from front to back) 

ReArrangedVar = [];
Separator = 0;
for i=1:length(Var)
    Var(i).value = Variable(Var(i).subset);
%     disp(Var(i).name);
%     disp(Var(i).value);
    ReArrangedVar = [ReArrangedVar Var(i).value];
    Separator = [Separator Separator(end)+numel(Var(i).subset)];
    Var(i).results = plotSummaryResults(Var(i).value, 0, Var(i).subset); 
    VarMean(i) = Var(i).results.mean;
    VarStd(i) = Var(i).results.std;
    VarStdErr(i) = Var(i).results.stderr;
end
 
separator=Separator(2:end);
resultsV = plotSummaryResults(ReArrangedVar,nPlot); 

%% fig(1) [MEAN]

if nPlot == 1 

    opt = 2;

    figure(resultsV.fig(1)); ax1 = gca; hold on;
    yR = get(ax1, 'YLim');

    for i=1:length(Var)
    %     SHAD_COL = [1-0.2*i 0 1];
        if mod(i,2)==0
            grey = grey2; 
        else
            grey = grey1; 
        end
        p= patch([0.5 separator(i)+0.5 separator(i)+0.5 0.5], [yR(1) yR(1) yR(2) yR(2)],...
            grey, 'EdgeColor', 'none');   

        if opt == 1
            plot([Separator(i)+0.5 Separator(i+1)+0.5],[Var(i).results.mean Var(i).results.mean],...
                'b-','LineWidth',1.5);   
            v1 = [Separator(i)+0.5 VarMean(i)+VarStdErr(i); Separator(i+1)+0.5 VarMean(i)+VarStdErr(i);...
                Separator(i+1)+0.5 VarMean(i)-VarStdErr(i); Separator(i)+0.5 VarMean(i)-VarStdErr(i)];
            v2 = [Separator(i)+0.5 VarMean(i)+VarStd(i); Separator(i+1)+0.5 VarMean(i)+VarStd(i);...
                Separator(i+1)+0.5 VarMean(i)-VarStd(i); Separator(i)+0.5 VarMean(i)-VarStd(i)];
            f = [1 2 3 4];
            patch('Faces',f,'Vertices',v2,'EdgeColor',SHAD_BLUE_2,'FaceColor','none');
            patch('Faces',f,'Vertices',v1,'EdgeColor','b','FaceColor','none');
        end

        TickPos(i) = Separator(i) + 0.5 + (Separator(i+1)-Separator(i))/2;
        TickLabel(i) = {Var(i).name};
        uistack(p,'bottom') 
    end

    if opt == 2
        errorbar(TickPos,VarMean,VarStd,'s','Color', 'b','MarkerFaceColor','b', ...
            'MarkerEdgeColor','b','Color',SHAD_BLUE_2);
        errorbar(TickPos,VarMean,VarStdErr,'sb');   
    end

    set(gca,'SortMethod','childorder');
    set(gca,'xTick', TickPos,'xticklabels',TickLabel,'TickLength',[0 0]);
    xlabel('Cortical Areas');
    ylabel(inputname(1));

end

%% mean_CorticalAreas with errorbar (Box Plots)
% 
% figure
% hold on
% 
% err1 = [resultsV.mean + resultsV.stderr resultsV.mean + resultsV.stderr...
%     resultsV.mean - resultsV.stderr resultsV.mean - resultsV.stderr];
% err2 = [resultsV.mean + resultsV.std resultsV.mean + resultsV.std...
%     resultsV.mean - resultsV.std resultsV.mean - resultsV.std];
% n=3;
% hlgn(n)=patch([0.5 5.5 5.5 0.5], err2, SHAD_RED_2, 'EdgeColor', 'none');
% slgn{n} = 'standard deviation';
% n=2;
% hlgn(n)=patch([0.5 5.5 5.5 0.5], err1, SHAD_RED_1, 'EdgeColor', 'none');
% slgn{n} = 'standard error';
% 
% errorbar([1:length(Var)],VarMean,VarStdErr,'sb', 'LineWidth', 1.5);
% set(gca,'xticklabels',[Var.name]);
% n=1;
% hlgn(n)=plot([0.5 5.5], [resultsV.mean resultsV.mean], 'r-', 'LineWidth', 1.5);
% slgn{n} = 'mean';

% xlabel('Cortical Area');
% ylabel('Mu2');
% 
% legend(hlgn, slgn, 'Location','southeast')
% 
% set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
% print('-deps2c', [AnalysisDir 'Summary/' 'Mu2Mean_CorticalAreas.eps']);

%var = [var sprintf('resultsMu2%c.mean ',area)]; 


%%
%%
% % % 
% % % xmax=33;
% % % if exist('subset') == 1
% % %     xmax=numel(subset)+1;
% % % end
% % % 
% % % % for i = 1:3
% % % %   eval(['val' num2str(i) '=' num2str(i * 10)]);
% % % % end
% % % %variable.name(1) = sprintf('mean%s',inputname(1))
% % % 
% % % meanVar = [];
% % % stderrVar = [];



% for area = CorticalArea(1:end)
%     if area == 'V'     % V = Visual Cortex
%     subset = [18 19 20 22 23 24 25 27 28 29 30 31 32];
%     Mu2V = Mu2(subset);
%     resultsMu2V = plotSummaryResults(Mu2V, subset);   
%     figure(resultsMu2V.fig(1));
%     ylabel('Mu2V');    
%     for i = 2:3; close(resultsMu2V.fig(i)); end;  
%     meanMu2 = [meanMu2 resultsMu2V.mean];
%     stderrMu2 = [stderrMu2 resultsMu2V.stderr];
% elseif area == 'R' % R = Retrosplenial Cortex 
%     subset = [13 17 21 26];
%     Mu2R = Mu2(subset);
%     resultsMu2R = plotSummaryResults(Mu2R, subset);   
%     figure(resultsMu2R.fig(1));
%     ylabel('Mu2R');
%     for i = 2:3; close(resultsMu2R.fig(i)); end; 
%     meanMu2 = [meanMu2 resultsMu2R.mean];
%     stderrMu2 = [stderrMu2 resultsMu2R.stderr];
% elseif area == 'P' % P = Parietal Association Area (PtA)
%     subset = [14 15 16];
%     Mu2P = Mu2(subset);
%     resultsMu2P = plotSummaryResults(Mu2P, subset);    
%     figure(resultsMu2P.fig(1));
%     ylabel('Mu2P');
%     for i = 2:3; close(resultsMu2P.fig(i)); end; 
%     meanMu2 = [meanMu2 resultsMu2P.mean];
%     stderrMu2 = [stderrMu2 resultsMu2P.stderr];
% elseif area == 'S' % S = Somatosensory Cortex
%     subset = [4 5 7 8 10 11 12];
%     Mu2S = Mu2(subset);
%     resultsMu2S = plotSummaryResults(Mu2S, subset);   
%     figure(resultsMu2S.fig(1));
%     ylabel('Mu2S');
%     for i = 2:3; close(resultsMu2S.fig(i)); end;    
%     meanMu2 = [meanMu2 resultsMu2S.mean];
%     stderrMu2 = [stderrMu2 resultsMu2S.stderr];
% elseif area == 'M' % M = Motor Cortex
%     subset = [1 2 3 6 9];
%     Mu2M = Mu2(subset);
%     resultsMu2M = plotSummaryResults(Mu2M, subset);    
%     figure(resultsMu2M.fig(1));
%     ylabel('Mu2M');
%     for i = 2:3; close(resultsMu2M.fig(i)); end;  
%     meanMu2 = [meanMu2 resultsMu2M.mean];
%     stderrMu2 = [stderrMu2 resultsMu2M.stderr];
% end
% end
% 
% %% Mu2 -- figure with subplots 
% %
% 
% figure(resultsMu2.fig(1)); ax1 = gca; 
% fig1 = get(ax1,'children'); xR1 = get(ax1, 'XLim'); 
% xL1 = get(ax1, 'XTickLabel'); xT1 = get(ax1, 'XTick');
% figure(resultsMu2V.fig(1)); ax2 = gca; 
% fig2 = get(ax2,'children'); xR2 = get(ax2, 'XLim'); yR2 = get(ax2, 'YLim'); 
% xL2 = get(ax2, 'XTickLabel'); xT2 = get(ax2, 'XTick');
% figure(resultsMu2R.fig(1)); ax3 = gca; 
% fig3 = get(ax3,'children'); xR3 = get(ax3, 'XLim'); yR3 = get(ax3, 'YLim'); 
% xL3 = get(ax3, 'XTickLabel'); xT3 = get(ax3, 'XTick');
% figure(resultsMu2P.fig(1)); ax4 = gca; 
% fig4 = get(ax4,'children'); xR4 = get(ax4, 'XLim'); yR4 = get(ax4, 'YLim'); 
% xL4 = get(ax4, 'XTickLabel'); xT4 = get(ax4, 'XTick');
% figure(resultsMu2S.fig(1)); ax5 = gca; 
% fig5 = get(ax5,'children'); xR5 = get(ax5, 'XLim'); yR5 = get(ax5, 'YLim'); 
% xL5 = get(ax5, 'XTickLabel'); xT5 = get(ax5, 'XTick');
% figure(resultsMu2M.fig(1)); ax6 = gca; 
% fig6 = get(ax6,'children'); xR6 = get(ax6, 'XLim'); yR6 = get(ax6, 'YLim'); 
% xL6 = get(ax6, 'XTickLabel'); xT6 = get(ax6, 'XTick');
% 
% figure; %create new figure
% s1 = subplot(3,3,[1 2 4 5]); copyobj(fig1,s1); title('Mu2');
% set(s1,'XLim',xR1,'XTickLabel',xL1,'XTick',xT1); 
% s2 = subplot(3,3,3); copyobj(fig2,s2); title('Mu2V');
% set(s2,'XLim',xR2,'YLim',yR2,'XTickLabel',xL2,'XTick',xT2);%,'FontSize', 4); 
% s3 = subplot(3,3,6); copyobj(fig3,s3); title('Mu2R');
% set(s3,'XLim',xR3,'YLim',yR3,'XTickLabel',xL3,'XTick',xT3); 
% s4 = subplot(3,3,7); copyobj(fig4,s4); title('Mu2P');
% set(s4,'XLim',xR4,'YLim',yR4,'XTickLabel',xL4,'XTick',xT4); 
% s5 = subplot(3,3,8); copyobj(fig5,s5); title('Mu2S');
% set(s5,'XLim',xR5,'YLim',yR5,'XTickLabel',xL5,'XTick',xT5); %,'FontSize', 8); 
% s6 = subplot(3,3,9); copyobj(fig6,s6); title('Mu2M');
% set(s6,'XLim',xR6,'YLim',yR6,'XTickLabel',xL6,'XTick',xT6); 
% 
% set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
% print('-deps2c', [AnalysisDir 'Summary/' 'Mu2_CorticalAreas.eps']);
% 
% %% Mu2Mean -- Cortical Areas
% %
% 
% figure
% hold on
% err1Mu2 = [resultsMu2.mean + resultsMu2.stderr resultsMu2.mean + resultsMu2.stderr...
%     resultsMu2.mean - resultsMu2.stderr resultsMu2.mean - resultsMu2.stderr];
% err2Mu2 = [resultsMu2.mean + resultsMu2.std resultsMu2.mean + resultsMu2.std...
%     resultsMu2.mean - resultsMu2.std resultsMu2.mean - resultsMu2.std];
% 
% n=3;
% hlgn(n)=patch([0.5 5.5 5.5 0.5], err2Mu2, SHAD_RED_2, 'EdgeColor', 'none');
% slgn{n} = 'standard deviation';
% n=2;
% hlgn(n)=patch([0.5 5.5 5.5 0.5], err1Mu2, SHAD_RED_1, 'EdgeColor', 'none');
% slgn{n} = 'standard error';
% errorbar([1:numel(CorticalArea)],meanMu2,stderrMu2,'ob', 'LineWidth', 1.);
% set(gca,'xticklabels',[' ' cellstr(CorticalArea')' ' ']);
% n=1;
% hlgn(n)=plot([0.5 5.5], [resultsMu2.mean resultsMu2.mean], 'r-', 'LineWidth', 1.5);
% slgn{n} = 'mean';
% 
% xlabel('Cortical Area');
% ylabel('Mu2');
% 
% legend(hlgn, slgn, 'Location','southeast')
% 
% set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
% print('-deps2c', [AnalysisDir 'Summary/' 'Mu2Mean_CorticalAreas.eps']);
% 
% %var = [var sprintf('resultsMu2%c.mean ',area)]; 