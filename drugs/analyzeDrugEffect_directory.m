function M = analyzeDrugEffect_directory
%function M = analyzeDrugEffect_directory

clear; 

cd \\VBOXSVR\mbreyes\ufabc\dados\AB

% baseline = 8;
% saline   = 1;
% apo_low  = 2;
% apo_high = 3;
% halo_low = 4;
% halo_high= 5;
% olan_low = 6;
% olan_high= 7;

apo_high = 1;
%halo_high= 2;
allSessions =  [ apo_high ];
sessLabels  =  {'apo_high'};

%allSessions =  [saline   apo_low   apo_high   halo_low   halo_high   olan_low   olan_high baseline   ];
%sessLabels  =  {'saline' 'apo_low' 'apo_high' 'halo_low' 'halo_high' 'olan_low' 'olan_high' 'baseline'};


%% --- 
var     = 1;
colCond = 3;
cond    = 1;

plotFlag = false;

% --- Parameters for peak detection ---
dt = 0.02;
rng = 0:dt:6;         % range of times for binning histogram
sigma = 0.2;
gauss = dt/sqrt(2*pi())/sigma*exp(-0.5*((rng-mean(rng))/sigma).^2);

%%
M=[];

countSess = 1;

for thisSession = allSessions
    cd(sessLabels{thisSession})
    
    thisData = analyzeDataFromThisSession
    lst = ls;                   % gets all file names
    subjects = []; sessions = [];
    % --- getting subject numbers ---
    for k=3:size(lst,1)
        subjects = [subjects str2num(lst(k,5:6))]; %#ok<ST2NM,AGROW>
        sessions = [sessions str2num(lst(k,8:10))]; %#ok<ST2NM,AGROW>
    end
    
    countSub  = 1;
    for l = 1:length(subjects)
        thisSubject = subjects(l);
        thisSessionNumber = sessions(l);
        %D = gatherDrrd(thisSubject,thisSession,plotFlag);
        D = drrd('AB1', thisSubject, thisSessionNumber, plotFlag,false);
        
        % --- Peak of the distribution ---
        n = (histc(D(:,1),rng))/length(D(:,1));
        C = conv(n,gauss,'same');
        ind = find(C == max(C),1,'last');
        peak = rng(ind) + dt/2;
        M.peak(countSub,countSess) = peak;
        
        % --- Median of primed trials ---
        ind = D(:,colCond) == cond;
        M.median_reinf(countSub,countSess) = median(D(ind,var));
        
        % --- Median time (all responses) ---
        M.median_all(countSub,countSess) = median(D(:,var));
        
        % --- Median iti (all responses) ---
        M.iti(countSub,countSess) = median(D(1:end-1,2));
        
        % --- ratio of reinforced trials ---
        M.reinf_ratio(countSub,countSess) = sum(D(:,3)==1)/length(D(:,3));
        
        % --- number of lever presses ---
        M.presses(countSub,countSess) = length(D(:,1));
        
        % --- peak normalized by opportunity ---
        %M.peakByOpp(countSub,countSess) = ...
        %    compareSessionDistributionsByOpp(thisSubject,thisSession,false);
        
        
        countSub = countSub + 1;                        % increments subject number
    end
    countSess = countSess + 1;      % increments the session number
    cd ..
end
close all;
subplot(3,3,1); hold on;
[~,p] = anova_rm(M.peak);close;
plotM(M.peak,allSessions,sessLabels,'time (s)',['Peak Pos. p='  num2str(p{2,6},'%.2f')]);

subplot(3,3,2); hold on;
[~,p] = anova_rm(M.median_reinf);close;
plotM(M.median_reinf,allSessions,sessLabels,'time (s)',['Median Reinf p=' num2str(p{2,6},'%.2f')]);

subplot(3,3,3); hold on;
[~,p] = anova_rm(M.median_all);close;
plotM(M.median_all,allSessions,sessLabels,'time (s)',['Median all resps p=' num2str(p{2,6},'%.2f')]);

subplot(3,3,4); hold on;
[~,p] = anova_rm(M.reinf_ratio);close;
plotM(M.reinf_ratio,allSessions,sessLabels,'ratio (%)',['Reinf. ratio(%) p=' num2str(p{2,6},'%.2f')]);

subplot(3,3,5); hold on;
[~,p] = anova_rm(M.presses);close;
plotM(M.presses,allSessions,sessLabels,'Number of events',['# presses p=' num2str(p{2,6},'%.2f')]);

subplot(3,3,6); hold on;
[~,p] = anova_rm(M.iti);close;
plotM(M.iti,allSessions,sessLabels,'time (s)',['Median ITI p=' num2str(p{2,6},'%.2f')]);

%subplot(3,3,7); hold on;
%[~,p] = anova_rm(M.peakByOpp);close;
%plotM(M.peakByOpp,allSessions,sessLabels,'time (s)',['Peak by opp. p=' num2str(p{2,6},'%.2f')]);


print -dpng -r300 lixo.png


function plotM(M,allSessions,sessLabels,x_axis_label,title_label)
%%
yl = 'auto';%[0 2.2];
bar(mean(M),'facecolor',[.5 .5 .5]); 
%plot(M','o-','linewidth',2,'markerfacecolor','w','markersize',8);
errorbar(mean(M),std(M)/3,'k.','linewidth',2)
set(gca,'xtick',1:1:length(allSessions),'XTickLabel',sessLabels);
title(title_label);
ylabel(x_axis_label);
ylim(yl);
