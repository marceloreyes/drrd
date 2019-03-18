function M = analyzeDrugEffect
%function analyzeDrugEffect

clear; 
%%sessions
SAL   = 7;
HAL1  = 8;
HAL2  = 10;
BAS1  = 11;
APO   = 12;
BAS2  = 14;
AMPH  = 15;
BAS3  = 16;
BAS4  = 17;
SalIP = 18;
HalIP = 19;
BAS5  = 20;
BAS6  = 21;
BAS7  = 22;
HalIPHigh = 23;
BAS8 = 24;
BAS9 = 25;
SAL2 = 26;
HAL3 = 27;
BAS10 = 28;
HALIPLow = 29;
BAS11 = 30;
SalIP2 = 31;
BAS12 = 32;
OLAN1m  = 33;
BAS13 = 34; 
OLAN2m = 35;
BAS14 = 36;
OLAN3m = 37;

%allSessions = [ SAL   HAL1   HAL2   BAS1   APO   BAS2   AMPH   SalIP   HalIP   BAS3   BAS4];
%sessLabels =  {'SAL' 'HAL1' 'HAL2' 'BAS1' 'APO' 'BAS2' 'AMPH' 'SalIP' 'HalIP' 'BAS3' 'BAS4'};

% allSessions = [ SalIP   HalIP    BAS3   BAS4   BAS5  ];
% sessLabels =  {'SalIP'  'HalIP' 'BAS3' 'BAS4' 'BAS5' };

% allSessions = [ BAS2    HalIP   BAS7 HalIP2];
% sessLabels =  { 'BAS2' 'HalIP' 'BAS7' 'HalIP2'};

allSessions =  [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
sessLabels  =  {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15'};

%allSessions =  [  OLAN2m   OLAN3m];
%sessLabels  =  { 'OLAN2m' 'OLAN3m'};

%allSessions = [ SAL2 HAL3];
%sessLabels =  { 'SAL2' 'HAL3'};



var     = 1;
colCond = 3;
cond    = 1;
%subjects = [1 2 3 4 5 6];
subjects = 7:12;

plotFlag = false;

% --- Parameters for peak detection ---
dt = 0.02;
rng = 0:dt:6;         % range of times for binning histogram
sigma = 0.2;
gauss = dt/sqrt(2*pi())/sigma*exp(-0.5*((rng-mean(rng))/sigma).^2);

%%
M=[];

thisSubject = 7;
thisSession = SAL;


countSess = 1;
for thisSession = allSessions
    countSub  = 1;
    for thisSubject = subjects
        %D = gatherDrrd(thisSubject,thisSession,plotFlag);
        D = drrd('AB1', thisSubject, thisSession, plotFlag,false);
        
        % --- Peak of the distribution ---
        n = (histc(D(:,1),rng))/length(D(:,1));
        C = conv(n,gauss,'same');
        ind = find(C == max(C),1,'last');
        peak = rng(ind) + dt/2;
        M.peak(countSub,countSess) = peak;
        
        % --- Median of primed trials ---
        ind = find(D(:,colCond) == cond);
        M.median_reinf(countSub,countSess) = median(D(ind,var));
        
        % --- Median time (all responses) ---
        M.median_all(countSub,countSess) = median(D(:,var));
        
        % --- Median iti (all responses) ---
        M.iti(countSub,countSess) = median(D(1:end-1,2));
        
        % --- ratio of reinforced trials ---
        M.reinf_ratio(countSub,countSess) = sum(D(:,3)==1)/length(D(:,3));
        
        % --- number of lever presses ---
        M.presses(countSub,countSess) = length(D(:,1));
        
        countSub = countSub + 1;                        % increments subject number
    end
    countSess = countSess + 1;      % increments the session number
end
close all;
subplot(2,3,1); hold on;
[~,p] = anova_rm(M.peak);close;
plotM(M.peak,allSessions,sessLabels,'time (s)',['Peak Pos. p='  num2str(p{2,6},'%.2f')]);

subplot(2,3,2); hold on;
[~,p] = anova_rm(M.median_reinf);close;
plotM(M.median_reinf,allSessions,sessLabels,'time (s)',['Median Reinf p=' num2str(p{2,6},'%.2f')]);

subplot(2,3,3); hold on;
[~,p] = anova_rm(M.median_all);close;
plotM(M.median_all,allSessions,sessLabels,'time (s)',['Median all resps p=' num2str(p{2,6},'%.2f')]);

subplot(2,3,4); hold on;
[~,p] = anova_rm(M.reinf_ratio);close;
plotM(M.reinf_ratio,allSessions,sessLabels,'ratio (%)',['Reinf. ratio(%) p=' num2str(p{2,6},'%.2f')]);

subplot(2,3,5); hold on;
[~,p] = anova_rm(M.presses);close;
plotM(M.presses,allSessions,sessLabels,'Number of events',['# presses p=' num2str(p{2,6},'%.2f')]);

subplot(2,3,6); hold on;
[~,p] = anova_rm(M.iti);close;
plotM(M.iti,allSessions,sessLabels,'time (s)',['Median ITI p=' num2str(p{2,6},'%.2f')]);

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
