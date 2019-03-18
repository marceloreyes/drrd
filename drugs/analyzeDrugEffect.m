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
SalIC  = 32;
basIC  = 34;
drug   = 35;
basIC2 = 36;
drugHigh = 37;
basIC3   = 38;
drugHigher=39;
ethanol = 42;

%allSessions = [ SAL   HAL1   HAL2   BAS1   APO   BAS2   AMPH   SalIP   HalIP   BAS3   BAS4];
%sessLabels =  {'SAL' 'HAL1' 'HAL2' 'BAS1' 'APO' 'BAS2' 'AMPH' 'SalIP' 'HalIP' 'BAS3' 'BAS4'};

% allSessions = [ SalIP   HalIP    BAS3   BAS4   BAS5  ];
% sessLabels =  {'SalIP'  'HalIP' 'BAS3' 'BAS4' 'BAS5' };

% allSessions = [ BAS2    HalIP   BAS7 HalIP2];
% sessLabels =  { 'BAS2' 'HalIP' 'BAS7' 'HalIP2'};

%allSessions =  [ HalIPHigh   HALIPLow  SalIP2   OLAN1m   OLAN2m OLAN3m];
%sessLabels  =  {'HAL1' 'HAL0.2' 'SAL' 'OLAN1' 'OLAN2' 'OLAN3'};

%allSessions =  [  OLAN2m   OLAN3m];
%sessLabels  =  { 'OLAN2m' 'OLAN3m'};

%allSessions = [ SAL2 HAL3];
%sessLabels =  { 'SAL2' 'HAL3'};

%allSessions =  [ BasIC    SalIC    Drug    basIC2  DrugHigh  DrugHigher];
%sessLabels  =  {'basIC'  'SalIC'   'Drug' 'basIC2' 'DrugHigh' 'DrugHigher'};

%allSessions =  [ SalIC   ethanol   DrugHigh   DrugHigher];
%sessLabels  =  {'SalIC' 'ethanol' 'DrugHigh' 'DrugHigher'};

allSessions =  [ ethanol drugHigher];
sessLabels  =  { 'ethanol' 'drugHigher' };

%allSessions =  [basIC basIC2  SalIC ethanol ];
%sessLabels  = {'basIC' 'basIC2' 'SalIC' 'ethanol'};

%allSessions =  [ ethanol drugHigher];
%sessLabels  =  { 'ethanol' 'drugHigher' };


var     = 1;
colCond = 3;
cond    = 1;
%subjects = [7 8 9 10 11 12 13 14 15];
%subjects = [42 52 53 54 61 63]; % olanzapine
%subjects = [45 50 59 60 62]; % apo
subjects = [36 48 55 56 57]; % haloperidol

%subjects = [42 52 53 54 61 63 45 50 59 60 62 36 48 55 56 57];

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
    countSub  = 1;
    allD{countSess} = []; %#ok<AGROW>
    for thisSubject = subjects
        %D = gatherDrrd(thisSubject,thisSession,plotFlag);
        D = drrd('AB1', thisSubject, thisSession, plotFlag,false);
        out = floor(size(D,1)/3);
        %out = 50;
        %out = aux;
        D = D(out+1:end,:);
        
        allD{countSess} = [allD{countSess}; D(:,1)]; %#ok<AGROW>
        
        % --- Peak of the distribution ---
        M.peak(countSub,countSess) = find_peak(D,gauss,rng,dt);
        
        % --- Peak disconsidering premature responses ---
        cutoff = 0.5; % time in seconds
        Dcut = D(D(:,1)>cutoff,:);
        M.peakCut(countSub,countSess) = find_peak(Dcut,gauss,rng,dt);
        
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
        M.peakByOpp(countSub,countSess) = ...
            compareSessionDistributionsByOpp(thisSubject,thisSession,false);

        % --- STD of the distribution ---
        M.std(countSub,countSess) = std(D(:,1));
        
        % --- cumulative function ---
        n = histc(D(:,1),rng)/length(D(:,1)); 
        M.cumsum{countSub,countSess} = cumsum(n);
        
        countSub = countSub + 1;                        % increments subject number
    end
    countSess = countSess + 1;      % increments the session number
end
close all;
subplot(3,3,1); hold on;
[~,p] = anova_rm(M.peak);close;
plotM(M.peak,allSessions,sessLabels,'time (s)',['Peak Pos. p='  num2str(p{2,6},'%.3f')]);

subplot(3,3,2); hold on;
[~,p] = anova_rm(M.peakCut);close;
plotM(M.peakCut,allSessions,sessLabels,'time (s)',['PeakCut Pos. p='  num2str(p{2,6},'%.3f')]);

subplot(3,3,3); hold on;
[~,p] = anova_rm(M.median_reinf);close;
plotM(M.median_reinf,allSessions,sessLabels,'time (s)',['Median Reinf p=' num2str(p{2,6},'%.3f')]);

subplot(3,3,4); hold on;
[~,p] = anova_rm(M.median_all);close;
plotM(M.median_all,allSessions,sessLabels,'time (s)',['Median all resps p=' num2str(p{2,6},'%.3f')]);

subplot(3,3,5); hold on;
[~,p] = anova_rm(M.reinf_ratio);close;
plotM(M.reinf_ratio,allSessions,sessLabels,'ratio (%)',['Reinf. ratio(%) p=' num2str(p{2,6},'%.3f')]);

subplot(3,3,6); hold on;
[~,p] = anova_rm(M.presses);close;
plotM(M.presses,allSessions,sessLabels,'Number of events',['# presses p=' num2str(p{2,6},'%.3f')]);

%subplot(3,3,7); hold on;
% [~,p] = anova_rm(M.iti);close;
% plotM(M.iti,allSessions,sessLabels,'time (s)',['Median ITI p=' num2str(p{2,6},'%.3f')]);

subplot(3,3,7); hold on;
[~,p] = anova_rm(M.peakByOpp);close;
plotM(M.peakByOpp,allSessions,sessLabels,'time (s)',['Peak by opp. p=' num2str(p{2,6},'%.3f')]);

subplot(3,3,8); hold on;
[~,p] = anova_rm(M.std);close;
plotM(M.std,allSessions,sessLabels,'time (s)',['STD  p='  num2str(p{2,6},'%.3f')]);

subplot(3,3,9); hold on;
for k = 1:length(allSessions)
    allcum = horzcat(M.cumsum{:,k})
    mcum = mean(allcum,2);
    stdcum = std(allcum,0,2)/sqrt(size(M.cumsum,1));
    %n = histc(x,rng)/length(x);
    color = [k k k]/3;
    errorbar(rng,mcum,stdcum,'color',color);
    %plot(rng,cumsum(n),'k-');
end
%plotM(M.std,allSessions,sessLabels,'time (s)',['STD  p='  num2str(p{2,6},'%.3f')]);



%print -dpng -r300 lixo.png


function plotM(M,allSessions,sessLabels,x_axis_label,title_label)
%%
yl = 'auto';%[0 2.2];
%bar(mean(M),'facecolor',[.5 .5 .5]); 
%plot(M','o-','linewidth',2,'markerfacecolor','w','markersize',8);
errorbar(mean(M),std(M)/sqrt(size(M,1)),'k.','linewidth',2);
hold on;
for k = 1: size(M,2)
    plot(k,M(:,k),'ko');
end
set(gca,'xtick',1:1:length(allSessions),'XTickLabel',sessLabels);
title(title_label);
ylabel(x_axis_label);
ylim(yl);

%%
function peak = find_peak(D,gauss,rng,dt)
n = (histc(D(:,1),rng))/length(D(:,1));
C = conv(n,gauss,'same');
ind = find(C == max(C),1,'last');
peak = rng(ind) + dt/2;
        