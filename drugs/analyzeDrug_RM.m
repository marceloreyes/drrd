function M = analyzeDrug_RM(allSessions,eliminatedSubjs)
%function M = analyzeDrug_RM(allSessions,eliminatedSubjs)
% Examples: 
% M = analyzeDrug_RM([5 19], [35 47 55]);
% In this case it will analyze the columns 5 (saline ip) and 19
% (haloperidol 0.07), all subjects but the 35, 47 and 55.

% Todo
% put a routine to find the session by the name

%clear;
%eliminateMissing = true;

%% Parameters for analysis ---

PREFIX = 'AB1';

% --- Columns for repeated measures
if ~exist('allSessions','var') 
    %allSessions = {2,5};      % Baseline and controls - all sessions to be analyzed 
    %allSessions = [8 32 33 34];       % APOMORPHINE IC - all sessions to be analyzed
    %allSessions = [8 34];             % APOMORPHINE IC - all subjects (more than before)
    %allSessions = [6 28 29 30 31];    % APOMORPHINE IP - all sessions to be analyzed
    %allSessions = {[5 6], 31, 28, 29};    % APOMORPHINE IP - all sessions to be analyzed
    % allSessions = [6 16 17 18 19];    % HALOPERIDOL IP - all sessions to be analyzed
    %allSessions =  {[5 6], 19};
    %M = analyzeDrug_RM({[5 6], 19, 18, [16 17]}) 
     %allSessions = [8 13 14 15];       % HALOPERIDOL IC - all sessions to be analyzed
    %allSessions = [8 15];             % HALOPERIDOL IC - all subjects (more than before)
    allSessions = {[5 6], 24, 23};       % OLANZAPINE IP  - all sessions to be analyzed
    %allSessions = [6 20 21];       % OLANZAPINE IP  - all sessions to be analyzed
    %allSessions = [6 20];             % OLANZAPINE IP  - biggest difference
    %allSessions = [8 25 26 27];       % OLANZAPINE IC  - all sessions to be analyzed
    %allSessions = [8 25];             % OLANZAPINE IC  - biggest difference
end

% --- Input option to eliminate specific subjects (int array).
if ~exist('eliminatedSubjs','var') 
    eliminatedSubjs = [];
end

if ~exist('mergeSessions','var') 
    mergeSessions = [];
end

if iscell(allSessions)
    allSessionsCell = allSessions;
    allSessions     = horzcat(allSessions{:}); % gets all sessions numbers in an array
else
    for k = 1:length(allSessions)
        allSessionsCell{k} = allSessions(k);
    end
end


% --- Parameters for peak detection ---
dt = 0.02;
rng = 0:dt:6;         % range of times for binning histogram
sigma = 0.2;
gauss = dt/sqrt(2*pi())/sigma*exp(-0.5*((rng-mean(rng))/sigma).^2);

% --- range for cumulative sum graphs ---
rng_cum = 0:0.1:4;

% ---  Identification of the columns for the data matrix D
colCond = 3;
cond    = 1;
var     = 1;


%%
% --- Loading information from sessions and animals ---
% --- this info is in the sessions.xlsx worksheet ---
% sessLables are the header, which are the session names
% S contains the number of the session, except in the first column, which
% contains the subject numbers
[conditions S] = loadSessionsList;

% --- sessions to be analyzed by a repeated measure (RM) statistics ---
for k = allSessions
    ind = find(S(:,k));                  % index for sessions for that manipulation
    if k == allSessions(1)               % if first iteraction:
        allInd = ind;                    % gets all the indexes from first manipulation
    else
        allInd = intersect(ind,allInd);  % if not: gets only the sessions present
                                         % in this manipulationand AND in the previous
    end
end

for k = eliminatedSubjs
    allInd = setdiff(allInd,find(S(:,1)==k));
end
    

ind = allInd;
for k = 1:length(allSessionsCell)
    data(k).sess  = S(ind,allSessionsCell{k});          % number of the sessions for that manipulation
    data(k).subj  = S(ind,1);                           % all the animals that received this manipulation
    for l = 1:length(allSessionsCell{k})
        sessLabels{k}    = strcat(conditions{allSessionsCell{k}});      %#ok<AGROW> % the label of that session
    end
    N = length(ind);
end

%% Initiating the analysis for all animals
M = [];

% ---  Testing if there were animals
if isempty(ind)
    disp('No animals were subjetct to all drug sessions, aborting');
    return;
end

for countSess = 1:length(allSessionsCell)
    countSub  = 1;
    allD{countSess} = [];                                       %#ok<AGROW>
    for k = 1:N                          % repeats for all animals in this condition
        thisSubject = data(countSess).subj(k);
        thisSession = data(countSess).sess(k,:);
        D = gatherDrrd(thisSubject,thisSession,false);
        %D = drrd(PREFIX, thisSubject, thisSession, false,false);
        %out = floor(size(D,1)/4);
        %out = 50;
        %out = aux;
        %D = D(out+1:end,:);
        
        % --- Gathering all the durations. Skipping in case ---
        % --- there was no data                             ---
        if ~isempty(D)
            allD{countSess} = [allD{countSess}; D(:,1)];            %#ok<AGROW>
        
        
        % --- Recording filenames just for verification ---
        M.filename{countSub,countSess} = [PREFIX '0' num2str(thisSubject) '.' num2str(thisSession)];
        
        % --- Peak of the distribution ---
        %M.peak(countSub,countSess) = find_peak(D,gauss,rng,dt);
        
        % --- peak with ksdensity ---
        f = ksdensity(D(:,1),rng);
        aux = find(f==max(f),1,'last');
        peak = rng(aux);
        M.peak(countSub,countSess) = peak;
        
        
        % --- Peak disconsidering premature responses ---
        cutoff = 0.5; % time in seconds
        Dcut = D(D(:,1)>cutoff,:);
        fcut = ksdensity(Dcut(:,1),rng);
        aux = find(fcut==max(fcut),1,'last');
        peakCut = rng(aux);
        M.peakCut(countSub,countSess) = peakCut;
        
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
        M.presses(countSub,countSess) = length(D(:,1))/length(thisSession);
        
        % --- peak normalized by opportunity ---
        M.peakByOpp(countSub,countSess) = ...
            compareSessionDistributionsByOppSingle(D,false);
        
        
              
        % --- STD of the distribution ---
        M.std(countSub,countSess) = std(D(:,1));
        
        % --- cumulative function ---
        %n = histc(D(:,1),rng_cum)/length(D(:,1));
        %M.cumsum{countSub,countSess} = cumsum(n);
        n = histc(Dcut(:,1),rng_cum)/length(Dcut(:,1));
        M.cumsum{countSub,countSess} = cumsum(n);

        % --- average ksdensity ---
        n = histc(D(:,1),rng_cum)/length(D(:,1));
        f = ksdensity(D(:,1),rng_cum);
        f = f(:);
        M.ksdensity{countSub,countSess} = f;
        
        
        else
            msg = ['Skipping subject ' num2str(data(1).subj(countSub))];
            warning(msg); 
        end
        
        countSub = countSub + 1;                        % increments subject number
    end
end
close all;
subplot(3,3,1); hold on;
[~,p] = anova_rm(M.peak);close;
df = size(M.peak);
df(2) = df(2)-1;
df(1) = df(2)*(df(1)-1);
df = ['(' num2str(df(2)) ',' num2str(df(1)) ')'];
plotM(M.peak,allSessions,sessLabels,'time (s)',{['Peak Pos. p='  num2str(p{2,6},'%.3f')],...
    ['f' df '= ' num2str(p{2,5},'%.3f')]});

subplot(3,3,2); hold on;
[~,p] = anova_rm(M.peakCut);close;
plotM(M.peakCut,allSessions,sessLabels,'time (s)',{['PeakCut Pos. p='  num2str(p{2,6},'%.3f')],...
    ['f' df '= ' num2str(p{2,5},'%.3f')]});

subplot(3,3,3); hold on;
[~,p] = anova_rm(M.median_reinf);close;
plotM(M.median_reinf,allSessions,sessLabels,'time (s)',{['Median Reinf p=' num2str(p{2,6},'%.3f')],...
    ['f' df '= ' num2str(p{2,5},'%.3f')]});

subplot(3,3,4); hold on;
[~,p] = anova_rm(M.median_all);close;
plotM(M.median_all,allSessions,sessLabels,'time (s)',{['Median all resps p=' num2str(p{2,6},'%.3f')],...
    ['f' df '= ' num2str(p{2,5},'%.3f')]});

subplot(3,3,5); hold on;
[~,p] = anova_rm(M.reinf_ratio);close;
plotM(M.reinf_ratio,allSessions,sessLabels,'ratio (%)',{['Reinf. ratio(%) p=' num2str(p{2,6},'%.3f')],...
    ['f' df '= ' num2str(p{2,5},'%.3f')]});

subplot(3,3,6); hold on;
[~,p] = anova_rm(M.presses);close;
plotM(M.presses,allSessions,sessLabels,'Number of events',{['# presses p=' num2str(p{2,6},'%.3f')],...
    ['f' df '= ' num2str(p{2,5},'%.3f')]});

subplot(3,3,7); hold on;
for k = 1:length(allSessionsCell)
    allks = horzcat(M.ksdensity{:,k});
    mks = mean(allks,2);
    stdks = std(allks,0,2)/sqrt(size(M.ksdensity,1));
    %n = histc(x,rng)/length(x);
    color = [k k k]/6;
    errorbar(rng_cum,mks,stdks,'color',color);
    %plot(rng,cumsum(n),'k-');
end


%subplot(3,3,7); hold on;
%[~,p] = anova_rm(M.peakByOpp);close;
%plotM(M.peakByOpp,allSessions,sessLabels,'time (s)',['Peak by opp. p=' num2str(p{2,6},'%.3f')]);

subplot(3,3,8); hold on;
[~,p] = anova_rm(M.std);close;
plotM(M.std,allSessions,sessLabels,'time (s)',{['STD  p='  num2str(p{2,6},'%.3f')],...
    ['f' df '= ' num2str(p{2,5},'%.3f')]});

subplot(3,3,9); hold on;
for k = 1:length(allSessionsCell)
    allcum = horzcat(M.cumsum{:,k});
    mcum = mean(allcum,2);
    stdcum = std(allcum,0,2)/sqrt(size(M.cumsum,1));
    %n = histc(x,rng)/length(x);
    color = [k k k]/6;
    errorbar(rng_cum,mcum,stdcum,'color',color);
    %plot(rng,cumsum(n),'k-');
end


%print -dpng -r300 lixo.png


function plotM(M,allSessions,sessLabels,x_axis_label,title_label)
%%
yl = 'auto';%[0 2.2];
bar(mean(M),'facecolor',[.5 .5 .5]);
%plot(M','o-','linewidth',2,'markerfacecolor','w','markersize',8);
errorbar(mean(M),std(M)/sqrt(size(M,1)),'k.','linewidth',2);
hold on;
% for k = 1: size(M,2)
%     plot(k,M(:,k),'ko');
% end
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
