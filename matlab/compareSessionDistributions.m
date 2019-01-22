function peak = compareSessionDistributions(animalID,sessions,useLastTrials,eliminFirst)
%function peak = compareSessionDistrutions(animalID,sessions)
% example: peak = compareSessionDistrutions(1,1:10);
% analyzes the sessions 1 to 10 from animal 1. It returns the maximum value
% of the distribution of responses
% example: D = compareSessionDistribution(1,1:9,30) % this will use only the 
% last 30 trials
% of the session. If the animal executed less than 30 it will use all
% trials
% peak = compareSessionDistribution(1,1:9,0,20); elimins the first 30 trials 
% of the session for analysis 

if ~exist('useLastTrials','var')
    useLastTrials = 0;              %here if useLastTrials == 0 means that all trials should be used
end

if ~exist('eliminFirst','var')
    eliminFirst = 0;              %here if useLastTrials == 0 means that all trials should be used
end

close all;
prefix = 'AB1';
dt = 0.02;
rng = 0:dt:6;         % range of times for binning histogram
sigma = 0.2;

lnClr = {'k' [.3 .3 .3] [.4 .4 .4] [.5 .5 .5] [.6 .6 .6] [.7 .7 .7] [.8 .8 .8] 'r' 'm' 'g' 'c' 'y' };
maxClr = 12; 

gauss = dt/sqrt(2*pi())/sigma*exp(-0.5*((rng-mean(rng))/sigma).^2);

count = 1;
hc = rng(:);

for k = sessions
    D = drrd(prefix,animalID,k,false);
    
    if useLastTrials ~= 0
        if useLastTrials < size(D,1) - 10
            D = D(end-useLastTrials+1:end,:);
            disp('Used last trials');
        else
            warning('animal produced less reponses than 10 responses');
        end
    elseif eliminFirst < size(D,1) -10 
            D = D((eliminFirst+1):end,:);
            disp('Eliminated first trials');
        else
            warning('Selected less than 10 responses: EliminFirst ignored');
    end     
    
    n = histc(D(:,1),rng);
    n = n/length(D(:,1))/dt;
    
    %stairs(rng,n,'-','color',lnClr{clrCount},'linewidth',2) ; hold on;
    
    hc(:,count) = n(:); 
    
    C = conv(n,gauss,'same');
    hold on;
    disp(mod(count,maxClr));
    plot(rng+(dt/2),C,'-','markerfacecolor','w',...
        'color',lnClr{mod(count-1,maxClr)+1},'linewidth',4,'markersize',5);
    lgnd{count} = ['session ' num2str(k,'%g')];
    
    ind = find(C == max(C),1,'last');
    peak(count,:) = [rng(ind) + dt/2 C(ind)];
    disp(peak);
    
    count = count+1;   
end

legend(lgnd,'location','NE');
xlim([min(rng) max(rng)]);
set(gca, 'box','on');
title(['Animal ' num2str(animalID,'%d')]);

%% plotting the peak positions
for k = 1:count-1
    plot(peak(k,1),peak(k,2),'o-','markerfacecolor','w',...
        'color',lnClr{mod(k-1,maxClr)+1},'linewidth',2);
end

figure; hold on;
plot(sessions,peak(:,1));

