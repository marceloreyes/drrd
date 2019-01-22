
AllSess   =  1:2;                           % number of sessions to be gathered
Nma     = 20;                           % number of points for moving average (20)
Ntrials = 100;                          % max number of trials to analyze performance
NN      = 50;                           % for analysis of the first NN and last NN trials
%subjs = [8:23 25:28 34:75];            % subjects to be analyze
subjs = 1:75;
d1 = [];
d5 = [];
d7 = [];
d9 = [];
d11= [];
d12= [];

for k = subjs
    close all;
    D = gatherdrrd(k,AllSess,false);          % put together the first sessions
    if size(D,1) < Ntrials
        warning('Too few trials');
    else
        %D = D(1:Ntrials,:);
        strt = D(1,5);                     % initial criteria
        if round(strt*1000) == 100 || round(strt*1000) == 102
            d1(:,end+1).ma = movingAverage(D(:,1),Nma);
            d1(:,end+1).beg = mean(D(1:NN,1));
            d1(:,end+1).end = mean(D(end-NN+1:end,1));
        elseif round(strt*1000) == 500
            d5(:,end+1).ma = movingAverage(D(:,1),Nma);
            d5(:,end+1).beg = mean(D(1:NN,1));
            d5(:,end+1).end = mean(D(end-NN+1:end,1));
        elseif round(strt*1000) == 700
            d7(:,end+1).ma = movingAverage(D(:,1),Nma);
            d7(:,end+1).beg = mean(D(1:NN,1));
            d7(:,end+1).end = mean(D(end-NN+1:end,1));
        elseif round(strt*1000) == 900
            d9(:,end+1).ma = movingAverage(D(:,1),Nma);
            d9(:,end+1).beg = mean(D(1:NN,1));
            d9(:,end+1).end = mean(D(end-NN+1:end,1));
        elseif round(strt*1000) == 1100
            d11(:,end+1).ma = movingAverage(D(:,1),Nma);
            d11(:,end+1).beg = mean(D(1:NN,1));
            d11(:,end+1).end = mean(D(end-NN+1:end,1));
        elseif round(strt*1000) == 1200
            d12(:,end+1).ma = movingAverage(D(:,1),Nma);
            d12(:,end+1).beg = mean(D(1:NN,1));
            d12(:,end+1).end = mean(D(end-NN+1:end,1));
        end
    end
end


clf; hold on
% plot(mean(horzcat(d1.ma),2),'k');
% plot(mean(horzcat(d5.ma),2),'r');
% plot(mean(d7,2),'k');
% plot(mean(d9,2),'b');
% plot(mean(d11,2),'m');
% plot(mean(d12,2),'c');

%%
clf; hold on;
M = [ mean([d1.beg])  mean([d1.end])
      mean([d5.beg])  mean([d5.end])];
S= [ std([d1.beg])  std([d1.end])
     std([d5.beg])  std([d5.end])];
S = S/sqrt(length([d1.beg]));

M = M';
S = S';

d = 0.14;
bar([1 2],M,1);
errorbar([1-d 2-d],M(:,1),S(:,1),'k.','markersize',0.01,'linewidth',2);
errorbar([1+d 2+d],M(:,2),S(:,2),'k.','markersize',0.01,'linewidth',2);


colormap([1 .5 .0; 0.2 0.2 1]);               % color adjustment
ylabel('Tempo (s)','fontsize',18);          
ylim([0 1.2]);
set(gca,'xtick',[1 2],'box','on');
set(gca,'xticklabels',{'in�cio', 'final'},'fontsize',16);
legend({'Com Timeout','Sem Timeout'},'Location','NorthWest');

%print -dpng -r300 trainingDifferences.png

%%
figure; hold on; clear M S;
M = [ mean([d12.beg])  mean([d12.end])
      mean( [d5.beg])  mean([d5.end])];
S= [ std([d12.beg])  std([d12.end])
     std([d5.beg])  std([d5.end])];
S = S/sqrt(length([d12.beg]));

M = M';
S = S';

d = 0.14;
bar([1 2],M,1);
errorbar([1-d 2-d],M(:,1),S(:,1),'k.','markersize',0.01,'linewidth',2);
errorbar([1+d 2+d],M(:,2),S(:,2),'k.','markersize',0.01,'linewidth',2);


colormap([1 .5 .0; 0.2 0.2 1]);               % color adjustment
ylabel('Tempo (s)','fontsize',18);          
ylim([0 1.5]);
set(gca,'xtick',[1 2],'box','on');
set(gca,'xticklabels',{'in�cio', 'final'},'fontsize',16);
legend({'Com Timeout','Sem Timeout'},'Location','NorthWest');
