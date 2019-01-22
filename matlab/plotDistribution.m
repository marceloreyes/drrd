function hc = compareSessionDistrutions(D)
%function D = compareSessionDistrutions(animalID,sessions)
% example: D = gatherDrrd(1,1:9)
% runs the sessions 1 through 9 for animal 1.

%close all;
hold on; 

dt = 0.02;
rng = 0:dt:6;         % range of times for binning histogram
sigma = 0.2;

lnClr = {'k' 'r' 'm' 'g' 'c' 'y' [.3 .3 .3] [.4 .4 .4] [.5 .5 .5] [.6 .6 .6] [.7 .7 .7] [.8 .8 .8]};
maxClr = 12;

gauss = dt/sqrt(2*pi())/sigma*exp(-0.5*((rng-mean(rng))/sigma).^2);

count = 1;
hc = rng(:);


n = histc(D(:,1),rng);
n = n/length(D(:,1))/dt;

hc(:,count) = n(:);

C = conv(n,gauss,'same');
hold on;
disp(mod(count,maxClr));
plot(rng+(dt/2),C,'-','markerfacecolor','w',...
    'color',lnClr{mod(count-1,maxClr)+1},'linewidth',2,'markersize',5);
lgnd{count} = [num2str(length(D),'%g') ' trials'];

ind = find(C == max(C),1,'last');
peak(count,:) = [rng(ind) + dt/2 C(ind)];
disp(peak);


legend(lgnd,'location','NE');
xlim([min(rng) max(rng)]);
set(gca, 'box','on');

%% plotting the peak positions
%for k = 1:count-1
%    plot(peak(k,1),peak(k,2),'o-','markerfacecolor','w',...
%        'color',lnClr{mod(k-1,maxClr)+1},'linewidth',2);
%end


