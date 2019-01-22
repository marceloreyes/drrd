function ret = fitQuality(animalID,sessions,dt)
%function hc = compareSessionDistributions(animalID,sessions)

%function D = compareSessionDistrutions(animalID,sessions)
% example: D = gatherDrrd(1,1:9)
% runs the sessions 1 through 9 for animal 1.

close all;
prefix = 'AB1';

if ~exist('dt','var')
    dt = 0.02;               % delta T for histogram calculation 0.02
end
rng = 0:dt:3;           % range of times for binning histogram

%sigma = 0.2;            % standard deviation for the gaussian for smoothing

%lnClr = {'k' 'r' 'm' 'g' 'c' 'y' [.3 .3 .3] [.4 .4 .4] [.5 .5 .5] [.6 .6 .6] [.7 .7 .7] [.8 .8 .8]};
%maxClr = 12;

%gauss = dt/sqrt(2*pi())/sigma*exp(-0.5*((rng-mean(rng))/sigma).^2);

%hc = rng(:);


D = gatherDrrd(animalID,sessions,false);
if isempty(D) || ~isempty(find(D(:,1)<0, 1))
    warning('No data was returned or negative values found');
    ret = NaN;
else
    n   = histc(D(:,1),rng);
    n   = n/length(D(:,1))/dt;
    n   = n(:);                     % transforms to a column vector
    rng = rng+dt/2;
    ksn = ksdensity(D(:,1),rng,'support','positive');
    ksn = ksn(:);
    close all; hold on;
    for k = 1:length(rng)
        plot([rng(k) rng(k)],[n(k) ksn(k)],'r-','linewidth',0.5);
    end
    %plot(rng,  n,'ko','linewidth',2,'markerfacecolor','w');
    %plot(rng,ksn,'k--','linewidth',2);
    
    ret = w2(n,ksn);
    %disp([mean(n-ksn) mean((n-ksn).^2) ret]);
end


function ret = w2(Observed, Predicted)

Mean_obs = mean(Observed);

Total = sum((Observed-Mean_obs).^2);

Unexplained = sum((Observed-Predicted).^2);

Explained = Total-Unexplained;

ret = Explained/Total;




%ret = omegaSquare(n,ksn);



%stairs(rng,n,'-','color',lnClr{clrCount},'linewidth',2) ; hold on;

%hc(:,count) = n(:);         % gets the histogram counts (hc) in a columnn vector
%C = conv(n,gauss,'same');   % convolves the histogram counts with a gaussian for smoothig
%hold on;
%disp(mod(count,maxClr));
% plot(rng+(dt/2),C,'o-','markerfacecolor','w',...    % plots the smooth function
%     'color',lnClr{mod(count-1,maxClr)+1},...
%     'linewidth',2,'markersize',5);
% lgnd{count} = ['session ' num2str(k,'%g')];
%
% stairs(rng,hc,'k-','linewidth',1);
%
% ind = find(C == max(C),1,'last');
% peak(count,:) = [rng(ind) + dt/2 C(ind)];
% disp(peak);
%
%
% legend(lgnd,'location','NE');
% xlim([min(rng) max(rng)]);
% set(gca, 'box','on');

%% plotting the peak positions
% for k = 1:count-1
%     plot(peak(k,1),peak(k,2),'o-','markerfacecolor','w',...
%         'color',lnClr{mod(k-1,maxClr)+1},'linewidth',2);
% end
%
% figure;
% plot(hc-C);
% disp(mean(hc-C));
% disp(mean(hc-C).^2);
%

%figure; hold on;
%plot(sessions,peak(:,1));

