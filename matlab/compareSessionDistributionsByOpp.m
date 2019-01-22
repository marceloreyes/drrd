function M = compareSessionDistributionsByOpp(animalID,sessions,plotFlag)
%function D = compareSessionDistrutions(animalID,sessions)
% example: D = gatherDrrd(1,1:9)
% runs the sessions 1 through 9 for animal 1.


MIN_NUM_TRIALS = 5;             % minimal number of trials for this analysis


if ~exist('plotFlag','var')
    plotFlag = true;
end

close all;
prefix = 'AB1';
dt = 0.02;
rng = 0:dt:15;          % range of times for binning the histogram
sigma = 0.2;            % starndard deviation of the gaussian for smoothing
% the histograms

perOppFlag = true;      % flag to see if wants to normalize per opportunity.

lnClr = {'k' 'r' 'm' 'g' 'c' 'y'};  % line colors for displayin multiple sessions

clrCount = 1;           % counter for number of sessions (also counts the
% color for lines)

% -- Build the gaussian function for future smootihg the histogram ---
gauss = dt/sqrt(2*pi())/sigma*exp(-0.5*((rng-mean(rng))/sigma).^2);

% --- Loop for all sessions ---
for k = sessions
    % --- Analyzing data from one particular session
    D = drrd(prefix,animalID,k,false);
    
    % --- Checking whether there were enough responses
    if size(D,1)< MIN_NUM_TRIALS 
        M(clrCount) = NaN;              % if too small, do not even try to analyze
    else
        % --- Counting the events for the bins
        n = histc(D(:,1),rng);
        if perOppFlag
            opp = sum(n) - [0 cumsum(n(1:end-1))'];
            opp = opp(:);           % make it column vector
            n = n./opp;             % devide the counts by the opportunities. Notice
            % that the divisions by zero will generates
            % NaN values, which will happen after the
            % time of the longest trial.
            n(isnan(n)) = 0;        % All NaN responses are replaced by zeros
        else
            n = n/dt/length(D(:,1));    % This normalizes by # of responses
        end
        
        C = conv(n,gauss,'same');       % convolves the histogram created with
        % gaussian, producing a smoother curve
        
        peakInd = findPeak(C);          % Look for all peaks in the series
        % If there are multiple peaks, have to
        % choose the highest
        if length(peakInd) > 1
            peakInd = peakInd(1:end-1);                 % eliminates position
            % of last peak (always present)
            goodPeak = C(peakInd) == max(C(peakInd));   % looks for highest peak
            peakInd = peakInd(goodPeak);                % eliminates other peaks
        end
        M(clrCount) = rng(peakInd)+dt/2;        % line vector with results from
        % each session in each column
        
        % --- Graphical part (only if flag is TRUE ---
        if plotFlag
            hold on;
            plot(rng(1:length(C))+(dt/2),C,'o-','markerfacecolor','w',...
                'color',lnClr{clrCount},'linewidth',2,'markersize',5);
            lgnd{clrCount} = ['session ' num2str(k,'%g')];
            plot(rng(peakInd)+dt/2,C(peakInd),'k+','markersize',20);
        end
        
        clrCount = clrCount+1;          % increments counter for session and
        % graph color
    end
end

if plotFlag
    legend(lgnd,'location','NW');
    xlim([min(rng) max(rng)]);
    plot([1.2 1.2],ylim,'k--');
    %xlim([0 4]);
end

function peakInd = findPeak(C)
x = diff(C);
ind = find((x(1:end-1).*x(2:end))<=0);
ind2 = x(ind)>0;
peakInd = ind(ind2);
peakInd = peakInd + 1;
