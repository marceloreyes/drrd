function ret = plotDrrd(D, title_label)
%function plotDrrd(D)
% Each line of the matrix D is a trial
% Collumn 1 is the duration of the lever press
% Collumn 2 is the the time between the lever release and the next lever press (ITI)
% Collumn 3 is 1 for the reinforced trials
% Collumn 4 is 1 for trials where the light was on (valid trials)
% Collumn 5 shows the criterion (prime time) for each trial

if nargin == 1
	title_label = [];
end

primed  = 3;		% collum
valid   = 4;
primeT  = 5;
session = 6; 	% column with the session number
N = size(D,1);

% --- looking for the specific trials ---
validPrimed 	= find(D(:,primed)==1 	& D(:,valid)==1);
validNonPrimed 	= find(D(:,primed)==0 	& D(:,valid)==1);
invalid			= find(D(:, valid)==0); 

clf; hold on;

% --- plotting the prime times ---
%plot(D(:,primeT),1:N,'r','linewidth', 1.5);

% --- alternative: patch ---
%patch([ D(:,5); D(end,5); 0.00; 0.00], [1:N N+5 N+5 0], [.7 .8 .7] ,'EdgeColor' ,'none');% % [.7 .8 .7]
%patch([ D(:,5); D(end,5); 0.00; 0.00], [1:N N+5 N+5 0], [0 110 144]/255 ,'EdgeColor' ,'none');% % [.7 .8 .7]
patch([ D(:,5); D(end,5); 0.00; 0.00], [1:N N+5 N+5 0], [0.8 0.8 0.8] ,'EdgeColor' ,'none');% % [.7 .8 .7]

%plot(D(:,5),1:N,'r--','linewidth',2);

% --- Plotting the moving average of the lever press durations ---
%plot(movingAverage(D(:,1),20),1:N,'linewidth',2);

% --- Plotting each trial in a different style ---
plot(D(validPrimed,1)   ,validPrimed,   'k.','markersize',8, 'markerfacecolor','w');
plot(D(validNonPrimed,1),validNonPrimed,'k.','markersize',8, 'markerfacecolor','w');
%plot(D(validNonPrimed,1),validNonPrimed,'ko','markersize',5, 'markerfacecolor','w','linewidth',1);
plot(D(invalid,1)		,invalid,		'k.','markersize',8);

% --- setting up the scale and title ---
xlim([0 3]);
ylim([0 N+5]);
title(title_label);
set(gca,'box','on');

% --- printing the lines dividing the sessions ---
div = find(diff(D(:,session)));
for k = 1:length(div)
	plot(xlim,[div(k) div(k)],'k--');
end

% --- mounting return variable ---
ret = [length(validPrimed)/N length(validNonPrimed)/N length(invalid)/N] *100;


