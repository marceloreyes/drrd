
% ___________________________________________________________________________________________
% File:             drrd.m
% File type:        Function
% Created on:       September 25, 2013
% Created by:       Marcelo Bussotti Reyes
% Last modified on: August 29, 2018
% Last modified by: Marcelo Bussotti Reyes
% Modification:     Save variable D as empty matrix when med2tec returns an
%                   empty matrix. This indicates that the file parsed was
%                   not found. Previously it generated an error.
%	
% Purpose:          A function that analyzes the performance of rats in the drrd
%                   procedure. The animals are supposed to press a lever for a
%					duration longer than the prime time in order to receive food.
%
% Input:            File name prefix as a string, animal id as a number and 
%					session as a number. Filene name are in the format:
%					[prefix00A.00S] where A is the animal and S is the sesson 
%
% Output:           D, a matrix with 6 columns containing data from a trial in 
%					each line.
%                   In case saveMatFlag is parsed as true, the program will
%                   save a matlab file (.mat) with the same name as the
%                   original data file.
% Coments:          Uses functions: med2tec.m, and plotDrrd.m
%
% Format:           drrd(prefix,animalID,session))
% Example:          D = drrd('AB1',1,2'); 
%					this will analyze the file AB1001.002, animal 1 and session 2 
% Previous Modifications:
% Modification:     Included the option of saving a matlab file with the
%                   matrix D. This helps to speed up the analysis when
%                   several sessions and animals are to be analyzed.
%                   Included session column (6th) and the input format
%                   (Marcelo B. Reyes 2013)

%Todo
% To make sure that the number of events is correct, e.g. if the number of trials
% when the lever was pressed is the same as the releases. There is a possible problem
% that the animal may be pressin the lever when the trial starts. This would produce
% a lever release before a lever press (maybe?)
%
% make a correction for the animal ID, in case it is larger than 9 (two or three
% digits.

function D = drrd(prefix, animalID, session, plotFlag,saveMatFlag)
%function D = drrd(prefix, animalID, session, plotFlag,saveMatFlag)


filename = makeFileName(prefix,animalID,session);

% --- Indexes for each of the output variable columns
dtCol 		= 1;
itiCol		= 2;
primedCol 	= 3;
validCol	= 4;
phaseCol	= 5;
sessionCol	= 6;			% session variable column index

dPh			= 0.1;			% initiates the variable just in case it cannot be 
							% obtained from the data, as for example when no 
							% trials were primed (reinforced) in one session

if ~exist('filename','var')
    disp('Missing input variable: name of the file to be analyzed');
end
if ~exist('plotFlag','var')
    plotFlag = true;
end
if ~exist('saveMatFlag','var');
    saveMatFlag = false;
end


data = med2tec(filename); 		% reads data from medpc format to time-event code
if isempty(data)
    fprintf('No data analyzed\n');
    D = [];
    return;
end


% small correction for a bug in the med-pc file ---
% if the animal presses the lever at the same cycle of the Start command in the
% box, the first time can be registered wrong, so this sets it to zero as it 
% should be
data = fix_clock_reset(data);

% --- look for indexes of temporal events ---
startIndex      = find(data(:,2)== 1);  
endIndex        = find(data(:,2)== 3);  
primeIndex      = find(data(:,2)==18); 
lightOnIndex    = find(data(:,2)==11);
lightOffIndex   = find(data(:,2)==21);
phaseAdvIndex   = find(data(:,2)==17);        % indexes of trials where phase was advanced
phaseBckIndex   = find(data(:,2)==27);        % indexes of trials where phase was retreated

% checking if trials start with an startIndex, otherwise eliminates the
% first endIndex (there is no way to measure the tive of the trial if we do
% not know when it started)
[startIndex,endIndex] = checkAndFixIndices(startIndex,endIndex);


% --- searching for trials in which the animals received food. We call these "primed" ----
primedTrials	= findTrial(startIndex,primeIndex);

% --- searching for trials in which animals progressed or retreated phase
phaseAdvTrials	= findTrial(startIndex,phaseAdvIndex);
phaseBckTrials  = findTrial(startIndex,phaseBckIndex);

% --- searching for trials in which the animals responded with the light on 
% (not in timeout). We'll call these trials "valid". 
validTrials     = findValidTrial(startIndex,lightOnIndex,lightOffIndex);

% --- search for valid trials in which animals were and were not reinforded ---
validPrimed 	= intersect(validTrials,primedTrials);	% looks for trials when the animals received food and 
validNonPrimed 	= setdiff(  validTrials,primedTrials);
invalid			= setdiff(1:length(startIndex), validTrials);

% --- gets the initial prime time ---
if length(primeIndex) >= 1
	iniPh 			= data(primeIndex(1),1) - data(primeIndex(1)-1,1);
else
	iniPh = 1.2;
end
% --- Organizing data in one single matriz: D --- 
D = zeros(length(startIndex),6);        	% Initiates the vector for speed

% --- Calculating the duration of the lever presses ---
D(:,dtCol)  				= data(endIndex,1) - data(startIndex,1);
D(1:end-1,itiCol) 			= data(startIndex(2:end),1) - data(endIndex(1:end-1),1);
D(end,itiCol) 				= NaN;
D(primedTrials,primedCol)   = 1;             % sets to 1 all the trials that were primed
D(validTrials ,validCol)	= 1;             % sets to 1 all the trials that were primed
D(phaseAdvTrials,phaseCol)	=  dPh;
D(phaseBckTrials,phaseCol) 	= -dPh;
D(:,phaseCol) 	= cumsum(D(:,phaseCol))+iniPh;
D(:,sessionCol)		= session;			% puts a mark (1) on the last trial showing that 
									% it was the end of session (eos)
									
% --- graphical part ---
if plotFlag
    hold on;
    plotDrrd(D,filename);
end
%__________________________________________________________________________

N  = length(D(:,1));
vp = length(validPrimed);
vnp= length(validNonPrimed);
iv = length(invalid);

fprintf('Rat%d Trials:%d vreinf:%d(%.1f%%) vnreinf:%d(%.1f%%) Inv:%d(%.1f%%)\n',animalID,N,vp,vp/N*100,vnp,vnp/N*100,iv,iv/N*100);

% --- saving matlab file in case it was requested ---
if saveMatFlag
    save([filename '.mat'],'D');
end
%__________________________________________________________________________


function ret = findTrial(st, v)
% --- Looks for the trial in which the events occurred
% v is a list of indexes of temporal events. For example, if you know
% that an event ocurred in the index 102, this function will look in the
% indexes of the starts of the trial (st) and count the number of trials
% that had started before that particular event, let's say N, and hence 
% return that the event belongs to that trial N. If v is a vector, the
% function returns all the trials the events belong to

v = v(:);
ret = nan(size(v));
for k = 1:length(v)
    ret(k) = length(st(st<v(k)));
end
%__________________________________________________________________________


function ret = findValidTrial(st,u,v)
% --- Looks for trials that occurred between the events u and v. The most typical
% example is in the drrd trials: only the lever presses that ocurred while the
% light was on are valid. Hence, this function is used to find the trials that
% started after the light was turned on, and before the light was turned off.
% If there were multiple events, e.g. light turns on and off more than once, the 
% junction will look for events that are in between u(i) and v(i) where i is the
% ith time the light was turned on.

u = u(:); v = v(:);				% makes sure the inputs are column vectors
if length(u) ~= length(v)		% makes sure that the events have the same size
	if length(u) == length(v)+1;
		v(end+1) = st(end)+1; 	% just makes the last event one index after the last data set
	else
		disp('Incompatible number of events');
		exit(0);
	end
end

ret = [];						% initializes the return variable as empty

for k = 1:length(u)							% loop for all the events
    Nu = length(st(st<u(k)));   			% select the starts that happened before the event u(k)     
    Nv = length(st(st<v(k)));        		% select the starts that happened before the event v(k),
    										% supposedly after u(k)    
    ret = vertcat(ret,(Nu+1:Nv)');  		%#ok<AGROW> % adds to the return variable the trials between u(k) and v(k) 
end
%__________________________________________________________________________

function filename = makeFileName(prefix, animalID, session)
    % --- Function to put the parts of the file name together ---

    animalID = num2str(animalID/1000,'%.3f'); 
    animalID = animalID(3:end);			% same with animalID to 3 digits
    session  = num2str(session/1000,'%.3f'); 
    session  = session(2:end);			% converts the session to .+3 digits
    filename=[int16(prefix) animalID session]; 	% put together the parts 
%__________________________________________________________________________

function ret = fix_clock_reset(data)
    %sometimes the clock reset of the hardware (medPC) happens after
    % the beginning of the experiment. Typically it happens in the second
    % trial, but we detected that is happened later. This small routine 
    % fix this bug by eliminating all the trials before the reset. 

    % x must be a nx2 array

    % --- detecting if the event times in the fist column are allways increasing
    % and keep the last time it happened
    ind = find(diff(data(:,1))<-10,1,'last');

    if ~isempty(ind)
        fprintf('Warning: eliminated first %d trial(s) due to clock reset issue\n\n',ind)
        ret = data((ind+1):end,:);
    else
        ret = data;
    end

%__________________________________________________________________________

function [startIndex,endIndex] = checkAndFixIndices(startIndex,endIndex)
    % Fist checking if the data begins in the right order: with a
    % trial-start. Hence if the fist entry of the endIndex is smaller than
    % the startIndex, there is something wrong
    if endIndex(1)<startIndex(1)
        ind = find(endIndex<startIndex(1),1,'last');   % finds all start smaller than
        if ~isempty(ind)                                % if finds any
            if length(endIndex)>ind                     % checks if there are more indices to keep
                endIndex = endIndex(ind+1:end);           % cuts the first trials
            end
        end
    end
        
    % now ckecking if there are trials whose start was recorded but not the end    
    startIndex      = startIndex(1:length(endIndex));	% eliminates the last trial in case it was incomplete

