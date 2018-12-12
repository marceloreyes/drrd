function D = gatherDrrd(prefix,animalID,sessions,plotFlag)
%function D = gatherDrrd(animalID,sessions)
% prefix = the code for the experiment
% example: D = gatherDrrd('AB1',1,1:9)
% runs the sessions 1 through 9 for animal 1 of the AB1 experiment.

if ~exist('plotFlag','var')
    plotFlag = true;
end

D = [];

for k = sessions
    D = [D ;drrd(prefix,animalID,k,plotFlag)];
end;

if plotFlag
    plotDrrd(D,['Animal: ' num2str(animalID) ' sessions: ' num2str(sessions)]);
end
