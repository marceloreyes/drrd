close all;
clear M;
addpath ~/ufabc/prg/CP/

allSubjs = 13:75;

%badSubjs = [31];

%allSubjs = setdiff(allSubjs,badSubjs);
%M = NaN(length(allSubjs),3);
MM = zeros(0,5);
k = 1;
for thisSubj = allSubjs
    
    %D = gatherDrrd(thisSubj,1,1,false);
    D = drrd('AB1',thisSubj,1,0);
    if length(D)>100
        cp = cp_wrapper(D(:,1),1,3,6);
        
        initDt = round(D(1,5)*100)/100;     % avoid round errors
        
        switch initDt
            case(0.5)
                col = 1;
            case(0.7)
                col = 2;
            case(0.9)
                col = 3;
            case(1.1)
                col = 4;
            case(1.2)
                col = 5;
        end
        
        if length(cp)>2                % no change point
            %    changePoint = NaN;
            %else
            ct = cp(2:3,1);
            meanBef = mean(D(1:ct(1),1));
            meanAft = mean(D(ct(1)+1:ct(2),1));
            if meanBef < meanAft
                changePoint = cp(2,1);
            end
            
            M(k,:) = [thisSubj initDt changePoint]; %#ok<SAGROW>
            MM(length(nonzeros(MM(:,col)))+1,col) = changePoint;
            k = k + 1;
        end
    end
end


plot(M(:,2),M(:,3),'ko');