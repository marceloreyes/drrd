function analyzeDrugEffect_organized
%function analyzeDrugEffect_organized


loadSessionsList;           % loads all infor about sessions and animals
                            % see loadSessionsList.m

% --- Analyzys for group Halo (goup 1)            
for thisGroup = halo
    for thisDose = [d1mg]
        
        
        peaks = get_peak_position(...
            group(thisGroup,thisDose).subjs,...
            group(thisGroup,thisDose).sessions);
        
        %M.peak(countSubject,countSession) = ...
        %    get_peak_position(group(thisGroup   
    
    end
end


function get_peak_position (subjs,sess)

% verifying if sessions and subjects are the same size
if size(subjs) ~= size(sess)
    return(NaN);
end




    
    
