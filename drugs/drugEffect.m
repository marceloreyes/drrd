sal = 6;
hal = 7;
group = 7:15;
M=[];
for s = group
    D = drrd('AB1',s,sal,false); 
    M(end+1,1) = median(D(:,1)); 
    D = drrd('AB1',s,hal,false); 
    M(end,2) = median(D(:,1));
end

