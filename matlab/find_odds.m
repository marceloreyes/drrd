function aveOdds = find_odds(D,N,Ncut)

count = 1;
allOdds = nan(Ncut, length(N));

for k = N
    odds = cp_reyes(D,k); 
    odds = odds(1:Ncut); 
    odds = odds/sum(odds); 
    allOdds(:,count) = odds(:);
    count = count+1;
    disp(k);
end

aveOdds = mean(allOdds,2);