function y = movingAverage(x,N)
% --- calculates a moving average of the input vector (x) using N values
% the input is a single column vector
% output: single column vector with N less entries than x 

x = x(:);
if N >= length(x)
    y = ones(length(x),1)*mean(x);
    return;
end

%y = NaN(size(x,1)-N, N);

% --- making a matrix with the "delayed" vector ---	
for i = 1:N
	y(:,i) = x(i:end-N+i);
end
y = mean(y,2);                                  % averaging the values	
y 	= vertcat(ones(int32(N/2),1)*y(1), y);		% padding the vector to have the same size as x
y(end+1:length(x)) 	= y(end);                   % same at the end 
y = y(:);
end
