function [ value ] = objectiveFunctionKHM( d, p )
    if ~exist('p', 'var') || isempty(p)
        p = 2;
    end

    n = size(d,2);
    k = size(d,1);

    value = 0;    
    for i=1:n
        sum = 0;
        for j=1:k
            sum = sum + 1.0/((d(j,i))^p);
        end        
        value = value + k/sum;
    end
end

