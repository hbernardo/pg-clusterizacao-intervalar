function [ C, objvalue ] = KHM( initialC, X, p, iterStop, minChanging)
    if ~exist('p', 'var') || isempty(p)
        p = 2;
    end
    if ~exist('iterStop', 'var') || isempty(iterStop)
        iterStop = 100;
    end
    if ~exist('minChanging', 'var') || isempty(minChanging)
        minChanging = 0.01;
    end
    
    C = initialC;

    n=size(X,1);
    na=size(X,2)/2;
    k = length(C)/(na*2);    

    d = euclideanDistance(C, X', 1);
    
    objvalue = objectiveFunctionKHM(d,p);
    
    changing = 1;
    h=0;
    while changing > minChanging && h < iterStop
        m = zeros(k,n);
        for i=1:n
            for j=1:k
                sum=0;
                for j2=1:k
                    sum = sum + d(j2,i)^(-p-2);
                end
                m(j,i) = (d(j,i)^(-p-2))/sum;
            end
        end
        
        w = zeros(1,n);
        for i=1:n
            sum1=0;
            for j=1:k
                sum1 = sum1 + d(j,i)^(-p-2);
            end
            sum2=0;
            for j=1:k
                sum2 = sum2 + d(j,i)^(-p);
            end
            w(i) = sum1/(sum2^2);
        end
        
        for j=1:k
            sum1 = zeros(1,na*2);
            sum2 = 0;
            for i=1:n
                mw = m(j,i)*w(i);
                sum1 = sum1 + mw*X(i,:);
                sum2 = sum2 + mw;
            end
            C((2*na*(j-1)+1):j*2*na) = sum1/sum2;
        end
        
        d = euclideanDistance(C, X', 1);
        
        last_objvalue = objvalue;
        
        objvalue = objectiveFunctionKHM(d,p);
        
        changing = (abs(objvalue-last_objvalue))/last_objvalue;
        
        h=h+1;
    end
end

