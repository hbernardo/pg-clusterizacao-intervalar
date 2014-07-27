function [Cindex] = c_index_interval(inputs, output)
    nd=size(inputs,1);
    k=length(unique(output));
    
    %%%%%%%%%%%%%%% B %%%%%%%%%%%%%%%%%
    B=0;
    ab_=0;
    for i3=1:nd
        ab = inputs(i3,:);
        ab_=ab_+ab;
    end
    ab_=ab_/nd;        
    a_ = ab_(1:2:end);
    b_ = ab_(2:2:end);    

    for i2=1:nd
        ab=inputs(i2,:);        
        a = ab(1:2:end);
        b = ab(2:2:end);
        B=B+(norm(a-a_)^2+norm(b-b_)^2);
    end
    B=B/nd;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    d = ones(nd,nd)*(-1);
    
    for j=1:(nd-1)
        ab=inputs(j,:);
        a = ab(1:2:end);
        b = ab(2:2:end);
        for i=(j+1):nd
            ab_=inputs(i,:);
            a_ = ab_(1:2:end);
            b_ = ab_(2:2:end); 
            d(j,i) = 1.0/((exp(-(norm(a-a_)^2+norm(b-b_)^2)/B)));
        end
    end
    
    clusters = zeros(k,nd);
    for i=1:k
        count=1;
        for j=1:nd
            if output(j)==i %j is in cluster i
                clusters(i,count) = j;
                count=count+1;
            end
        end
    end
    
    %Testado até aqui
    
    ele = 0;
    S=0;
    for h=1:k
        j=1;
        while (j<nd && clusters(h,j)~=0)
            ab=inputs(clusters(h,j),:);
            a = ab(1:2:end);
            b = ab(2:2:end);
            i=(j+1);
            while (i<=nd && clusters(h,i)~=0)
                ele=ele+1;
                ab_=inputs(clusters(h,i),:);
                a_ = ab_(1:2:end);
                b_ = ab_(2:2:end); 
                dista = 1.0/((exp(-(norm(a-a_)^2+norm(b-b_)^2)/B)));
                S=S+dista;
                i=i+1;
            end
            j=j+1;
        end
    end
    
    [sortedValues,sortIndex] = sort(d(:),'descend');
    maxValues = sortedValues(1:ele);
    Smax = sum(maxValues);
    
    [sortedValues,sortIndex] = sort(d(:),'ascend');
    sortedValues(sortedValues==-1) = [];
    minValues = sortedValues(1:ele);
    Smin = sum(minValues);
    
    Cindex = (S-Smin)/(Smax-Smin);        % C index, Hubert and Levin
end