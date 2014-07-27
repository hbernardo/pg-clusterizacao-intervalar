function fposition=Live_fn(x,inputs)
    nd=size(inputs,1);
    na=size(inputs,2);
    nc=size(x,1)/na;

    sum1=0.0;
    w=10;

    for i=1:nd
        distances=zeros(1,nc);
        for j=1:nc
            for k=1:na
                distances(j)=distances(j)+(inputs(i,k)-x((j-1)*na+k))^2;
            end
            %distances(j)=sqrt(distances(j));
        end
        sum2=sum(distances);
        [m ind]=min(distances);
        not_min=distances;
        not_min(ind)=[];
        sum_not_min=sum(not_min);
        %sum1=sum1+(sum2+(w-1)*m);
        %sum1=sum1+(sum_not_min+w*m);
        %sum1=sum1+(sum2*m);
        sum1=sum1+(sum_not_min*m);
        %sum1=sum1+m;
        %sum1=sum1+((1-w)*sum_not_min+w*m);
        %sum1=sum1+((1-w)*sum2+w*m);
    end
    
fposition=sum1;%/nd;