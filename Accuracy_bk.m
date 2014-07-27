function acc=Accuracy(centroids,inputs,outputs)
    nd=size(inputs,1);
    na=size(inputs,2);
    nc=size(centroids,1)/na;
    
    hits=0.0;
    premap=zeros(nc,nc);
    distances=zeros(nd,nc);
    for i=1:nd
        for j=1:nc
            for k=1:na
                distances(i,j)=distances(i,j)+(inputs(i,k)-centroids((j-1)*na+k))^2;
            end
        end
        [m ind] = min(distances(i,:));
%         [ind outputs(i)]
        premap(ind,outputs(i)) = premap(ind,outputs(i))+1;
    end
    
    map=zeros(1,nc);
    for i=1:nc
        [m ind] = max(premap(i,:));
        map(i)=ind;
    end
    
    for i=1:nd
        [m ind] = min(distances(i,:));
        if (map(ind)==outputs(i))
             hits=hits+1.0;
         end
    end

acc=hits/nd;