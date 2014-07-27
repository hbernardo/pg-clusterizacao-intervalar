function fposition=hclusteringni(x,inputs,gama)
    nd=size(inputs,1);
    na=size(inputs,2);
    x=x';
    nc=size(x,1)/na;
    
    distances=zeros(nd,nc);
    for i2=1:nd
        ab = inputs(i2,:);
        for j=1:nc
            g = x((na*(j-1)+1):j*na)';
            distances(i2,j) = norm(g-ab);
        end            
    end

    inds = zeros(nd,1);
    for i3=1:nd
        [m,ind] = min(distances(i3,:));
        inds(i3)=ind;
    end
    
    
    
    Js=0;
    
    for k=1:nc
        distances = zeros(1,nd);
        
        g = x((na*(k-1)+1):k*na)';              
        
        np = 0;
        for i=1:nd
            if inds(i)==k
                np = np+1;
                ab = inputs(i,:);

                distances(i) = distances(i) + norm(g-ab);
            end
        end
        if np == 0
            Js=Inf;
        else
            Js=Js+sum(distances)/np;
        end
    end
    
    fposition=Js/nc;

% d = intervalDistance(x, inputs', gama);
% fposition = (objectiveFunctionKHM(d,2));
end