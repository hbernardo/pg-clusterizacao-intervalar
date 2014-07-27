function fposition=hclustering(x,inputs,gama)
    nd=size(inputs,1);
    na=size(inputs,2)/2;
    x=x';
    nc=size(x,1)/(na*2);
    
    %adjusting centroids positions (beta greater than alfa) 
    x2 = [];
    for k=1:nc
        g = x((2*na*(k-1)+1):k*2*na)';

        for h=1:na
            if g(h*2-1)>g(h*2)
                aux = g(h*2-1);
                g(h*2-1) = g(h*2);
                g(h*2) = aux;
            end
        end

        x2 = vertcat(x2,g');
    end
    
    x = x2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    
    
    distances=zeros(nd,nc);
    for i2=1:nd
        ab = inputs(i2,:);
        a = ab(1:2:end);
        b = ab(2:2:end);
        for j=1:nc
            g = x((2*na*(j-1)+1):j*2*na)';
            alpha = g(1:2:end);
            beta = g(2:2:end);
            distances(i2,j) = 1.0/((exp(-(norm(a-alpha)^2+norm(b-beta)^2)/B)));
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
        
        g = x((2*na*(k-1)+1):k*2*na)';
        alpha = g(1:2:end);
        beta = g(2:2:end);               
        
        np = 0;
        for i=1:nd
            if inds(i)==k
                np = np+1;
                ab = inputs(i,:);
                a = ab(1:2:end);
                b = ab(2:2:end); 

                S = exp(-(norm(a-alpha)^2+norm(b-beta)^2)/B);
                f=S^gama;
                distances(i) = distances(i) + (1.0)/f;
%                 distances(i) = distances(i) + f;
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