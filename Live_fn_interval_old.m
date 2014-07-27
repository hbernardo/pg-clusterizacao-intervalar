function fposition=Live_fn_interval(x,inputs)
    nd=size(inputs,1);
    na=size(inputs,2)/2;
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
    %
    
    Js=0;
    
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
%         B=B+(pdist(vertcat(ab,ab_),'euclidean'));
        
        a = ab(1:2:end);
        b = ab(2:2:end);
        B=B+(norm(a-a_)^2+norm(b-b_)^2);
    end
    B=B/nd;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Ss = zeros(nd,nc);
    %calculate all S(xi,gi)
    for k=1:nc
        g = x((2*na*(k-1)+1):k*2*na)';
        alpha = g(1:2:end);
        beta = g(2:2:end);
        for i=1:nd
            ab = inputs(i,:);
%             S=exp(-((pdist(vertcat(ab,g),'euclidean')))/B);

            a = ab(1:2:end);
            b = ab(2:2:end);
            S = exp(-(norm(a-alpha)^2+norm(b-beta)^2)/B);
            
            Ss(i,k) = S;
        end
    end
    
    gama = Estimate_gama(Ss,inputs);
    for k=1:nc        
        for i=1:nd
            S=Ss(i,k);
            f=S^gama;
            Js=Js+f;
        end
    end
    
fposition=((1.0)/Js);