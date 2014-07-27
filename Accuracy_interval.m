function [acc, best_inds]=Accuracy_interval(centroids,inputs,outputs)
    nd=size(inputs,1);
    na=size(inputs,2)/2;
    nc=size(centroids,1)/(na*2);
    
    %adjusting centroids positions (beta greater than alfa)
    centroids2 = [];
    for k=1:nc
        g = centroids((2*na*(k-1)+1):k*2*na)';

        for h=1:na
            if g(h*2-1)>g(h*2)
                aux = g(h*2-1);
                g(h*2-1) = g(h*2);
                g(h*2) = aux;
            end
        end

        centroids2 = vertcat(centroids2,g');
    end

    centroids = centroids2;
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
    
    
    best_acc = 0;
    best_inds = zeros(nd,1);
    
    P = perms(1:nc);
    
    ite = length(P);
    
    if (ite>720)
        ite = 720;
    end
    
    for i = 1 : ite
        p_centroids = zeros(1,nc*na*2);
        for j = 1:length(P(i,:))
            for k=1:na
                p_centroids((j-1)*na*2+k*2-1) = centroids((P(i,j)-1)*na*2+k*2-1);
                p_centroids((j-1)*na*2+k*2) = centroids((P(i,j)-1)*na*2+k*2);
            end
        end
        p_centroids = p_centroids';
        
        hits=0.0;

        distances=zeros(nd,nc);
        for i2=1:nd
            ab = inputs(i2,:);
            a = ab(1:2:end);
            b = ab(2:2:end);
            for j=1:nc
                g = p_centroids((2*na*(j-1)+1):j*2*na)';
                alpha = g(1:2:end);
                beta = g(2:2:end);
                distances(i2,j) = 1.0/((exp(-(norm(a-alpha)^2+norm(b-beta)^2)/B)));%(norm(a-alpha)^2+norm(b-beta)^2);
%                 distances(i,j) = pdist(vertcat(ab,g),'euclidean');
            end            
        end

        inds = zeros(nd,1);
        for i3=1:nd
            [m,ind] = min(distances(i3,:));
            inds(i3)=ind;
            if (ind==outputs(i3))
                 hits=hits+1.0;
            end
        end
        
        ac = hits/nd;
        if (ac>best_acc)
            best_acc=ac;
            best_inds = inds;
        end
        
    end

acc=best_acc;