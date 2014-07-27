function [acc, best_inds]=Accuracy(centroids,inputs,outputs)
    nd=size(inputs,1);
    na=size(inputs,2);
    nc=size(centroids,1)/na;
    
    best_acc = 0;
    best_inds = zeros(nd,1);
    
    P = perms(1:nc);
    
    ite = length(P);
    
    if (ite>720)
        ite = 720;
    end
    
    for i = 1 : ite
        p_centroids = [];
        for j = 1:length(P(i,:))
            for k=1:na
                p_centroids((j-1)*na+k) = centroids((P(i,j)-1)*na+k);
            end
        end
        p_centroids = p_centroids';
        
        hits=0.0;
        distances=zeros(nd,nc);
        for i2=1:nd
            for j=1:nc
                for k=1:na
                    distances(i2,j)=distances(i2,j)+(inputs(i2,k)-p_centroids((j-1)*na+k))^2;
                end
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