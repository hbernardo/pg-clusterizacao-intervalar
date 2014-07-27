function d = euclideanDistance(Xx, Y, gama)
    inputs = Y';
    x = (reshape(Xx, 1, numel(Xx)))';

    nd=size(inputs,1);
    na=size(inputs,2);
    nc=size(x,1)/na;
    
    Ss = zeros(nd,nc);
    %calculate all S(xi,gi)
    for k=1:nc
        g = x((na*(k-1)+1):k*na)';
        for i=1:nd
            ab = inputs(i,:);
%             S=exp(-((pdist(vertcat(ab,g),'euclidean')))/B);
            S = norm(g-ab);            
            Ss(i,k) = S;
        end
    end
    
d=Ss';