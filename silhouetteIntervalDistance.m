function [ d ] = silhouetteIntervalDistance( dp, inputs )
    nd=size(inputs,1);
    
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
    
    d = zeros(nd,1);
    
    a = dp(1:2:end);
    b = dp(2:2:end);
    for i=1:nd
        ab_=inputs(i,:);
        a_ = ab_(1:2:end);
        b_ = ab_(2:2:end); 
        d(i) = 1.0/((exp(-(norm(a-a_)^2+norm(b-b_)^2)/B)));
    end
end

