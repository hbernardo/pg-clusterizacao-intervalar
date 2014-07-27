function gama=Estimate_gama(inputs)
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
%         B=B+(pdist(vertcat(ab,ab_),'euclidean'));

        a = ab(1:2:end);
        b = ab(2:2:end);
        B=B+(norm(a-a_)^2+norm(b-b_)^2);
    end
    B=B/nd;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    delta=0.9999;
    
    m=1;

    Jv_=zeros(nd,1);
    Jv_1=zeros(nd,1);

    gama_m=5*m;
    gama_m1=5*(m+1);

    for i=1:nd
        ab = inputs(i,:);
        a = ab(1:2:end);
        b = ab(2:2:end);
        for j=1:nd
            ab2 = inputs(j,:);
            a2 = ab2(1:2:end);
            b2 = ab2(2:2:end);
            S = exp(-(norm(a2-a)^2+norm(b2-b)^2)/B);
        
            Jv_(i)=Jv_(i)+S^gama_m;
            Jv_1(i)=Jv_1(i)+S^gama_m1;  
        end               
    end
    
    figure;
    plot(Jv_);
    hold on;
    pause;
    plot(Jv_1);    

    m=m+1;

    while corr(Jv_,Jv_1)<delta
        Jv_=zeros(nd,1);
        Jv_1=zeros(nd,1);

        gama_m=5*m;
        gama_m1=5*(m+1);

        for i=1:nd
            ab = inputs(i,:);
            a = ab(1:2:end);
            b = ab(2:2:end);
            for j=1:nd
                ab2 = inputs(j,:);
                a2 = ab2(1:2:end);
                b2 = ab2(2:2:end);
                S = exp(-(norm(a-a2)^2+norm(b-b2)^2)/B);

                Jv_(i)=Jv_(i)+S^gama_m;
                Jv_1(i)=Jv_1(i)+S^gama_m1;  
            end               
        end
        
        pause;
        plot(Jv_1);

        m=m+1;
    end
    
    hold off;
    pause;
    figure;
    plot(Jv_);
    gama_m
    pause;
       
gama=gama_m;