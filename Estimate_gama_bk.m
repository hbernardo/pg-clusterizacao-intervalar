function gama=Estimate_gama(alpha,beta,inputs)
    nd=size(inputs,1);
    na=size(inputs,2)/2;
    
    m=1;
    delta=0.97;
    
    J_=0;
    
    while J_<delta
        J_=0;
        gama_m=5*m;
        for i=1:nd
            for j=1:na                
                %%%%%%%%%%%%%%% B %%%%%%%%%%%%%%%%%
                B=0;
                a_=0;
                b_=0;
                for i3=1:nd
                    for j3=1:na
                        a=inputs(i3,(j3*2-1)); 
                        b=inputs(i3,j3*2);

                        a_=a_+a;
                        b_=b_+b;
                    end
                end
                a_=a_/(na*nd);
                b_=b_/(na*nd);

                for i2=1:nd
                    for j2=1:na
                        a=inputs(i2,(j2*2-1));
                        b=inputs(i2,j2*2);

                        B=B+((a-a_)^2+(b-b_)^2);
                    end
                end
                B=B/(na*nd);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                a=inputs(i,(j*2-1));
                b=inputs(i,j*2);

                f=(exp(-((a-alpha)^2+(b-beta)^2)/B))^gama_m;
                J_=J_+f;
            end
        end
        m=m+1;
        J_=1/J_;
    end   
    
gama=gama_m;