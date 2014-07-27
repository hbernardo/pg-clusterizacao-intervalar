function fposition=Live_fn_interval(x,inputs)
    nd=size(inputs,1);
    na=size(inputs,2)/2;
    nc=size(x,1)/2;

    Js=0;
    for k=1:nc
        alpha=x(k*2-1);
        beta=x(k*2);
        gama=Estimate_gama(alpha,beta,inputs);
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
                
                S=exp(-((a-alpha)^2+(b-beta)^2)/B);
                f=S^gama;
                Js=Js+f;
            end
        end
    end
    
fposition=(1/Js);