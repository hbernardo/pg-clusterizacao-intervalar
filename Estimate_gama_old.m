function gama=Estimate_gama(Ss,inputs)
    nd=size(inputs,1);
    nc=size(Ss,2);
    
    m=1;
    delta=0.97;
    
    Jv_=zeros(nc*nd,1);
    Jv_1=zeros(nc*nd,1);
    
    gama_m=5*m;
    gama_m1=5*(m+1);
    
    for j=1:nc
        for i=1:nd
            S=Ss(i,j);
            Jv_((nd*(j-1))+i)=S^gama_m;
            Jv_1((nd*(j-1))+i)=S^gama_m1;          
        end        
    end
    m=m+1;
    
    while corr(Jv_,Jv_1)<delta
        Jv_=zeros(nc*nd,1);
        Jv_1=zeros(nc*nd,1);   
        x=zeros(nc*nd,1);

        gama_m=5*m;
        gama_m1=5*(m+1);

        for j=1:nc
            for i=1:nd
                S=Ss(i,j);
                Jv_((nd*(j-1))+i)=S^gama_m;
                Jv_1((nd*(j-1))+i)=S^gama_m1;          
            end        
        end
        m=m+1;
    end   
    
gama=gama_m;