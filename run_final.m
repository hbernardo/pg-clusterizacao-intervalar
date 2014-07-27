function [ results ] = run_final( inputs, output, name)
    results = zeros(4,16);

    disp(strcat('!!! ',name, '!!!'));
    
    nd = length(inputs);
    na = (length(inputs(1,:))/2); % número de atributos intevalares
    nc = length(unique(output));

    % Soluções do k-harmonic means para diferentes números de clusters
    n_obv = 100; %número de observações para se tirar a média
    obvalues = zeros(6,n_obv);
    obv_min = zeros(1,6);
    obv_avg = zeros(1,6);

    inferior = min(inputs(:));
    superior = max(inputs(:));
    for k=2:6
        dim = 2*k*na;
        VRmin=ones(dim,1)*superior/2;
        VRmax=ones(dim,1)*superior;
        C = zeros(dim,1);
        for iter=1:n_obv
            for d=1:dim
                C(d,1) = unifrnd(VRmin(d), VRmax(d));
            end
            [centroids, value] = KHMinterval( C, inputs, 2);
            [acc cidx] = Accuracy_interval(centroids,inputs,output);
            cindex = c_index_interval(inputs, cidx);
            obvalues(k,iter) = cindex;
        end
        obv_avg(k) = sum(obvalues(k,:))/n_obv;
    end

    % Plotando c-index para diferentes números de clusters
    figure;
    c_index_values = obv_avg(2:6)
    plot(2:6,c_index_values);
    title(name);
    xlabel('Número de clusters');
    ylabel('C-index');

   
    % comparação de resultados
    
    steps = 10;
    steps2 = 10;
    
    % Tirando a média
    med_inputs = zeros(nd,na);

    for i=1:nd
        for j=1:na
           new = (inputs(i,(j*2-1))+inputs(i,j*2))/2;
           med_inputs(i,j)=new;
        end
    end
    %%%
    
%     % k-means
%     disp('==> k-means');
%     col = 1;
%     
%     accs = zeros(1,steps); 
%     sils = zeros(1,steps); 
%     cahs = zeros(1,steps); 
%     dabs = zeros(1,steps); 
%     solutions = zeros(nd, steps); 
% 
%     for i=1:steps
%         [Cluster Codebook] = cvKmeans(med_inputs',nc,.05,@cvEucdist);
%         centroids = (reshape(Codebook, 1, numel(Codebook)))';
%         
%         [accs(1,i) cidx] = Accuracy(centroids,med_inputs,output);
%         solutions(:,i) = cidx;
%         eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
%         cahs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
%         dabs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'silhouette');
%         sils(1,i) = eva.CriterionValues;
%     end
%     
%     cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
%     [acc_max, amax_id] = max(accs);
%     
%     [cah_min, cmin_id] = min(cahs);
%     [cah_max, cmax_id] = max(cahs);
%     
%     [dab_min, dmin_id] = min(dabs);
%     [dab_max, dmax_id] = max(dabs);
%     
%     [sil_min, smin_id] = min(sils);
%     [sil_max, smax_id] = max(sils);
%     
% %     worst_solution = solutions(:,amin_id)
% %     best_solution = solutions(:,amax_id)
%     
% %     acc_min
% %     acc_max
%     acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
%     
% %     cah_min
% %     cah_max
%     cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
%     
% %     dab_min
% %     dab_max
%     dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
%     
% %     sil_min
% %     sil_max
%     sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
    
    % k-means intervalar
    disp('==> k-means intervalar');
    col=2;
    
    accs = zeros(1,steps); 
    sils = zeros(1,steps); 
    cahs = zeros(1,steps); 
    dabs = zeros(1,steps); 
    solutions = zeros(nd, steps); 
    
    for i=1:steps
        [Cluster Codebook] = cvKmeans(inputs',nc,.05,@intervalDistance);
        centroids = (reshape(Codebook, 1, numel(Codebook)))';
        
        [accs(1,i) cidx] = Accuracy_interval(centroids,inputs,output);
        solutions(:,i) = cidx;
        eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
        cahs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
        dabs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'silhouette');
        sils(1,i) = eva.CriterionValues;
    end

    cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
    [acc_max, amax_id] = max(accs);
    
    [cah_min, cmin_id] = min(cahs);
    [cah_max, cmax_id] = max(cahs);
    
    [dab_min, dmin_id] = min(dabs);
    [dab_max, dmax_id] = max(dabs);
    
    [sil_min, smin_id] = min(sils);
    [sil_max, smax_id] = max(sils);
    
%     worst_solution = solutions(:,amin_id)
%     best_solution = solutions(:,amax_id)
    
%     acc_min
%     acc_max
    acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
    
%     cah_min
%     cah_max
    cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
    
%     dab_min
%     dab_max
    dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
    
%     sil_min
%     sil_max
    sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
    
%     % rcpso
%     disp('==> rcpso');
%     col=15;
%     
%     accs = zeros(1,steps2); 
%     sils = zeros(1,steps2); 
%     cahs = zeros(1,steps2); 
%     dabs = zeros(1,steps2); 
%     solutions = zeros(nd, steps2); 
%     
%     dim = nc*na;
%     for i=1:steps2
%         [out, rate, centroids] = chaoticpso('hclusteringni',dim,0,'pso','yes',0.007,0,med_inputs);
%         
%         [accs(1,i) cidx] = Accuracy(centroids',med_inputs,output);
%         solutions(:,i) = cidx;
%         eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
%         cahs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
%         dabs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'silhouette');
%         sils(1,i) = eva.CriterionValues;
%     end
% 
%     cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
%     [acc_max, amax_id] = max(accs);
%     
%     [cah_min, cmin_id] = min(cahs);
%     [cah_max, cmax_id] = max(cahs);
%     
%     [dab_min, dmin_id] = min(dabs);
%     [dab_max, dmax_id] = max(dabs);
%     
%     [sil_min, smin_id] = min(sils);
%     [sil_max, smax_id] = max(sils);
%     
% %     worst_solution = solutions(:,amin_id)
% %     best_solution = solutions(:,amax_id)
%     
% %     acc_min
% %     acc_max
%     acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
%     
% %     cah_min
% %     cah_max
%     cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
%     
% %     dab_min
% %     dab_max
%     dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
%     
% %     sil_min
% %     sil_max
%     sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
    
    % rcpso intervalar
    disp('==> rcpso intervalar');
    col=16;
    
    accs = zeros(1,steps2); 
    sils = zeros(1,steps2); 
    cahs = zeros(1,steps2); 
    dabs = zeros(1,steps2); 
    solutions = zeros(nd, steps2); 
    
    dim = 2*nc*na;
    for i=1:steps2
        [out, rate, centroids] = chaoticpso('hclustering',dim,0,'pso','yes',0.007,0,inputs);
        
        [accs(1,i) cidx] = Accuracy_interval(centroids',inputs,output);
        solutions(:,i) = cidx;
        eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
        cahs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
        dabs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'silhouette');
        sils(1,i) = eva.CriterionValues;
    end

    cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
    [acc_max, amax_id] = max(accs);
    
    [cah_min, cmin_id] = min(cahs);
    [cah_max, cmax_id] = max(cahs);
    
    [dab_min, dmin_id] = min(dabs);
    [dab_max, dmax_id] = max(dabs);
    
    [sil_min, smin_id] = min(sils);
    [sil_max, smax_id] = max(sils);
    
%     worst_solution = solutions(:,amin_id)
%     best_solution = solutions(:,amax_id)
    
%     acc_min
%     acc_max
    acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
    
%     cah_min
%     cah_max
    cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
    
%     dab_min
%     dab_max
    dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
    
%     sil_min
%     sil_max
    sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
    
    %%% KHMS %%%
    
    pe = 2;
    
%     % k-harmonic means    
%     disp(strcat('==> k-harmonic means (p=',num2str(pe),')'));
%     col=3;
%     
%     accs = zeros(1,steps); 
%     sils = zeros(1,steps); 
%     cahs = zeros(1,steps); 
%     dabs = zeros(1,steps); 
%     solutions = zeros(nd, steps); 
%     
%     inferior = min(med_inputs(:));
%     superior = max(med_inputs(:));
%     dim = nc*na;
%     VRmin=ones(dim,1)*superior/2;
%     VRmax=ones(dim,1)*superior;
%     C = zeros(dim,1);
%     for i=1:steps
%         for d=1:dim
%             C(d,1) = unifrnd(VRmin(d), VRmax(d));
%         end
%         [centroids, value] = KHM( C, med_inputs, pe);
%         
%         [accs(1,i) cidx] = Accuracy(centroids,med_inputs,output);
%         solutions(:,i) = cidx;
%         eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
%         cahs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
%         dabs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'silhouette');
%         sils(1,i) = eva.CriterionValues;
%     end
% 
%     cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
%     [acc_max, amax_id] = max(accs);
%     
%     [cah_min, cmin_id] = min(cahs);
%     [cah_max, cmax_id] = max(cahs);
%     
%     [dab_min, dmin_id] = min(dabs);
%     [dab_max, dmax_id] = max(dabs);
%     
%     [sil_min, smin_id] = min(sils);
%     [sil_max, smax_id] = max(sils);
%     
% %     worst_solution = solutions(:,amin_id)
% %     best_solution = solutions(:,amax_id)
%     
% %     acc_min
% %     acc_max
%     acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
%     
% %     cah_min
% %     cah_max
%     cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
%     
% %     dab_min
% %     dab_max
%     dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
%     
% %     sil_min
% %     sil_max
%     sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
    
    % k-harmonic means intervalar
    disp(strcat('==> k-harmonic means intervalar (p=',num2str(pe),')'));
    col=6;
    
    accs = zeros(1,steps); 
    sils = zeros(1,steps); 
    cahs = zeros(1,steps); 
    dabs = zeros(1,steps); 
    solutions = zeros(nd, steps); 
    
    inferior = min(inputs(:));
    superior = max(inputs(:));
    dim = 2*nc*na;
    VRmin=ones(dim,1)*superior/2;
    VRmax=ones(dim,1)*superior;
    C = zeros(dim,1);
    for i=1:steps
        for d=1:dim
            C(d,1) = unifrnd(VRmin(d), VRmax(d));
        end
        [centroids, value] = KHMinterval( C, inputs, pe);
        
        [accs(1,i) cidx] = Accuracy_interval(centroids,inputs,output);
        solutions(:,i) = cidx;
        eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
        cahs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
        dabs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'silhouette');
        sils(1,i) = eva.CriterionValues;
    end

    cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
    [acc_max, amax_id] = max(accs);
    
    [cah_min, cmin_id] = min(cahs);
    [cah_max, cmax_id] = max(cahs);
    
    [dab_min, dmin_id] = min(dabs);
    [dab_max, dmax_id] = max(dabs);
    
    [sil_min, smin_id] = min(sils);
    [sil_max, smax_id] = max(sils);
    
%     worst_solution = solutions(:,amin_id)
%     best_solution = solutions(:,amax_id)
    
%     acc_min
%     acc_max
    acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
    
%     cah_min
%     cah_max
    cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
    
%     dab_min
%     dab_max
    dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
    
%     sil_min
%     sil_max
    sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
    
%     % psokhm
%     disp(strcat('==> psokhm (p=',num2str(pe),')'));
%     col=9;
%     
%     accs = zeros(1,steps2); 
%     sils = zeros(1,steps2); 
%     cahs = zeros(1,steps2); 
%     dabs = zeros(1,steps2); 
%     solutions = zeros(nd, steps2); 
%     
%     dim = nc*na;
%     for i=1:steps2
%         [out, rate, centroids] = chaoticpsoKHM(dim,'pso','yes',0.007,0,med_inputs,pe,0);
% 
%         [accs(1,i) cidx] = Accuracy(centroids',med_inputs,output);
%         solutions(:,i) = cidx;
%         eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
%         cahs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
%         dabs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'silhouette');
%         sils(1,i) = eva.CriterionValues;
%     end
% 
%     cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
%     [acc_max, amax_id] = max(accs);
%     
%     [cah_min, cmin_id] = min(cahs);
%     [cah_max, cmax_id] = max(cahs);
%     
%     [dab_min, dmin_id] = min(dabs);
%     [dab_max, dmax_id] = max(dabs);
%     
%     [sil_min, smin_id] = min(sils);
%     [sil_max, smax_id] = max(sils);
%     
% %     worst_solution = solutions(:,amin_id)
% %     best_solution = solutions(:,amax_id)
%     
% %     acc_min
% %     acc_max
%     acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
%     
% %     cah_min
% %     cah_max
%     cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
%     
% %     dab_min
% %     dab_max
%     dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
%     
% %     sil_min
% %     sil_max
%     sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
    
    % psokhm intervalar
    disp(strcat('==> psokhm intervalar (p=',num2str(pe),')'));
    col=12;
    
    accs = zeros(1,steps2); 
    sils = zeros(1,steps2); 
    cahs = zeros(1,steps2); 
    dabs = zeros(1,steps2); 
    solutions = zeros(nd, steps2); 
    
    dim = 2*nc*na;
    for i=1:steps2
        [out, rate, centroids] = chaoticpsoKHM(dim,'pso','yes',0.007,0,inputs,pe,1);
        
        [accs(1,i) cidx] = Accuracy_interval(centroids',inputs,output);
        solutions(:,i) = cidx;
        eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
        cahs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
        dabs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'silhouette');
        sils(1,i) = eva.CriterionValues;
    end

    cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
    [acc_max, amax_id] = max(accs);
    
    [cah_min, cmin_id] = min(cahs);
    [cah_max, cmax_id] = max(cahs);
    
    [dab_min, dmin_id] = min(dabs);
    [dab_max, dmax_id] = max(dabs);
    
    [sil_min, smin_id] = min(sils);
    [sil_max, smax_id] = max(sils);
    
%     worst_solution = solutions(:,amin_id)
%     best_solution = solutions(:,amax_id)
    
%     acc_min
%     acc_max
    acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
    
%     cah_min
%     cah_max
    cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
    
%     dab_min
%     dab_max
    dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
    
%     sil_min
%     sil_max
    sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
    
    %%%%%
    pe = 2.5;
    %%%%%
    
%     % k-harmonic means    
%     disp(strcat('==> k-harmonic means (p=',num2str(pe),')'));
%     col=4;
%     
%     accs = zeros(1,steps); 
%     sils = zeros(1,steps); 
%     cahs = zeros(1,steps); 
%     dabs = zeros(1,steps); 
%     solutions = zeros(nd, steps); 
%     
%     inferior = min(med_inputs(:));
%     superior = max(med_inputs(:));
%     dim = nc*na;
%     VRmin=ones(dim,1)*superior/2;
%     VRmax=ones(dim,1)*superior;
%     C = zeros(dim,1);
%     for i=1:steps
%         for d=1:dim
%             C(d,1) = unifrnd(VRmin(d), VRmax(d));
%         end
%         [centroids, value] = KHM( C, med_inputs, pe);
%         
%         [accs(1,i) cidx] = Accuracy(centroids,med_inputs,output);
%         solutions(:,i) = cidx;
%         eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
%         cahs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
%         dabs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'silhouette');
%         sils(1,i) = eva.CriterionValues;
%     end
% 
%     cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
%     [acc_max, amax_id] = max(accs);
%     
%     [cah_min, cmin_id] = min(cahs);
%     [cah_max, cmax_id] = max(cahs);
%     
%     [dab_min, dmin_id] = min(dabs);
%     [dab_max, dmax_id] = max(dabs);
%     
%     [sil_min, smin_id] = min(sils);
%     [sil_max, smax_id] = max(sils);
%     
% %     worst_solution = solutions(:,amin_id)
% %     best_solution = solutions(:,amax_id)
%     
% %     acc_min
% %     acc_max
%     acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
%     
% %     cah_min
% %     cah_max
%     cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
%     
% %     dab_min
% %     dab_max
%     dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
%     
% %     sil_min
% %     sil_max
%     sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
    
    % k-harmonic means intervalar
    disp(strcat('==> k-harmonic means intervalar (p=',num2str(pe),')'));
    col=7;
    
    accs = zeros(1,steps); 
    sils = zeros(1,steps); 
    cahs = zeros(1,steps); 
    dabs = zeros(1,steps); 
    solutions = zeros(nd, steps); 
    
    inferior = min(inputs(:));
    superior = max(inputs(:));
    dim = 2*nc*na;
    VRmin=ones(dim,1)*superior/2;
    VRmax=ones(dim,1)*superior;
    C = zeros(dim,1);
    for i=1:steps
        for d=1:dim
            C(d,1) = unifrnd(VRmin(d), VRmax(d));
        end
        [centroids, value] = KHMinterval( C, inputs, pe);
        
        [accs(1,i) cidx] = Accuracy_interval(centroids,inputs,output);
        solutions(:,i) = cidx;
        eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
        cahs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
        dabs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'silhouette');
        sils(1,i) = eva.CriterionValues;
    end

    cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
    [acc_max, amax_id] = max(accs);
    
    [cah_min, cmin_id] = min(cahs);
    [cah_max, cmax_id] = max(cahs);
    
    [dab_min, dmin_id] = min(dabs);
    [dab_max, dmax_id] = max(dabs);
    
    [sil_min, smin_id] = min(sils);
    [sil_max, smax_id] = max(sils);
    
%     worst_solution = solutions(:,amin_id)
%     best_solution = solutions(:,amax_id)
    
%     acc_min
%     acc_max
    acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
    
%     cah_min
%     cah_max
    cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
    
%     dab_min
%     dab_max
    dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
    
%     sil_min
%     sil_max
    sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
    
%     % psokhm
%     disp(strcat('==> psokhm (p=',num2str(pe),')'));
%     col=10;
%     
%     accs = zeros(1,steps2); 
%     sils = zeros(1,steps2); 
%     cahs = zeros(1,steps2); 
%     dabs = zeros(1,steps2); 
%     solutions = zeros(nd, steps2); 
%     
%     dim = nc*na;
%     for i=1:steps2
%         [out, rate, centroids] = chaoticpsoKHM(dim,'pso','yes',0.007,0,med_inputs,pe,0);
% 
%         [accs(1,i) cidx] = Accuracy(centroids',med_inputs,output);
%         solutions(:,i) = cidx;
%         eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
%         cahs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
%         dabs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'silhouette');
%         sils(1,i) = eva.CriterionValues;
%     end
% 
%     cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
%     [acc_max, amax_id] = max(accs);
%     
%     [cah_min, cmin_id] = min(cahs);
%     [cah_max, cmax_id] = max(cahs);
%     
%     [dab_min, dmin_id] = min(dabs);
%     [dab_max, dmax_id] = max(dabs);
%     
%     [sil_min, smin_id] = min(sils);
%     [sil_max, smax_id] = max(sils);
%     
% %     worst_solution = solutions(:,amin_id)
% %     best_solution = solutions(:,amax_id)
%     
% %     acc_min
% %     acc_max
%     acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
%     
% %     cah_min
% %     cah_max
%     cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
%     
% %     dab_min
% %     dab_max
%     dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
%     
% %     sil_min
% %     sil_max
%     sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
    
    % psokhm intervalar
    disp(strcat('==> psokhm intervalar (p=',num2str(pe),')'));
    col=13;
    
    accs = zeros(1,steps2); 
    sils = zeros(1,steps2); 
    cahs = zeros(1,steps2); 
    dabs = zeros(1,steps2); 
    solutions = zeros(nd, steps2); 
    
    dim = 2*nc*na;
    for i=1:steps2
        [out, rate, centroids] = chaoticpsoKHM(dim,'pso','yes',0.007,0,inputs,pe,1);
        
        [accs(1,i) cidx] = Accuracy_interval(centroids',inputs,output);
        solutions(:,i) = cidx;
        eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
        cahs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
        dabs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'silhouette');
        sils(1,i) = eva.CriterionValues;
    end

    cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
    [acc_max, amax_id] = max(accs);
    
    [cah_min, cmin_id] = min(cahs);
    [cah_max, cmax_id] = max(cahs);
    
    [dab_min, dmin_id] = min(dabs);
    [dab_max, dmax_id] = max(dabs);
    
    [sil_min, smin_id] = min(sils);
    [sil_max, smax_id] = max(sils);
    
%     worst_solution = solutions(:,amin_id)
%     best_solution = solutions(:,amax_id)
    
%     acc_min
%     acc_max
    acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
    
%     cah_min
%     cah_max
    cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
    
%     dab_min
%     dab_max
    dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
    
%     sil_min
%     sil_max
    sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
    
    %%%%%
    pe = 3;
    %%%%%
    
%     % k-harmonic means    
%     disp(strcat('==> k-harmonic means (p=',num2str(pe),')'));
%     col=5;
%     
%     accs = zeros(1,steps); 
%     sils = zeros(1,steps); 
%     cahs = zeros(1,steps); 
%     dabs = zeros(1,steps); 
%     solutions = zeros(nd, steps); 
%     
%     inferior = min(med_inputs(:));
%     superior = max(med_inputs(:));
%     dim = nc*na;
%     VRmin=ones(dim,1)*superior/2;
%     VRmax=ones(dim,1)*superior;
%     C = zeros(dim,1);
%     for i=1:steps
%         for d=1:dim
%             C(d,1) = unifrnd(VRmin(d), VRmax(d));
%         end
%         [centroids, value] = KHM( C, med_inputs, pe);
%         
%         [accs(1,i) cidx] = Accuracy(centroids,med_inputs,output);
%         solutions(:,i) = cidx;
%         eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
%         cahs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
%         dabs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'silhouette');
%         sils(1,i) = eva.CriterionValues;
%     end
% 
%     cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
%     [acc_max, amax_id] = max(accs);
%     
%     [cah_min, cmin_id] = min(cahs);
%     [cah_max, cmax_id] = max(cahs);
%     
%     [dab_min, dmin_id] = min(dabs);
%     [dab_max, dmax_id] = max(dabs);
%     
%     [sil_min, smin_id] = min(sils);
%     [sil_max, smax_id] = max(sils);
%     
% %     worst_solution = solutions(:,amin_id)
% %     best_solution = solutions(:,amax_id)
%     
% %     acc_min
% %     acc_max
%     acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
%     
% %     cah_min
% %     cah_max
%     cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
%     
% %     dab_min
% %     dab_max
%     dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
%     
% %     sil_min
% %     sil_max
%     sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
    
    % k-harmonic means intervalar
    disp(strcat('==> k-harmonic means intervalar (p=',num2str(pe),')'));
    col=8;
    
    accs = zeros(1,steps); 
    sils = zeros(1,steps); 
    cahs = zeros(1,steps); 
    dabs = zeros(1,steps); 
    solutions = zeros(nd, steps); 
    
    inferior = min(inputs(:));
    superior = max(inputs(:));
    dim = 2*nc*na;
    VRmin=ones(dim,1)*superior/2;
    VRmax=ones(dim,1)*superior;
    C = zeros(dim,1);
    for i=1:steps
        for d=1:dim
            C(d,1) = unifrnd(VRmin(d), VRmax(d));
        end
        [centroids, value] = KHMinterval( C, inputs, pe);
        
        [accs(1,i) cidx] = Accuracy_interval(centroids,inputs,output);
        solutions(:,i) = cidx;
        eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
        cahs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
        dabs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'silhouette');
        sils(1,i) = eva.CriterionValues;
    end

    cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
    [acc_max, amax_id] = max(accs);
    
    [cah_min, cmin_id] = min(cahs);
    [cah_max, cmax_id] = max(cahs);
    
    [dab_min, dmin_id] = min(dabs);
    [dab_max, dmax_id] = max(dabs);
    
    [sil_min, smin_id] = min(sils);
    [sil_max, smax_id] = max(sils);
    
%     worst_solution = solutions(:,amin_id)
%     best_solution = solutions(:,amax_id)
    
%     acc_min
%     acc_max
    acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
    
%     cah_min
%     cah_max
    cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
    
%     dab_min
%     dab_max
    dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
    
%     sil_min
%     sil_max
    sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
    
%     % psokhm
%     disp(strcat('==> psokhm (p=',num2str(pe),')'));
%     col=11;
%     
%     accs = zeros(1,steps2); 
%     sils = zeros(1,steps2); 
%     cahs = zeros(1,steps2); 
%     dabs = zeros(1,steps2); 
%     solutions = zeros(nd, steps2); 
%     
%     dim = nc*na;
%     for i=1:steps2
%         [out, rate, centroids] = chaoticpsoKHM(dim,'pso','yes',0.007,0,med_inputs,pe,0);
% 
%         [accs(1,i) cidx] = Accuracy(centroids',med_inputs,output);
%         solutions(:,i) = cidx;
%         eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
%         cahs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
%         dabs(1,i) = eva.CriterionValues;
%         eva = evalclusters(med_inputs,cidx,'silhouette');
%         sils(1,i) = eva.CriterionValues;
%     end
% 
%     cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
%     [acc_max, amax_id] = max(accs);
%     
%     [cah_min, cmin_id] = min(cahs);
%     [cah_max, cmax_id] = max(cahs);
%     
%     [dab_min, dmin_id] = min(dabs);
%     [dab_max, dmax_id] = max(dabs);
%     
%     [sil_min, smin_id] = min(sils);
%     [sil_max, smax_id] = max(sils);
%     
% %     worst_solution = solutions(:,amin_id)
% %     best_solution = solutions(:,amax_id)
%     
% %     acc_min
% %     acc_max
%     acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
%     
% %     cah_min
% %     cah_max
%     cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
%     
% %     dab_min
% %     dab_max
%     dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
%     
% %     sil_min
% %     sil_max
%     sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
    
    % psokhm intervalar
    disp(strcat('==> psokhm intervalar (p=',num2str(pe),')'));
    col=14;
    
    accs = zeros(1,steps2); 
    sils = zeros(1,steps2); 
    cahs = zeros(1,steps2); 
    dabs = zeros(1,steps2); 
    solutions = zeros(nd, steps2); 
    
    dim = 2*nc*na;
    for i=1:steps2
        [out, rate, centroids] = chaoticpsoKHM(dim,'pso','yes',0.007,0,inputs,pe,1);
        
        [accs(1,i) cidx] = Accuracy_interval(centroids',inputs,output);
        solutions(:,i) = cidx;
        eva = evalclusters(med_inputs,cidx,'CalinskiHarabasz');
        cahs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'DaviesBouldin');
        dabs(1,i) = eva.CriterionValues;
        eva = evalclusters(med_inputs,cidx,'silhouette');
        sils(1,i) = eva.CriterionValues;
    end

    cahs(isnan(cahs)) = []; dabs(isnan(dabs)) = []; sils(isnan(sils)) = []; [acc_min, amin_id] = min(accs);
    [acc_max, amax_id] = max(accs);
    
    [cah_min, cmin_id] = min(cahs);
    [cah_max, cmax_id] = max(cahs);
    
    [dab_min, dmin_id] = min(dabs);
    [dab_max, dmax_id] = max(dabs);
    
    [sil_min, smin_id] = min(sils);
    [sil_max, smax_id] = max(sils);
    
%     worst_solution = solutions(:,amin_id)
%     best_solution = solutions(:,amax_id)
    
%     acc_min
%     acc_max
    acc_avg = 100*mean(accs); results(1,col) = acc_avg; acc_avg
    
%     cah_min
%     cah_max
    cah_avg = mean(cahs); results(2,col) = cah_avg; cah_avg
    
%     dab_min
%     dab_max
    dab_avg = mean(dabs); results(3,col) = dab_avg; dab_avg
    
%     sil_min
%     sil_max
    sil_avg = mean(sils); results(4,col) = sil_avg; sil_avg
    
end

