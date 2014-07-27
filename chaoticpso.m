function [out, rate, posbest] = chaoticpso(functname,D,noise,flag_method,flag_jump,eta_percentage,display_process,inputs)
    if ~exist('inputs', 'var') || isempty(inputs)
        inputs = [];
    end
   
   rand('twister',sum(100*clock));

%   - functname (string): name of the benchmark

%   - D (int): dimension number

%   - noise (real): noise level to corrupt the function landscape

%   - flag_method (string): which method to use ?
%         'pso': canonical particle swarm optimization
%         'lbest': standard pso (local topology)
%         'fips': fully informed pso
%         'bbpso': bare bones pso

%   - flag_jump (string):
%         'yes': to use the jump strategy
%         'no': otherwise

%   - eta_percentage (real):
%         -1: to consider default eta
%         otherwise: contains the size of the jump as percentage of the search space.

%   - display_process (int):
%         >0: extremely verbose mode


   
   % - eta_p: default eta value given as percentage of the search space.
   % - inferior and superior: limits of the search space (the same for all dimension).
   if strcmp(functname,'rosenbrock')
      inferior = -50;
      superior =  50;
      eta_p = 1;
   elseif strcmp(functname,'rastrigin')
      inferior = -5.12;
      superior =  5.12;
      eta_p = 10;
   elseif strcmp(functname,'ackley')
      inferior = -32;
      superior =  32;
      eta_p = 1.6;%1;
   elseif strcmp(functname,'griewank')
      inferior = -600;
      superior =  600;
      eta_p = 0.1;
   elseif strcmp(functname,'schaffer')
      D = 2;
      inferior = -100;
      superior =  100;
      eta_p = 1;
   elseif strcmp(functname,'sphere')
      inferior = -100;
      superior =  100;
      eta_p = 0.5;
   elseif strcmp(functname,'schwefel')
      inferior = -500;
      superior =  500;
      eta_p = 3.5;%1.8;
   elseif strcmp(functname,'hclustering') || strcmp(functname,'hclusteringni')
%       gama = Estimate_gama(inputs);
      gama = 1;
      inferior = min(inputs(:));
      superior = max(inputs(:));
      eta_p = 0.007;
   end


   %asimetrically initialization
   if strcmp(functname,'schwefel')
      VRmin=ones(D,1)*inferior;
      VRmax=ones(D,1)*inferior/2;
   else
      VRmin=ones(D,1)*superior/2;
      VRmax=ones(D,1)*superior;
   end


   % common pso parameters
   me = 1500;
   ps = 20;
   mv = 6;
   epsilon = 1e-4;


   % use default eta percentage (-1: yes)
   if eta_percentage == -1
      eta_percentage = eta_p;
   end

   % not default
   % eta = percentage of the search space
   if eta_percentage > 0
      eta = (eta_percentage/100.0) * (superior - inferior);
   end


   message = sprintf('PSO: %%g/%%g iterations, GBest = %%g\n');


   %------------------------------------------------------------------------------------------------------------------- 
   % initialize population of particles and their velocities at time zero
   pos = zeros(ps,D);
   vel = zeros(ps,D);
   
   
   if strcmp(functname,'hclustering') || strcmp(functname,'hclusteringni')
       na=size(inputs,2)/2;
       nc=D/(na*2);
       pos = zeros(ps, D);
       
%        [idx, C] = cvKmeans(inputs',nc,.05,@intervalDistance,false,1.0);
%        C = (reshape(C, 1, numel(C)));
%        pos(1,:) = C;
       
       for i=1:ps %each particle
           possibles = inputs;    
           particle = zeros(1,D);   
           for j=1:nc
               possibles_size = size(possibles,1);
               random_dataset = randi([1 possibles_size],1,1);        
               centroid = possibles(random_dataset,:);
               possibles(random_dataset,:) = [];
               particle(2*na*(j-1)+1:2*na*j)=centroid;
           end
           pos(i,:) = particle;
       end
       
       e = zeros(ps,1);
       for j=1:ps
          e(j) = fitness(functname, pos(j,1:D), noise, inputs, gama);
       end
       pbest = pos;
       pbestval = e;
       [gbestval,idx] = min(pbestval);
       
       size_seq_gbest = round(0.5*me);
       seq_gbest = zeros(1,size_seq_gbest);
       
       pos_seq = 1;
       seq_gbest(pos_seq) = gbestval;
       
       h=1;
       while h<me && length(unique(seq_gbest))~=1
           pos_seq=pos_seq+1;
           if (pos_seq > size_seq_gbest)
             pos_seq = 1;  
           end
           
%            [idx, C] = cvKmeans(inputs',nc,.05,@intervalDistance,false,1.0);
%            C = (reshape(C, 1, numel(C)));
%            pos(1,:) = C;
%            e(1) = fitness(functname, C, noise, inputs, gama);
%            
%            if (e(1) < pbestval(1))
%                pbestval(1) = e(1);
%                pbest(1,:) = pos(1,:);               
%            end
           
           for i=1:ps %each particle
               possibles = inputs;    
               particle = zeros(1,D);   
               for j=1:nc
                   possibles_size = size(possibles,1);
                   random_dataset = randi([1 possibles_size],1,1);        
                   centroid = possibles(random_dataset,:);
                   possibles(random_dataset,:) = [];
                   particle(2*na*(j-1)+1:2*na*j)=centroid;
               end
               pos(i,:) = particle;               
               e(i) = fitness(functname, particle, noise, inputs, gama);
           end
           
           [wval,idwv] = max(pbestval);
           for i=1:ps              
              if isempty(find(pbestval==e(i))) && e(i) < wval
                pbestval(idwv) = e(i);
                pbest(idwv,:) = pos(i,:);
                [wval,idwv] = max(pbestval);
              end
           end
           
           [gbestval,idx] = min(pbestval);
           seq_gbest(pos_seq) = gbestval;
%            seq_gbest(pos_seq) = sum(pbestval);
           
%            [gbv,idx] = min(pbestval);
%            if (gbv<gbestval)
%              gbestval = gbv;
%              pbest = pos;
%            end
%            seq_gbest(pos_seq) = gbestval;
            if display_process > 0
                pbestval'
                gbestval
            end
           h=h+1;
       end
       
       pos = pbest;
   else
       for d=1:D
          pos(1:ps,d) = unifrnd(VRmin(d)*ones(ps,1), VRmax(d)*ones(ps,1));
          %vel(1:ps,d) = unifrnd(-mv*ones(ps,1), mv*ones(ps,1)); % no need
       end
   end


   % added to test FIPS topology, take the best 3 neighbors
   if strcmp(flag_method,'fips') == 1
      load('fips');
   end


   fail_count = zeros(ps,1); % contador de falhas
   local_search_count = 5; % numero maximo de falhas antes de chavear para salto


   pulou=0; % flag para controle de salto: sim(1) ou nao(0).
   success_counter = zeros(me,2); % contador: jumps com sucesso por iteracao

   e = zeros(ps,1);
   for j=1:ps
      e(j) = fitness(functname, pos(j,1:D), noise, inputs, gama);
   end
   

   % initial pbest positions values
   pbest = pos;

   pbestval = e;
   [gbestval,idx] = min(pbestval);
   gbest = pbest(idx,:);  % this is gbest position


   % inicializar sequencia caotica com valor diferente de: 0,0.25,0.5,0.75,1.0
   cj=rand;
   while (cj==0) | (cj==0.25) | (cj==0.5) | (cj==0.75) | (cj==1.0)
      cj=rand;
   end

   size_seq_gbest = 20;
   seq_gbest = zeros(1,size_seq_gbest);
   pos_seq = 1;
   seq_gbest(pos_seq) = gbestval;
   
   if display_process > 0
    figure;
   end
   
   last_pbestval = pbestval;
   improvement = gbestval;
   
   min_improvement = 0.0001;
   
   i=1;
   while i<=me && (i<size_seq_gbest || improvement > min_improvement)%length(unique(seq_gbest))~=1 % start epoch loop (iterations)
       pos_seq=pos_seq+1;
       if (pos_seq > size_seq_gbest)
         pos_seq = 1;  
       end 
           
      for j=1:ps  % start particle loop

         % organizacao de estrutura em anel para metodo lbest
         if strcmp(flag_method,'lbest')
            n = j+1;
            p = j-1;
            if (j==1) p=ps; end
            if (j==ps) n=1; end

            m = pbestval(p); lbest=p;
            if (pbestval(j) < m)  m=pbestval(j); lbest=j; end
            if (pbestval(n) < m)  lbest=n; end
         end

         pulou=0;
         for dimcnt=1:D
            if strcmp(flag_jump, 'no') | (fail_count(j) <= local_search_count)

               if strcmp(flag_method,'pso')
                  vel(j,dimcnt) = 0.729*vel(j,dimcnt)...
                              + 1.49*rand*(pbest(j,dimcnt)-pos(j,dimcnt))...
                              + 1.49*rand*(gbest(1,dimcnt)-pos(j,dimcnt));
                  pos(j,dimcnt) = pos(j,dimcnt) + vel(j,dimcnt);

               elseif strcmp(flag_method,'lbest')
                  vel(j,dimcnt) = 0.729*(vel(j,dimcnt)...
                              + 2.05*rand*(pbest(j,dimcnt)-pos(j,dimcnt))...
                              + 2.05*rand*(pbest(lbest,dimcnt)-pos(j,dimcnt)));
                  pos(j,dimcnt) = pos(j,dimcnt) + vel(j,dimcnt);
                  
               elseif strcmp(flag_method,'fips')
                  xsum=0;
                  for degree=1:3
                     xsum=xsum + rand*4.1*(pbest(fips(j,degree),dimcnt)-pos(j,dimcnt));
                  end
                  vel(j,dimcnt) = 0.729*(vel(j,dimcnt) + xsum/3);
                  pos(j,dimcnt) = pos(j,dimcnt) + vel(j,dimcnt);
                  
               elseif strcmp(flag_method,'bbpso')
                  u = (pbest(j,dimcnt) + gbest(1,dimcnt))/2;
                  s = abs(pbest(j,dimcnt) - gbest(1,dimcnt));
                  pos(j,dimcnt) = u + s*randn;
               end
            else %jump
               pulou = 1;

               cj = 4*cj*(1-cj);
               pos(j,dimcnt) = pbest(j,dimcnt)*(1 + eta*(2*cj-1));
            end
         end  % end of the for dimcnt loop   

         % limit position components to maximums
         for dimcnt=1:D
            if pos(j,dimcnt)>superior
               pos(j,dimcnt)=pbest(j,dimcnt);
            end
            if pos(j,dimcnt)<inferior
               pos(j,dimcnt)=pbest(j,dimcnt);
            end
            if vel(j,dimcnt)>mv
               vel(j,dimcnt)=mv;
            end
            if vel(j,dimcnt)<-mv
               vel(j,dimcnt)=-mv;
            end    
         end


         % contador: total jumps (iteration i)
         if pulou==1
            success_counter(i,2) = success_counter(i,2) + 1;
         end

         % contador: reseta contador de falha apos um salto.
         if fail_count(j) > local_search_count
            fail_count(j) = 0;
         end
            

         % atualizar pbest se necessario
         e(j) = fitness(functname, pos(j,1:D), noise, inputs, gama);
         if pbestval(j) > e(j)
               pbestval(j) = e(j);
               pbest(j,:) = pos(j,:);

               % contador: successfull jump
               if pulou==1
                  success_counter(i,1) = success_counter(i,1) + 1;
               end
         else
               % contador: usado para decidir quando saltar
               fail_count(j) = fail_count(j) + 1;
         end
   
         % assign gbest by finding minimum of all particle pbests 
         [iterbestval,idx] = min(pbestval);
         if gbestval > iterbestval
               gbestval = iterbestval;
               gbest = pbest(idx,:);
         end

      end  % end particle loop
      

%       seq_gbest(pos_seq) = gbestval;
%       seq_gbest(pos_seq) = sum(pbestval);

      best_improvement = max(rdivide((last_pbestval - pbestval),last_pbestval));
      last_pbestval = pbestval;
      seq_gbest(pos_seq) = best_improvement;
      
      % mostrar mensangem de acompanhamento do processo
      if display_process > 0
         fprintf(message,i,me,gbestval);
         pbestval'
         if i > size_seq_gbest
            improvement = (sum(seq_gbest)/size_seq_gbest)
         end 
         x=pos(:,1);
         y=pos(:,2); 
         clf    
         plot(x, y , 'h')   
         axis([inferior superior inferior superior]);
         pause(.01)
      end
      
%          if strcmp(functname,'hclustering') || strcmp(functname,'hclusteringni')
%            if sum(seq_gbest)~=0 && length(unique(seq_gbest))==1
%                na=size(inputs,2)/2;
%                nc=D/(na*2);
%                for i2=1:ps %each particle
%                    possibles = inputs;    
%                    particle = zeros(1,D);   
%                    for j2=1:nc
%                        possibles_size = size(possibles,1);
%                        random_dataset = randi([1 possibles_size],1,1);        
%                        centroid = possibles(random_dataset,:);
%                        possibles(random_dataset,:) = [];
%                        particle(2*na*(j2-1)+1:2*na*j2)=centroid;
%                    end
%                    pos(i2,:) = particle;
%                end
%            end          
%          end
   
      % cumulative sucess sum
      success_counter(i+1,:) = success_counter(i,:);
      
      i=i+1;
   end  % end epoch loop
   

   % function output
   out = gbestval;
   posbest = gbest;

   % salvar o contador de jump apenas quando estiver em modo verbose.
   if display_process > 0
      save_sucess(success_counter, functname, flag_method);
   end

   % taxa de sucesso com jump.
   rate = success_counter(end,1)*100/success_counter(end,2);
end


function f = fitness(functname, pos, noise, inputs, gama)
   if strcmp(functname,'hclustering') || strcmp(functname,'hclusteringni')
      f = feval(functname, pos, inputs, gama);
   else
      f = feval(functname, pos);
   end
   f = f + noise*randn;
end


function save_sucess(success_counter, bench, method)

   sucess_rate = success_counter(:,1)*100 ./ success_counter(:,2);
   
   sucess_rate(isnan(sucess_rate)) = 100;

   filename = strcat(bench,'_');
   filename = strcat(filename,method);
   filename = strcat(filename,datestr(fix(clock),'_yyyy-mm-dd_HH:MM:SS'));
   filename = strcat(filename,'_success.pts');

   %save(filename,'sucess_rate','-ascii');
   
end

