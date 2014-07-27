clear
clc
close all

a1 = normrnd(25,4,50,1);
b1 = normrnd(22,3,50,1);
input1 = horzcat(a1,b1);
output1 = ones(1,50) * 1;

scatter(a1,b1);
hold on;

a2 = normrnd(60,3,50,1);
b2 = normrnd(20,5,50,1);
input2 = horzcat(a2,b2);
output2 = ones(1,50) * 2;

scatter(a2,b2);

a3 = normrnd(30,4,50,1);
b3 = normrnd(60,5,50,1);
input3 = horzcat(a3,b3);
output3 = ones(1,50) * 3;

scatter(a3,b3);

a4 = normrnd(60,5,50,1);
b4 = normrnd(50,5,50,1);
input4 = horzcat(a4,b4);
output4 = ones(1,50) * 4;

scatter(a4,b4);
hold off;

pause;

inputs = vertcat(input1,input2,input3,input4);
output = horzcat(output1,output2,output3,output4);

na = length(inputs(1,:))/2;
nc = 4;
no_data = length(inputs);
nd = no_data;

min_v = min(inputs(:));
max_v = max(inputs(:));
gama = Estimate_gama(inputs);

% Parameters
n = 1000;          % Size of the swarm " no of birds "
bird_setp  = 200; % Maximum number of "birds steps"
dim = 2*nc*na;          % Dimension of the problem

c2 =1.2;          % PSO parameter C1 
c1 = 0.12;        % PSO parameter C2 
w =0.9;           % pso momentum or inertia  
fitness=0*ones(n,bird_setp);

                                       %-----------------------------%
                                       %    initialize the parameter %
                                       %-----------------------------%
                                       
R1 = rand(dim, n);
R2 = rand(dim, n);
current_fitness =0*ones(n,1);

                                 %------------------------------------------------%
                                 % Initializing swarm and velocities and position %
                                 %------------------------------------------------%

current_position = zeros(dim, n);
for i=1:n %each particle
    possibles = inputs;    
    particle = zeros(dim,1);   
    for j=1:nc
        possibles_size = size(possibles,1);
        random_dataset = randi([1 possibles_size],1,1);        
        centroid = possibles(random_dataset,:);
        possibles(random_dataset,:) = [];
        particle(2*na*(j-1)+1:2*na*j)=centroid;
    end
    current_position(:,i) = particle;
end
                                 
% current_position = (max_v-min_v)*(rand(dim, n))+min_v;
velocity = .3*randn(dim, n) ;
local_best_position  = current_position ;


                                 %-------------------------------------------%
                                 %     Evaluate initial population           %           
                                 %-------------------------------------------%

for i = 1:n
    current_fitness(i) = Live_fn_interval(current_position(:,i),inputs,gama);    
end


local_best_fitness  = current_fitness ;
[global_best_fitness,g] = min(local_best_fitness) ;

for i=1:n
    globl_best_position(:,i) = local_best_position(:,g) ;
end
                                               %-------------------%
                                               %  VELOCITY UPDATE  %
                                               %-------------------%

velocity = w *velocity + c1*(R1.*(local_best_position-current_position)) + c2*(R2.*(globl_best_position-current_position));

                                               %------------------%
                                               %   SWARMUPDATE    %
                                               %------------------%
                                               
            
current_position = current_position + velocity ;

                                               %------------------------%
                                               %  evaluate anew swarm   %
                                               %------------------------%
                                               

figure;
%% Main Loop
iter = 0 ;        % Iterationscounter
while  ( iter < bird_setp )
iter = iter + 1;

for i = 1:n,
current_fitness(i) = Live_fn_interval(current_position(:,i), inputs, gama) ;    

end


for i = 1 : n
        if current_fitness(i) < local_best_fitness(i)
           local_best_fitness(i)  = current_fitness(i);  
           local_best_position(:,i) = current_position(:,i)   ;
        end   
 end

  
 [current_global_best_fitness,g] = min(local_best_fitness);
  
    
if current_global_best_fitness < global_best_fitness
   global_best_fitness = current_global_best_fitness;
   
    for i=1:n
        globl_best_position(:,i) = local_best_position(:,g);
    end
   
end


 velocity = w *velocity + c1*(R1.*(local_best_position-current_position)) + c2*(R2.*(globl_best_position-current_position));
 current_position = current_position + velocity; 
  
 

 
x=current_position(1,:);
y=current_position(2,:);



% if (iter>0.9*bird_setp)
% [Jbest_min,I] = min(current_fitness); % minimum fitness
% centroids = current_position(:,I); % best solution
% 
% Live_fn_interval_test(centroids, inputs, gama);
% 
% else 
clf    
    plot(x, y , 'h')   
    axis([0 80 0 80]);
    
pause(.05)
% end

end % end of while loop its mean the end of all step that the birds move it 
                      

              [Jbest_min,I] = min(current_fitness) % minimum fitness
              centroids = current_position(:,I) % best solution               

              Accuracy_interval(centroids,inputs,output)
    

%


