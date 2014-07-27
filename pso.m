%% Particle Swarm Optimization Simulation
% Animiation of birds movement of a swarm to get the global minimum solution 
%
% Author: Wael Mansour (wael192@yahoo.com)
%
% MSc Student, Electrical Enginering Dept, 
% Faculty of Engineering Cairo University, Egypt

%% Initialization
clear
clc

%Initialization of the input file
% fileName = 'iris.data';
% fid = fopen(fileName);
% pcursor = 0;

cities_temp = [
[[-4,4], [-5,3], [2,12], [5,15], [7,17], [10,20], [10,20], [12,23], [10,20], [5,15], [1,10], [-1,4]];
[[6,12], [6,12], [8,16], [11,19], [16,25], [19,29], [22,32], [22,32], [19,28], [16,23], [11,18], [8,14]];
[[13,19], [14,19], [17,23], [21,27], [25,32], [28,34], [29,36], [30,36], [28,34], [24,31], [20,26], [15,21]];
[[19,28], [19,28], [22,30], [24,32], [27,33], [26,32], [25,30], [25,30], [24,30], [24,32], [23,32], [20,30]];
[[8,20], [9,22], [11,25], [14,29], [17,33], [20,35], [22,36], [22,35], [20,33], [18,31], [14,26], [10,20]];
[[13,27], [16,29], [21,34], [24,36], [26,36], [26,33], [26,32], [26,32], [26,32], [24,32], [23,32], [20,30]];
[[22,30], [22,30], [23,31], [24,31], [25,31], [25,30], [25,29], [25,29], [25,30], [24,29], [23,29], [22,30]];
[[-2,2], [-3,2], [-1,5], [3,10], [8,16], [11,20], [14,22], [14,21], [11,18], [7,12], [3,7], [1,4]];
[[13,23], [14,24], [17,28], [19,31], [22,34], [25,36], [28,39], [28,39], [25,37], [21,34], [17,30], [14,26]];
[[-10,9], [-8,10], [-4,17], [0,24], [3,27], [7,30], [8,32], [8,31], [5,27], [0,22], [-3,14], [-8,10]];
[[-3,5], [-6,6], [3,9], [7,13], [10,17], [15,17], [16,24], [16,23], [11,19], [6,13], [3,8], [-2,6]];
[[13,17], [12,16], [15,19], [19,23], [22,27], [25,29], [25,30], [25,30], [25,29], [22,27], [18,23], [14,19]];
[[22,31], [23,32], [23,33], [23,33], [23,32], [23,32], [23,31], [23,32], [23,32], [23,31], [23,31], [23,31]];
[[8,13], [8,14], [9,16], [11,18], [13,21], [16,24], [17,26], [18,27], [17,24], [14,21], [11,17], [8,14]];
[[2,6], [2,7], [3,10], [5,13], [8,17], [11,20], [13,22], [13,21], [11,19], [8,14], [5,10], [3,7]];
[[20,30], [20,31], [22,33], [26,35], [28,39], [27,38], [26,36], [26,35], [25,34], [24,32], [22,30], [21,29]];
[[1,9], [1,12], [3,16], [6,19], [9,24], [13,29], [16,34], [16,33], [13,28], [8,20], [4,14], [1,9]];
[[21,27], [22,27], [24,29], [24,31], [25,31], [25,31], [23,29], [24,28], [25,28], [24,29], [22,28], [22,27]];
[[22,28], [22,29], [22,29], [21,28], [19,25], [18,24], [17,23], [17,23], [17,24], [18,25], [19,27], [21,28]];
[[6,22], [15,23], [17,25], [18,27], [18,27], [18,27], [18,27], [18,26], [18,26], [16,25], [19,27], [21,28]];
[[-13,-6], [-12,-5], [-8,0], [0,8], [7,18], [11,23], [13,24], [11,22], [6,16], [1,8], [-5,0], [-11,-5]];
[[-6,1], [-5,3], [-2,9], [3,14], [7,18], [10,21], [12,23], [11,23], [8,20], [4,13], [0,7], [-4,3]];
[[12,25], [13,26], [14,25], [14,24], [13,22], [12,21], [11,21], [11,21], [11,24], [13,24], [13,23], [13,23]];
[[-2,4], [-3,4], [1,9], [6,15], [12,22], [17,27], [21,29], [20,28], [16,24], [11,19], [5,12], [-2,6]];
[[-2,4], [-3,4], [1,9], [6,15], [12,22], [17,27], [21,29], [20,28], [16,24], [11,19], [5,12], [-2,6]];
[[1,7], [1,7], [2,12], [5,16], [8,19], [12,22], [14,24], [13,24], [11,21], [7,16], [4,10], [1,6]];
[[4,11], [5,13], [7,16], [10,19], [13,23], [17,18], [20,31], [20,31], [17,27], [13,21], [9,16], [5,12]];
[[6,13], [6,14], [7,17], [8,18], [10,19], [11,21], [12,22], [12,22], [12,23], [11,22], [8,18], [6,14]];
[[0,7], [1,6], [1,8], [6,16], [12,22], [16,25], [18,31], [16,30], [9,28], [3,24], [7,19], [1,8]];
[[23,30] [23,30], [24,31], [24,31], [24,30], [25,30], [25,30], [25,30], [24,30], [24,30], [24,30], [23,30]];
[[-9,-5], [-9,-6], [-4,2], [1,8], [6,15], [11,19], [14,22], [13,20], [9,15], [5,9], [1,4], [-2,2]];
[[20,30], [20,30], [18,26], [16,23], [12,20], [5,17], [8,16], [9,17], [11,20], [13,22], [16,26], [20,30]];
[[0,5], [5,8], [10,15], [15,18], [20,25], [28,30], [36,38], [38,40], [29,30], [18,20], [9,12], [-5,0]];
[[0,9], [0,10], [3,13], [9,18], [14,23], [18,25], [22,29], [13,31], [20,27], [13,21], [8,16], [2,21]];
[[-8,-1], [-8,-1], [-4,4], [-2,11], [-8,18], [13,24], [16,27], [16,26], [12,22], [6,14], [-1,17], [-5,1]];
[[-2,1], [-1,3], [1,8], [5,14], [10,19], [13,22], [15,24], [14,23], [11,19], [7,13], [2,7], [1,3]];
[[-11,9], [-8,15], [-7,18], [-1,21], [2,27], [6,30], [10,31], [8,25], [5,23], [3,22], [0,19], [-11,8]]
];

na = length(cities_temp(1,:))/2;
%four kinds of temperature
nc = 4;
%Number of data points in the dataset
no_data = length(cities_temp);
nd = no_data;
 
%Flower categories saved as category numbers in output array
% str_out1 = 'Iris-setosa';
% str_out2 = 'Iris-versicolor';
% str_out3 = 'Iris-virginica';

% %The first 4 columns of data
% inputs = zeros(no_data,na);
% %Flower categories 1, 2 or 3 is saved in this variable.
% outputs = zeros(no_data,1);

for i=1:nd
    for j=1:na
       new = (cities_temp(i,(j*2-1))+cities_temp(i,j*2))/2;
       med_cities_temp(i,j)=new;
    end
end

inputs = cities_temp;
output = zeros(1,nd);
output(1)=1;output(2)=1;output(8)=1;output(10)=1;output(11)=1;output(14)=1;output(15)=1;output(17)=1;output(21)=1;output(22)=1;output(25)=1;output(26)=1;output(27)=1;output(28)=1;output(29)=1;output(31)=1;output(34)=1;output(35)=1;output(36)=1;output(37)=1;
output(3)=2;output(4)=2;output(5)=2;output(6)=2;output(7)=2;output(9)=2;output(12)=2;output(13)=2;output(16)=2;output(18)=2;output(19)=2;output(24)=2;output(30)=2;
output(20)=3;output(23)=3;output(32)=3;
output(33)=4;

min_v = min(inputs(:));
max_v = max(inputs(:));
gama = Estimate_gama(inputs);
%gama = 1;

% while(1)
%     pcursor = pcursor + 1;
%     %Read next line in the file
%     tline = fgetl(fid);
%     %Check for EOF
%     if (length(tline) < 2)
%         break;
%     end
%     %Find location of commas
%     commaLocs=findstr(',',tline);
%     %Extract data from the data files
%     if (strcmp(fileName,'iris.data'))
%         %Input assignment
%         inputs(pcursor,1) = str2double(tline(1:(commaLocs(1)-1)));
%         inputs(pcursor,2) = str2double(tline((commaLocs(1)+1):(commaLocs(2)-1)));
%         inputs(pcursor,3) = str2double(tline((commaLocs(2)+1):(commaLocs(3)-1)));
%         inputs(pcursor,4) = str2double(tline((commaLocs(3)+1):(commaLocs(4)-1)));
%         %Output assignment
%         str_out = tline((commaLocs(4)+1):length(tline));
%         switch str_out
%             case str_out1
%                 outputs(pcursor,:) = 1;
%             case str_out2
%                 outputs(pcursor,:) = 2;
%             case str_out3
%                 outputs(pcursor,:) = 3;
%         end
%     end
%     if (strcmp(fileName,'wine.data'))
%         %Output assignment
%         outputs(pcursor,:) = str2double(tline(1:(commaLocs(1)-1)));
%         %Input assignment
%         for i=1:(na-1)
%             inputs(pcursor,i) = str2double(tline((commaLocs(i)+1):(commaLocs(i+1)-1)));
%         end
%         inputs(pcursor,na) = str2double(tline((commaLocs(na)+1):length(tline)));
%     end
% end
% 
% fid = fclose(fid);

% Parameters
n = 1000;          % Size of the swarm " no of birds "
bird_setp  = 1000; % Maximum number of "birds steps"
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
    axis([-50 50 -50 50]);
    
pause(.2)
% end

end % end of while loop its mean the end of all step that the birds move it 
                      

              [Jbest_min,I] = min(current_fitness) % minimum fitness
              centroids = current_position(:,I) % best solution               

              Accuracy_interval(centroids,inputs,output)
    

%