clear
clc
close all

% Data 1

x = normrnd(25,4,50,1);
y = normrnd(22,3,50,1);
a1 = x - 5*rand(1,1);
b1 = x + 5*rand(1,1);
a2 = y - 5*rand(1,1);
b2 = y + 5*rand(1,1);
input1 = horzcat(a1,b1,a2,b2);
output1 = ones(1,50) * 1;

scatter(x,y);
hold on;

x = normrnd(60,3,50,1);
y = normrnd(20,5,50,1);
a1 = x - 5*rand(1,1);
b1 = x + 5*rand(1,1);
a2 = y - 5*rand(1,1);
b2 = y + 5*rand(1,1);
input2 = horzcat(a1,b1,a2,b2);
output2 = ones(1,50) * 2;

scatter(x,y);

x = normrnd(30,4,50,1);
y = normrnd(60,5,50,1);
a1 = x - 5*rand(1,1);
b1 = x + 5*rand(1,1);
a2 = y - 5*rand(1,1);
b2 = y + 5*rand(1,1);
input3 = horzcat(a1,b1,a2,b2);
output3 = ones(1,50) * 3;

scatter(x,y);

x = normrnd(60,5,50,1);
y = normrnd(50,5,50,1);
a1 = x - 5*rand(1,1);
b1 = x + 5*rand(1,1);
a2 = y - 5*rand(1,1);
b2 = y + 5*rand(1,1);
input4 = horzcat(a1,b1,a2,b2);
output4 = ones(1,50) * 4;

scatter(x,y);
hold off;

pause;

inputs = vertcat(input1,input2,input3,input4)
output = horzcat(output1,output2,output3,output4);

%[s,h]=silhouette(inputs,output,@silhouetteIntervalDistance)

[idx centroids] = KHMinterval( 4, inputs, 2 );
[acc cidx] = Accuracy_interval(centroids,inputs,output)
