clc; clear variables; close all;

%number of users
NoU = 15;

%data size
N = 10^3;

%emergency data generation time 
%based on poisson distribution 
emergGenerationTime = rand(1,N);

%waiting times of each device
waitingTimes = rand(1,N);

%number of codes related to each square tile in network area
numberofloccodes = 100;

%collision occur if the codes and resources are similar 
%probability of collision
%devicesSameCode = randi([1 NoU]);%variable
devicesSameCode =2;

a = prod(randperm(devicesSameCode))
b = nchoosek(NoU,devicesSameCode)%each device selecting a code 
c= NoU^devicesSameCode

%probability of collision due to selecting the same code based on the number of codes
%available 
probCollision = 1 -  b/c*a
