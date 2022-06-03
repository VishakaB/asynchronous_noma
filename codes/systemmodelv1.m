clc; clear variables; close all;

%number of users
NoU = 15000;

%data size
N = 10^3;

%emergency data generation time 
%based on poisson distribution 
emergGenerationTime = rand(1,N);

%waiting times of each device
waitingTimes = rand(1,N);

%location of devices

%number of codes related to each square tile in network area
numberofloccodes = 100;

%probability of devices having the same resources in terms of time and code
%assuming carrier is same


%collision occur if the codes and resources are similar 
%number of collisions 
%collisionCount = 
