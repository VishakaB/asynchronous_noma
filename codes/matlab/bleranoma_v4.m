clc;
clear all;
close all;

%% Block Error Rate Simulation

%% Simulation Configuration
numTrBlks = 5;              % Number of simulated transport blocks 
SNRdB = [-20 -18 -15 -12.5 -10 -6.4 -3.5 0.7];% Range of SNR values in dB
simReps = [2 16 64];        % Repetitions to simulate
N  = 10;% nb bits transmitted
d1 = 10; d2 = 500; % Distances of users from base station (BS)
d3 = 5;
eta = 4;           % Path loss exponent

%% Propagation Channel Model Configuration

for repIdx = 1:numel(simReps)
   for snrIdx = 1:numel(SNRdB)
   
   end
end