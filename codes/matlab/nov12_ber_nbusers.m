clc;
clear all;
close all;
%% ber vs total number of users
% K = 5;
% K = 10;
%% for  K = 15 in total 
% intermediate K = 5; K = 4; K = 4; K = 2;

% sumsym_dur_vec = [0.1,0.1,0.1;0.2,0.3,0.4;0.5,0.4,0.3;]*0.01;
% power_vec = [0.1;0.01;0.001;0.0001;0.00001];
k_vec = linspace(1,15,15);
u = 1;%initiate u
nbiter = 100;
iter_K_vec = [4,2,2,2;];
rratio_vec = linspace(1,5,20);