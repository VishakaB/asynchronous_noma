%start date: Sep 3 2022
%last update: 07/09/2022

%Title: A-NOMA D2D URSLLC analysis
%goal1: average BLER metric analysis analytical
%goal2: comparison with simulation results

%goal2: latency metric analysis
%goal3: secrecy capacity metric analysis
%goal4: optimization of ursllc metrics

clc;
clear all;
close all;

%% input data: environmnet
N = 10^4;%total number of bits
m = 200;%block length
rc = N/m;%channel coding rate

%transmission power
%--------------------------------------------------------------------------
%SINR at the receiver 

ui = 1;%maximum sinr limit
vi = 2;%minimum sinr limit
ck = 3;%transmit snr

chi = (2*pi*(2^(2*ri)-1)/m)^(-0.5);

%% average bler calculation
avg_bler = chi*(ui - vi + ck^2*(exp(-ui/ck)/ui^2 ...
- exp(-vi/ck)/vi^2));

fprintf('avg bler %f\n',avg_bler)

%% average latency calculation
avg_bler = chi*(ui - vi + ck^2*(exp(-ui/ck)/ui^2 ...
- exp(-vi/ck)/vi^2));

%% average secrecy capacity calculation
avg_bler = chi*(ui - vi + ck^2*(exp(-ui/ck)/ui^2 ...
- exp(-vi/ck)/vi^2));

%% optimization
%minimize the average block error rate 
%under bler, latency, secrecy rate constraints
cvx_begin quiet
   variable m_k(K,1)
   dual variables var1 var2 var3 
   minimize(avg_bler)
   subject to
      var1 : -sum(decision_uk)+ sum(decision_uk.^2)<=0;

      var2: decision_uk'*sumsym_dur_vec-1/(priority+0.001)...
          *priority_max/timeslot <= 0
           
      var3: decision_uk.*((noisepower^2 + interf_vec +...
          power_vec(1:K).*mean(g_vec(1:K,:),2))...
                -power_vec(1:K).*mean(g_vec(1:K,:),2)*...
                (1+1/sinr_th)) <= 0                  
cvx_end