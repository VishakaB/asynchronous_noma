clc;
close all;
clear all;

%initialization
k1=5; %Rician factor %ref: https://www.researchgate.net/publication/263669548_Probability_Distribution_of_Rician_K-Factor_in_Urban_Suburban_and_Rural_Areas_Using_Real-World_Captured_Data
mean1=sqrt(k1/(k1+1));% mean
sigma=sqrt(1/(2*(k1+1)));%variance
N=3;  % Number of Bits for data_user1
d1 = 0.8; d2 = 500;    %Distances of users from base station (BS)
d3 = 5;d4 = 5;d5 = 5;d6 = 5;d7 = 5;d8 = 5;d9 = 5;dA=0.8;
eta = 4;            %Path loss exponent
n = 3;

% Number of Bits
 
random_iterations=3;
communication_radius = 30;%change this 
max_dist     = 100;%meters
max_eta      = 15;
etath        = 4;%change this 
noisepower   = 0.1;
max_tx_power = 20;%change this
B            = 1;%channel bandwidth

pth          = max_tx_power.*communication_radius^-etath;
h_th         = sqrt(communication_radius^-etath)*sqrt(pth/2)*(randn(1,N)+...
1i*randn(1,N))/sqrt(2);
g_th         = (abs(h_th)).^2;

rate_th      = log2( 1 + sqrt(pth/2)*g_th/noisepower);

eth          = 1;
timeslot     = 1;

K = 3;%number of superimposed data
%%vectors
%Distances of users from rx
dist_k = max_dist*abs(randn(K,1));

%sorted transmit power vector %descending 
dist_vec = sort(dist_k,'ascend'); 

% Path loss exponent
eta_k   = max_eta*abs(randn(K,1));
eta_vec = sort(eta_k,'ascend'); 

% unsorted transmit power vector
transmitpow_k = max_tx_power*abs(randn(K,1));

%sorted transmit power vector %descending 
power_vec = sort(transmitpow_k,'descend'); 

%tx power vec
%power_vec = tx_power_percentage_vec.*tx_pow_k;
pathloss_exp = sqrt(dist_vec.^-eta_vec);

%channel coefficients of each user vec
h_vec =  pathloss_exp.*sqrt(power_vec/2).*(randn(1,N)+1i*randn(1,N))/sqrt(2);

%channel gains of each user vec
g_vec = (abs(h_vec)).^2;

%symbol interference vec
symdur_k = 1*abs(randn(K,1));
sym_dur_vec = sort(symdur_k,'descend');%change here

for k = 1:K
    K_vec(k,1) = K-(k-1);
end

total_sym = sum(K_vec);
est_sym = zeros(total_sym,1)

userdata_vec = rand(K,N)>0.5;          % Generation of data for user1

moddata_vec = 2.*userdata_vec(:,:)-1;    % BPSK modulation 0 -> -1; 1 -> 0 

[sim_delay0] = sim_delayfunc(K, h_vec, userdata_vec, random_iterations,K_vec)

function [sim_delay] = sim_delayfunc(K, h_vec, userdata_vec, random_iterations,K_vec)

for i = 1: random_iterations%random iterations
proptstart(i) = tic; 

%superimposed data
super_signal = userdata_vec.*h_vec;

%AWGN noise 

noise= 1/sqrt(2)*[randn(1,length(super_signal)) + j*randn(1,length(super_signal))];

%received signal 

y = super_signal + noise'; %Addition of Noise

%equalization 
eq_vec = y./h_vec;


est_sym = zeros(3,3);

%sic decoding for each user symbol 
for nbsym = 1:length(K_vec)
    if nbsym ==1
        est_sym(nbsym,:) = eq_vec(K_vec(nbsym),:)>0;
    else
        sub_vec = sum(eq_vec(K_vec(1:nbsym),:),1) - ...
            sum(est_sym(1:nbsym-1,:),1)
        est_sym(nbsym,:) = sub_vec>0
    end
end

propend(i) = toc(proptstart(i));
end
sim_delay = mean(propend);
end
