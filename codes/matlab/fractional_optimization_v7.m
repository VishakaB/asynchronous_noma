clc;
clear all;
close all;

%number of users 
K = 4;
%% environmnet
N=10^3;  % Number of Bits
d4 = 10; d3 = 9; d2 = 4; d1 = 3;    %Distances of users from rx
communication_radius = d4;
eta1 = 1;            % Path loss exponent
eta2 = 2;
eta3 = 3;
eta4 = 7;
etath = 4;

noisepower   = 0.1;
max_tx_power = 3;
B           = 10^6;%channel bandwidth
%p1 = max_tx_power ;
%p2 = max_tx_power ;
%p3 = max_tx_power ;
%p4 = max_tx_power ;

%transmission power of each user
transmitpow_k = max_tx_power*abs(randn(K,1));% unsorted transmit power vector

tx_pow_k = sort(transmitpow_k,'descend'); %sorted transmit power vector %descending 

p1 =tx_pow_k(1);%nearest user
p2 =tx_pow_k(2);
p3 =tx_pow_k(3);
p4 =tx_pow_k(4);%farthest user

%channel coefficients of each user 
h1 = sqrt(d1^-eta1)*sqrt(p1/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);%near user
h2 = sqrt(d2^-eta2)*sqrt(p2/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);%second near user
h3 = sqrt(d3^-eta3)*sqrt(p3/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);%third near user
h4 = sqrt(d4^-eta4)*sqrt(p4/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);%far user

h_th = sqrt(communication_radius^-etath)*sqrt(max_tx_power/2)*(randn(1,N)+...
1i*randn(1,N))/sqrt(2);

%channel gains
g1 = (abs(h1)).^2;
g2 = (abs(h2)).^2;
g3 = (abs(h3)).^2;
g4 = (abs(h4)).^2;
g_th = (abs(h_th)).^2;

%% throughput of each user
%considering synchronous uplink noma
SNR = 20:2:40;
snr = db2pow(SNR);

%add the symbol duration????????????????
C1 = B*(log2(1 +g1*p1./(g2.*p2+g3.*p3+g4.*p4+1)));
C2 = B*(log2(1 +g2.*p2./(g3.*p3+g4.*p4+1)));
C3 = B*(log2(1 +g3.*p3./(g4.*p4+1)));
C4 = B*(log2(1 +g4.*p4));

%% optimization problem

%symbols
s1   = 0.6;
s2   = 0.5;
s3   = 0.2;
s4   = 0.1;

%vectors
interf_vec     = zeros(K,1);
factorialk_vec = zeros(K,1);
sumsym_dur_vec = zeros(K,1);
sym_dur_vec = [s1;s2;s3;s4];
power_vec   = [p1;p2;p3;p4];
g_vec       = [g1;g2;g3;g4];%channel gain vector

rate_th     = log2( 1 + sqrt(max_tx_power/2)*g_th/noisepower);
desired_id  = 1;
eth         = 1;
timeslot    = 1;

%factorial k vector
for k = 1:K
    factorialk_vec(k) = K-(k-1);
    k = k-1;
end

for j = 1:K
for k =1:K%interference %only from the next neighbor user
    if k ~= desired_id & k == desired_id+1
        interf_vec(desired_id,1) = power_vec(k)*1*sym_dur_vec(k)...
            + interf_vec(desired_id,1);
        sumsym_dur_vec(desired_id,1)= interf_vec(desired_id,1)...
            + sumsym_dur_vec(desired_id,1);
    end   
   
end
 desired_id =desired_id+1;
end
       
K_vec = [4;3;2;1];
k=1;
for i  =1:10
cvx_begin quiet
    
    variable decision_uk(K,1) binary
   
    %objective
     maximize(decision_uk'*K_vec+(...
     1*sum(decision_uk)-...
     1*sum(decision_uk.^2)))
    %+ 0.5*sum(decision_uk -...
     %  decision_uk.^2)-1)
    
    %constraints
    subject to
    
    %decision_uk <= ones(K,1) 
    %zeros(K,1)  <= decision_uk 
    %sum(decision_uk - decision_uk.^2) <= 0
    %1-decision_uk(k,1) <= 0 
    %decision_uk(k+1,1) <= decision_uk(k,1)
    %decision_uk(k+3,1)<= decision_uk(k+2,1)<= decision_uk(k+1,1)<= decision_uk(k,1)
    %1-sum(decision_uk) <= 0
    %check how to use uvec
    %decision_uk(K,1) <= decision_uk(K-1,1) %use uvec vector
    
    decision_uk'*(max_tx_power*mean(g_th,2)*...
                  (noisepower^2+(interf_vec))... 
                  -power_vec.*mean(g_vec,2)*noisepower^2) <= 0
              
    decision_uk'*sumsym_dur_vec <=  sum(sumsym_dur_vec)/timeslot
    
cvx_end

decision_uk
end