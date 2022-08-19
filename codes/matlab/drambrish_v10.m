%date: 13 july 2022
%goal: optimize the secrecy 
%close all
clc
clear all

%% initialization 
K = 2;
max_dist     = 100;%meters
max_eta      = 10;
etath        = 4;%change this 
noisepower   = 0.1;
max_tx_power = 1000;%change this
B            = 1;%channel bandwidth
timeslot     = 1;
N            = 10^3;

transmit_snrdb_vec = linspace(1,100,20);
total = length(transmit_snrdb_vec);
t = 1;%sinr index 
cs = zeros(total,1);
for transmit_snrdb = transmit_snrdb_vec
  
time_offset = 0.15;

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

power_vec(1) = max_tx_power;
%receive_pow_ratio = receive_pow_ratio_vec(receive_pow_ratioi);
%change here
power_vec(1) =  10^(transmit_snrdb/10)*noisepower;

%tx power vec
%power_vec = tx_power_percentage_vec.*tx_pow_k;
eta_vec = 1.5;
pathloss_exp = sqrt(dist_vec.^-eta_vec);

%channel coefficients of each user vec
h_vec =  pathloss_exp.*sqrt(power_vec/2).*(randn(1,N)+1i*randn(1,N))/sqrt(2);

%channel gains of each user vec
g_vec = (abs(h_vec)).^2;

for k = 1:K
    nsym(k,1) = K-(k-1);
    user_strength(k,1) = k;
end

delta_mat = zeros(K,K);
delta_mat(1,:) = zeros(K,1);
delta_mat(2:K,:) = time_offset*ones(K-1,K);%B1,... Bn, C1....,Cn, ....... %Z1,....Zn

reverse_delta_mat(K,:)  = zeros(K,1);
reverse_delta_mat(1:K-1,:) = time_offset*ones(K-1,K);%A1, A2

desired_id = 1;
for j = 1:K%interference vector loop
for k = 1:K
    %interference vec %only from the next neighbor user
    if k ~= desired_id & k == desired_id+1 & desired_id == 1        
        sumsym_dur_vec(desired_id,1) = sum(delta_mat(desired_id+1,:))
        
    elseif k ~= desired_id & k == desired_id+1 & desired_id > 1 & desired_id <K
        sumsym_dur_vec(desired_id,1) = sum(delta_mat(desired_id+1,:))+...
            sum(reverse_delta_mat(desired_id-1,:))
    elseif desired_id==K
        sumsym_dur_vec(desired_id,1) = sum(reverse_delta_mat(desired_id-1,:));       
    end
 desired_id = desired_id+1;
end
end

interf_vec = zeros(K,1);
desired_id = 1;
for j = 1:K%interference vector loop
for k = 1:K%neighbor users index
     if k ~= desired_id  %strongest
        interf_vec(desired_id) = interf_vec(desired_id)+sum(power_vec(k).*...
            mean(g_vec(k,:),2).*delta_mat(desired_id,:));
     elseif desired_id == 1
        if(k<K)
        interf_vec(desired_id) = sum(power_vec(k+1).*...
            mean(g_vec(k+1,:),2).*delta_mat(desired_id+1,:));
        end
     end
end
desired_id = desired_id+1;
end

%% optimization
%capacity of legitimate user 
alpha1 = 0.9;%power allocation coefficients
alpha2 = 0.1;
rho    = power_vec(1);%transmission power
omega  = 0.8;
nev    = -0.1;

num = alpha1.*rho.*omega;
den = (1-rho).*nev+1;
SINRd = num/den;

num2 = alpha2.*rho;
den2 = (1-rho).*nev+1;
SINRe = num2/den2;

gamma_leg = log((1+SINRd)/(1+SINRe));

sinr_thmin = 0.1;
sinr_thmax = 1;
 
cvx_begin quiet
   variable decision_uk
   dual variables var1 var2 var3
   minimize(-decision_uk*abs(log((1+SINRd)/(1+SINRe))))%maximize secrecy
   subject to
      var1:(log(1+SINRd)-sinr_thmin)...
                 -(sinr_thmax-log(1+SINRe))      
      var2: decision_uk*((log(1+SINRd)-sinr_thmin)...
                 -(sinr_thmax-log(1+SINRe)))>=0    
      var3: decision_uk >=0  
cvx_end

decision_uk = decision_uk>0.5
var3;
cs(t) = decision_uk*(log((1+SINRd)/(1+SINRe)));
t = t+1;
end

%% performance analysis
%secrecy capacity vs SINRd
figure
plot(transmit_snrdb_vec,cs,'r--')
xlabel('SNR legitimate')
ylabel('Secrecy capacity')