%date: 13 july 2022
%goal: optimize the secrecy 
%for irc
close all
clc
clear all

%% initialization 
K = 2;
max_dist     = 100;%meters
max_eta      = 10;
etath        = 4;%change this 
noisepower   = 0.1;
max_tx_power = 500;%change this
B            = 1;%channel bandwidth
timeslot     = 1;
N            = 10^3;

transmit_snrdb_vec = linspace(1,100,20);
total = length(transmit_snrdb_vec);
t = 1;%sinr index 
cs = zeros(3,total);
alpha1_vec = [0.0001,0.001,0.1,0.5];
alpha2_vec = [0.4,0.3,0.15,0.1];

for cases = 1:4
t = 1; 
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
power_vec_e = sort(transmitpow_k,'descend'); 
power_vec(1) = max_tx_power;
%receive_pow_ratio = receive_pow_ratio_vec(receive_pow_ratioi);
%change here
power_vec(1) =  10^(transmit_snrdb/10)*noisepower;
power_vec_e(1) = (alpha1_vec(cases))*10^(transmit_snrdb/10)*noisepower;
%tx power vec
%power_vec = tx_power_percentage_vec.*tx_pow_k;
eta_vec = 1.5;
pathloss_exp = sqrt(dist_vec.^-eta_vec);

%channel coefficients of each user vec
h_vec =  pathloss_exp.*sqrt(power_vec/2).*(randn(1,N)+1i*randn(1,N))/sqrt(2);
h_vec_e =  pathloss_exp.*sqrt(power_vec_e/2).*(randn(1,N)+1i*randn(1,N))/sqrt(2);

%channel gains of each user vec
g_vec = (abs(h_vec)).^2;
g_vec_e = (abs(h_vec_e)).^2;

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

sinr_thmin = 1;
sinr_thmax = 0.001;
 
SINRd = sum(power_vec(1:K).*mean(g_vec(1:K,:),2)./(noisepower^2 + (interf_vec)));
SINRe = sum(power_vec_e(1:K).*mean(g_vec_e(1:K,:),2)./(noisepower^2 + (interf_vec)));

cvx_begin quiet
   variable decision_uk(K,1)
   dual variables var1 var2 var3
   minimize (-((((1+sum(decision_uk.*power_vec(1:K).*mean(g_vec(1:K,:),2)./...
(noisepower^2 + (interf_vec))))))...
-(1+sum(decision_uk.*power_vec_e(1:K).*mean(g_vec_e(1:K,:),2)./...
(noisepower^2 + (interf_vec))))))%maximize secrecy
   subject to
      var1: decision_uk.*((noisepower^2 + interf_vec +...
          power_vec(1:K).*mean(g_vec(1:K,:),2))...
                -power_vec(1:K).*mean(g_vec(1:K,:),2)*...
                (1+1/sinr_thmin)- (decision_uk.*((noisepower^2 + interf_vec +...
          power_vec_e(1:K).*mean(g_vec_e(1:K,:),2))...
                -power_vec_e(1:K).*mean(g_vec_e(1:K,:),2)*...
                (1+1/sinr_thmax)))) <= 0  
      var3: decision_uk >=0  
cvx_end

decision_uk = decision_uk>0.5
var3;
cs(cases,t) = ((log((1+sum(decision_uk.*power_vec(1:K).*mean(g_vec(1:K,:),2)./(noisepower^2 + (interf_vec)))))));
t = t+1;
end
C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colros.
%% performance analysis
%secrecy capacity vs SINRd
hold on;
grid on;
plot(transmit_snrdb_vec,cs(cases,:),'color',C{cases},'marker','o');
xlabel('SNR legitimate');
ylabel('Secrecy capacity');

end
legend({'\alpha_{{leg}_1} = 0.001';'\alpha_{{leg}_2} = 0.01';'\alpha_{{leg}_3} = 0.1';'\alpha_{{leg}_4} = 0.5';})