%v16
%initialized date: 23 may 2022
%last updated: 25 may 2022
%energy efficiency of NOMA asynchronous D2D SIC decoding
%Output: Energy efficiency based ...
%on number of users in ...
%proposed optimized sic traingle decoding method
%results: ber vs nbsuperimposed_data
clc;
clear all;
close all;


%% input data: environmnet
%--------------------------------------------------------------------------
%%scalars
%number of users 
alldatadecoded = false;
mpriority = 20;
x = zeros(mpriority-1,1);
y = zeros(mpriority-1,1);
z = zeros(mpriority-1,1);
zz = zeros(mpriority-1,1);
r = zeros(mpriority-1,1);
s = zeros(mpriority-1,1);

%% 
%Number of Bits
N=10^4;  

communication_radius = 30;%change this 
max_dist     = 100;%meters
max_eta      = 15;%change here
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

%random iterations
%--------------------------------------------------------------------------
K = 3;%number of superimposed data
for indx = 2:1:mpriority 
initialK = indx;
K  = indx;
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
symdur_k = abs(randn(K,1));%change here
sym_dur_vec = sort(symdur_k,'descend');%change here

pr_vec = [0.5;1;1.5;2;2.5;3;3.5;4;4.5;5;5.5;6;6.5;7.5;8;8.5;10;12;15;20];
   
%fprintf("indx pr  %i %f\n",initialK,pr_vec(8));%priority: 4%change here
K = initialK;
%[x(initialK-1),y(initialK-1),z(initialK-1),zz(initialK-1),r(initialK-1),s(initialK-1)] = seqsic(initialK,alldatadecoded,K,...
    %pr_vec(4),power_vec,sym_dur_vec,g_vec,max_tx_power,timeslot);
r;
s;

end

save x.mat;
save y.mat;
save z.mat;
save zz.mat;


for i = 1:100 %random iterations 
v =1;
while (alldatadecoded == false & K>1) 

%nsymbols vector of each user: K vec #loop
clear K_vec;
for k = 1:K
    K_vec(k,1) = K-(k-1);
end
initialK_vec = K_vec;

miter =10;
priority_max = 200;
priority = 0.5;
lambda1 = priority ;%change this%energy saving priority %left energy is low
learn_rate = 0.4;
tolerance2 = 0.5;%lambda
tolerance = 0.02;%uk
Rmin = 0.000001;
sinr_th = 1e-6;
%--------------------------------------------------------------------------

%vectors
interf_vec     = zeros(K,1);
factorialk_vec = zeros(K,1);
sumsym_dur_vec = zeros(K,1);
desired_id   = 1;

clear nsym;
clear user_strength;
for k = 1:K
    nsym(k,1) = K-(k-1);
    user_strength(k,1) = k;
end

clear delta_mat;
clear reverse_delta_mat;
%time offsets between users 
%delta_mat: rows -> user index, columns-> symbol index %time offset with
delta_mat = zeros(K,K);
delta_mat(1,:) = zeros(K,1);
delta_mat(2:K,:) = 0.5*rand(K-1,K);%B1,... Bn, C1....,Cn, ....... %Z1,....Zn

reverse_delta_mat(K,:)  = zeros(K,1);
reverse_delta_mat(1:K-1,:) = 0.5*rand(K-1,K);%A1, A2

for j = 1:K%interference vector loop
for k =1:K
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

%% optimization problem
nbiter = 10;
noisepower   = 0.1;
if(K>1)
opt_decision_uk = ones(K,1);

%decision_uk = ones(K,1);%initialize uk
learn_rateuk = 0.3;
n_iter = 1000;
toleranceuk = 0.5;
n = K;
%% 
cvx_begin
   variable decision_uk(n,1) 
   dual variables lam mu vo
   minimize(-decision_uk'*K_vec)
   subject to
      lam : -sum(decision_uk)+ sum(decision_uk.^2)<=0;
      %mu : sum(decision_uk) > 0;
      mu: decision_uk'*sumsym_dur_vec- ...
               1/(priority+0.001)*priority_max/timeslot <=0
      vo:  decision_uk.*((noisepower^2 + interf_vec + power_vec(1:K).*mean(g_vec(1:K,:),2))...
                -power_vec(1:K).*mean(g_vec(1:K,:),2)*(1+1/sinr_th) ) <= 0   
cvx_end

echo off

decision_uk>0.8
lam
mu
vo

if K<=1 
    alldatadecoded=true;
    
    disp('break');
    %break;
end%end if 

end

%% throughput of each user
%considering synchronous uplink noma
E_max = 20;

SINR_k = power_vec(1:K).*mean(g_vec(1:K,:),2)./(interf_vec(1:K)+noisepower^2);

throughput_vec = log(1+SINR_k);

total_throughput = sum(throughput_vec);
total_throughput = 2;%fix here????
%% energy efficiency 
%proposed optimal sic
for k = 1:length(opt_decision_uk)
    K_vec(k,1) = length(opt_decision_uk)-(k-1);
end
total_energ_consump = E_max - E_max^(exp(-log(2)/1000*opt_decision_uk'*K_vec));
if(total_energ_consump<0)
    fprintf("opt_decision_uk %i %f %f %f %f\n",K,length(opt_decision_uk),K_vec,E_max,E_max^(exp(-log(2)/1000*opt_decision_uk'*K_vec)));
end
energy_eff(v) = total_throughput/(total_energ_consump +0.01);
energy_eff;
%%conv sic

total_energ_consump_conv = E_max - E_max^(exp(-log(2)/1000*sum(initialK_vec)));
energy_eff_conv(v) = total_throughput/(total_energ_consump_conv+0.01);
energy_eff_conv;
if(energy_eff_conv<0)
    pause on;
end
%% complexity analysis 

%only sic decoding 
sic_complextiyprop(v) = sum(opt_decision_uk)^2*log(1/0.01);
sic_complextiyconv(v) = sum(initialK)^2;

v = v+1;
K = K-sum(opt_decision_uk);%update K
end%end while

avgenergy_effconv(i) = abs(mean(energy_eff_conv));
avgenergy_eff(i) = abs(mean(energy_eff));

avgcomplexity_conv(i) = mean(sic_complextiyconv);
avgcomplexity_prop(i) = mean(sic_complextiyprop);


end

a = abs(mean(energy_eff_conv));
b = abs(mean(energy_eff));

c = mean(sic_complextiyconv);
d = mean(sic_complextiyprop);




fprintf("avg energy eff proposed %f\n",mean(energy_eff));
fprintf("avg energy eff conv %f\n",mean(energy_eff_conv));
fprintf("complexity conv %f\n",mean(sic_complextiyconv));
fprintf("complexity prop %f\n",mean(sic_complextiyprop));


