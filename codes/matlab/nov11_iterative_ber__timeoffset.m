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
iter_K_vec = [4;2;2;2;10];
timeoffset_vec = 0.1*linspace(1,5,20);

for indx = 1: length(timeoffset_vec)

for j = 1:nbiter
rng(j)%random seed
for i = 1: length(iter_K_vec)
initialK = iter_K_vec(i);
K = initialK;
n = K;    
kopt_idx = 1;
timeoff_min = 0.1;
timeoff_max = 0.5;
mod_order = 4;
noise = 0.1;
N = 10^4;   % Number of Bits

% time_offset = 0.5;
max_tx_power = 2;
% receive_pow_ratio = 0;
noisepower = 0.1;
time_offset = 0.5;
%receive_pow_ratio = 5;
% unsorted transmit power vector
transmitpow_k = max_tx_power*abs(randn(K,1));

%sorted transmit power vector %descending 
power_vec = sort(transmitpow_k,'descend'); 
transmit_snrdb = 23;
power_vec(1) =  10^(transmit_snrdb/10)*noisepower/1000;

receive_pow_ratio = 2;
%change here
for d = 2: K
    power_vec(d) = power_vec(d-1)/10^(receive_pow_ratio);
end
clear delta_mat;
%time offsets between users 
%delta_mat: rows -> user index, columns-> symbol index %time offset with
delta_mat = zeros(K,K);
delta_mat(1:K,:) = abs(timeoffset_vec(indx)*randn(K,K));%B1,... Bn, C1....,Cn, ....... %Z1,....Zn

sumsym_dur_vec = tril(delta_mat);

% Path loss exponent
max_eta  = 10;
max_dist = 100;
noisepower = 0.1;
eta_k   = max_eta*abs(randn(K,1));
eta_vec = sort(eta_k,'ascend'); 
eta_vec = 1.5;

% Distances of users from rx
dist_k = max_dist*abs(randn(K,1));

% sorted transmit power vector %descending 
dist_vec = sort(dist_k,'ascend'); 
pathloss_exp = sqrt(dist_vec.^-eta_vec);

%channel coefficients of each user vec
h_vec =  pathloss_exp.*sqrt(power_vec(1:K)/2).*(randn(1,N)+1i*randn(1,N))/sqrt(2);
g_vec = (abs(h_vec)).^2;%channel gains of each user vec
interf_vec = sum(mean(1-g_vec(1:K,:),2).*sumsym_dur_vec,2);
priority = 1;
sinr_th = 1e-6;

clear K_vec;
for k = 1:K
    K_vec(k,1) = K-(k-1);
end

if(K==initialK)

%% strongest user
    delta_i = sumsym_dur_vec(1,1);%known
    p_d     = power_vec(1); %desired power
    p_iw    = power_vec(2:end);%interferes power 
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*p_iw)/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
    p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
    p_bita1_1    = p_err_sym1/log(mod_order);

    delta_i = sum(sumsym_dur_vec(1,1:initialK-1));%known
    p_d     = power_vec(1); %desired power
    p_iw    = power_vec(2:end);%interferes power 
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*p_iw)/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
    p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
    p_bita1_2    = p_err_sym1/log(mod_order);

    delta_i = sum(sumsym_dur_vec(1,1:initialK));%known
    p_d     = power_vec(1); %desired power
    p_iw    = power_vec(2:end);%interferes power 
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*p_iw)/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
    p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
    p_bita1_3    = p_err_sym1/log(mod_order);

    ber_1(u) = 0.3*(p_bita1_1 + p_bita1_2 + p_bita1_3);

    %% second strongest user
    delta_i = sum(sumsym_dur_vec(2,1:initialK-1));%known
    p_d     = power_vec(round(initialK/2)); %desired power
    p_iw    = power_vec(round(initialK/2)+1);%interferes power 
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (delta_i.*p_iw/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
    p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
    p_bita2_1    = p_err_sym1/log(mod_order);

    delta_i = sum(sumsym_dur_vec(2,1:initialK));%known
    p_d     = power_vec(round(initialK/2)); %desired power
    p_iw    = power_vec(round(initialK/2)+1);%interferes power 
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (delta_i.*p_iw/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
    p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
    p_bita2_2    = p_err_sym1/log(mod_order);

    ber_2(u) = 0.5*(p_bita2_1 + p_bita2_2);
%     *[p_bita1_1*p_bita1_2*p_bita1_3 +...
%      (1-p_bita1_1)*(1-p_bita1_2)*(1-p_bita1_3)];

    %% weakest user
    delta_w = [(sum(sumsym_dur_vec(end,1:initialK-1)))];%known
    p_d     = power_vec(end); %desired power
    p_iw    = power_vec(1);   %interferes power 
    fun = @(delta_i) (1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        ((6/(mod_order -1))*sum(delta_i.*[p_iw;])/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
    pp_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_w));
    p_bita31    = pp_err_sym1/log(mod_order);
    
    delta_w = [sum(sumsym_dur_vec(end,1:initialK))];%known
    p_d     = power_vec(end); %desired power
    p_iw    = power_vec(2);  %interferes power 
    fun = @(delta_i) (1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        ((6/(mod_order -1))*sum(delta_i.*[p_iw;])/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
    pp_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_w));
    p_bita32    = pp_err_sym1/log(mod_order);
    
    ber_3(u) = 0.5*(p_bita31+p_bita32);

end 

if i == 5
    %% strongest user
    delta_i = sumsym_dur_vec(1,1);%known
    p_d     = power_vec(1); %desired power
    p_iw    = power_vec(2:end);%interferes power 
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*p_iw)/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
    p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
    p_bita1_1    = p_err_sym1/log(mod_order);

    delta_i = sum(sumsym_dur_vec(1,1:initialK-1));%known
    p_d     = power_vec(1); %desired power
    p_iw    = power_vec(2:end);%interferes power 
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*p_iw)/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
    p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
    p_bita1_2    = p_err_sym1/log(mod_order);

    delta_i = sum(sumsym_dur_vec(1,1:initialK));%known
    p_d     = power_vec(1); %desired power
    p_iw    = power_vec(2:end);%interferes power 
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*p_iw)/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
    p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
    p_bita1_3    = p_err_sym1/log(mod_order);

    ber_1(u) = 1/3*(p_bita1_1 + p_bita1_2 + p_bita1_3);
    
end

if i ==1
    ber1_1 = ber_1(1);
elseif i == 2
    ber_2_1 = ber_2(1);
elseif i == 4 
    ber_3_end = ber_3(end);
elseif i==5
    ber_1_initial = ber_1;
    ber_2_initial = ber_2;
    ber_3_initial = ber_3;
end

avg_ber_1(i) = mean(ber_1);%mean ber of each step
avg_ber_2(i) = mean(ber_2);
avg_ber_3(i) = mean(ber_3);

end

avgavgavg_ber_1(j) = mean(avg_ber_1);%nearest
avgavgavg_ber_2(j) = mean(avg_ber_2);%middle
avgavgavg_ber_3(j) = mean(avg_ber_3);%farthest

user_near(j) = ber1_1;
user_middle(j) = ber_2_1;
user_far(j) = ber_3_end;

user_conv_near(j) = ber_1_initial;
user_conv_middle(j) = ber_2_initial;
user_conv_far(j) = ber_3_initial;
end

propber_1(indx)  = mean(avgavgavg_ber_1);
propber_2(indx)  = mean(avgavgavg_ber_2);
propber_3(indx)  = mean(avgavgavg_ber_3);

user_1(indx)  = mean(user_near);
user_2(indx)  = mean(user_middle);
user_3(indx)  = mean(user_far);

user_1c(indx)  = mean(user_conv_near);
user_2c(indx)  = mean(user_conv_middle);
user_3c(indx)  = mean(user_conv_far);

end

% figure (1)
% semilogy(1:5,smooth(propber_1),'bo-','LineWidth',2)
% hold on 
% semilogy(1:5,smooth(propber_2),'rs-','LineWidth',2)
% hold on
% semilogy(1:5,smooth(propber_3),'k*-','LineWidth',2)

%each user BER vs received power ratio
figure (3)
semilogy(timeoffset_vec ,smooth(user_1),'bo-','LineWidth',1)
hold on 
semilogy(timeoffset_vec ,smooth(user_2),'rs-','LineWidth',1)
hold on
semilogy(timeoffset_vec ,smooth(user_3),'k*-','LineWidth',1)
hold on 
semilogy(timeoffset_vec ,smooth(user_1c),'mo--','LineWidth',1)
hold on 
semilogy(timeoffset_vec,smooth(user_2c),'rs--','LineWidth',1)
hold on
semilogy(timeoffset_vec,smooth(user_3c),'k*--','LineWidth',1)
grid on
