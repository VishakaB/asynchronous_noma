clc;
clear all;
close all;
%% ber vs total number of users
% K = 5;
% K = 10;

% sumsym_dur_vec = [0.1,0.1,0.1;0.2,0.3,0.4;0.5,0.4,0.3;]*0.01;
% power_vec = [0.1;0.01;0.001;0.0001;0.00001];
k_vec = linspace(1,15,15);
u = 1;%initiate u
nbiter = 3;

for kindx = 2: length(k_vec)
K = k_vec(kindx);

initialK = K;
for j = 1:nbiter
rng(j)%random seed
alldatadecoded = false;
% ber_1 = zeros(K,1);
% ber_2 = zeros(K,1);
% ber_3 = zeros(K,1);
% for v = 1:K
% ber_1 = [];
% ber_2 = [];
% ber_3 = [];
kopt_idx = 1;
while (alldatadecoded == false && K> 1) 
timeoff_min = 0.1;
timeoff_max = 0.5;
mod_order = 4;
noise = 0.1;
N=10^4;   % Number of Bits
n = K;
alldatadecoded = false;
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

power_vec(1) =  max_tx_power;
receive_pow_ratio = 5;
%change here
for d = 2: K
    power_vec(d) = power_vec(d-1)/10^(receive_pow_ratio);
end
clear delta_mat;
%time offsets between users 
%delta_mat: rows -> user index, columns-> symbol index %time offset with
delta_mat = zeros(K,K);
delta_mat(1:K,:) = abs(time_offset*randn(K,K));%B1,... Bn, C1....,Cn, ....... %Z1,....Zn

sumsym_dur_vec = tril(delta_mat);

% Path loss exponent
max_eta      = 10;
max_dist = 100;
noisepower = 0.1;
eta_k   = max_eta*abs(randn(K,1));
eta_vec = sort(eta_k,'ascend'); 
eta_vec = 1.5;

%Distances of users from rx
dist_k = max_dist*abs(randn(K,1));

%sorted transmit power vector %descending 
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
    %% strongest user conv

%% strongest user
    delta_i = sumsym_dur_vec(1,1);%known
    p_d     = power_vec(1); %desired power
    p_iw    = power_vec(2);%interferes power 
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (delta_i.*p_iw/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max);
    p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
    p_bita1_1    = p_err_sym1/log(mod_order);

    delta_i = sum(sumsym_dur_vec(1,1:initialK-1));%known
    p_d     = power_vec(1); %desired power
    p_iw    = power_vec(2);%interferes power 
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (delta_i.*p_iw/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max);
    p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
    p_bita1_2    = p_err_sym1/log(mod_order);

    delta_i = sum(sumsym_dur_vec(1,1:initialK));%known
    p_d     = power_vec(1); %desired power
    p_iw    = power_vec(2);%interferes power 
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (delta_i.*p_iw/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max);
    p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
    p_bita1_3    = p_err_sym1/log(mod_order);

    ber_1(u) = 1/3*(p_bita1_1 + p_bita1_2 + p_bita1_3);

    %% second strongest user
    delta_i = sum(sumsym_dur_vec(2,1:initialK-1));%known
    p_d     = power_vec(2); %desired power
    p_iw    = power_vec(1);%interferes power 
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        ((6/(mod_order -1))*delta_i.*p_iw/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max);
    p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
    p_bita2_1    = p_err_sym1/log(mod_order);

    delta_i = sum(sumsym_dur_vec(2,1:initialK));%known
    p_d     = power_vec(2); %desired power
    p_iw    = power_vec(1);%interferes power 
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        ((6/(mod_order -1))*delta_i.*p_iw/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max);
    p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
    p_bita2_2    = p_err_sym1/log(mod_order);

    ber_2(u) = 0.5*(p_bita2_1 + p_bita2_2);
    % *[pp_bita1_1*pp_bita1_2*pp_bita1_3 +...
    %   (1-pp_bita1_1)*(1-pp_bita1_2)*(1-pp_bita1_3)]

    %% weakest user
    delta_i = [(sum(sumsym_dur_vec(end,1:initialK-1)));...
        (sum(sumsym_dur_vec(end,1:initialK)))];%known
    p_d     = power_vec(end); %desired power
    p_iw    = power_vec(2);   %interferes power 
    p_iw1   = power_vec(1);
    fun = @(delta_i) 1-(1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        ((6/(mod_order -1))*sum(delta_i.*[p_iw;p_iw1;])/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
    pp_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
    p_bita3    = pp_err_sym1/log(mod_order);

    ber_3(u) = p_bita3;

end    
%% 
cvx_begin quiet
   variable decision_uk(n,1)
   dual variables lam gan ha
   minimize(-decision_uk'*K_vec)
   subject to
      lam: -sum(decision_uk)+ sum(decision_uk.^2)<=0;

      gan: decision_uk'*sum(sumsym_dur_vec,2)-1.2 <= 0
           
      ha: decision_uk.*((noisepower^2 + interf_vec(1:K) + power_vec(1:K).*mean(g_vec(1:K,:),2))...
                -power_vec(1:K).*mean(g_vec(1:K,:),2)*(1+1/sinr_th)) <= 0  
cvx_end

decision_uk = round(decision_uk);
kopt = sum(decision_uk);
sindx = length(sumsym_dur_vec);

%% strongest user
delta_i = sumsym_dur_vec(1,1);%known
p_d     = power_vec(1); %desired power
p_iw    = power_vec(1:end);%interferes power 
fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
    (delta_i.*p_iw/2+noise))))).^2;
q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
pp_bita1_1    = p_err_sym1/log(mod_order);

delta_i = sum(sumsym_dur_vec(1,1:kopt-1));%known
p_d     = power_vec(1); %desired power
p_iw    = power_vec(1:end);%interferes power 
fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
    (delta_i.*p_iw/2+noise))))).^2;
q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
pp_bita1_2    = p_err_sym1/log(mod_order);

delta_i = sum(sumsym_dur_vec(1,1:kopt));%known
p_d     = power_vec(1); %desired power
p_iw    = power_vec(1:end);%interferes power 
fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
    (delta_i.*p_iw/2+noise))))).^2;
q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
pp_bita1_3    = p_err_sym1/log(mod_order);

pber_1(u) = 1/3*(pp_bita1_1 + pp_bita1_2 + pp_bita1_3)

%% second strongest user
delta_i = sum(sumsym_dur_vec(2,1:kopt-1));%known
p_d     = power_vec(2); %desired power
p_iw    = power_vec(end);%interferes power 
fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
    ((6/(mod_order -1))*delta_i.*p_iw/2+noise))))).^2;
q   = integral(fun,timeoff_min,timeoff_max);
p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
pp_bita2_1    = p_err_sym1/log(mod_order);

delta_i = sum(sumsym_dur_vec(2,1:kopt));%known
p_d     = power_vec(2); %desired power
p_iw    = power_vec(end);%interferes power 
fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
    ((6/(mod_order -1))*delta_i.*p_iw/2+noise))))).^2;
q   = integral(fun,timeoff_min,timeoff_max);
p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
pp_bita2_2    = p_err_sym1/log(mod_order);

pber_2(u) = 0.5*(pp_bita2_1 + pp_bita2_2)
% *[pp_bita1_1*pp_bita1_2*pp_bita1_3 +...
%   (1-pp_bita1_1)*(1-pp_bita1_2)*(1-pp_bita1_3)]

%% weakest user
delta_i = [(sum(sumsym_dur_vec(end,1:kopt-1)));...
   (sum(sumsym_dur_vec(end,kopt:end)))];%known
p_d     = power_vec(end); %desired power
p_iw    = power_vec(2);   %interferes power 
p_iw1   = power_vec(1);
p_iw_vec = power_vec(1:end-1,1);
fun = @(delta_i) 1-(1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
    ((6/(mod_order -1))*sum(delta_i.*p_iw_vec)/2+noise))))).^2;
q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
pp_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
pp_bita3    = pp_err_sym1/log(mod_order);

pber_3(u) = pp_bita3

%*[ber_1(u)*ber_2(u)+ber_1(u)*(1 - ber_2(u))+...
 %   (1-ber_1(u))*ber_2(u)+(1-ber_1(u))*(1-ber_2(u))]
if K<=1 
    alldatadecoded = true;
    kopt_idx = 1;
end

K = K-sum(decision_uk);
%fprintf('K: %d\n',K);
u = u+1;%number of times iteratively decoded
% %update K  

% end
end%end

avg_pber_1(j) = mean(pber_1);
avg_pber_2(j) = mean(pber_2);
avg_pber_3(j) = mean(pber_3);

avg_ber_1(j) = mean(ber_1);
avg_ber_2(j) = mean(ber_2);
avg_ber_3(j) = mean(ber_3);

end

avgavgpber_1(kindx) = mean(avg_pber_1);
avgavgpber_2(kindx) = mean(avg_pber_2);
avgavgpber_3(kindx) = mean(avg_pber_3);

avgavgber_1(kindx) = mean(avg_ber_1);
avgavgber_2(kindx) = mean(avg_ber_2);
avgavgber_3(kindx) = mean(avg_ber_3);

end

avgavgavg_ber_1 = mean(avgavgber_1);%nearest
avgavgavg_ber_2 = mean(avgavgber_2);%middle
avgavgavg_ber_3 = mean(avgavgber_3);%farthest

avgavgavg_convber_1 = mean(avgavgber_1);%nearest
avgavgavg_convber_2 = mean(avgavgber_2);%middle
avgavgavg_convber_3 = mean(avgavgber_3);%farthest

figure (4)
semilogy(k_vec,smooth(avgavgpber_1),'bo-','LineWidth',1)
hold on 
semilogy(k_vec,smooth(avgavgber_1),'bo--','LineWidth',2)
hold on
semilogy(k_vec,smooth(avgavgpber_2),'rs-','LineWidth',1)
hold on 
semilogy(k_vec,smooth(avgavgber_2),'rs--','LineWidth',2)
hold on 
semilogy(k_vec,smooth(avgavgpber_3),'k+-','LineWidth',1)
hold on 
semilogy(k_vec,smooth(avgavgber_3),'k+--','LineWidth',2)
grid on;