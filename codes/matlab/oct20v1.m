a = zeros(5,1);
b = zeros(5,1);

for iter = 1: 5
transmit_snrdb = 10;
sim_delay_prop = 0;
sim_delay_conv = 0;
nbiterations = 1;
nbrandom_iterations = 3;
receive_pow_ratio = 0.5;
initialK = 4;
N = 10^4;
priority = 1;

for receive_pow_ratioi = 1: 1%number of superimposed data loop
alldatadecoded = false;
%fprintf('receive_pow_ratio %i\n',receive_pow_ratio);

for i = 1:nbrandom_iterations %random iterations 
v =1;
K = initialK;
max_dist     = 100;%meters
max_eta      = 10;
etath        = 4;%change this 
noisepower   = 0.1;
max_tx_power = 2;%change this
B            = 1;%channel bandwidth
timeslot     = 1;

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
for d = 2: K
    power_vec(d) = power_vec(d-1)/10^(receive_pow_ratio);
end

%tx power vec
%power_vec = tx_power_percentage_vec.*tx_pow_k;
eta_vec = 1.5;
pathloss_exp = sqrt(dist_vec.^-eta_vec);

%channel coefficients of each user vec
h_vec =  pathloss_exp.*sqrt(power_vec/2).*(randn(1,N)+1i*randn(1,N))/sqrt(2);

%channel gains of each user vec
g_vec = (abs(h_vec)).^2;

%symbol interference vec
symdur_k = 0.5*abs(randn(K,1));
sym_dur_vec = sort(symdur_k,'descend');%change here

%fprintf('iteration count %f\n',i)
%fprintf('K %i\n',K);
while (alldatadecoded == false & K> 1) 
%fprintf('v %f\n',v)
%nsymbols vector of each user: K vec #loop
clear K_vec;
for k = 1:K
    K_vec(k,1) = K-(k-1);
end

initialK_vec = K_vec;
miter = 10;
priority_max = 120;
%lambda1 = priority;%change this%energy saving priority %left energy is low
learn_rate = 0.4;
tolerance2 = 0.5;%lambda
tolerance = 0.02;%uk
Rmin = 1e-6;
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
time_offset = 0.15;
delta_mat = zeros(K,K);
delta_mat(1,:) = zeros(K,1);
delta_mat(1:K,:) = abs(time_offset*randn(K,K));%B1,... Bn, C1....,Cn, ....... %Z1,....Zn

reverse_delta_mat(K,:)  = zeros(K,1);
reverse_delta_mat(1:K,:) = abs(time_offset*randn(K,K));%A1, A2

desired_id = 1;
desired_id = 1; k= 1;
for desired_id = 1:K%interference vector loop
for k = 1:K
    %interference vec %only from the next neighbor user
    if k ~= desired_id & desired_id == 1        
        sumsym_dur_vec(desired_id,1) = sum(sum(delta_mat(desired_id+1:end,:)))
        
    elseif k ~= desired_id & desired_id <K
        sumsym_dur_vec(desired_id,1) = sum(sum(delta_mat(desired_id+1:end,:)))+...
            sum(reverse_delta_mat(desired_id-1,:))
    elseif desired_id==K
        sumsym_dur_vec(desired_id,1) = sum(sum(delta_mat(desired_id-1:end,:))); 
       
    end
end
end

sumsym_dur_vec = [0.01;0.1;0.3;0.8];

desired_id = 1;
for desired_id = 1:K%interference vector loop
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

interf_vec = [0.03;0.5;0.7;0.9];
%% optimization problem
pastdelay = 0;
nbiter = 10;
noisepower   = 0.1;

if(K>1)
opt_decision_uk = 0.1*ones(K,1);

proptstart(v)=tic;

n = K;
%% 
cvx_begin quiet
   variable decision_uk(n,1)
   dual variables var1 var2 var3 var4
   minimize(-decision_uk'*K_vec)
   subject to
      var1 : -sum(decision_uk)+ sum(decision_uk.^2)<=0;

      var2: decision_uk'*sumsym_dur_vec-1/(priority+0.001)*priority_max/timeslot <= 0
           
      var3: decision_uk.*((noisepower^2 + interf_vec + power_vec(1:K).*mean(g_vec(1:K,:),2))...
                -power_vec(1:K).*mean(g_vec(1:K,:),2)*(1+1/sinr_th)) <= 0  
      var4: decision_uk(1:K-1)>= decision_uk(2:K)     
cvx_end

echo off
%fprintf('cvx_slvtol %f\n',cvx_iterations);
diary on;
decision_uk = decision_uk>0.8;

%complexity prop
sic_complextiyprop(v) = sum(decision_uk)^2*log(1/tolerance)*log(1/tolerance);

proptend(v)    = toc(proptstart(v));

%complexity analysis
%proptend(v) = toc(proptstart(v));

random_iterations = 10;
N=10^4;  

userdata_vec = rand(initialK,N)>0.5;          % Generation of data for user1
clear K_vec;
for k = 1:K
    K_vec(k,1) = length(opt_decision_uk)-(k-1);
end

%% ber analysis
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
delta_mat(2:K,:) = 0.5*randn(K-1,K);%B1,... Bn, C1....,Cn, ....... %Z1,....Zn

reverse_delta_mat(K,:)  = zeros(K,1);
reverse_delta_mat(1:K-1,:) = 0.5*randn(K-1,K);%A1, A2

mod_order = 4;
timeoff_min = 0.01;
timeoff_max = 0.5;

% if K==initialK
%     [berfinalconv] = berfunc(power_vec, noisepower, nsym, user_strength,...
%     delta_mat,reverse_delta_mat,mod_order,timeoff_min,timeoff_max);
% end
% 
% if K>=2
% [berfinal0] = berfunc(power_vec, noisepower, nsym, user_strength,...
%     delta_mat,reverse_delta_mat,mod_order,timeoff_min,timeoff_max);
% end

% ber_propuser1(v) = berfinal0(1);
% ber_convuser1(v) = berfinalconv(1);
% ber_propuserw(v) = berfinal0(K);
% ber_convuserw(v) = berfinalconv(K);
% ber_propuseri(v) = berfinal0(round(K/2)+1);
% ber_convuseri(v) = berfinalconv(round(K/2)+1);
%% sim delay

% if(K>1)
%     [sim_delay_prop(v),ber_prop(v)] = sim_delayfunc(K, h_vec(1:K,:), userdata_vec(1:K,:), random_iterations,K_vec);
% end
% if(K>1)
%     [sim_delay_conv(v),ber_conv(v)] = sim_delayfunc(initialK, h_vec(1:initialK,:), userdata_vec, random_iterations,initialK_vec);
% end

if K<=1 
    alldatadecoded = true;
    %proptend(v)    = toc(proptstart(v));
    %disp('break');
    nbiterations  = nbiterations+1;
    iterations(v) = nbiterations;
else
    nbiterations  = nbiterations+1;
    iterations(v) = nbiterations;
    %fprintf('nbiterations %i\n',nbiterations);
    decision_uk(1) = 1;%make sure decision_uk is non zero
    %fprintf('recalc sic %f\n',v)
    %break;
end%end if 
decision_uk ;
K = K-sum(decision_uk);%update K
end%end if 

%% throughput of each user
%considering synchronous uplink noma
E_max = 10;

SINR_k = power_vec(1:K).*mean(g_vec(1:K,:),2)./(interf_vec(1:K)+noisepower^2);

throughput_vec = log(1+SINR_k);

total_throughput = sum(throughput_vec);
total_throughput = 2;%fix here????

%% energy efficiency 
%proposed optimal sic
for k = 1:length(decision_uk)
    K_vec(k,1) = length(decision_uk)-(k-1);
end

total_energ_consump = E_max - E_max^(exp(-log(2)/1000*decision_uk'*K_vec));
energy_eff(v) = total_throughput/(total_energ_consump +0.01);
energy_eff;

%%conv sic
total_energ_consump_conv = E_max - E_max^(exp(-log(2)/1000*sum(initialK_vec)));
energy_eff_conv(v) = total_throughput/(total_energ_consump_conv+0.01);
energy_eff_conv;
if(energy_eff_conv<0)
    pause on;
end

%complexity analysis
convtstart(v) = tic;
for l = 1: sum(initialK_vec)
   for g = 1:sum(initialK_vec) 
   end
end
convtend(v) = toc(convtstart(v));

%% complexity analysis 

%only sic decoding 
sic_complextiyconv(v) = sum(initialK)^2;

v = v+1;%all user decoding index v

end %end while

%i: random iteration index
avgenergy_eff(i) = abs(mean(energy_eff));
avgenergy_effconv(i) = abs(mean(energy_eff_conv));

avgcomplexity_prop(i) = mean(sic_complextiyprop);
avgcomplexity_conv(i) = mean(sic_complextiyconv);

avgdelay_conv(i) = mean(sim_delay_conv);
avgdelay_prop(i) = mean(sim_delay_prop);

%avgiterations = mean(propend);
totaldelay_prop(i) = mean(sim_delay_prop+iterations*0.01);
avggradientdelay(i) = mean(proptend);

% avgberconv(i) = mean(ber_conv);
% avgberprop(i) = mean(ber_prop);
% 
% avgberconvuser1(i) = mean(ber_convuser1);
% avgberuser1(i) = mean(ber_propuser1);
% 
% avgberweakprop(i) = mean(ber_propuserw);
% avgberweakconv(i) = mean(ber_convuserw);
% 
% avgberinterprop(i) = mean(ber_propuseri);
% avgberinterconv(i) = mean(ber_convuseri);

end

end

a(iter) = abs(mean(energy_eff_conv));
b(iter) = abs(mean(energy_eff));

c = mean(sic_complextiyconv);
d = mean(avgcomplexity_prop);

e = mean(avgdelay_conv);
f = mean(avgdelay_prop);

g = mean(totaldelay_prop);
h = mean(avggradientdelay);

% i = mean(avgberconv);
% j = mean(avgberprop);

% strongconv = mean(avgberconvuser1); 
% strongprop = mean(avgberuser1); 
% 
% weakconv = mean(avgberweakconv);
% weakprop = mean(avgberweakprop);
% 
% intermconv = mean(avgberinterconv); 
% intermprop = mean(avgberinterprop); 

end