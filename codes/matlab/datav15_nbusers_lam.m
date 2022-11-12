%v16
%initialized date: 23 may 2022
%last updated: 25 may 2022
%energy efficiency of NOMA asynchronous D2D SIC decoding
%Output: Energy efficiency based ...
%on number of users in ...
%proposed optimized sic traingle decoding method
%results testing: ee vs total interference in sic triangle
clc;
clear all;
close all;

%% input data: environment
%--------------------------------------------------------------------------
%%scalars
%number of users 
alldatadecoded = false;
mpriority = 20;
x = zeros(mpriority,1);
y = zeros(mpriority,1);
z = zeros(mpriority,1);
zz = zeros(mpriority,1);
%% 
% Number of Bits
N=10^4;  

communication_radius = 30;%change this 
max_dist     = 100;%meters
max_eta      = 15;
etath        = 4;%change this 
noisepower   = 0.1;
max_tx_power = 2;%change this
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
userK_vec = [3,5,8,15,20];
K = userK_vec(5);%number of superimposed data
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

pr_vec = [0.5;1;1.5;2;2.5;3;3.5;4;4.5;5;5.5;6;6.5;7.5;8;8.5;10;12;15;20];

for indx = 1:1:mpriority 
    %rng(1);%same random seed
    fprintf("indx pr  %i %f\n",indx,pr_vec(indx));
    
    initialK = K;
    [x(indx),y(indx),z(indx),zz(indx)] = seqsic(initialK,alldatadecoded,K,...
        pr_vec(indx),power_vec,sym_dur_vec,g_vec,max_tx_power,timeslot);
    x
    y
end

save x.mat;
save y.mat;
save z.mat;
save zz.mat;


function [a,b,c,d] = seqsic(initialK,alldatadecoded,K,priority,power_vec,sym_dur_vec,...
g_vec,max_tx_power,timeslot)

for nbusers = initialK: initialK%number of superimposed data loop
for i = 1:100 %random iterations 
v =1;
while (alldatadecoded == false) 

%nsymbols vector of each user: K vec #loop
clear K_vec;
for k = 1:K
    K_vec(k,1) = K-(k-1);
end
initialK_vec = K_vec;

miter =10;
priority_max = 200;%change here

lambda1 = priority;%change this%energy saving priority %left energy is low
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
for j = 1:K%interference vector loop
for k =1:K
    %interference vec %only from the next neighbor user
    if k ~= desired_id & k == desired_id+1
        
        interf_vec(desired_id,1) = power_vec(k)*1*sym_dur_vec(k)...
            + interf_vec(desired_id,1);
        sumsym_dur_vec(desired_id,1)= interf_vec(desired_id,1)...
            + sumsym_dur_vec(desired_id,1);
    end      
end
 desired_id = desired_id+1;
end
       
%% optimization problem
nbiter = 10;
noisepower   = 0.1;
if(K>1)
opt_decision_uk = ones(K,1);

for j = 1:3%avoid null decision_uk loop
for m = 1:nbiter%lambda converge until loop
    
decision_uk = ones(K,1);%initialize uk

grad_lam = -sum(decision_uk)...
        +sum(decision_uk.^2);
diff  = -learn_rate*grad_lam;

    if (abs(diff) <= tolerance)
        convergedlam  = true;
        kk=1;
        cvx_begin quiet
    
            variable decision_uk(K,1) 
            
            %objective
            minimize(-decision_uk'*K_vec -lambda1*sum(decision_uk)...
                +lambda1*sum(decision_uk.^2))

            %constraints
            subject to

            decision_uk <= ones(K,1) 
            zeros(K,1)  <= decision_uk 
   
            for kk = 1: K-1
                decision_uk(kk+1,1)<=decision_uk(kk,1)
            end

            sum(decision_uk)>=1
            1<= sum(decision_uk)<= K
            decision_uk.*((noisepower^2 + interf_vec + power_vec(1:K).*mean(g_vec(1:K,:),2))...
                -power_vec(1:K).*mean(g_vec(1:K,:),2)*(1+1/sinr_th) ) <= 0
           
            0.00001*max_tx_power/timeslot <= decision_uk'*sumsym_dur_vec <= ...
               1/(priority+0.001)*priority_max/timeslot%change interference threshold limits
           %energy saving priority =lambda1 %if lambda high interference
           %energy saving priority =lambda1 %if lambda high interference
           %threshold limit that a device can handle goes low
        cvx_end
        
        %sumsym_dur_vec
        %disp('1');
        sumsym_dur_vec;
        decision_uk = round(decision_uk);
        %fprintf("K decision_uk %i %f\n",K,decision_uk);
        %fprintf("m i %f %f\n",m,i);
        if (isnan(decision_uk))
            %disp('NaN uk');
            %decision_u= ones(K,1);%re-initialize uk
            opt_decision_uk = ones(K,1);%previous answer
            break;

        else
            %stochastic gradient descent
            grad_uk = -K_vec-lambda1...
                            +2*lambda1'*decision_uk;
            %grad_uk = decision_uk'*K_vec-lambda1...
               %             +lambda1'*decision_uk;
            diffuk  = -learn_rate*(1/K)*(grad_uk);
            %fprintf("abs diff uk %f\n",1/K*abs(mean(diffuk)));

            if (1/K*abs(mean(diffuk)) <= tolerance2)
                convergeduk  = true;
                %fprintf("actual decision uk: %f\n",decision_uk);
                %fprintf("yes converge uk %i\n",m);
                opt_decision_uk = decision_uk>0.8;
                break;
            else
                decision_uk    = decision_uk + diffuk;
                if ((decision_uk) <zeros(K,1))
                     decision_uk = zeros(K,1);
                end
            end 
        end
    else
        %disp('lambda')
        lambda1    = (lambda1 + diff);
        if (lambda1 <0)
                 lambda1 = 0;
        end
    end 
    opt_decision_uk = decision_uk>0.8;
    if (isnan(opt_decision_uk))
        opt_decision_uk = ones(K,1);
    end
end%end lambda converge
    %opt_decision_uk = decision_uk>0.8;
    %disp(decision_uk);
end%end null uk
K = K-sum(opt_decision_uk);%update K
if K<=1 
    alldatadecoded=true;
    
    disp('break');
    %break;
end%end if 
end

%% throughput of each user
%considering synchronous uplink noma
E_max = 10;

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
end%end while

avgenergy_eff(i) = abs(mean(energy_eff));
avgenergy_effconv(i) = abs(mean(energy_eff_conv));

avgcomplexity_prop(i) = mean(sic_complextiyprop);
avgcomplexity_conv(i) = mean(sic_complextiyconv);

end
a = abs(mean(energy_eff_conv));
b = abs(mean(energy_eff));

c = mean(sic_complextiyconv);
d = mean(sic_complextiyprop);

fprintf("nbusers %i\n",nbusers);
fprintf("avg energy eff proposed %f\n",mean(energy_eff));
fprintf("avg energy eff conv %f\n",mean(energy_eff_conv));
fprintf("complexity conv %f\n",mean(sic_complextiyconv));
fprintf("complexity prop %f\n",mean(sic_complextiyprop));
end
end