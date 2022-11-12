%v16
%initialized date: 23 may 2022
%last updated: 25 may 2022
%energy efficiency of NOMA asynchronous D2D SIC decoding
%Output: Energy efficiency based ...
%on number of users in ...
%proposed optimized sic traingle decoding method
clc;
clear all;
close all;

%% input data: environmnet
%--------------------------------------------------------------------------
%%scalars
%number of users 
alldatadecoded = false;
K = 4;%number of superimposed data
initialK=K;
for nbusers = 5: 5%number of superimposed data loop
for i = 1:4 %random iterations 

while (alldatadecoded == false) 
% Number of Bits
N=10^4;  

communication_radius = 30;%change this 
max_dist     = 100;
max_eta      = 15;
etath        = 4;%change this 
noisepower   = 0.1;
max_tx_power = 10;%change this
B            = 10^6;%channel bandwidth

pth          = max_tx_power.*communication_radius^-etath;
h_th         = sqrt(communication_radius^-etath)*sqrt(pth/2)*(randn(1,N)+...
1i*randn(1,N))/sqrt(2);
g_th         = (abs(h_th)).^2;

rate_th      = log2( 1 + sqrt(pth/2)*g_th/noisepower);

desired_id   = 1;
eth          = 1;
timeslot     = 1;

%random iterations
%--------------------------------------------------------------------------

%%vectors
%Distances of users from rx
dist_k = max_dist*abs(randn(K,1));

%sorted transmit power vector %descending 
dist_vec = sort(dist_k,'ascend'); 

% Path loss exponent
eta_k = max_eta*abs(randn(K,1));
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

%nsymbols vector of each user: K vec #loop
clear K_vec;
for k = 1:K
    K_vec(k,1) = K-(k-1);
end

miter =10;
priority_max = 200;
priority = 1.5;
lambda1 = priority;%change this%energy saving priority %left energy is low
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
if(K>0)
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
            minimize(-decision_uk'*K_vec)

            %constraints
            subject to
   
            for kk = 1: K-1
                decision_uk(kk+1,1)<=decision_uk(kk,1)
            end

            sum(decision_uk)>=1
            1<= sum(decision_uk)<= K
            decision_uk.*((noisepower^2 + interf_vec + power_vec.*mean(g_vec,2))...
                -power_vec.*mean(g_vec,2)*(1+1/sinr_th) ) <= 0
           
            0.001*max_tx_power/timeslot <= decision_uk'*sumsym_dur_vec <= ...
               1/(priority+0.001)*priority_max/timeslot%change interference threshold limits
           %energy saving priority =lambda1 %if lambda high interference
           %threshold limit that a device can handle goes low
        cvx_end
        
        %sumsym_dur_vec
        sumsym_dur_vec;
        decision_uk = round(decision_uk);
        fprintf("K decision_uk %i %f\n",K,decision_uk);
        fprintf("m i %f %f\n",m,i);
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
if K==0 
    alldatadecoded=true;
    disp('break')
    break;
end%end if 
end


%% throughput of each user
%considering synchronous uplink noma
E_max = 10;

SINR_k = power_vec.*mean(g_vec,2)./(interf_vec+noisepower^2);

throughput_vec = log(1+SINR_k);

total_throughput = sum(throughput_vec);

%% energy efficiency 
%proposed optimal sic
total_energ_consump = E_max - E_max^(exp(-log(2)/1000*opt_decision_uk'*K_vec));
energy_eff(i) = total_throughput/(total_energ_consump +0.01);

%%conv sic
total_energ_consump_conv = E_max - E_max^(exp(-(log(2)/1000)*sum(K_vec)));
energy_eff_conv(i) = total_throughput/(total_energ_consump_conv+0.01);

%% complexity analysis 
L = 3;%fixed number of iterations
ProbConverged = 0.5; %depends on the BER of decoded data

%conv SIC
bigOorderconv = 0;

if(K>L)
    bigOorderconv(i)  = K^2;%conv SIC over L iterations %if K>L
else 
    bigOorderconv(i)  = K*L;
end

%optimized SIC
bigOorderprop = 0;
if(ProbConverged == 0.5)
 iterationsconverge = 1;
 bigOorderoptprob = K^3*(K^2+3);
 bigOorderprop = sum(decision_uk).^2+bigOorderoptprob;
 %only sic decoding 
 sic_complextiyprop(i) = sum(opt_decision_uk)^2*log(1/0.01);
 sic_complextiyconv(i) = sum(initialK)^2;
end
end%end while
end
fprintf("nbusers %i\n",nbusers);
fprintf("avg energy eff proposed %f\n",mean(energy_eff));
fprintf("avg energy eff conv %f\n",mean(energy_eff_conv));
fprintf("complexity conv %f\n",mean(sic_complextiyconv));
fprintf("complexity prop %f\n",mean(sic_complextiyprop));
end
