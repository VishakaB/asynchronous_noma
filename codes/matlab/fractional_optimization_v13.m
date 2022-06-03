clc;
clear all;
close all;

%number of users 
K = 5;
%% environmnet
N=10^6;  % Number of Bits
d5 =50; d4 = 10; d3 = 9; d2 = 4; d1 = 3;%Distances of users from rx
communication_radius = 30;%change this
eta1 = 0.01;            % Path loss exponent
eta2 = 0.2;
eta3 = 1;
eta4 = 5;
eta5 = 12;
etath = 4;%change this %chang

noisepower   = 5;
max_tx_power = 1;%change this
B            = 10^6;%channel bandwidth

for i = 1:1%random iterations
%transmission power of each user
transmitpow_k = max_tx_power*abs(randn(K,1));% unsorted transmit power vector

tx_pow_k = sort(transmitpow_k,'descend'); %sorted transmit power vector %descending 

p1 =tx_pow_k(1);%nearest user
p2 =tx_pow_k(2);
p3 =tx_pow_k(3);
p4 =tx_pow_k(4);%farthest user
p5 =tx_pow_k(5);

%channel coefficients of each user 
h1 = sqrt(d1^-eta1)*sqrt(p1/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);%near user
h2 = sqrt(d2^-eta2)*sqrt(p2/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);%second near user
h3 = sqrt(d3^-eta3)*sqrt(p3/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);%third near user
h4 = sqrt(d4^-eta4)*sqrt(p4/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);%far user
h5 = sqrt(d5^-eta5)*sqrt(p4/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);%far user

pth = max_tx_power.*communication_radius^-etath;
h_th = sqrt(communication_radius^-etath)*sqrt(pth/2)*(randn(1,N)+...
1i*randn(1,N))/sqrt(2);

%channel gains
g1 = (abs(h1)).^2;
g2 = (abs(h2)).^2;
g3 = (abs(h3)).^2;
g4 = (abs(h4)).^2;
g5 = (abs(h5)).^2;
g_th = (abs(h_th)).^2;

%% throughput of each user
%considering synchronous uplink noma
SNR = 20:2:40;
snr = db2pow(SNR);

%add the symbol duration
C1 = B*(log2(1 +g1*p1./(g2.*p2+g3.*p3+g4.*p4+1)));
C2 = B*(log2(1 +g2.*p2./(g3.*p3+g4.*p4+1)));
C3 = B*(log2(1 +g3.*p3./(g4.*p4+1)));
C4 = B*(log2(1 +g4.*p4));

%% optimization problem

%symbols
s1   = 0.8;%change this
s2   = 0.7;
s3   = 0.1;
s4   = 0.08;
s5   = 0.09;

%vectors
interf_vec     = zeros(K,1);
factorialk_vec = zeros(K,1);
sumsym_dur_vec = zeros(K,1);
sym_dur_vec = max_tx_power*[s1;s2;s3;s4;s5];%change this
power_vec   = [p1;p2;p3;p4;p5];
g_vec       = [g1;g2;g3;g4;g5];%channel gain vector

rate_th     = log2( 1 + sqrt(p5/2)*g_th/noisepower);
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
    if k ~= desired_id & k == desired_id+1 %channel gain is too small
        interf_vec(desired_id,1) = power_vec(k)*1*sym_dur_vec(k)...
            + interf_vec(desired_id,1);
        sumsym_dur_vec(desired_id,1)= interf_vec(desired_id,1)...
            + sumsym_dur_vec(desired_id,1);
    end      
end
 desired_id = desired_id+1;
end
       
K_vec = [5;4;3;2;1];
k=1;
energy_priority = 0.1;%change this
lambda1 = energy_priority*ones(K,1);%change this%energy saving priority
lambda2 = 0.75;%change this
lambda3 = 2.5;%change thi

miter  = 10;
nbiter = 10;
convergedlam = false;
learn_rate= 0.1;
tolerance = 0.1;
tolerance2 = 1;
%%testing constraints
y1 =pth*mean(g_th,2)*...%check this
                  (noisepower^2+(interf_vec))... 
                  -power_vec.*mean(g_vec,2)*noisepower^2; 

y2 = power_vec.*mean(g_vec,2)*noisepower^2;
y2 - y1;

%0.1*max_tx_power/timeslot

%%stochastic gradient descent to converge

for n = 1:3
     
    cvx_begin quiet

                variable decision_uk(K,1) 
                %objective
                minimize(-decision_uk'*K_vec -lambda1'*decision_uk...
                    +lambda1'*decision_uk.^2)

                %constraints
                subject to

                decision_uk <= ones(K,1) 
                zeros(K,1)  <= decision_uk 

                %remove this constraint 
                decision_uk(k+4,1)< decision_uk(k+3,1)< decision_uk(k+2,1)...
                    <decision_uk(k+1,1)< decision_uk(k,1) %selecting a range of users


                decision_uk'*(pth*mean(g_th,2)*...%check this
                              (noisepower^2+(interf_vec))... 
                              -power_vec.*mean(g_vec,2)*noisepower^2) <= 0

                0.01*max_tx_power/timeslot <= decision_uk'*sumsym_dur_vec <= ...
                   max_tx_power/timeslot%change interference threshold limits

    cvx_end

    
    %sum(decision_uk)
    %if(lambda1<energy_left )
    grad_lam = -sum(decision_uk)...
            +sum(decision_uk.^2);
    diff  = -learn_rate*grad_lam;
    %fprintf("difflam %f %i\n",diff, n);
    if (abs(diff) <= tolerance)
        convergedlam  = true;
        %fprintf("%lambda f\n",lambda1);

        for m =1:10
             %% cvx tool
             cvx_begin quiet

                variable decision_uk(K,1) 
                %objective
                minimize(-decision_uk'*K_vec -lambda1'*decision_uk...
                    +lambda1'*decision_uk.^2)

                %constraints
                subject to

                decision_uk <= ones(K,1) 
                zeros(K,1)  <= decision_uk 

                %remove this constraint 
                decision_uk(k+4,1)< decision_uk(k+3,1)< decision_uk(k+2,1)...
                    <decision_uk(k+1,1)< decision_uk(k,1) %selecting a range of users


                decision_uk'*(pth*mean(g_th,2)*...%check this
                              (noisepower^2+(interf_vec))... 
                              -power_vec.*mean(g_vec,2)*noisepower^2) <= 0

                0.01*max_tx_power/timeslot <= decision_uk'*sumsym_dur_vec <= ...
                   max_tx_power/timeslot%change interference threshold limits

            cvx_end

            grad_uk = decision_uk'*K_vec-lambda1...
                +lambda1'*decision_uk;
            diffuk  = -learn_rate*grad_uk;

            if (abs(mean(diffuk)) <= tolerance2 & convergedlam ==true)
                convergedlam  = true;
                fprintf("%f\n",decision_uk);
                fprintf("yes converge uk %i\n",m);
                break;
            end 
        end %miter end

    else
        lambda1    = abs(lambda1 + diff);

    end  

end%m iter   
        
end%uk iter    
