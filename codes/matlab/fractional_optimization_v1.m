clc; clear variables; 
close all;

% the sic triangle optimization 
% by selecting the optimal number of user messages to decode
% considering the energy efficiency objective function

%% initialization of parameters
%number of users 
K = 4;

%use stochastic gradient descent to converge
nbiter       = 1000;%maximum number of iterations to converge

%initialize the optimizing variables
lam_opt(1)  = 0;%lambda 
up_opt(1)   = 0;%upsilon
beta(1)     = 2;%dinkelbach parameter
user_opt(1) = 0;%optimizing variable

E_circuit   = 5.0; %fixed circuit power consumption
B           = 10^6;%channel bandwidth
E_dec       = 0.1;% decoding energy per symbol

ep_userk    = 1; %epsilon 1
ep_lamk     = 1;%epsilon 2
ep_upk      = 1;%epsilon 2

grad_u      = zeros(K,nbiter);%gradient of lagrangian function
t           = 0.4; %fixed learning rate

%% environment 
k1=5; %Rician factor %ref: https://www.researchgate.net/publication/263669548_Probability_Distribution_of_Rician_K-Factor_in_Urban_Suburban_and_Rural_Areas_Using_Real-World_Captured_Data
mean=sqrt(k1/(k1+1));% mean
sigma=sqrt(1/(2*(k1+1)));%variance

N=10^3;  % Number of Bits
d4 = 10; d3 = 9; d2 = 4; d1 = 3;    %Distances of users from rx

eta = 4;            % Path loss exponent

noisepower   = 1.0;
max_tx_power = 100;

%transmission power of each user
transmitpow_k = max_tx_power*abs(randn(K,1));% unsorted transmit power vector

tx_pow_k = sort(transmitpow_k,'descend'); %sorted transmit power vector %descending 

p1 =tx_pow_k(1);%nearest user
p2 =tx_pow_k(2);
p3 =tx_pow_k(3);
p4 =tx_pow_k(4);%farthest user

%channel coefficients of each user 
h1 = sqrt(d1^-eta)*sqrt(p1/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
h2 = sqrt(d2^-eta)*sqrt(p2/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
h3 = sqrt(d3^-eta)*sqrt(p3/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
h4 = sqrt(d4^-eta)*sqrt(p4/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);

%channel gains
g1 = (abs(h1)).^2;
g2 = (abs(h2)).^2;
g3 = (abs(h3)).^2;
g4 = (abs(h4)).^2;

%% throughput of each user
%considering synchronous uplink noma
SNR = 20:2:40;
snr = db2pow(SNR);

%add the symbol duration????????????????
C1 = B*(log2(1 +g1*p1./(g2.*p2+g3.*p3+g4.*p4+1)));
C2 = B*(log2(1 +g2.*p2./(g3.*p3+g4.*p4+1)));
C3 = B*(log2(1 +g3.*p3./(g4.*p4+1)));
C4 = B*(log2(1 +g4.*p4));

% lambda values 
l1 = 1;l2 = 2;l3 = 3;l4 = 4;
v1 = 1;v2 = 2;v3 = 3;v4 = 4;

%% grad uk

%calculate the grad of lagrangian objective function
%optimizing variable 

U = triu( ones(K,K) );

rate_vec = [C1;C2;C3;C4];
up_vec = [1;2;3;4];
rate_th = B*(log2(1 +g4.*p4));
C = E_dec*(sum(U));
lambda = [l1;l2;l3;l4];

%grad function of u %optimizing variable
grad_u(:,1)   = sum(rate_vec,2)- beta*C'-...
                lambda+up_vec'*(sum(rate_th - rate_vec,2)) ;

%grad function of upsilon and lambda
grad_up(:,1)  = sum(rate_th)*ones(K,1); 
grad_lm(:,1)  = zeros(K,1); 

%% iterations
%consider the grad of the objective function 

n = 1;
u(:,n) =ones(K,1);

for iter =1:nbiter
    disp(n);
    %define u vector and give an initial value to it
    u(:,n+1) = u(:,n) + t*grad_u(n);
    %u(:,n+1)  = u(:,n+1) >0;%convert to binary
    if(u(:,n+1)==u(:,n))
        disp('yes');
        disp(u(:,n+1));
        disp(u(:,n));
        break;
    end
    n = n+1;
end

%first determine the lamda and mew optimal
%use such to determine the u optimal

%check convergence

%% optimization of the binary vector 
%converting it into a positive semidefinite matrix

%define the constraints


%binary vector is not positive semidefinite
%create a posi defnite dme
%{
u = randn(4,1)
disp(u)
U = u*u';
disp(U)
disp(diag(U));%should be a binary vector
u'*sum_rate
%}
%optimizing variable: diag of U

%{
%while (uk(n-1) ~= u_k(n))%run until convergence %approx equal
    for k = 1: K%calculate throughput and energy consumption for each user
        %store in SINR_k and Ec_k in vectors
        RecSignalpower_k = channelgain_k(k).^2*transmitPower_k(k) ; 
        Inteferencepower_k = ;
        SINR_k(k) = RecSignalpower_k/(Inteferencepower_k+noisepower); 

        A_k(k) = B*log(1+SINR_k(k));%per user throughput
        C_k(k) = E_dec*sum(sum(U));%per user energy consumption
    end
    user_opt(n) = user_opt(n-1) + t*grad_uk(n-1);%stochastic gradient descent to predict the current user_opt
   
    %calculate the grad of lagrangian objective function
    grad_uk(n) = sum (A_k)-beta*E_circuit-sum(lamda_opt)+...
    sum(upsilon_opt.*(rate_th - A_k));%store value of grad uk 
    
    for i = 1: snr %for each snr
       m = 1;%initialize n
       for j = 1: nbiter
        
           
        
       end
    end 
    
end
%}
