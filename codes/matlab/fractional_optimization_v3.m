%date: 10/05/2022
%version 7 of optimization problem
%maximize the number of users
%subject to sinr and channel gain constraints
%constraint to constrict the optimal variable to a binary variable
clc;
clear all;
close all;

%% randomization
randnumber = randi(1000);%rand seed
rng(randnumber);

%% initialization
%number of users 
K = 4;

%use stochastic gradient descent to converge
nbiter      = 3;%maximum number of iterations to converge

%initialize the optimizing variables
lam_opt(1)  = 0;%lambda 
up_opt(1)   = 0;%upsilon
beta(1)     = 2;%dinkelbach parameter
user_opt(1) = 0;%optimizing variable

E_circuit   = 5.0; %fixed circuit power consumption
B           = 10^6;%channel bandwidth
E_dec       = 0.1;% decoding energy per symbol

ep_userk    = 1;%epsilon 1
ep_lamk     = 1;%epsilon 2
ep_upk      = 1;%epsilon 2

grad_u      = zeros(K,nbiter);%gradient of lagrangian function
t           = 0.4; %fixed learning rate

%% environmnet
N=10^3;  % Number of Bits
d4 = 10; d3 = 9; d2 = 4; d1 = 3;    %Distances of users from rx
communication_radius = d4;
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

h_th = sqrt(communication_radius^-eta)*sqrt(max_tx_power/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);

%channel gains
g1 = (abs(h1)).^2;
g2 = (abs(h2)).^2;
g3 = (abs(h3)).^2;
g4 = (abs(h4)).^2;
g_th = (abs(h_th)).^2;

%% throughput of each user
%considering synchronous uplink noma
SNR = 20:2:40;
snr = db2pow(SNR);

%add the symbol duration????????????????
C1 = B*(log2(1 +g1*p1./(g2.*p2+g3.*p3+g4.*p4+1)));
C2 = B*(log2(1 +g2.*p2./(g3.*p3+g4.*p4+1)));
C3 = B*(log2(1 +g3.*p3./(g4.*p4+1)));
C4 = B*(log2(1 +g4.*p4));

%% optimization problem
n = 3;

lam1(1)     = 1;
lam2(1)     = 1;
lam3(1)     = 1;
lam4(1)     = 1;
lam5(1)     = 1;

%symbols
s1   = 0.5;
s2   = 0.2;
s3   = 0.01;
s4   = 0.005;

%vectors
interf_vec  = zeros(K,1);
factorialk_vec = zeros(K,1);
sumsym_dur_vec = zeros(K,1);
sym_dur_vec = [s1;s2;s3;s4];
power_vec   = [p1;p2;p3;p4];
g_vec       = [g1;g2;g3;g4];%channel gain vector
g_th        = communication_radius.^g_th;
rate_th     = log2( 1 + max_tx_power*g_th/noisepower);
desired_id  = 1;

%factorial k vector
for k = 1:K
    factorialk_vec(k) = K-(k-1);
    k = k-1;
end

for k =1:K%interference
    if k ~= desired_id
        interf_vec(k) = power_vec(k)*g_vec(k)*sym_dur_vec(k) + interf_vec(k);
        sumsym_dur_vec(k)= sym_dur_vec(k) + sumsym_dur_vec(k);
    end   
end

n = nbiter;
grad_uk(:,n) = zeros(K,1);

uk    = ones(K,nbiter+1);
lam1k = ones(K,nbiter+1);
lam2k = ones(K,nbiter+1);
lam3k = ones(K,nbiter+1);
lam4k = ones(K,nbiter+1);
lam5k = ones(K,nbiter+1);

grad_lam1k = ones(K,1);
grad_lam2k = ones(K,1);
grad_lam3k = ones(K,1);
grad_lam4k = ones(K,1);
grad_lam5k = ones(K,1);

%stop criterian %check convergence
euk   = 0.01;%uk
elam1k= 0.01;
elam2k= 100;
elam3k= 100;
elam4k= 100;
elam5k= 100;

%% iterations
%consider the grad of the objective function 
%check for different snr values
tic
for j = 1:1
for iter = 1:nbiter %until lambda converge
    n = iter;
    
    %define u vector and give an initial value to it
    uk(:,n+1)      = uk(:,n) + t*grad_uk(:,n);
    
    %u(:,n+1) = u(:,n+1) >0;%convert to binary
    if(n >2 & abs(sum(uk(:,n+1))-sum(uk(:,n)))< 0.01)%uk converge 
        disp('yes');
        disp(n);
        disp(uk(:,n+1)>0);
        disp(uk(:,n)>0);
        break;
    else
        %check convergence of lamda and if not update using gradient descent
        if (n >2 & lam1k(:,n+1)-lam1k(:,n)> elam1k)
            lam1k(:,n+1)  = lam1k(:,n) + t*grad_lam1k(:,n);
            grad_lam1k(:,n+1) = - sum(uk(:,n)-1);
        elseif (n==1)
            grad_lam1k(:,n+1) = - sum(uk(:,n)-1);
            lam1k(:,n+1)  = lam1k(:,n) + t*grad_lam1k(:,n);
        elseif (n>2 & lam1k(:,n+1)-lam1k(:,n)< elam1k)
            disp(n)
            disp(abs(lam1k(:,n+1)-lam1k(:,n))< elam1k)
            disp('lam1k');
            break
        end   
        
        grad_lam2k(:,n+1)= - sum(uk(:,n)-uk(:,n).^2);
        if (n >2 & lam2k(:,n+1)-lam2k(:,n)> elam2k)
            lam2k(:,n+1) = lam2k(:,n) + t*grad_lam2k(:,n);
            grad_lam2k(:,n+1)= - sum(uk(:,n)-uk(:,n).^2);
        else
            grad_lam2k(:,n+1) = grad_lam2k(:,n);
        end
        
        uvec = uk(:,n);
        lam3k(:,n+1) = lam3k(:,n) + t*grad_lam3k(:,n);
        for k = 1: k-1
            uvec(k,n) = uk(k+1,n) - uk(k,n);
        end
        
        if (n >2 & lam3k(:,n+1)-lam3k(:,n)> elam3k)
            %update here
            grad_lam3k(:,n+1)= mean(-lam3k.*uvec,2);
        else
            grad_lam3k(:,n+1) = grad_lam3k(:,n);
        end
        
        lam4k(:,n+1) = lam4k(:,n) + t*grad_lam4k(:,n);
        if (n >2 & lam4k(:,n+1)-lam4k(:,n)> elam4k)
            grad_lam4k(:,n+1)= - sum(uk(:,n).*ones(K,1).*max_tx_power...
            *mean(g_th,2).*(noisepower^2+(interf_vec)) ...
            - power_vec.*mean(g_vec,2)*noisepower^2);
        else
            grad_lam4k(:,n+1) = grad_lam4k(:,n);
        end
          
        lam5k(:,n+1) = lam5k(:,n) + t*grad_lam5k(:,n);
        if (n >2 & lam4k(:,n+1)-lam4k(:,n)> elam4k)        
            grad_lam5k(:,n+1)= - sum(uk(:,n).*sumsym_dur_vec);
        else
            grad_lam5k(:,n+1) = grad_lam5k(:,n);
        end
    end
    
    grad_uk(:,n+1) = factorialk_vec.*ones(K,1) - lam1k(:,n) - lam2k(:,n) ...
              -2*lam2k(:,n).*uk(:,n) ...
              -lam3k(:,n).*(uvec(:,n)) ...
              -lam4k(:,n).*(max_tx_power*mean(g_th,2)*(noisepower^2+(interf_vec))... 
              -power_vec.*mean(g_vec,2)*noisepower^2)...
              -lam5k(:,n).*sumsym_dur_vec;
          
    %include kkt consitions as additional constraints to narrow the uk
    %solution   
    %all lamda values must be equal zero or positive  
    %the grad of the lagrangian wrt to each variable must be equal to zero
    %the 
end
end
toc