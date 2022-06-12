%date: 10/05/2022
%version 7 of optimization problem
%maximize the number of users
%subject to sinr and energy threshold constraints
%constraints to constrict the optimal variable to a binary variable
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
nbiter      = 50;%maximum number of iterations to converge

%initialize the optimizing variables
lam_opt(1)  = 0;%lambda 
up_opt(1)   = 0;%upsilon
beta(1)     = 2;%dinkelbach parameter
user_opt(1) = 0;%optimizing variable

E_circuit   = 5.0; %fixed circuit power consumption
B           = 10^6;%channel bandwidth
E_dec       = 0.1;% decoding energy per symbol

%% environmnet
N=10^3;  % Number of Bits
d4 = 10; d3 = 9; d2 = 4; d1 = 3;    %Distances of users from rx
communication_radius = d4;
eta1 = 1;            % Path loss exponent
eta2 = 2;
eta3 = 3;
eta4 = 7;
etath = 4;

noisepower   = 0.1;
max_tx_power = 3;
%p1 = max_tx_power ;
%p2 = max_tx_power ;
%p3 = max_tx_power ;
%p4 = max_tx_power ;

%transmission power of each user
transmitpow_k = max_tx_power*abs(randn(K,1));% unsorted transmit power vector

tx_pow_k = sort(transmitpow_k,'descend'); %sorted transmit power vector %descending 

p1 =tx_pow_k(1);%nearest user
p2 =tx_pow_k(2);
p3 =tx_pow_k(3);
p4 =tx_pow_k(4);%farthest user

%channel coefficients of each user 
h1 = sqrt(d1^-eta1)*sqrt(p1/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);%near user
h2 = sqrt(d2^-eta2)*sqrt(p2/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);%second near user
h3 = sqrt(d3^-eta3)*sqrt(p3/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);%third near user
h4 = sqrt(d4^-eta4)*sqrt(p4/2)*(randn(1,N)+1i*randn(1,N))/sqrt(2);%far user

h_th = sqrt(communication_radius^-etath)*sqrt(max_tx_power/2)*(randn(1,N)+...
1i*randn(1,N))/sqrt(2);

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

%symbols
s1   = 0.6;
s2   = 0.5;
s3   = 0.2;
s4   = 0.1;

%vectors
interf_vec     = zeros(K,1);
factorialk_vec = zeros(K,1);
sumsym_dur_vec = zeros(K,1);
sym_dur_vec = [s1;s2;s3;s4];
power_vec   = [p1;p2;p3;p4];
g_vec       = [g1;g2;g3;g4];%channel gain vector

rate_th     = log2( 1 + sqrt(max_tx_power/2)*g_th/noisepower);
desired_id  = 1;
eth         = 10;
timeslot    = 1;

%factorial k vector
for k = 1:K
    factorialk_vec(k) = K-(k-1);
    k = k-1;
end

for j = 1:K
for k =1:K%interference %only from the next neighbor user
    if k ~= desired_id & k == desired_id+1
        interf_vec(desired_id,1) = power_vec(k)*1*sym_dur_vec(k) + interf_vec(desired_id,1);
        sumsym_dur_vec(desired_id,1)= interf_vec(desired_id,1)  + sumsym_dur_vec(desired_id,1);
    end   
   
end
 desired_id =desired_id+1;
end

max = nbiter+1;

k0 = 0.4;
uk    = zeros(K,max);
uvec  = zeros(K,max);
lam1k = k0*ones(K,max);
lam2k = k0*ones(K,max);
lam3k = k0*ones(K,max);
lam4k = k0*ones(K,max);
lam5k = k0*ones(K,max);

k0 = 1;
grad_uk = k0*ones(K,max);
grad_lam1k = k0*ones(K,max);
grad_lam2k = k0*ones(K,max);
grad_lam3k = k0*ones(K,max);
grad_lam4k = k0*ones(K,max);
grad_lam5k = k0*ones(K,max);

%stop criterian %check convergence
tolerance    = 1;
learn_rate   = 0.1; %fixed learning rate
euk   = tolerance ;%uk
elam1k= tolerance ;
elam2k= tolerance ;
elam3k= tolerance ;
elam4k= tolerance ;
elam5k= tolerance ;

%% iterations
%consider the grad of the objective function 
%check for different snr values

for j = 1:20%number of random iterations
   %m = 1; 
   convergeduk = false;
   ukconverged = false;
   converged1  = false;
   converged2  = false;
   converged3  = false;
   converged4  = false;
   converged5  = false;
   
   for m = 1:nbiter %until lambda converge
        %% uk converge using stockastic gradient descent
        for n= 1:nbiter%until uk converge
            uvec(:,n) = uk(:,n);
            %uvec(1,n) = 1-uk(1,n);
            for k = 2: K
                if (k<K)
                    uvec(k,n) =  (uk(k,n)-uk(k+1,n)) ;
                else
                    uvec(k,n) = -uk(k,n)+1;%last element+1 to prevent...
                    %it satifying when its zero
                end
            end
            
            grad_uk(:,n+1) = (factorialk_vec.*ones(K,1) -...
                lam1k(:,n) - lam2k(:,n) ...
                  -2*lam2k(:,n).*uk(:,n) ...
                  -lam3k(:,n).*(uvec(:,n)) ...
                  -lam4k(:,n).*(max_tx_power*mean(g_th,2)*...
                  (noisepower^2+(interf_vec))... 
                  -power_vec.*mean(g_vec,2)*noisepower^2)...
                  -lam5k(:,n).*sumsym_dur_vec);

            diffuk  = learn_rate*grad_uk(:,n+1);

            if (mean(abs(diffuk)) <= tolerance)
                convergeduk = true; 
                uk(:,n+1)    = (uk(:,n) + diffuk)>0;
            else
                uk(:,n+1)    = (uk(:,n) + diffuk)>0;
            end  
 
        end
          
        %% lambda converge using stochastic gradient descent
        %check convergence of lamda and if not update using gradient descent
        grad_lam1k(:,m) = - (uk(:,m)-1);
        diff1  = learn_rate*grad_lam1k;
        if (abs(diff1(m)) <= tolerance)
            converged1 = true;
        else
            lam1k(:,m+1)    = lam1k(:,m) + diff1(m);
        end   

        grad_lam2k(:,m) = - (uk(:,m)-uk(:,m).^2);
        diff2  = learn_rate*grad_lam2k;
        if (abs(diff2(m)) <= tolerance)
            converged2 = true;
        else
            lam2k(:,m+1)    = lam2k(:,m) + diff2(m);
        end  

        grad_lam3k(:,m) =  -uvec(:,n);
        diff3  = learn_rate*grad_lam3k;
        if (abs(diff3(m)) <= tolerance)
            converged3 = true;
        else
            lam3k(:,m+1)    = lam3k(:,m) + diff3(m);
        end  

        grad_lam4k(:,m) =  -(uk(:,m).*ones(K,1).*max_tx_power...
            *mean(g_th,2).*(noisepower^2+(interf_vec)) ...
            - power_vec.*mean(g_vec,2)*noisepower^2);
        diff4  = learn_rate*grad_lam4k;
        if (abs(diff4(m)) <= tolerance)
            converged4 = true;
        else
            lam4k(:,m+1)    = lam4k(:,m) + diff4(m);
        end  

        grad_lam5k(:,m) = - (uk(:,m).*sumsym_dur_vec - (eth/timeslot).*ones(K,1));
        diff5  = learn_rate*grad_lam5k;
        if (abs(diff5(m)) <= tolerance)
            converged5 = true;
        else
            lam5k(:,m+1)    = lam5k(:,m) + diff5(m);
        end  

        if (m>9 & convergeduk == true & converged1==true & converged2==true & converged3==true & ...
                    converged4 == true & converged5 == true)
                fprintf('iter: %i\n',m);
                disp('stop yes')
                fprintf('uk: %g\n',uk(:,n));
                break;
        end
   end 
   if (m ==1000 & convergeduk == false)
         fprintf('do not converge\n');
         fprintf('uk: %g\n',uk(:,n));
         disp(uk);
   end
end
 
%check uvec again
%check the constraints again