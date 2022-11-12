%v16
%last update: 14 june 2022

% %energy efficiency of NOMA asynchronous D2D SIC decoding
%Output: Energy efficiency based ...
%on number of users in ...
%proposed optimized sic traingle decoding method
%testing complity results: ber vs received power ratio
clc;
clear all;
close all;

%% input data: environmnet

%--------------------------------------------------------------------------
%%scalars
%number of users 
alldatadecoded = false;
receive_pow_ratio_vec = linspace(1,10,10);%change here
mpriority = 20;

z = zeros(mpriority,1);
zz = zeros(mpriority,1);
e = zeros(mpriority,1);
f = zeros(mpriority,1);
g = zeros(mpriority,1);
h = zeros(mpriority,1);
i = zeros(mpriority,1);
j = zeros(mpriority,1);

strongconv = zeros(mpriority,1);
strongprop = zeros(mpriority,1);
weakconv = zeros(mpriority,1);
weakprop = zeros(mpriority,1);
intermconv = zeros(mpriority,1);
intermprop = zeros(mpriority,1);

%% 
N=10^4;   % Number of Bits

communication_radius = 30;%change this 
max_dist     = 100;%meters
max_eta      = 15;
etath        = 4;%change this 
noisepower   = 0.1;
max_tx_power = 0.2;%change this
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
K = 10;%number of superimposed data
transmit_snrdb_vec  =linspace(10,20,3);
timeoffset_vec = 0.01*linspace(1,15,mpriority);

for indx = 1:length(transmit_snrdb_vec)
    
initialK = 3;
K = 3;%number of superimposed data
transmit_snrdb = transmit_snrdb_vec(indx);
receive_pow_ratio = 5;%????????
time_offset = 0.15;
pr_vec = [0.5;1;1.5;2;2.5;3;3.5;4;4.5;5;5.5;6;6.5;7.5;8;8.5;10;12;15;20];

[strongconv(indx),strongprop(indx),weakconv(indx),...
    weakprop(indx),intermconv(indx),intermprop(indx)] =...
    seqsic(initialK,alldatadecoded,K,...
    pr_vec(2),N,receive_pow_ratio_vec,receive_pow_ratio,time_offset,transmit_snrdb);

strongconv;
strongprop;
intermconv;
intermprop;
weakconv;
weakprop;
end 
 
save strongconv.mat;
save strongprop.mat;
save weakconv.mat;
save weakprop.mat;
save weakconv.mat;
save weakprop.mat;

function [strongconv,strongprop,weakconv,weakprop,intermconv,intermprop]...
    = seqsic(initialK,alldatadecoded,K,priority,...
N,receive_pow_ratio_vec,receive_pow_ratio,time_offset,transmit_snrdb)

sim_delay_prop = 0;
sim_delay_conv = 0;
nbiterations = 1;
nbrandom_iterations = 3;

for receive_pow_ratioi = 1: 1%number of superimposed data loop
alldatadecoded = false;
fprintf('transmit_snrdb %i\n',transmit_snrdb);

for i = 1:nbrandom_iterations %random iterations 
v =1;
while (alldatadecoded == false & K> 1) 
v =1;
K = initialK;
max_dist     = 100;%meters
max_eta      = 10;
etath        = 4;%change this 
noisepower   = 0.1;
max_tx_power = 0.2;%change this
B            = 1;%channel bandwidth
timeslot     = 1;

%Distances of users from rx
dist_k = max_dist*abs(randn(K,1));

%sorted transmit power vector %descending 
dist_vec = sort(dist_k,'ascend'); 

% Path loss exponent
eta_k   = max_eta*abs(randn(K,1));
eta_vec = sort(eta_k,'ascend'); 

%unsorted transmit power vector
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
save power_vec.mat;
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
%nsymbols vector of each user: K vec #loop
clear K_vec;
for k = 1:K
    K_vec(k,1) = K-(k-1);
end

initialK_vec = K_vec;
miter = 10;
priority_max = 120;
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
delta_mat(1:K,:) = 0.5*rand(K,K);%B1,... Bn, C1....,Cn, ....... %Z1,....Zn

reverse_delta_mat(K,:)  = zeros(K,1);
reverse_delta_mat(1:K,:) = 0.5*rand(K,K);%A1, A2

for j = 1:K%interference vector loop
for k = 1:K
    %interference vec %only from the next neighbor user
    if k ~= desired_id & k == desired_id+1 & desired_id == 1        
        sumsym_dur_vec(desired_id,1) = sum(delta_mat(desired_id+1,:))
        
    elseif k ~= desired_id & k == desired_id+1 & desired_id > 1 & desired_id <K
        sumsym_dur_vec(desired_id,1) = sum(delta_mat(desired_id+1,:))+...
            sum(delta_mat(desired_id-1,:))
    elseif desired_id==K
        sumsym_dur_vec(desired_id,1) = sum(delta_mat(desired_id-1,:)); 
       
    end
 desired_id = desired_id+1;
end
end
sumsym_dur_vec = tril(delta_mat);

desired_id = 1;
for j = 1:K%interference vector loop
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
interf_vec = tril(delta_mat);
%% optimization problem
pastdelay = 0;
nbiter = 10;
noisepower = 0.1;

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
           
      var3: decision_uk.*((noisepower^2 + sum(interf_vec,2) + power_vec(1:K).*mean(g_vec(1:K,:),2))...
                -power_vec(1:K).*mean(g_vec(1:K,:),2)*(1+1/sinr_th)) <= 0  
      var4: decision_uk(1:K-1)>= decision_uk(2:K)     
cvx_end

echo off
diary off;
decision_uk = round(decision_uk);

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

mod_order = 4;
timeoff_min = 0.01;
timeoff_max = 0.5;

if K==initialK
    [berfinalconv] = berfunc(power_vec, noisepower, nsym, user_strength,...
    delta_mat,reverse_delta_mat,mod_order,timeoff_min,timeoff_max,sumsym_dur_vec);
end

if K>=2
[berfinal0] = berfunc(power_vec, noisepower, nsym, user_strength,...
    delta_mat,reverse_delta_mat,mod_order,timeoff_min,timeoff_max,sumsym_dur_vec);
end

ber_propuser1(v) = berfinal0(1);
ber_convuser1(v) = berfinalconv(1);
ber_propuserw(v) = berfinal0(K);
ber_convuserw(v) = berfinalconv(K);
ber_propuseri(v) = berfinal0(round(K/2)+1);
ber_convuseri(v) = berfinalconv(round(K/2)+1);

K = K-sum(decision_uk);%update K
end%end if 

v = v+1;%all user decoding index v

end %end while

avgberconv(i) = mean(ber_convuser1);
avgberprop(i) = mean(ber_propuser1);

avgberweakprop(i) = mean(ber_propuserw);
avgberweakconv(i) = mean(ber_convuserw);

avgberinterprop(i) = mean(ber_propuseri);
avgberinterconv(i) = mean(ber_convuseri);

end

end

strongconv = mean(avgberconv); 
strongprop = mean(avgberprop); 

weakconv = mean(avgberweakconv);
weakprop = mean(avgberweakprop);

intermconv = mean(avgberinterconv); 
intermprop = mean(avgberinterprop); 

end


function [berfinal] = berfunc(power_vec, noise, nsym, user_strength,...
    delta_mat,reverse_delta_mat,mod_order,timeoff_min,timeoff_max,sumsym_dur_vec)

for i = 1: length(user_strength)
for j = 1: nsym(i)
if user_strength(i) == 1 && j == 1%A1 strongest user decoding
    %interf by B1, s+1
    delta_i = sum(sumsym_dur_vec(1,:),2);%known
    p_d     = power_vec(i); %desired power
    p_iw    = power_vec(i+1);%interferes power 
    %disp('yea')
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
    (delta_i.*p_iw/2+noise))))).^2;
    %fprintf('length(user_strength) %f\n',length(user_strength))  
    q   = integral(fun,timeoff_min,timeoff_max);
        
    p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
    p_bita1    = p_err_sym1/log(mod_order);
    ber_vec(i,1) = p_bita1;

elseif (user_strength(i) == 1 && j > 1)%A2.... An
    %disp('no')
    %interf->B1 and B2 %s-1 and s+1
    delta_i = sum(sumsym_dur_vec(user_strength(i),:),2);
    p_d     = power_vec(i); %desired power
    p_iw     = power_vec(i+1);%interferes power
    power_v  = [p_iw/2;p_iw/2];
    fun = @(delta_i) (1 - qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v)+noise))))).^2;
    q   = 1 - integral(fun,timeoff_min,timeoff_max);
    p_err_sym2 = (1/(timeoff_max -timeoff_min))*mean(q*sum(delta_i));
    p_bita2 = p_err_sym2/log(mod_order); 
    ber_vec(i,j) = p_bita2;
    
elseif (user_strength(i) > 1 && j == 1 && i~=length(user_strength))%B1%intermediate
    delta_i = sum(sumsym_dur_vec(user_strength(i),1),2);
    detected_vec = [(6/(mod_order-1));(6/(mod_order-1));1;];
    p_d     = power_vec(i); %desired power
    p_is    = power_vec(i-1);%interferes power
    p_iw    = power_vec(i+1);%interferes power
    power_v  = [p_is/2;p_is/2;p_iw/2];
    
    %change here %error vec1
    error_vec1 = [1;1;1];%one error one correct
    %include the delta i into the func expression???????
    %change here
    fun1 = @(delta_i) (1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v.*detected_vec.*(error_vec1))+noise))))).^2;
    q1   = 1- integral(fun1,timeoff_min,timeoff_max);
    p_interb1 = (1/(timeoff_max -timeoff_min))*mean(q1*sum(delta_i.*error_vec1));
    p_err_sym_ec1 = p_interb1*(ber_vec(i-1,j))*(ber_vec(i-1,j+1))*log(mod_order);
    
    %change here %error vec1
    error_vec1 = [1;0;1];%one error one correct
    %include the delta i into the func expression??????)
    fun1 = @(delta_i) (1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v.*detected_vec.*(error_vec1))+noise))))).^2;
    q2   = integral(fun1,timeoff_min,timeoff_max);
    p_interb2 = (1/(timeoff_max -timeoff_min))*mean(q2*sum(delta_i.*error_vec1));
    p_err_sym_ec2 = p_interb2*(ber_vec(i-1,j))*(1-ber_vec(i-1,j+1))*log(mod_order);
    
    %change here %error vec2
    error_vec1 = [0;1;1];%one error one correct
    %include the delta i into the func expression???????
    fun1 = @(delta_i) (1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v.*detected_vec.*(error_vec1))+noise))))).^2;
    q3   = 1-integral(fun1,timeoff_min,timeoff_max);
    p_interb3 = (1/(timeoff_max -timeoff_min))*mean(q3*sum(delta_i.*error_vec1));
    p_err_sym_ec3 = p_interb3*(1-ber_vec(i-1,j))*(ber_vec(i-1,j+1))*log(mod_order);
    
    %change here %error vec3
    error_vec1 = [0;0;1];%one error one correct
    %include the delta i into the func expression???????
    fun1 = @(delta_i) (1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v.*detected_vec.*(error_vec1))+noise))))).^2;
    q4   = 1-integral(fun1,timeoff_min,timeoff_max);
    p_interb4 = (1/(timeoff_max -timeoff_min))*mean(q4*sum(delta_i.*error_vec1));
    p_err_sym_ec4 = p_interb4*(1-ber_vec(i-1,j))*(1-ber_vec(i-1,j+1))*log(mod_order);
    
    avgperr = 1/4*(p_err_sym_ec1 +p_err_sym_ec2 +p_err_sym_ec3 +p_err_sym_ec4);

    p_bitb1 = avgperr/log(mod_order); 
    ber_vec(i,j) = p_bitb1;
    
elseif (user_strength(i) > 1 && j > 1 && i~=round(length(user_strength)))%B2.....Bn & not last user 
    %disp('here1 ')
    %inter -> A2, A3 and C1, C2
    nbinterf = 4;
    delta_i = sum(sumsym_dur_vec(user_strength(i),:),2);;
   
    detected_vec = [(6/(mod_order-1));(6/(mod_order-1));1;1;];
    p_d     = power_vec(i); %desired power
    p_is    = power_vec(i-1);%interferes power
    p_iw    = power_vec(i+1);%interferes power
    power_v  = [p_is/2;p_is/2;p_iw/2;p_iw/2];
    
    %change here %error vec1
    error_vec1 = [1;1;1;1];%one error one correct
    %include the delta i into the func expression???????
    fun1 = @(delta_i) (1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v.*detected_vec.*(error_vec1))+noise))))).^2;
    q1   = integral(fun1,timeoff_min,timeoff_max);
    p_interb1 = (1/(timeoff_max -timeoff_min))*mean(q1*sum(delta_i.*error_vec1));
    p_err_sym_ec1 = p_interb1*(ber_vec(i-1,j))*(ber_vec(i-1,j+1))*log(mod_order);
    
    %change here %error vec1
    error_vec2 = [1;0;1;1];%one error one correct
    %include the delta i into the func expression???????
    fun1 = @(delta_i) (1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v.*detected_vec.*(error_vec1))+noise))))).^2;
    q2   = integral(fun1,timeoff_min,timeoff_max);
    p_interb2 = (1/(timeoff_max -timeoff_min))*mean(q2*sum(delta_i.*error_vec2));
    p_err_sym_ec2 = p_interb2*(ber_vec(i-1,j))*(1-ber_vec(i-1,j+1))*log(mod_order);
    
    %change here %error vec1
    error_vec3 = [0;1;1;1];%one error one correct
    %include the delta i into the func expression???????
    fun1 = @(delta_i) (1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v.*detected_vec.*(error_vec1))+noise))))).^2;
    q3   = integral(fun1,timeoff_min,timeoff_max);
    p_interb3 = (1/(timeoff_max -timeoff_min))*mean(q3*sum(delta_i.*error_vec3));
    p_err_sym_ec3 = p_interb3*(1-ber_vec(i-1,j))*(ber_vec(i-1,j+1))*log(mod_order);
    
    %change here %error vec1
    error_vec4 = [0;0;1;1];%one error one correct
    %include the delta i into the func expression???????
    fun1 = @(delta_i) (1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v.*detected_vec.*(error_vec1))+noise))))).^2;
    q4   = integral(fun1,timeoff_min,timeoff_max);
    p_interb4 = (1/(timeoff_max -timeoff_min))*mean(q4*sum(delta_i.*error_vec4));
    p_err_sym_ec4 = p_interb4*(1-round(ber_vec(i-1,j)))*log(mod_order)*(1-round(ber_vec(i-1,j+1)))*log(mod_order);
    %p_err_sym_ec4 = p_err_sym_ec3;
    avgperr2 = 1/4*(p_err_sym_ec1 +p_err_sym_ec2 +p_err_sym_ec3 +p_err_sym_ec4);

    p_bitb2 = avgperr2/log(mod_order); 
    ber_vec(i,j) = p_bitb2;
    
elseif (user_strength(i) > 1 && j == 1 && i==round(length(user_strength)))%weakest user
    %disp('here2 ')
    %interf-> B1 and B2
    delta_i =[(delta_mat(user_strength(i-1),j)); ...
        0.5-delta_mat(user_strength(i-1),j+1)];
    detected_vec = [(6/(mod_order-1));(6/(mod_order-1))];
    p_d     = power_vec(i); %desired power
    p_is    = power_vec(i-1);%interferes power
    power_v  = [p_is/2;p_is/2;];
    
    error_vec1 = [1;1;];%one error one correct
    %include the delta i into the func expression???????
    fun1 = @(delta_i) (1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v.*detected_vec.*(error_vec1))+noise))))).^2;
    q1   = integral(fun1,timeoff_min,timeoff_max);
    p_interb1 = (1/(timeoff_max -timeoff_min))*mean(q1*sum(delta_i.*error_vec1));
    p_err_sym_ec1 = p_interb1*(ber_vec(i-1,j))*(ber_vec(i-1,j+1))*...
        log(mod_order);
    
    error_vec2 = [1;0;];%one error one correct
    %include the delta i into the func expression???????
    fun1 = @(delta_i) (1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v.*detected_vec.*(error_vec2))+noise))))).^2;
    q2   = integral(fun1,timeoff_min,timeoff_max);
    p_interb2 = (1/(timeoff_max -timeoff_min))*mean(q2*sum(delta_i.*error_vec2));
    p_err_sym_ec2 = p_interb2*(ber_vec(i-1,j))*(1-ber_vec(i-1,j+1))*...
        log(mod_order);
    
    error_vec3 = [0;1;];%one error one correct
    %include the delta i into the func expression???????
    fun1 = @(delta_i) (1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v.*detected_vec.*(error_vec3))+noise))))).^2;
    q3   = integral(fun1,timeoff_min,timeoff_max);
    p_interb3 = (1/(timeoff_max -timeoff_min))*mean(q3*sum(delta_i.*error_vec3));
    p_err_sym_ec3 = p_interb3*(1-ber_vec(i-1,j))*(ber_vec(i-1,j+1))*...
        log(mod_order);
   
    error_vec4 = [0;0;];%one error one correct
    %include the delta i into the func expression???????
    fun1 = @(delta_i) (1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v.*detected_vec.*(error_vec4))+noise))))).^2;
    q4   = integral(fun1,timeoff_min,timeoff_max);
    p_interb4 = (1/(timeoff_max -timeoff_min))*mean(q4*sum(delta_i.*error_vec4));
    p_err_sym_ec4 = p_interb4*(1-ber_vec(i-1,j))*(1-ber_vec(i-1,j+1))*...
        log(mod_order);
    %p_err_sym_ec4 = p_err_sym_ec3
    avgperr2 = 1/4*(p_err_sym_ec1 +p_err_sym_ec2 +p_err_sym_ec3 +p_err_sym_ec4);

    p_bitb2 = avgperr2/log(mod_order); 
    ber_vec(i,j) = p_bitb2;
   
else
    disp('noo')
end
end
end
berfinal = sum(ber_vec,2);
end
