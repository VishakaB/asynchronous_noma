
clear all; 
clc;
rng(0)
%% initial input data 
mod_order = 4;
timeoff_min = 0.01;
timeoff_max = 0.5;
max_tx_power = 1000;

K = 20;

% unsorted transmit power vector
transmitpow_k = max_tx_power*abs(randn(K,1));

%sorted transmit power vector %descending 
power_vec = sort(transmitpow_k,'descend'); 

noise = 0.1;

for k = 1:K
    nsym(k,1) = K-(k-1);
    user_strength(k,1) = k;
end

%time offsets between users 
%delta_mat: rows -> user index, columns-> symbol index %time offset with
delta_mat = zeros(K,K);
delta_mat(1,:) = zeros(K,1);
delta_mat(2:K,:) = 0.5*rand(K-1,K);%B1,... Bn, C1....,Cn, ....... %Z1,....Zn

reverse_delta_mat(K,:)  = zeros(K,1);
reverse_delta_mat(1:K-1,:) = 0.5*rand(K-1,K);%A1, A2

%only for one iteration
                                                                                     
%function for ber of k users in T-SIC method 
%input data
ber_vec = zeros(K,K);

[berfinal0] = berfunc(power_vec, noise, nsym, user_strength,...
    delta_mat,reverse_delta_mat,mod_order,timeoff_min,timeoff_max)

function [berfinal] = berfunc(power_vec, noise, nsym, user_strength,...
    delta_mat,reverse_delta_mat,mod_order,timeoff_min,timeoff_max)
for i = 1: length(user_strength)
for j = 1: nsym(i)
    
if user_strength(i) == 1 & j == 1%A1
    %interf by B1, s+1
    delta_i = delta_mat(user_strength(i+1),j);%known
    p_d     = power_vec(i); %desired power
    p_iw    = power_vec(i+1);%interferes power 
    %disp('yea')
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
    (delta_i.*p_iw/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max);
    p_err_sym1 = (1/(timeoff_max -timeoff_min))*q*delta_i;
    p_bita1    = p_err_sym1/log(mod_order);
    ber_vec(i,1) = p_bita1;
    
elseif (user_strength(i) == 1 & j > 1)%A2.... An
    %disp('no')
    %interf->B1 and B2 %s-1 and s+1
    delta_i = [(1-delta_mat(user_strength(i+1),j-1)); ...
        delta_mat(user_strength(i+1),j)];
    p_d     = power_vec(i); %desired power
    p_iw     = power_vec(i+1);%interferes power
    power_v  = [p_iw/2;p_iw/2];
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v)+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max);
    p_err_sym2 = (1/(timeoff_max -timeoff_min))*q*sum(delta_i);
    p_bita2 = p_err_sym2/log(mod_order); 
    ber_vec(i,j) = p_bita2;
    
elseif (user_strength(i) > 1 & j == 1 & i~=length(user_strength))%B1
    %inter -> A1, A2 and C1
    %disp('here0')
    delta_i = [(reverse_delta_mat(user_strength(i-1),j));...
        (1-reverse_delta_mat(user_strength(i-1),j+1)); ...
        delta_mat(user_strength(i+1),j);];
    detected_vec = [(6/(mod_order-1));(6/(mod_order-1));1;];
    p_d     = power_vec(i); %desired power
    p_is    = power_vec(i-1);%interferes power
    p_iw    = power_vec(i+1);%interferes power
    power_v  = [p_is/2;p_is/2;p_iw/2];
    %include the delta i into the func expression???????
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v)+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max);
    p_interb1 = (1/(timeoff_max -timeoff_min))*q*sum(delta_i);
    p_err_sym2 = p_interb1*(ber_vec(i-1,j)*ber_vec(i-1,j+1)+...
        (1-ber_vec(i-1,j))*ber_vec(i-1,j+1)+(ber_vec(i-1,j))*(1-ber_vec(i-1,j+1))...
        +(1-ber_vec(i-1,j))*(1-ber_vec(i-1,j+1)));
    p_bita2 = p_err_sym2/log(mod_order); 
    ber_vec(i,j) = p_bita2;
    
elseif (user_strength(i) > 1 & j > 1 & i~=length(user_strength))%B2.....Bn & not last user 
    %disp('here1 ')
    %inter -> A2, A3 and C1, C2
    nbinterf = 4;
    delta_i = zeros(nbinterf,1);
    for h = 1: nbinterf
        if h == 1
            delta_i(h) = (reverse_delta_mat(user_strength(i-1),j));
        elseif h == 2
            delta_i(h) = 1-reverse_delta_mat(user_strength(i-1),j+1);
        elseif h == 3
            delta_i(h) = (1-delta_mat(user_strength(i+1),j-1));
        else
            delta_i(h) = delta_mat(user_strength(i+1),j);
        end
    end
    detected_vec = [(6/(mod_order-1));(6/(mod_order-1));1;1;];
    p_d     = power_vec(i); %desired power
    p_is    = power_vec(i-1);%interferes power
    p_iw    = power_vec(i+1);%interferes power
    power_v  = [p_is/2;p_is/2;p_iw/2;p_iw/2];
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(detected_vec.*delta_i.*power_v) +noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max);
    p_interb2 = (1/(timeoff_max - timeoff_min))*q*sum(delta_i);
    p_err_sym2 = p_interb2*(ber_vec(i-1,j)*ber_vec(i-1,j+1)+...
        (1-ber_vec(i-1,j))*ber_vec(i-1,j+1)+(ber_vec(i-1,j))*...
        (1-ber_vec(i-1,j+1))...
        +(1-ber_vec(i-1,j))*(1-ber_vec(i-1,j+1)));
    p_bita3 = p_err_sym2/log(mod_order);
    ber_vec(i,j) = p_bita3;
    
elseif (user_strength(i) > 1 & j >= 1 & i==length(user_strength))
    %disp('here2 ')
    %interf-> B1 and B2
    delta_i =[(reverse_delta_mat(user_strength(i-1),j)); ...
        1-reverse_delta_mat(user_strength(i-1),j+1)];
    detected_vec = [(6/(mod_order-1));(6/(mod_order-1))];
    p_d     = power_vec(i); %desired power
    p_is    = power_vec(i-1);%interferes power
    power_v  = [p_is/2;p_is/2;];
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(detected_vec.*delta_i.*power_v)+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max);
    p_interb3 = (1/(timeoff_max -timeoff_min))*q*sum(delta_i);
    p_err_sym3 = p_interb3*(ber_vec(i-1,j)*ber_vec(i-1,j+1)+...
        (1-ber_vec(i-1,j))*ber_vec(i-1,j+1)+(ber_vec(i-1,j))*...
        (1-ber_vec(i-1,j+1))...
        +(1-ber_vec(i-1,j))*(1-ber_vec(i-1,j+1)));
    p_bita4 = p_err_sym3/log(mod_order);
    ber_vec(i,j) = p_bita4;
else
    %disp('oops')
end
end
end
berfinal = sum(ber_vec,2);
end
