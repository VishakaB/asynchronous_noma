
clear all; 
clc;

%% initial input data 

mod_order = 4;
timeoff_min = 0.01;
timeoff_max = 0.5;
p_a = 1000; %40 dB
p_b = 500;
p_c = 100;
noise = 0.1;
nsym = [3,2,1];
user_strength = [1,2,3];
%time offsets between users 
%A1->B2
delta_a1_b2 = 0.1;
%A2->B3
delta_a2_b3 = 0.2;

%delta_mat: rows -> user index, columns-> symbol index %time offset with
%other users from left to right
delta_mat = [0,0,0;0.5,0.2,0.4;0.4,0.5,0.8;]
reverse_delta_mat = [0.5,0.2,0.4;0.4,0.5,0.8;0,0,0]

%pseudo code 
%only for one iteration
%function for ber of k users in T-SIC method 
%input data
%change power of users??????????
for i = 1: length(user_strength)
for j = 1: nsym(i)
    
if user_strength(i) == 1 & j == 1%A1
    %interf by B1, s+1
    delta_i = delta_mat(user_strength(i+1),j);%known
    p_d     = p_a; %desired power
    p_i     = p_b;%interferes power 
    disp('yea')
    fun = @(delta_i) (qfunc(sqrt(3*p_a./(2.*(mod_order-1).*...
    (delta_i.*p_b/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max);
    p_err_sym1 = (1/(timeoff_max -timeoff_min))*q*delta_i;
    p_bita1    = p_err_sym1/log(mod_order)
elseif (user_strength(i) == 1 & j > 1)%A2.... An
    disp('no')
    %interf->B1 and B2 %s-1 and s+1
    delta_i = (1-delta_mat(user_strength(i+1),j-1)) +...
        delta_mat(user_strength(i+1),j);
    fun = @(delta_i) (qfunc(sqrt(3*p_a./(2.*(mod_order-1).*...
        (delta_i.*p_b/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max);
    p_err_sym2 = (1/(timeoff_max -timeoff_min))*q*delta_i;
    p_bita2 = p_err_sym2/log(mod_order) 
elseif (user_strength(i) > 1 & j == 1 & i~=length(user_strength))%B1
    %inter -> A1, A2 and C1
    delta_i = (reverse_delta_mat(user_strength(i-1),j))+...
        (1-reverse_delta_mat(user_strength(i-1),j+1)) +...
        delta_mat(user_strength(i+1),j);
    fun = @(delta_i) (qfunc(sqrt(3*p_a./(2.*(mod_order-1).*...
        (delta_i.*p_b/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max);
    p_err_sym2 = (1/(timeoff_max -timeoff_min))*q*delta_i;
    p_bita2 = p_err_sym2/log(mod_order) 
elseif (user_strength(i) > 1 & j > 1 & i~=length(user_strength))%B2.....Bn & not last user 
    disp('here1 ')
    %inter -> A2, A3 and C1, C2
    delta_i =(reverse_delta_mat(user_strength(i-1),j)) +...
        1-delta_mat(user_strength(i-1),j+1)+...
        (1-delta_mat(user_strength(i+1),j-1)) +...
        delta_mat(user_strength(i+1),j);
    fun = @(delta_i) (qfunc(sqrt(3*p_a./(2.*(mod_order-1).*...
        (delta_i.*p_b/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max);
    p_err_sym2 = (1/(timeoff_max -timeoff_min))*q*delta_i;
    p_bita2 = p_err_sym2/log(mod_order) 
elseif (user_strength(i) > 1 & j >= 1 & i==length(user_strength))
    disp('here2 ')
    %interf-> B1 and B2
    delta_i =(reverse_delta_mat(user_strength(i-1),j)) +...
        1-delta_mat(user_strength(i-1),j+1);
    fun = @(delta_i) (qfunc(sqrt(3*p_a./(2.*(mod_order-1).*...
        (delta_i.*p_b/2+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max);
    p_err_sym2 = (1/(timeoff_max -timeoff_min))*q*delta_i;
    p_bita2 = p_err_sym2/log(mod_order) 
else
    disp('oops')
end
end
end