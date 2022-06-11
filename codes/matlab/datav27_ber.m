%june 9 2022
%goal: ber vs snr in anoma
%ref:https://in.mathworks.com/matlabcentral/...
%answers/383872-matrix-dimensions-must-agree-using-integral
clear all; 
clc;

%initial input data 
mod_order = 4;
timeoff_min = 0.01;
timeoff_max = 0.5;
p_a = 1000; %40 dB
p_b = 200;
p_c = 100;
noise = 0.1;

% bit error rate general 
delta_i = 0.003;
fun = @(delta_i) 1-(1 - qfunc(sqrt(3*p_a./(2.*(mod_order-1).*(delta_i.*p_b/2+noise))))).^2;
q   = integral(fun,0.1,1);
p_err_sym = (1/(timeoff_max -timeoff_min))*q*delta_i;
p_bit = p_err_sym/log(mod_order);

%% first user BER
% bit error rate symbol wise 
% first user first symbol %s = 11 %A1
%A1 -> interf by B1
delta_21 = 0.1;
delta_22 = 0.05;
s = 1;
delta_i = delta_21;
fun = @(delta_i) 1-(1 - qfunc(sqrt(3*p_a./(2.*(mod_order-1).*(delta_i.*p_b/2+noise))))).^2;
p_err_sym1 = (1/(timeoff_max -timeoff_min))*q*delta_i;
p_bita1 = p_err_sym1/log(mod_order)

% first user second symbol s = 12 %A2
s = 2;
%A2-> interf by B1 and B2
delta_i = delta_21 + delta_22;
fun = @(delta_i) 1-(1 - qfunc(sqrt(3*p_a./(2.*(mod_order-1).*(delta_i.*p_b/2+noise))))).^2;
p_err_sym2 = (1/(timeoff_max -timeoff_min))*q*delta_i;
p_bita2 = p_err_sym2/log(mod_order)

%% second user BER
% second user second symbol s = 12 %B1
%B1-> interf by A1 and A2, and C1
delta_11 = 0.2;
delta_12 = 0.1;
delta_31 = 0.04;
delta_i = delta_11 + delta_12 ; %interference from first user detected symbols

fun = @(delta_i) 1-(1 - qfunc(sqrt(3*pk./(2.*(mod_order-1).*...
    (delta_i.(6/(mod_order-1))*p_a/2+ delta_31*p_c/2+noise))))).^2;
p_interb1 = (1/(timeoff_max -timeoff_min))*q*delta_i;
p_err_sym = p_interb1*p_err_sym1*p_err_sym2;
p_bit1 = p_err_sym/log(mod_order)


%check from here 
delta_i = delta_11; %interference from first user detected symbols
fun = @(delta_i) 1-(1 - qfunc(sqrt(3*pk./(2.*(mod_order-1).*...
    (delta_i.(6/(mod_order-1))*p_a/2+delta_31*p_c/2+noise))))).^2;
p_interb2 =(1/(timeoff_max -timeoff_min))*q*delta_i;
p_err_sym = p_interb2*(1-p_err_sym1)*p_err_sym2;
p_bit2 = p_err_sym/log(mod_order)

delta_i = delta_12; %interference from first user detected symbols
fun = @(delta_i) 1-(1 - qfunc(sqrt(3*pk./(2.*(mod_order-1).*...
    (delta_i.(6/(mod_order-1))*p_a/2+delta_31*p_c/2+noise))))).^2;
p_interb3 =(1/(timeoff_max -timeoff_min))*q*delta_i;
p_err_sym = p_interb3*p_err_sym1*(1-p_err_sym2);
p_bit3 = p_err_sym/log(mod_order)

delta_i = 0;%interference from first user detected symbols
fun = @(delta_i) 1-(1 - qfunc(sqrt(3*pk./(2.*(mod_order-1).*...
    (delta_i.(6/(mod_order-1))*p_a/2+delta_31*p_c+noise))))).^2;
p_interb4 = (1/(timeoff_max -timeoff_min))*q*delta_i;
p_err_sym = p_interb4*(1-p_err_sym1)*(1-p_err_sym2);
p_bit4 = p_err_sym/log(mod_order)

pbit = p_bit1 + p_bit2 + p_bit3 + p_bit4 

%B2-> interf by A2 and A3, C1 and C2
delta_13 = 0.02;
delta_32 = 0.01;
delta_i = delta_12 + delta_13 + delta_31 + delta_32; %interference from first user detected symbols

%% third user BER
%C1-> interf by B1 and B2
delta_21 = 0.2;
delta_22 = 0.1;
delta_i = delta_21 + delta_22; %interference from first user detected symbols
