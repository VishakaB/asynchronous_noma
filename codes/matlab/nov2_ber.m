sumsym_dur_vec = [0.1,0.1,0.1;0.2,0.3,0.4;0.5,0.4,0.3;]*0.01;
power_vec = [0.1;0.01;0.001];
timeoff_min = 0.1;
timeoff_max = 0.5;
mod_order = 4;
noise = 0.1;

%% strongest user
delta_i = sumsym_dur_vec(1,1);%known
p_d     = power_vec(1); %desired power
p_iw    = power_vec(2);%interferes power 
fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
    (delta_i.*p_iw/2+noise))))).^2;
q   = integral(fun,timeoff_min,timeoff_max);
p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
p_bita1_1    = p_err_sym1/log(mod_order);

delta_i = sumsym_dur_vec(1,1)+sumsym_dur_vec(1,2);%known
p_d     = power_vec(1); %desired power
p_iw    = power_vec(2);%interferes power 
fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
    (delta_i.*p_iw/2+noise))))).^2;
q   = integral(fun,timeoff_min,timeoff_max);
p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
p_bita1_2    = p_err_sym1/log(mod_order);

delta_i = sumsym_dur_vec(1,1)+sumsym_dur_vec(1,2)+sumsym_dur_vec(1,3);%known
p_d     = power_vec(1); %desired power
p_iw    = power_vec(2);%interferes power 
fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
    (delta_i.*p_iw/2+noise))))).^2;
q   = integral(fun,timeoff_min,timeoff_max);
p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
p_bita1_3    = p_err_sym1/log(mod_order);

ber_1 = 1/3*(p_bita1_1 + p_bita1_2 + p_bita1_3)

%% second strongest user
delta_i = sumsym_dur_vec(2,1)+sumsym_dur_vec(2,2)+sumsym_dur_vec(2,3);%known
p_d     = power_vec(2); %desired power
p_iw    = power_vec(1);%interferes power 
fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
    ((6/(mod_order -1))*delta_i.*p_iw/2+noise))))).^2;
q   = integral(fun,timeoff_min,timeoff_max);
p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
p_bita2_1    = p_err_sym1/log(mod_order);

delta_i = sumsym_dur_vec(2,1)+sumsym_dur_vec(2,2)+sumsym_dur_vec(2,3);%known
p_d     = power_vec(2); %desired power
p_iw    = power_vec(1);%interferes power 
fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
    ((6/(mod_order -1))*delta_i.*p_iw/2+noise))))).^2;
q   = integral(fun,timeoff_min,timeoff_max);
p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
p_bita2_2    = p_err_sym1/log(mod_order);

ber_2 = 0.5*(p_bita2_1 + p_bita2_2)*[p_bita1_1*p_bita1_2*p_bita1_3 +...
    (1-p_bita1_1)*(1-p_bita1_2)*(1-p_bita1_3)]

%% weakest user

delta_i = [0.5*(sumsym_dur_vec(3,1)+sumsym_dur_vec(3,2)+sumsym_dur_vec(3,3));...
    0.5*(sumsym_dur_vec(3,1)+sumsym_dur_vec(3,2)+sumsym_dur_vec(3,3))];%known
p_d     = power_vec(3); %desired power
p_iw    = power_vec(2);%interferes power 
p_iw1   = power_vec(1);
fun = @(delta_i) 1-(1-qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
    ((6/(mod_order -1))*sum(delta_i.*[p_iw;p_iw1;])/2+noise))))).^2;
q   = integral(fun,timeoff_min,timeoff_max,'ArrayValued', true);
p_err_sym1 = (1/(timeoff_max - timeoff_min))*mean(q.*(delta_i));
p_bita3    = p_err_sym1/log(mod_order);

ber_3 = p_bita3*[ber_1*ber_2+ber_1*(1 - ber_2)+...
    (1-ber_1)*ber_2+(1-ber_1)*(1-ber_2)]
