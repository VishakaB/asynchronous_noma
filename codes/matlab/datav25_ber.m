%june 9 2022
%goal: ber vs snr in anoma

%input data 
mod_order = 4;
timeoff_min = 0.1;
timeoff_max = 1;

user_no = 1;
symb_no = 1;

%prob_error = 1/(timeoff_max - timeoff_min)*
prob_error =1;
ber = prob_error/log(mod_order);
%delta = 2
%term1 = qfunc(sqrt(delta^-1))
pk = 1;
noise = 1;
p1  = 1;
delta_i = 1.2;

%b =3;
%test
%ref:https://in.mathworks.com/matlabcentral/...
%answers/383872-matrix-dimensions-must-agree-using-integral
%fcn = @(x)1./((x.^2)+(2.*b.*x)+1)
%fcn_b = integral(fcn,0,Inf)
fun = @(delta_i) 1-(1 - qfunc(sqrt(3*pk./(2.*(mod_order-1).*(delta_i.*pi/2+noise))))).^2;
q   = integral(fun,0.1,1)

p_err_sym = (1/(timeoff_max -timeoff_min))*q

p_bit = p_err_sym/log(mod_order)

fun = @(delta) 1 - ( 1 - qfunc(sqrt(3*received_sinr/(2*(mod_order - 1)*delta))));
q = integral(fun,timeoff_min,timeoff_max);


