%define delta_i as a vec 
%include in integral

delta_i = [(reverse_delta_mat(user_strength(i-1),j));...
        (1-reverse_delta_mat(user_strength(i-1),j+1)); ...
        delta_mat(user_strength(i+1),j);];
 
    p_d     = power_vec(i); %desired power
    p_is    = power_vec(i-1);%interferes power
    p_iw    = power_vec(i+1);%interferes power
    
    power_v  = [p_d/2;p_is/2;p_is/2];
    
    %continue here..................
    %include the delta i into the func expression???????
    fun = @(delta_i) (qfunc(sqrt(3*p_d./(2.*(mod_order-1).*...
        (sum(delta_i.*power_v)+noise))))).^2;
    q   = integral(fun,timeoff_min,timeoff_max)