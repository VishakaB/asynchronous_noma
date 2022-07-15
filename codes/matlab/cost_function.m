

0< omeg < 1
-1< mu <1
alpha1  = 0.4
alpha2  = 0.6
rho_s = 0.9;

lb = [0,-1];
ub = [1,1];

opt_var = [omeg,mu];
x0 = [0.1,0.1]

x = fmincon(cost_function,x0,[],[],[],[],lb,ub,constraint)

function [c,eq] = constraint(opt_var)
    c = opt_var;
    ceq = [];
end

function f = cost_function(omeg,mu,alpha1,alpha2,rho_s)

  num = alpha1.*rho_s.*omeg
  den = (1-rho_s).*mu+1
  s1 = num/den
  
  num2 = alpha2.*rho_s
  den2 = (1-rho_s).*mu+1
  s2 = num2/den2
 
  f = log2((1+s1)/(1+s2))
  
end