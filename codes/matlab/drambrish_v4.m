
alpha1  = 0.4
alpha2  = 0.6
rho_s = 0.9;

lb = [0,-1];
ub = [1,1];

%opt_var = [omeg,mu];
x0 = [0.1,0.1]%initial solution

z = @(x)log2((1+x)/(1+x));

x = fmincon(z,x0,[],[],[],[],lb,ub,constraint)

function [c,eq] = constraint(x)
    c = x;
    ceq = [];
end

 
function z = cost_function(x,alpha1,alpha2,rho_s)
  %omeg = x(1)
  %mu   = x(2)
  num = alpha1.*rho_s.*x(1)
  den = (1-rho_s).*x(2)+1
  s1 = num/den
  
  num2 = alpha2.*rho_s
  den2 = (1-rho_s).*mu+1
  s2 = num2/den2
 
  z = log2((1+s1)/(1+s2))
  z = @(z)log2((1+s1)/(1+s2));
  
end