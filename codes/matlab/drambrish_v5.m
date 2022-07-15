clear all;
close all;
clc;

x = optimvar('x');
y = optimvar('y');
 
obj = objfunx(x,y);
prob = optimproblem('Objective',obj);

initx =  0.1;
inity = 1;

for i = 1: 4
        if (i ==1)
            con1 = x <= 1; 
            prob.Constraints.constr = con1;
        elseif (i ==2)
            con2 = x >= 0; 
            prob.Constraints.constr = con2;
            
        elseif (i==3)
            con3 = y>=-1; 
            prob.Constraints.constr = con3;
        else
            con4 = y <= 1; 
            prob.Constraints.constr = con4;
        end
x0.x = initx;
x0.y = inity;
show(prob)
[sol,fval] = solve(prob,x0)
initx = sol.x;
inity = sol.y;
end

function f = objfunx(x,y)
rho = 0.9;
h = 1;
f= - log((1+((0.6.*rho.*x.*...
    abs(sum(h.^2)))/((1-rho)*y.*abs(sum(h.^2))+1)))...
./(1+(0.4.*rho.*abs(sum(h.^2))/((1-rho).*y.*abs(sum(h.^2))+1))))
end
 
%{
function f = objfunx(x,y)
f = log(1+x/y) - log(1+1/y);
end
%}
%{
x = optimvar('x','type','integer','LowerBound',0,'UpperBound',48);
y = optimvar('y','type','integer','LowerBound',0,'UpperBound',48);
obj = fcn2optimexpr(@objfunx,x,y);
prob = optimproblem('Objective',obj);
con1 = x <= 48; prob.Constraints.constr = con1;
con2 = y <= 48; prob.Constraints.constr = con2;
con3 = (1-((0.2.^x))) + ((0.2.^x)).*(1-((0.3.^y))) >= 0.99;
prob.Constraints.constr = con3;
x0.x = 0; x0.y = 0;
show(prob)
%}

 

