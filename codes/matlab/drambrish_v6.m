%start date: june 14 2022
%last update: 8 july 2022

%goal: maximize security capacity
clear all
close all
clc

%initalization 
h = 1;%channel state information
rho = 5;%transmission power
alpha1 = 0.9999;%power allocation coefficients
alpha2 = 0.0001;
  
%% optimization

%objective function
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
fun = @(x) -log2(((1+(alpha1.*rho.*x(1))/((1-rho)*x(2).*abs(sum(h.^2))+1)))...
./(1+(alpha2.*rho)/((1-rho).*x(2))+1));

x0 = [0.5,-0.4];%w,v,rho %initial solution
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0,-1];
ub = [1,1];

%constraints
nonlcon8 = @(x) constr5(x);
nonlcon9 = @(x) constr6(x,rho,h);
nonlcon12 = @(x)constr9(x,rho,alpha1);
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,...
    {nonlcon9,nonlcon8,nonlcon12},options);

%optimal CS value
num = alpha1.*rho.*x(1);
den = (1-rho).*x(2)+1;
s1 = num/den;

num2 = alpha2.*rho;
den2 = (1-rho).*x(2)+1;
s2 = num2/den2;

opt_omeg = x(1);
opt_mu   = x(2);
z = abs(log2((1+s1)/(1+s2)));

fprintf('omega:%f \n',opt_omeg)
fprintf('nev :%f\n',opt_mu)
fprintf('CS values:%f \n',z)

%% constraints

function [c,ceq] = constr5(x)
c = x(2)-x(1);
ceq = [];
end

function [c,ceq] = constr6(x,rho,h)%omeg = 0.3, new = 0.1
c = ((0.6.*rho.*0.3.*...
    abs(sum(h.^2)))/((1-rho)*0.1.*abs(sum(h.^2))+1))-((0.6.*rho.*x(1).*...
    abs(sum(h.^2)))/((1-rho)*x(2).*abs(sum(h.^2))+1));
ceq = [];
end

function [c,ceq] = constr9(x,rho,alpha1,alpha2)
num = alpha1.*rho.*x(1);
den = (1-rho).*x(2)+1;
s1 = num/den;
num2 = alpha2.*rho;
den2 = (1-rho).*x(2)+1;
s2 = num2/den2;
c  = s2 - s1;
ceq = [];
end