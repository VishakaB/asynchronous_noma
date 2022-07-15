%start date: june 14 2022
%last update: 8 july 2022

%goal: maximize security capacity
clear all
close all
clc

for i = 1: 10
%rng(2)%seed
%initalization 
N = 10^3;
rho = 6;%transmission power
h = 1/sqrt(rho)*abs(rand(N,1)+i*rand(N,1));%channel information
alpha1 = 0.999;%power allocation coefficients
alpha2 = 0.001;
  
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
    {nonlcon8,nonlcon12},options);

%optimal CS value
num = alpha1.*rho.*x(1);
den = (1-rho).*x(2)+1;
s1(i) = num/den;

num2 = alpha2.*rho;
den2 = (1-rho).*x(2)+1;
s2(i) = num2/den2;

opt_omeg(i) = x(1);
opt_mu(i)   = x(2);
z(i) = (log2((1+s1)/(1+s2)));
%{
fprintf('omega:%f \n',opt_omeg)
fprintf('nev :%f\n',opt_mu)
fprintf('CS values:%f \n',z)
%}
end

figure (1)
plot(s1,z,'b')

%% constraints
function [c,ceq] = constr5(x)%omega > nev
c = x(2)-x(1);
ceq = [];
end

function [c,ceq] = constr6(x,rho,h)%CS > 0
c = ((0.6.*rho.*0.3.*...
    abs(sum(h.^2)))/((1-rho)*0.1.*abs(sum(h.^2))+1))-((0.6.*rho.*x(1).*...
    abs(sum(h.^2)))/((1-rho)*x(2).*abs(sum(h.^2))+1));
ceq = [];
end

function [c,ceq] = constr9(x,rho,alpha1,alpha2) %sinr desired > sinr eavesdropper
num = alpha1.*rho.*x(1);
den = (1-rho).*x(2)+1;
s1 = num/den;
num2 = alpha2.*rho;
den2 = (1-rho).*x(2)+1;
s2 = num2/den2;
c  = s2 - s1;
ceq = [];
end

%{
function [c,ceq] = constr10(x,rho,alpha1,alpha2) %sinr desired > sinr eavesdropper
num = alpha1.*rho.*x(1);
den = (1-rho).*x(2)+1;
s1 = num/den;
num2 = alpha2.*rho;
den2 = (1-rho).*x(2)+1;
s2 = num2/den2;
c  = s1;
ceq = [];
end
%}