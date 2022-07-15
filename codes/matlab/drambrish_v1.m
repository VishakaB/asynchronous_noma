%date: june 14 2022

%goal: maximize security capacity
clear all
close all
clc
N = 10;
h = 1;
rho = 100;
alpha1 = 0.9;
alpha2 = 0.1;

for i = 1: 1
    
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
fun = @(x) -log2((1+((0.6.*rho.*x(1).*...
    abs(sum(h.^2)))/((1-rho)*x(2).*abs(sum(h.^2))+1)))...
./(1+(0.4.*rho.*abs(sum(h.^2))/((1-rho).*x(2).*abs(sum(h.^2))+1))));%alpha1=0.6, alpha2=0.4

A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];

nonlcon2 = @(x) constr1(x);
nonlcon3 = @(x) constr11(x);
nonlcon4 = @(x) constr2(x);
nonlcon5 = @(x) constr22(x);
nonlcon6 = @(x) constr3(x,rho);
nonlcon7 = @(x) constr4(x,rho);
nonlcon8 = @(x) constr5(x,rho);
nonlcon9 = @(x) constr6(x,rho,h);

x0 = [0.01,-0.01];%w,v,rho %initial solution
x  = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,{nonlcon2,...
    nonlcon3,nonlcon4,nonlcon5,nonlcon6,nonlcon7,nonlcon8,nonlcon9},options);

opt_omeg(i) = x(1)
opt_mu(i)   = x(2)

num = alpha1.*rho.*x(1);
den = (1-rho).*x(2)+1;
s1(i) = num/den;

num2 = alpha2.*rho;
den2 = (1-rho).*x(2)+1;
s2(i)= num2/den2;

z(i) = abs(log2((1+s1)/(1+s2)))

end

plot(opt_omeg,z,'--r')
hold on;
plot(opt_mu,z,'--k')
hold on
plot(s1,z,'g')
set(gca,'yscale','log')

function [c,ceq] = constr1(x)
c = x(1)-1;
ceq = [];
end

function [c,ceq] = constr11(x)
c = -x(1);
ceq = [];
end

function [c,ceq] = constr2(x)
c = -1 - x(2);
ceq = [];
end

function [c,ceq] = constr22(x)
c = x(2)-1;
ceq = [];
end

function [c,ceq] = constr3(x,rho)
c = -x(1)*rho/(x(2)+1)+1e6;
ceq = [];
end

function [c,ceq] = constr4(x,rho)
c = -rho/(x(2)+1)+1e6;
ceq = [];
end

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