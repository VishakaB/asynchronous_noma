%date: june 14 2022
%goal: maximize security capacity

N = 10^4;%number of bits
h = (randn(1,N)+1i*randn(1,N))/sqrt(2);
rho = 0.58;%power 

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');%tool

fun = @(x) -log((1+((0.6.*rho.*x(1).*abs(sum(h.^2)))/((1-rho)*x(2).*abs(sum(h.^2))+1)))...
./(1+(0.4.*rho.*abs(sum(h.^2))/((1-rho).*x(2).*abs(sum(h.^2))+1))));%alpha1=0.6, alpha2=0.4

A = [];
b = [];
Aeq = [];
beq = [];
lb = [0,-1];
ub = [1,1];
nonlcon=[]
x = [0.1,0.1];%w,v,rho
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
-log((1+((0.6.*rho.*x(1).*abs(sum(h.^2)))/((1-rho)*x(2).*abs(sum(h.^2))+1)))...
./(1+(0.4.*rho.*abs(sum(h.^2))/((1-rho).*x(2).*abs(sum(h.^2))+1))))
