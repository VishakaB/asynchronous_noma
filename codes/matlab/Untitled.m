close all
clear all
clc
%channel gain
N = 100;
communication_radius = 30;%change this 
max_dist     = 100;%meters
max_eta      = 15;
etath        = 4;%change this 

h_th         = sqrt(communication_radius^-etath)*sqrt(1/2)*(randn(1,N)+...
1i*randn(1,N))/sqrt(2);
g_th         = (abs(h_th)).^2;

alpha1 = 0.4
alpha2 = 1-alpha1
omeg = (linspace(0,1,100)) %variable

rho_s = 0.9
mu = (linspace(-1,1,100))%variable

num = alpha1.*rho_s.*omeg
den = (1-rho_s).*mu+1

s1 = num./den

num2 = alpha2.*rho_s
den2 = (1-rho_s).*mu+1

s2 = num2./den2

x1 = zeros(length(s1),1);
for i = 1:length(s1)
    if (s1(i)<0.0200)
        x1(i) = 0
    else
        x1(i) = s1(i)
    end
end


y1 = zeros(length(s2),1);
for j = 1:length(s2)
    if (s2(j)>0.5933)
        y1(j) = s2(j);
    else
        y1(j) = 0;
    end
end

x1
y1 

x2 = zeros(length(s1),1);
for i = 1:length(s1)
    if (s1(i)<0.0278)
        x2(i) = 0
    else
        x2(i) = s1(i)
    end
end

y2 = zeros(length(s2),1);
for j = 1:length(s2)
    if (s2(j)>0.5907)
        y2(j) = s2(j);
    else
        y2(j) = 0;
    end
end

z1 = log2(([1+x1])./([1+y1]))%objective
plot(x1,z1)

hold on;
z2 = log2(([1+x2])./([1+y2]))%objective
plot(x2,z2,'--r')
