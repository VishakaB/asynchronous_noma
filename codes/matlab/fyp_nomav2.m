%dte: july 26
%goal: maximize the receiver sum rate wrt transmit SNR
clc; clear variables; close all;

N  = 10;
Pt = 30;                        %Max BS Tx power (dBm)
pt = (10^-3)*db2pow(Pt);        %Max BS Tx power (Linear scale)
No = -114;                      %Noise power (dBm)
no = (10^-3)*db2pow(No);        %Noise power (Linear scale)

snr = 1:1:40;                 %Far user target rate range (R*)

df = 1000; dn = 500;            %Distances

eta = 4;                        %Path loss exponent

p1 = zeros(1,length(snr));
p2 = zeros(1,length(snr));
pa1 = zeros(1,length(snr));
pa2 = zeros(1,length(snr));

af = 0.75; an = 0.25;       %Fixed PA (for comparison)

hf = sqrt(df^-eta)*(randn(1,N) + 1i*randn(1,N))/sqrt(2);
hn = sqrt(dn^-eta)*(randn(1,N) + 1i*randn(1,N))/sqrt(2);

g1 = (abs(hf)).^2;
g2 = (abs(hn)).^2;

gamma_f = zeros(length(snr),N);
gamma_n = zeros(length(snr),N);
Cf = zeros(length(snr),N);
Cn = zeros(length(snr),N);

for u = 1:length(snr)
epsilon(u) = (2^(snr(u)))-1;         %Target SINR for far user 

n0 = (sqrt(1)*10^(-snr(u)/20))*no; %Addition of Noise
 
%rate of far user 
gamma_f(u,:) = pt*af*g1./(pt*g1*an + no);

%rate of near user 
gamma_n(u,:) = pt*g2*an/no;

Cf(u,:) = log2(1 + gamma_f(u,:));

Cn(u,:) = log2(1 + gamma_n(u,:));

sum_rate(u) = mean(Cn(u,:)+Cf(u,:));
end

figure;
plot(snr,sum_rate,'--+r','linewidth',2); hold on; grid on;

xlabel('SNR');
ylabel('Sum rate');
%legend('Far user','Near user','Far user (fair PA)','Near user (fair PA)');

















