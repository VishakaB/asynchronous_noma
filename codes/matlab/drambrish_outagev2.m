% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5G Wireless Communications
% Outage Analysis in MIMO systems
%
% Example MatLab script for Outage probability analysis in MIMO
%
% MSc in Wireless Communication
% Office of School of Postgraduate Studies & Research 
% Sri Lanka Technological Campus, Padukka
% Sri Lanka | SLTC | www.sltc.lk
% 
% October 2022
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inicializations
clc;
close all;
clear all;

mT = 4;
mR = 3;
N = 10^5;%number of iterations
SNRdB = 1:20;%SNR range in db scale
SNR = 10.^(SNRdB/10);

%capacity 
C_SISO = zeros(1,length(SNR));
C_SIMO = zeros(1,length(SNR));
C_MISO = zeros(1,length(SNR));
C_MIMO = zeros(1,length(SNR));

rateth = 5;
p1 = zeros(1,length(SNR));
p2 = zeros(1,length(SNR));

%% Capacity expressions

for ite = 1:N
h_SISO = (randn +1i*randn)/sqrt(2);
h_SIMO = (randn(mR,1)+1i*randn(mR,1))/sqrt(2);
h_MISO = (randn(1,mT)+1i*randn(1,mT))/sqrt(2);
h_MIMO = (randn(mR,mT)+1i*randn(mR,mT))/sqrt(2);

for K = 1:length(SNR)
C_SISO(K) = C_SISO(K) + log2(1+ SNR(K)*norm(h_SISO)^2);
C_SIMO(K) = C_SIMO(K) + log2(1+ SNR(K)*norm(h_SIMO)^2);
C_MISO(K) = C_MISO(K) + log2(1+ SNR(K)*norm(h_MISO)^2/mT);
C_MIMO(K) = C_MIMO(K) + log2(abs(det(eye(mR)+SNR(K)*h_MIMO*h_MIMO'/mT)));
end

end

%average rate
C_SISO = (C_SISO/N);
C_SIMO = (C_SIMO/N);
C_MISO = (C_MISO/N);
C_MIMO = (C_MIMO/N);

%% Outage calculations
%Check for outage
for k = 1:length(SNR)
    if C_SISO(k) < rateth
        p1(k) = p1(k)+1;
    end 
end

for k = 1:length(SNR)
    if C_SIMO(k) < rateth
        p2(k) = p2(k)+1;
    end 
end

for k = 1:length(SNR)
    if C_MISO(k) < rateth
        p3(k) = p3(k)+1;
    end 
end

for k = 1:length(SNR)
    if C_MIMO(k) < rateth
        p4(k) = p4(k)+1;
    end 
end

figure(1)
grid on;
plot(SNRdB,C_SISO,'r',SNRdB,C_SIMO,'b',SNRdB,C_MISO,...
    'm',SNRdB,C_MIMO,'k')
legend('SISO','SIMO','MISO','MIMO');
xlabel('SNR in dB');
ylabel('Capacity (b/s/Hz)');
title('Capacity Vs. SNR');

figure(2)
grid on;
plot(SNRdB,smooth(smooth(p1)),'r',SNRdB,smooth(smooth(p2)),'b--')
legend('SIMO','MISO');
xlabel('SNR in dB');
ylabel('Outage (b/s/Hz)');
title('Outage Vs. SNR');