clc;
close all;
clear all;

mT = 2;
mR = 2;
N = 10000; %number of bits
ITER = 1000;%number of trials
SNRdB = [0:25];%snr range in dB
SNR = 10.^(SNRdB/10);%SNR in linear scale
C_SISO = zeros(1,length(SNR));
C_SIMO = zeros(1,length(SNR));
C_MISO = zeros(1,length(SNR));
C_MIMO = zeros(1,length(SNR));
p_simo = zeros(1,length(SNR));
rateth = 1; snrth = 2;
p1 = zeros(1,length(SNR));
p2 = zeros(1,length(SNR));
count = 0;
count2 =0;

for ite = 1:ITER
    
    h_SISO = (randn +1i*randn)/sqrt(2);
    h_SIMO = (randn(mR,1)+1i*randn(mR,1))/sqrt(2);
    h_MISO = (randn(1,mT)+1i*randn(1,mT))/sqrt(2);
    h_MIMO = (randn(mR,mT)+1i*randn(mR,mT))/sqrt(2);

    for K = 1:length(SNR)
        C_SISO(K) = C_SISO(K) + log2(1+ SNR(K)*norm(h_SISO)^2);
        C_SIMO(K) = C_SIMO(K) + log2(1+ SNR(K)*norm(h_SIMO)^2);
        C_MISO(K) = C_MISO(K) + log2(1+ SNR(K)*norm(h_MISO)^2/mT);
        C_MIMO(K) = C_MIMO(K) + log2(abs(det(eye(mR)+SNR(K)*h_MIMO*h_MIMO'/mT)));

        for b = 1:N        
            s(b) = SNR(K)*norm(h_SIMO)^2;
            g(b) = SNR(K)*norm(h_MISO)^2/mT;
            if(s(b)<snrth)
                count = count+1;
            end    
            if(g(b)<snrth)
                count2 = count2+1;
            end   
        end
    p_simo(K)  = count/N;   
    p_miso(K)  = count2/N;  
    end

end

%average caapcity
C_SISO = (C_SISO/ITER);
C_SIMO = (C_SIMO/ITER);
C_MISO = (C_MISO/ITER);
C_MIMO = (C_MIMO/ITER);

%average outage
p_simo  = p_simo/ITER; 
p_miso  = p_miso/ITER; 

figure(1)
plot(SNRdB,C_SISO,'r',SNRdB,C_SIMO,'b',SNRdB,C_MISO,'m',SNRdB,C_MIMO,'k')
legend('SISO','SIMO','MISO','MIMO');
xlabel('SNR in dB');
ylabel('Capacity (b/s/Hz)');
title('Capacity Vs. SNR');
grid;

figure(2)
grid on;
plot(SNRdB,p_simo,'r')
hold on 
plot(SNRdB,p_miso,'b--')
legend('SIMO','MISO');
xlabel('SNR in dB');
ylabel('Outage (b/s/Hz)');
title('Outage Vs. SNR');
