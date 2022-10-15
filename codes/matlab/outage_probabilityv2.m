clc;
close all;
clear all;

% goal: outage probability analysis in MIMO
%% initial parameters
SNR = 10 ^ (SNR/10);    %SNR in linear scale
n_t = ;      % number of Tx antennas
n_r = ;      % number of Rx antennas
rateth = ;             % Target rate of users in bps/Hz

Pt = 0:2:40;                %Transmit power in dBm
pt = (10^-3)*10.^(Pt/10);   %Transmit power in linear scale
BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

p = length(Pt);

for u = 1:p
    %Calculate SNRs
    gamma_1 = ;
    gamma_12 = a1*pt(u)*g2./(a2*pt(u)*g2+no);
    gamma_2 = a2*pt(u)*g2/no;
    
    %Calculate achievable rates
    R1 = log2(1+gamma_1);
    R12 = log2(1+gamma_12);
    R2 = log2(1+gamma_2);
    
    %Find average of achievable rates
    R1_av(u) = mean(R1);
    R12_av(u) = mean(R12);
    R2_av(u) = mean(R2);
    
    %Check for outage
    for k = 1:N
        if R1(k) < rate1
            p1(u) = p1(u)+1;
        end
        if (R12(k) < rate1)||(R2(k) < rate2)
            p2(u) = p2(u)+1;
        end
    end
end

pout1 = p1/N; 
pout2 = p2/N;

figure;
semilogy(Pt, pout1, 'linewidth', 1.5); hold on; grid on;
semilogy(Pt, pout2, 'linewidth', 1.5);
xlabel('Transmit power (dBm)');
ylabel('Outage probability');
legend('User 1 (far user)','User 2 (near user)');