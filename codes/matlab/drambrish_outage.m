clc;clear all;close all;
%MISO
N = 10^5;
nt = 2; snrth = 1;
s =  zeros(1,N); snrdb = 2:2:12;
p_sim  = zeros(1,length(snrdb));
p1_ana = zeros(1,length(snrdb));
h1 = 1/sqrt(2)*(randn(1,N)+1i*randn(1,N));
h2 = 1/sqrt(2)*(randn(1,N)+1i*randn(1,N));
for jj = 1:length(snrdb)
    count = 0;
    snr = 10.^(snrdb(jj)/10);
    for i = 1:N
        s(i) = (abs(h1(i)).^2+abs(h2(i)).^2).*snr;
        if(s(i)<snrth)
            count = count+1;
        end    
    end
    p_sim(jj)  = count/N;
p1_ana(jj) = 1 - exp(-snrth/snr).*(snrth./snr+1);
end

figure (1);
semilogy(snrdb, p_sim,'r--');
hold on;
semilogy(snrdb,p1_ana,'b');
hold on;

%% SIMO

%MISO
N = 10^5;
nr = 2; snrth = 1;
s =  zeros(1,N); snrdb = 2:2:12;
p_sim  = zeros(1,length(snrdb));
p1_ana = zeros(1,length(snrdb));
h1 = 1/sqrt(2)*(randn(1,N)+1i*randn(1,N));

for jj = 1:length(snrdb)
    count = 0;
    snr = 10.^(snrdb(jj)/10);
    for i = 1:N
        s(i) = (abs(h1(i)).^2).*snr;
        if(s(i)<snrth)
            count = count+1;
        end    
    end
    p_sim2(jj)  = count/N;
%p1_ana(jj) = 1 - exp(-snrth/snr).*(snrth./snr+1);
end
hold on
semilogy(snrdb, p_sim2,'m');
legend('Simulation','Analytical','SIMO')