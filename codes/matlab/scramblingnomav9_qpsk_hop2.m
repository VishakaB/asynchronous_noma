%ref:https://in.mathworks.com/matlabcentral/fileexchange/23498-mc-cdma
clc;
close all;
clear all;
k1=3; %Rician factor %ref: https://www.researchgate.net/publication/263669548_Probability_Distribution_of_Rician_K-Factor_in_Urban_Suburban_and_Rural_Areas_Using_Real-World_Captured_Data
mean=sqrt(k1/(k1+1));% mean
sigma=sqrt(1/(2*(k1+1)));%variance
N=10^4;  % Number of Bits for data_user1
d1 = 5; d2 = 3;    %Distances of users from base station (BS)
d3 = 3;
eta = 1;                %Path loss exponent

%% ------------------Hop 1 first user ------------------------------
users=1;            % Total Number of relays 
%------------------Generation of Walsh code--------------------------------
n =8;                               %Number of  Data Sub-Carriers used
walsh=hadamard(n);              
code1=walsh(2,:);                   %Taking 2nd row of walsh code for User1

%------------------Generating data for User1-------------------------------                     
data_i = (2*(rand(1,N)>0.5)-1);
data_q = (2*(rand(1,N)>0.5)-1);
data_source_qpsk = data_i + j*data_q; %QPSK modulation

%------------------Spreading & IFFT for User1------------------------------
data_s11=data_source_qpsk';
spdata_s11=data_s11*code1;    % Spreading 
spdatas11t=(spdata_s11)';
ifftdata_s11=ifft(spdatas11t);      % Taking the IFFT
ifftdatas11=ifftdata_s11';

%------------------Append Cyclic Prefix1 for User1-------------------------
y1=[ifftdatas11(:,[(n-2):n]) ifftdatas11];
tx_user1=y1';                % Transmitting data for user1

%----------------------Adding data for Transmission of All User------------
xs11=tx_user1;

%----------------------Creating Rician Channel---------------------------
Taps=4;                                   % Number of Taps
p1=1;                                     % Power of Tap1
p2=1;                                     % Power of Tap2
p3=1;                                     % Power of Tap3
p4=1;
gain1=sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
gain2=sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain3=sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3
gain4=sqrt(p4/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap4
x11=xs11(:);
x12=reshape(x11,1,length(x11));
i=1:length(x12);        
delay1=1; 
for i=delay1+1:length(x12) % Producing one sample delay in Tap2 w.r.t. Tap1
   x13(i)=xs11(i-delay1);
end
delay2=2;
for i=delay2+1:length(x12) % Producing two sample delay in Tap2 w.r.t. Tap1
   x14(i)=xs11(i-delay2);
end
delay3=3;
for i=delay3+1:length(x12) % Producing three sample delay in Tap2 w.r.t. Tap1
   x15(i)=xs11(i-delay3);
end
x10=reshape(x13,(n+3),length(x13)/(n+3));
x20=reshape(x14,(n+3),length(x14)/(n+3));
x30=reshape(x15,(n+3),length(x15)/(n+3));
ch1=repmat(gain1,(n+3),1);     
ch2=repmat(gain2,(n+3),1);
ch3=repmat(gain3,(n+3),1);
ch4=repmat(gain4,(n+3),1);
data_channel1=xs11.*ch1+x10.*ch2+x20.*ch3+x30.*ch4;  % Passing data through channel 

%------------------------Addition of AWGN noise ---------------------------
data_channel =data_channel1;
data_noise=data_channel(:);
data_noise=reshape(data_noise,1,length(data_noise));
noise = 1/sqrt(2)*[randn(1,length(data_noise)) + j*randn(1,length(data_noise))]; 

%% ------ receiver ---------------------
snr = [0:20];                 % multiple Eb/N0 values

for i = 1:length(snr)
y= data_noise + (sqrt(1)*10^(-snr(i)/20))*noise; %Addition of Noise
  
%--------------------------Receiver ---------------------------------------
data_received = y;           %fadded data received with awgn noise
%---------------------Removing Cyclic Prefix-------------------------------
rx1=reshape(data_received,(n+3),length(data_received)/(n+3));
rx12=rx1';
rx13 = rx12(:,[(4:(n+3))]); 
rx14=rx13';

%-----------------Taking FFT ----------------------------------------------
fft_data_received =fft(rx14);

%----------------equilization of the channel-------------------------------
channel_response1=fft([gain1;gain2;gain3;gain4],n);
data_equilized1=fft_data_received.*conj(channel_response1);

%----------------BER of Data User1-----------------------------------------
recdatas11=(data_equilized1'*code1')';
r_s11_i =2*(real(recdatas11)>0)-1;
r_s11_q =2*(imag(recdatas11)>0)-1;
errors_user1i(i) = size(find([data_i- r_s11_i]),2); %Errors for User1
errors_user1q(i) = size(find([data_q- r_s11_q]),2); %Errors for User1
SBer1 = (errors_user1i+errors_user1q)/(2*N);                              % simulated ber user1
                           
end

%% ------------------Hop 2------------------------------

users=2;            % Number of Users
%------------------Generation of Walsh code--------------------------------
code2=walsh(6,:);                   %Taking 3rd row of walsh code for User2
%------------------est for User1-------------------------------

%------------------Spreading & IFFT for User1------------------------------
estdata_i = r_s11_i;
estdata_q = r_s11_q;
estdata_source_qpsk = estdata_i + j*estdata_q; %QPSK modulation

estdata_user11=estdata_source_qpsk';
spdata1_user11=estdata_user11*code1;    % Spreading 
spdata120=(spdata1_user11)';
ifftdata_user11=ifft(spdata120);      % Taking the IFFT
ifftdata12=ifftdata_user11';

%------------------Append Cyclic Prefix1 for User1-------------------------
y1=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata1=y1';
esttx_user1=transdata1;                % Transmitting data for user1
%------------------Generating data for User2-------------------------------
data2_i = (2*(rand(1,N)>0.5)-1);
data2_q = (2*(rand(1,N)>0.5)-1);
data_relay2_qpsk = data2_i + j*data2_q; %QPSK modulation

%-----------------Spreading & IFFT for User2-------------------------------
data_user21=data_relay2_qpsk';
spdata2=data_user21*code2;          % Spreading 
spdata22=(spdata2)';
ifftdata_user2=ifft(spdata22);      % Taking the IFFT
ifftdata22=ifftdata_user2';
%-----------------Append Cyclic Prefix1 for User2--------------------------
y2=[ifftdata22(:,[(n-2):n]) ifftdata22];
transdata2=y2';
tx_user2=transdata2;                % Transmitting data for user2
%----------------------Adding data for Transmission of All User------------
x2=esttx_user1++tx_user2;

%----------------------Creating a new Rician Channel---------------------------
gain1=sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
gain2=sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain3=sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3
gain4=sqrt(p4/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap4
x11=x2(:);
x12=reshape(x11,1,length(x11));
i=1:length(x12);        
delay1=1; 
for i=delay1+1:length(x12) % Producing one sample delay in Tap2 w.r.t. Tap1
   x13(i)=x2(i-delay1);
end
delay2=2;
for i=delay2+1:length(x12) % Producing two sample delay in Tap2 w.r.t. Tap1
   x14(i)=x2(i-delay2);
end
delay3=3;
for i=delay3+1:length(x12) % Producing three sample delay in Tap2 w.r.t. Tap1
   x15(i)=x2(i-delay3);
end
x11=reshape(x13,(n+3),length(x13)/(n+3));
x22=reshape(x14,(n+3),length(x14)/(n+3));
x33=reshape(x15,(n+3),length(x15)/(n+3));
ch1=repmat(gain1,(n+3),1);     
ch2=repmat(gain2,(n+3),1);
ch3=repmat(gain3,(n+3),1);
ch4=repmat(gain4,(n+3),1);
data_channel2=x2.*ch1+x11.*ch2+x22.*ch3+x33.*ch4;  % Passing data through channel 

%------------------------Addition of AWGN noise ---------------------------
data_noise2=data_channel2(:);
data_noise22=reshape(data_noise2,1,length(data_noise2));
noise22 = 1/sqrt(2)*[randn(1,length(data_noise22)) + j*randn(1,length(data_noise22))]; 
               
for i = 1:length(snr)
y2 = data_noise22 + (sqrt(1)*10^(-snr(i)/20))*noise22; %Addition of Noise
  
%--------------------------Receiver ---------------------------------------
data_received2 = y2;           %fadded data received with awgn noise

%---------------------Removing Cyclic Prefix-------------------------------
rx1=reshape(data_received2,(n+3),length(data_received2)/(n+3));
rx12=rx1';
rx13 = rx12(:,[(4:(n+3))]); 
rx14=rx13';

%-----------------Taking FFT ----------------------------------------------
fft_data_received =fft(rx14);

%----------------equilization of the channel-------------------------------
channel_response=fft([gain1;gain2;gain3;gain4],n);
data_equilized2=fft_data_received.*conj(channel_response);

%----------------BER of Data User1-----------------------------------------

recdata110=(data_equilized2'*code1')';
r_s11hop2_i =2*(real(recdata110)>0)-1;
r_s11hop2_q =2*(imag(recdata110)>0)-1;
errors_1hop2i(i) = size(find([data_i- r_s11hop2_i]),2); %Errors for User1
errors_1hop2q(i) = size(find([data_q- r_s11hop2_q]),2); %Errors for User1
SBer11 = (errors_1hop2i+errors_1hop2q)/(2*N);   

%----------------BER of Data User2-----------------------------------------

recdata21=(data_equilized2'*code2')';
r_s22hop2_i =2*(real(recdata21)>0)-1;
r_s22hop2_q =2*(imag(recdata21)>0)-1;
errors_2hop2i(i) = size(find([data2_i- r_s22hop2_i]),2); %Errors for User1
errors_2hop2q(i) = size(find([data2_q- r_s22hop2_q]),2); %Errors for User2
SBer22 = (errors_2hop2i+errors_2hop2q)/(2*N);                            % simulated ber user2

end

%% ------------------Hop 3------------------------------

users=3;            % Number of Users
%------------------Generation of Walsh code--------------------------------
code3=walsh(7,:);                   %Taking 2nd row of walsh code for User1

%------------------est for User1-------------------------------
estdatahop3_i = r_s11hop2_i;
estdatahop3_q = r_s11hop2_q;
estdata_sourcehop3_qpsk = r_s11hop2_i + j*r_s11hop2_q; %QPSK modulation

%------------------Spreading & IFFT for User1------------------------------
estdata_user111=estdata_sourcehop3_qpsk';
spdata1_user111=estdata_user111*code1;    % Spreading 
spdata121=(spdata1_user111)';
ifftdata_user111=ifft(spdata121);      % Taking the IFFT
ifftdata12=ifftdata_user111';

%------------------Append Cyclic Prefix1 for User1-------------------------
y11=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata11=y11';
esttx_user11=transdata11;                % Transmitting data for user1

%------------------est for User2-------------------------------

%------------------Spreading & IFFT for User2------------------------------
estdatahop2s2_i = r_s22hop2_i;
estdatahop2s2_q = r_s22hop2_q;
estdata_source2hop3_qpsk = r_s22hop2_i + j*r_s22hop2_q; %QPSK modulation

estdata_user222=estdata_source2hop3_qpsk';
spdata1_user222=estdata_user222*code2;    % Spreading 
spdata222=(spdata1_user222)';
ifftdata_user22=ifft(spdata222);      % Taking the IFFT
ifftdata122=ifftdata_user22';

%------------------Append Cyclic Prefix1 for User1-------------------------
y222=[ifftdata122(:,[(n-2):n]) ifftdata122];
transdata222=y222';
esttx_user22=transdata222;                % Transmitting data for user1

%------------------Generating data for User3-------------------------------
data3_i = (2*(rand(1,N)>0.5)-1);
data3_q = (2*(rand(1,N)>0.5)-1);
data_relay3_qpsk = data3_i + j*data3_q; %QPSK modulation

%-----------------Spreading & IFFT for User2-------------------------------
data_user31=data_relay3_qpsk';
spdata3=data_user31*code3;          % Spreading 
spdata33=(spdata3)';
ifftdata_user3=ifft(spdata33);      % Taking the IFFT
ifftdata33=ifftdata_user3';

%-----------------Append Cyclic Prefix1 for User2--------------------------
y3=[ifftdata33(:,[(n-2):n]) ifftdata33];
transdata3=y3';
tx_user3=transdata3;                % Transmitting data for user2
%----------------------Adding data for Transmission of All User------------
x3=esttx_user11+esttx_user22+tx_user3;
%----------------------Creating a new Rayleigh Channel---------------------------

gain1=sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
gain2=sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain3=sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3
gain4=sqrt(p4/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap4
x11=x3(:);
x12=reshape(x11,1,length(x11));
i=1:length(x12);        
delay1=1; 
for i=delay1+1:length(x12) % Producing one sample delay in Tap2 w.r.t. Tap1
   x13(i)=x3(i-delay1);
end
delay2=2;
for i=delay2+1:length(x12) % Producing two sample delay in Tap2 w.r.t. Tap1
   x14(i)=x3(i-delay2);
end
delay3=3;
for i=delay3+1:length(x12) % Producing three sample delay in Tap2 w.r.t. Tap1
   x15(i)=x3(i-delay3);
end
x111=reshape(x13,(n+3),length(x13)/(n+3));
x222=reshape(x14,(n+3),length(x14)/(n+3));
x333=reshape(x15,(n+3),length(x15)/(n+3));
ch1=repmat(gain1,(n+3),1);     
ch2=repmat(gain2,(n+3),1);
ch3=repmat(gain3,(n+3),1);
ch4=repmat(gain4,(n+3),1);
data_channel3=x3.*ch1+x111.*ch2+x222.*ch3+x333.*ch4;  % Passing data through channel 

%------------------------Addition of AWGN noise ---------------------------
data_noise3=data_channel3(:);
data_noise33=reshape(data_noise3,1,length(data_noise3));
noise33 = 1/sqrt(2)*[randn(1,length(data_noise33)) + j*randn(1,length(data_noise33))]; 

for i = 1:length(snr)
y3 = data_noise33 + (sqrt(1)*10^(-snr(i)/20))*noise33; %Addition of Noise
  
%--------------------------Receiver ---------------------------------------
data_received3 = y3;           %fadded data received with awgn noise
%---------------------Removing Cyclic Prefix-------------------------------
rx1=reshape(data_received3,(n+3),length(data_received3)/(n+3));
rx12=rx1';
rx13 = rx12(:,[(4:(n+3))]); 
rx14=rx13';
%-----------------Taking FFT ----------------------------------------------
fft_data_received =fft(rx14);

%----------------equilization of the channel-------------------------------
channel_response=fft([gain1;gain2;gain3;gain4],n);
data_equilized3=fft_data_received.*conj(channel_response);

%----------------BER of Data User1-----------------------------------------

recdata110=(data_equilized3'*code1')';
r_s11hop3_i =2*(real(recdata110)>0)-1;
r_s11hop3_q =2*(imag(recdata110)>0)-1;
errors_1hop3i(i) = size(find([data_i- r_s11hop3_i]),2); %Errors for User1
errors_1hop3q(i) = size(find([data_q- r_s11hop3_q]),2); %Errors for User1
SBer111 = (errors_1hop3i+errors_1hop3q)/(2*N);                              % simulated ber user1

%----------------BER of Data User2-----------------------------------------
recdata21=(data_equilized3'*code2')';
r_s22hop3_i =2*(real(recdata21)>0)-1;
r_s22hop3_q =2*(imag(recdata21)>0)-1;
errors_2hop3i(i) = size(find([data2_i- r_s22hop3_i]),2); %Errors for User1
errors_2hop3q(i) = size(find([data2_q- r_s22hop3_q]),2); %Errors for User1
SBer222= (errors_2hop3i+errors_2hop3q)/(2*N);                                  % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata31=(data_equilized3'*code3')';
r_s33hop3_i =2*(real(recdata31)>0)-1;
r_s33hop3_q =2*(imag(recdata31)>0)-1;
errors_3hop3i(i) = size(find([data3_i- r_s33hop3_i]),2); %Errors for User1
errors_3hop3q(i) = size(find([data3_q- r_s33hop3_q]),2); %Errors for User1
SBer333= (errors_3hop3i+errors_3hop3q)/(2*N);                               % simulated ber user2
end

% ------------------------Theoretical Result-------------------------------
snrlnr=10.^(snr/10);
TBerf = .5*erfc(sqrt(k1*snrlnr./(snrlnr+k1)));% theoretical BER fro Flat fadding

%-------------------Displaying Result--------------------------------------       
figure
grid on;
semilogy(snr,TBerf,'r-','LineWidth',2);
hold on;
semilogy(snr,SBer1,'kd-','LineWidth',2);
hold on;
semilogy(snr,SBer11,'kd--','LineWidth',2);
hold on;
semilogy(snr,SBer22,'go--','LineWidth',2);
hold on;
semilogy(snr,SBer111,'kd-','LineWidth',2);
hold on;
semilogy(snr,SBer222,'go-.','LineWidth',2);
hold on;
semilogy(snr,SBer333,'mo-.','LineWidth',2);
axis([0 20 10^-5 1]);
set(gca,'fontsize',14) % say...default here is 10
grid on
legend('Theoratical BER for BPSK on Rician Channel (urban)' ,'Simulated BER for User1 - hop 1','Simulated BER for User1 - hop 2','Simulated BER for User2 - hop 2','Simulated BER for User1 - hop 3','Simulated BER for User2 - hop 3','Simulated BER for User3 - hop 3','FontSize', 14);
xlabel('SNR, dB','FontSize', 18);
ylabel('Bit Error Rate','FontSize', 18);
title('BER Vs SNR on Rayleigh Channel','FontSize', 18)
