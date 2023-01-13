%ref:https://in.mathworks.com/matlabcentral/fileexchange/23498-mc-cdma
clc;
close all;
clear all;

k1 = 5;% Rician factor % ref: https://www.researchgate.net/publication/263669548_Probability_Distribution_of_Rician_K-Factor_in_Urban_Suburban_and_Rural_Areas_Using_Real-World_Captured_Data
mean = sqrt(k1/(k1+1));% mean
sigma = sqrt(1/(2*(k1+1)));% variance
N=10;  % Number of Bits for data_user1
d1 = 10;  d2 = 500;  % Distances of users from base station (BS)
d3 = 5;
eta = 4;            % Path loss exponent

%% ------------------Hop 1------------------------------
users=1;            % Number of Users

%------------------Generation of Walsh code--------------------------------
n =4;                              %Number of  Data Sub-Carriers
walsh=hadamard(n);              
code1=walsh(2,:);                   %Taking 2nd row of walsh code for User1

%------------------Generating data for User1-------------------------------                     
data_user1 = rand(1,N)>0.5;         % Generation of data for user1
data_user1bpsk = 2*data_user1-1;    % BPSK modulation 0 -> -1; 1 -> 0 

%------------------Spreading & IFFT for User1------------------------------
data_user11=data_user1bpsk';
spdata1_user1=data_user11*code1;    % Spreading 
spdata12=(spdata1_user1)';
ifftdata_user1=ifft(spdata12);      % Taking the IFFT
ifftdata12=ifftdata_user1';

%------------------Append Cyclic Prefix1 for User1-------------------------
y1=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata1=y1';
tx_user1=transdata1;                % Transmitting data for user1

%----------------------Adding data for Transmission of All User------------
x=tx_user1;
%----------------------Creating Rayleigh Channel---------------------------
Taps=4;                                   % Number of Taps
p1=1;                                     % Power of Tap1
p2=1;                                     % Power of Tap2
p3=1;                                     % Power of Tap3
p4=1;
gain1=sqrt(d1^-eta)*sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];% Gain for Tap1

ch1=repmat(gain1,(n+3),1);     

data_channel1=x.*ch1;  % Passing data through channel 

%% ------user 2 using same code but different channel
                  %Taking 2nd row of walsh code for User1
%------------------Generating data for User1-------------------------------                          
data_user2= rand(1,N)>0.5;          % Generation of data for user1
data_user2bpsk = 2*data_user2-1;    % BPSK modulation 0 -> -1; 1 -> 0 

%------------------Spreading & IFFT for User1------------------------------
data_user22=data_user2bpsk';
spdata1_user2=data_user22*code1;    % Spreading 
spdata122=(spdata1_user2)';
ifftdata_user2=ifft(spdata122);      % Taking the IFFT
ifftdata122=ifftdata_user2';

%------------------Append Cyclic Prefix1 for User1-------------------------
y2=[ifftdata122(:,[(n-2):n]) ifftdata122];
transdata2=y2';
tx_user2=transdata2;                % Transmitting data for user1

%----------------------Adding data for Transmission of All User------------
x2=tx_user2;

%----------------------Creating Rayleigh Channel---------------------------
Taps=4;                                        % Number of Taps
%p1=1;                                       % Power of Tap1
%p2=1;                                     % Power of Tap2
%p3=1;                                     % Power of Tap3
%p4=1;
gain1s=sqrt(d2^-eta)*sqrt(p1/2)*[(randn(1,N)*(sigma+mean)).^2 + j*(randn(1,N)*sigma).^2];   % Gain for Tap1
gain2s=sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain3s=sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3
gain4s=sqrt(p4/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap4
x110=x2(:);
x120=reshape(x110,1,length(x110));
i=1:length(x120);        
delay1=1; 
for i=delay1+1:length(x120) % Producing one sample delay in Tap2 w.r.t. Tap1
   x13(i)=x2(i-delay1);
end
delay2=2;
for i=delay2+1:length(x120) % Producing two sample delay in Tap2 w.r.t. Tap1
   x14(i)=x2(i-delay2);
end
delay3=3;
for i=delay3+1:length(x120) % Producing three sample delay in Tap2 w.r.t. Tap1
   x15(i)=x2(i-delay3);
end
x111=reshape(x13,(n+3),length(x13)/(n+3));
x222=reshape(x14,(n+3),length(x14)/(n+3));
x333=reshape(x15,(n+3),length(x15)/(n+3));
ch11=repmat(gain1s,(n+3),1);     
ch22=repmat(gain2s,(n+3),1);
ch33=repmat(gain3s,(n+3),1);
ch44=repmat(gain4s,(n+3),1);
data_channel2=x2.*ch11+x111.*ch22+x222.*ch33+x333.*ch44;  % Passing data through channel 
%------------------------Addition of AWGN noise ---------------------------
data_channel=data_channel1+ data_channel2;
data_noise=data_channel(:);
data_noise=reshape(data_noise,1,length(data_noise));
noise = 1/sqrt(2)*[randn(1,length(data_noise)) + j*randn(1,length(data_noise))]; 

%% ------ receiver ---------------------
snr = [0:20];                 % multiple Eb/N0 values

for i = 1:length(snr)
y2= data_noise + (sqrt(1)*10^(-snr(i)/20))*noise; %Addition of Noise
  
%--------------------------Receiver ---------------------------------------
data_received2 = y2;           %fadded data received with awgn noise
%---------------------Removing Cyclic Prefix-------------------------------
rx1=reshape(data_received2,(n+3),length(data_received2)/(n+3));
rx12=rx1';
rx13 = rx12(:,[(4:(n+3))]); 
rx14=rx13';

%-----------------Taking FFT ----------------------------------------------
fft_data_received2 =fft(rx14);

%----------------equilization of the channel-------------------------------

data_equilized1=fft_data_received2/sqrt(d1^-eta)*sqrt(p1/2);
channel_response2=fft([gain1s;gain2s;gain3s;gain4s],n);
data_equilized2=fft_data_received2.*conj(channel_response2);

%----------------BER of Data User1-----------------------------------------
recdata112=(data_equilized2'*code1')';
x11_hat = ones(1,N);
x11_hat(recdata112<0) = -1;
recdata11=(data_equilized1'*code1'-x11_hat')';
recdata12=real(recdata11)>0;
errors_user1(i) = size(find([data_user1- recdata12]),2); %Errors for User1
SBer1 = errors_user1/N;                              % simulated ber user1
                             
% simulated ber user2
% sic decoding
x12_hat = ones(1,N);
x12_hat(recdata11<0) = -1;
%recdata22=((data_equilized2'*code1') - x12_hat')';
recdata22=(data_equilized2'*code1' - x12_hat')';
recdata122=real(recdata22)>0;
errors_user2(i) = size(find([data_user2- recdata122]),2); %Errors for User1
SBer2 = errors_user2/N;                              % simulated ber user1

end

%% ------------------Hop 2------------------------------

users=2;            % Number of Users
%------------------Generation of Walsh code--------------------------------
code2=walsh(4,:);   
code3=walsh(3,:); %Taking 3rd row of walsh code for User2
%------------------est for User11-------------------------------

%------------------Spreading & IFFT for User1------------------------------
estdata_user11bpsk = 2*recdata12-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user11=estdata_user11bpsk';
spdata1_user11=estdata_user11*code1;    % Spreading 
spdata11=(spdata1_user11)';
ifftdata_user11=ifft(spdata11);      % Taking the IFFT
ifftdata11=ifftdata_user11';
%------------------Append Cyclic Prefix1 for User1-------------------------
y11=[ifftdata11(:,[(n-2):n]) ifftdata11];
transdata11=y11';
esttx_user11=transdata11;                % Transmitting data for user1

%------------------est for User12-------------------------------

%------------------Spreading & IFFT for User1------------------------------
estdata_user12bpsk = 2*recdata122-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user12=estdata_user12bpsk';
spdata1_user12=estdata_user12*code3;    % Spreading 
spdata12=(spdata1_user12)';
ifftdata_user12=ifft(spdata12);      % Taking the IFFT
ifftdata12=ifftdata_user12';
%------------------Append Cyclic Prefix1 for User1-------------------------
y12=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata12=y12';
esttx_user12=transdata12;                % Transmitting data for user1

%------------------Generating data for User2-------------------------------
M=N;                             % Number of Bits for  data_user2
data_user3= rand(1,M)>0.5;          % Generation of data for user2
data_user2bpsk = 2*data_user3-1;    % BPSK modulation 0 -> -1; 1 -> 0 
%-----------------Spreading & IFFT for User2-------------------------------
data_user21=data_user2bpsk';
spdata2=data_user21*code2;          % Spreading 
spdata22=(spdata2)';
ifftdata_user2=ifft(spdata22);      % Taking the IFFT
ifftdata22=ifftdata_user2';
%-----------------Append Cyclic Prefix1 for User2--------------------------
y2=[ifftdata22(:,[(n-2):n]) ifftdata22];
transdata2=y2';
tx_user2=transdata2;                % Transmitting data for user2
%----------------------Adding data for Transmission of All User------------
x3=esttx_user11+esttx_user12+tx_user2;
%----------------------Creating a new Rayleigh Channel---------------------------
p1=1;                                       % Power of Tap1
p2=1;                                     % Power of Tap2
p3=1;                                     % Power of Tap3
p4=1;
gain1ss=sqrt(d3^-eta)*sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
gain2ss=sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain3ss=sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3
gain4ss=sqrt(p4/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap4
x111=x3(:);
x122=reshape(x111,1,length(x111));
i=1:length(x122);        
delay1=1; 
for i=delay1+1:length(x122) % Producing one sample delay in Tap2 w.r.t. Tap1
   x133(i)=x3(i-delay1);
end
delay2=2;
for i=delay2+1:length(x122) % Producing two sample delay in Tap2 w.r.t. Tap1
   x144(i)=x3(i-delay2);
end
delay3=3;
for i=delay3+1:length(x122) % Producing three sample delay in Tap2 w.r.t. Tap1
   x155(i)=x3(i-delay3);
end
x11=reshape(x133,(n+3),length(x133)/(n+3));
x22=reshape(x144,(n+3),length(x144)/(n+3));
x33=reshape(x155,(n+3),length(x155)/(n+3));
ch11=repmat(gain1ss,(n+3),1);     
ch22=repmat(gain2ss,(n+3),1);
ch33=repmat(gain3ss,(n+3),1);
ch44=repmat(gain4ss,(n+3),1);
data_channel2=x3.*ch11+x11.*ch22+x22.*ch33+x33.*ch44;  % Passing data through channel 

%------------------------Addition of AWGN noise ---------------------------
data_noise2=data_channel2(:);
data_noise22=reshape(data_noise2,1,length(data_noise2));
noise22 = 1/sqrt(2)*[randn(1,length(data_noise22)) + j*randn(1,length(data_noise22))]; 
               
for i = 1:length(snr)
y2 = data_noise22 + (sqrt(1)*10^(-snr(i)/20))*noise22; %Addition of Noise
  
%--------------------------Receiver ---------------------------------------
data_received2 =y2;           %fadded data received with awgn noise
%---------------------Removing Cyclic Prefix-------------------------------
rx1=reshape(data_received2,(n+3),length(data_received2)/(n+3));
rx12=rx1';
rx13 = rx12(:,[(4:(n+3))]); 
rx14=rx13';
%-----------------Taking FFT ----------------------------------------------
fft_data_received =fft(rx14);
%----------------equilization of the channel-------------------------------
channel_response22=fft([gain1ss;gain2ss;gain3ss;gain4ss],n);
data_equilized22=fft_data_received.*conj(channel_response22);
%----------------BER of Data User1-----------------------------------------

recdata11=(data_equilized22'*code1')';
recdata120=real(recdata11)>0;
errors_user11(i) = size(find([data_user1- recdata120]),2); %Errors for User1
SBer11 = errors_user11/N;                              % simulated ber user1
%----------------BER of Data source User2-----------------------------------------

x12_hat = ones(1,N);
x12_hat(recdata11<0) = -1;
recdata12=(data_equilized22'*code3')';
recdata122=real(recdata12)>0;
errors_user12(i) = size(find([data_user2- recdata122]),2); %Errors for User1
SBer12 = errors_user12/N;                              % simulated ber user1

%----------------BER of Data relay User2-----------------------------------------
recdata21=(data_equilized22'*code2')';
recdata22=real(recdata21)>0;
errors_user22(i) = size(find([data_user3- recdata22]),2); %Errors for User1
SBer22 = errors_user22/M;     

end

%% ------------------Hop 3------------------------------

users=4;            % Number of Users
%------------------Generation of Walsh code--------------------------------
code4=walsh(1,:);                   %Taking 2nd row of walsh code for User1

%------------------est for User1-------------------------------

%------------------Spreading & IFFT for User1------------------------------
estdata_user1bpskhop3 = 2*recdata120-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user111=estdata_user1bpskhop3';
spdata1_user111=estdata_user111*code1;    % Spreading 
spdata121=(spdata1_user111)';
ifftdata_user111=ifft(spdata121);      % Taking the IFFT
ifftdata12=ifftdata_user111';

%------------------Append Cyclic Prefix1 for User1-------------------------
y11=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata11=y11';
esttx_user11=transdata11;                % Transmitting data for user1

%------------------est for User2-------------------------------

%------------------Spreading & IFFT for User1------------------------------
estdata_user2bpskhop3 = 2*recdata122-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user222=estdata_user2bpskhop3';
spdata1_user222=estdata_user222*code3;    % Spreading 
spdata222=(spdata1_user222)';
ifftdata_user22=ifft(spdata222);      % Taking the IFFT
ifftdata122=ifftdata_user22';

%------------------Append Cyclic Prefix1 for User1-------------------------
y222=[ifftdata122(:,[(n-2):n]) ifftdata122];
transdata222=y222';
esttx_user22=transdata222;                % Transmitting data for user1

%------------------est for User3-------------------------------

%------------------Spreading & IFFT for User1------------------------------
estdata_user3bpskhop3 = 2*recdata22-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user333=estdata_user3bpskhop3';
spdata1_user333=estdata_user333*code2;    % Spreading 
spdata333=(spdata1_user333)';
ifftdata_user33=ifft(spdata333);      % Taking the IFFT
ifftdata33=ifftdata_user33';

%------------------Append Cyclic Prefix1 for User1-------------------------
y333=[ifftdata33(:,[(n-2):n]) ifftdata33];
transdata333=y333';
esttx_user33=transdata333;                % Transmitting data for user1

%------------------Generating data for User4-------------------------------
M=N;                             % Number of Bits for  data_user2
data_user4= rand(1,M)>0.5;          % Generation of data for user2
data_user4bpsk = 2*data_user4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
%-----------------Spreading & IFFT for User2-------------------------------
data_user44=data_user4bpsk';
spdata4=data_user44*code4;          % Spreading 
spdata44=(spdata4)';
ifftdata_user4=ifft(spdata44);      % Taking the IFFT
ifftdata44=ifftdata_user4';
%-----------------Append Cyclic Prefix1 for User2--------------------------
y4=[ifftdata44(:,[(n-2):n]) ifftdata44];
transdata4=y4';
tx_user4=transdata4;                % Transmitting data for user2
%----------------------Adding data for Transmission of All User------------
x4=esttx_user11+esttx_user22+esttx_user33+tx_user4;
%----------------------Creating a new Rayleigh Channel---------------------------

gain1sss=sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
gain2sss=sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain3sss=sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3
gain4sss=sqrt(p4/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap4
x41=x4(:);
x42=reshape(x41,1,length(x41));
i=1:length(x42);        
delay1=1; 
for i=delay1+1:length(x42) % Producing one sample delay in Tap2 w.r.t. Tap1
   x43(i)=x4(i-delay1);
end
delay2=2;
for i=delay2+1:length(x42) % Producing two sample delay in Tap2 w.r.t. Tap1
   x44(i)=x4(i-delay2);
end
delay3=3;
for i=delay3+1:length(x42) % Producing three sample delay in Tap2 w.r.t. Tap1
   x45(i)=x4(i-delay3);
end
x1111=reshape(x43,(n+3),length(x43)/(n+3));
x2222=reshape(x44,(n+3),length(x44)/(n+3));
x3333=reshape(x45,(n+3),length(x45)/(n+3));
ch41=repmat(gain1sss,(n+3),1);     
ch42=repmat(gain2sss,(n+3),1);
ch43=repmat(gain3sss,(n+3),1);
ch44=repmat(gain4sss,(n+3),1);
data_channel4=x4.*ch41+x1111.*ch42+x2222.*ch43+x3333.*ch44;  % Passing data through channel 

%------------------------Addition of AWGN noise ---------------------------
data_noise4=data_channel4(:);
data_noise44=reshape(data_noise4,1,length(data_noise4));
noise44 = 1/sqrt(2)*[randn(1,length(data_noise44)) + j*randn(1,length(data_noise44))]; 

for i = 1:length(snr)
y4 = data_noise44 + (sqrt(1)*10^(-snr(i)/20))*noise44; %Addition of Noise
  
%--------------------------Receiver ---------------------------------------
data_received4 = y4;           %fadded data received with awgn noise
%---------------------Removing Cyclic Prefix-------------------------------
rx1=reshape(data_received4,(n+3),length(data_received4)/(n+3));
rx12=rx1';
rx13 = rx12(:,[(4:(n+3))]); 
rx14=rx13';
%-----------------Taking FFT ----------------------------------------------
fft_data_received4 =fft(rx14);

%----------------equilization of the channel-------------------------------
channel_response4=fft([gain1sss;gain2sss;gain3sss;gain4sss],n);
data_equilized4=fft_data_received4.*conj(channel_response4);

%----------------BER of Data User1-----------------------------------------
recdata110=(data_equilized4'*code1')';
recdata1200=real(recdata110)>0;
errors_user111(i) = size(find([data_user1- recdata1200]),2); %Errors for User1
SBer111 = errors_user111/N;                              % simulated ber user1

%----------------BER of Data User2-----------------------------------------
recdata2221=(data_equilized4'*code3')';
recdata222=real(recdata2221)>0;
errors_user222(i) = size(find([data_user2- recdata222]),2); %Errors for User1
SBer222 = errors_user222/M;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata3=(data_equilized4'*code2')';
recdata333=real(recdata3)>0;
errors_user333(i) = size(find([data_user3- recdata333]),2); %Errors for User1
SBer333 = errors_user333/M;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata4=(data_equilized4'*code4')';
recdata444=real(recdata4)>0;
errors_user444(i) = size(find([data_user4- recdata444]),2); %Errors for User1
SBer444 = errors_user444/M;                               % simulated ber user2
end

%% ------------------------Theoretical Result-------------------------------
snrlnr=10.^(snr/10);
TBerf = .5*erfc(sqrt(k1*snrlnr./(snrlnr+k1)));% theoretical BER fro Flat fadding
%% ------------------------Figures ------------------------
figure
grid on;
semilogy(snr,TBerf,'r-','LineWidth',2);
hold on;
semilogy(snr,SBer1,'kd-','LineWidth',2);%far 
hold on;
semilogy(snr,SBer2,'bd-','LineWidth',2);%near
hold on;
semilogy(snr,SBer11,'kd--','LineWidth',2);%near2
hold on;
semilogy(snr,SBer12','bd--','LineWidth',2);%far2
hold on;
semilogy(snr,SBer22,'gd-','LineWidth',2);%%relay 1 hop2
hold on;
semilogy(snr,SBer111 ,'kd-.','LineWidth',2);%near 1 hop3
hold on;
semilogy(snr,SBer222 ,'bd-.','LineWidth',2);%far 1 hop3
hold on;
semilogy(snr,SBer333 ,'gd--','LineWidth',2);%relay 1 hop3
hold on;
semilogy(snr,SBer444 ,'md-','LineWidth',2);%relay 2 hop 3
axis([0 20 10^-4 1]);
set(gca,'fontsize',14) % say...default here is 10
grid on
legend('Theoratical BER for BPSK on Rician Channel (Urban)' ,'Simulated BER for near User - hop 1','Simulated BER for far User - hop 1','Simulated BER for near User - hop 2','Simulated BER for far User - hop 2','Simulated BER for relay1 User - hop 2','Simulated BER for near User - hop 3','Simulated BER for far User - hop 3','Simulated BER for relay1 User - hop 3','Simulated BER for relay2 User - hop 3','FontSize', 12);
xlabel('SNR, dB','FontSize', 18);
ylabel('Bit Error Rate','FontSize', 18);
title('BER Vs SNR on Rayleigh Channel','FontSize', 18)
