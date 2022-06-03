%ref:https://in.mathworks.com/matlabcentral/fileexchange/23498-mc-cdma
clc;close all;
clear all;
k1=5; %Rician factor %ref: https://www.researchgate.net/publication/263669548_Probability_Distribution_of_Rician_K-Factor_in_Urban_Suburban_and_Rural_Areas_Using_Real-World_Captured_Data
mean=sqrt(k1/(k1+1));% mean
sigma=sqrt(1/(2*(k1+1)));%variance
N=10^4;  % Number of Bits for  data_user1
d1 = 1; d2 = 3;    %Distances of users from base station (BS)
d3 = 3;
eta = 1;                %Path loss exponent
n =8;                               %Number of  Data Sub-Carriers
walsh=hadamard(n);  
code1=walsh(2,:);   
code3=walsh(3,:);                   %Taking 2nd row of walsh code for User1

%% ------------------Hop 1------------------------------
users=1;            % Number of Users
%------------------Generation of Walsh code--------------------------------
%------------------Generating data for User1-------------------------------                     
data_user1= rand(1,N)>0.5;          % Generation of data for user1
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
Taps=4;                                        % Number of Taps
p1=1;                                       % Power of Tap1
p2=1;                                     % Power of Tap2
p3=1;                                     % Power of Tap3
p4=1;
gain1=sqrt(d1^-eta)*sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
gain2=sqrt(d1^-eta)*sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain3=sqrt(d1^-eta)*sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3
gain4=sqrt(d1^-eta)*sqrt(p4/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap4
x11=x(:);
x12=reshape(x11,1,length(x11));
i=1:length(x12);        
delay1=1; 
for i=delay1+1:length(x12) % Producing one sample delay in Tap2 w.r.t. Tap1
   x13(i)=x(i-delay1);
end
delay2=2;
for i=delay2+1:length(x12) % Producing two sample delay in Tap2 w.r.t. Tap1
   x14(i)=x(i-delay2);
end
delay3=3;
for i=delay3+1:length(x12) % Producing three sample delay in Tap2 w.r.t. Tap1
   x15(i)=x(i-delay3);
end
x10=reshape(x13,(n+3),length(x13)/(n+3));
x20=reshape(x14,(n+3),length(x14)/(n+3));
x30=reshape(x15,(n+3),length(x15)/(n+3));
ch1=repmat(gain1,(n+3),1);     
ch2=repmat(gain2,(n+3),1);
ch3=repmat(gain3,(n+3),1);
ch4=repmat(gain4,(n+3),1);
data_channel1=x.*ch1+x10.*ch2+x20.*ch3+x30.*ch4;  % Passing data through channel 

%------------------------Addition of AWGN noise ---------------------------

%% ------user 2 using same code but different channel
                  %Taking 2nd row of walsh code for User1
%------------------Generating data for User1-------------------------------                          
data_user2 = rand(1,N)>0.5;          % Generation of data for user1
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

%----------------------pathloss------------
x2=2*sqrt(d2^-eta)*tx_user2;

%----------------------Creating Rayleigh Channel---------------------------

data_channel2=x2;  % Passing data through channel 
%------------------------Addition of AWGN noise ---------------------------
data_channel =data_channel1+ data_channel2;
data_noise=data_channel(:);
data_noise=reshape(data_noise,1,length(data_noise));
noise = 1/sqrt(2)*[randn(1,length(data_noise)) + j*randn(1,length(data_noise))]; 

%% ------ receiver ---------------------
snr = [0:100];                 % multiple Eb/N0 values

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
channel_response1=fft([gain1;gain2;gain3;gain4],n);
data_equilized1=fft_data_received2.*conj(channel_response1);

data_equilized2=fft_data_received2;

%----------------BER of Data User1-----------------------------------------
recdata112=(data_equilized2'*code1')';
x11_hat = ones(1,N);
x11_hat(recdata112<0) =-1;
recdata11=(data_equilized1'*code1')';
recdata12=real(recdata11)>0;
errors_user1(i) = size(find([data_user1- recdata12]),2); %Errors for User1
SBer1 = errors_user1/N;                              % simulated ber user1
                             
% simulated ber user2
% sic decoding
x12_hat = ones(1,N);
x12_hat(recdata11<0) =-1;
%remodulate decoded symbol x1 

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
data_user3 = rand(1,N)>0.5;          % Generation of data for user1
data_user3bpsk = 2*data_user3-1;    % BPSK modulation 0 -> -1; 1 -> 0 

%-----------------Spreading & IFFT for User2-------------------------------
data_user3bpskt=data_user3bpsk';
spdata3=data_user3bpskt*code2;          % Spreading 
spdata33=(spdata3)';
ifftdata_user3=ifft(spdata33);      % Taking the IFFT
ifftdata33=ifftdata_user3';

%-----------------Append Cyclic Prefix1 for User2--------------------------
y3=[ifftdata33(:,[(n-2):n]) ifftdata33];
transdata3=y3';

tx_user3=transdata3;                % Transmitting data for user2
%----------------------Adding data for Transmission of All User------------
x3=(esttx_user11+esttx_user12+tx_user3);
%----------------------Creating a new Rayleigh Channel---------------------------
p1=1;                                       % Power of Tap1
p2=1;                                     % Power of Tap2
p3=1;                                     % Power of Tap3
p4=1;
gain1ss=sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
gain2ss=sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain3ss=sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3
gain4ss=sqrt(p4/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap4
x111=x3(:);
x122=reshape(x111,1,length(x111));
i=1:length(x12);        
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
data_channel3=x3.*ch11+x11.*ch22+x22.*ch33+x33.*ch44;  % Passing data through channel 

%------------------------Addition of AWGN noise ---------------------------
data_noise3=data_channel3(:);
data_noise33=reshape(data_noise3,1,length(data_noise3));
noise33 = 1/sqrt(2)*[randn(1,length(data_noise33)) + j*randn(1,length(data_noise33))]; 
               
for i = 1:length(snr)
y3 = data_noise33 + (sqrt(1)*10^(-snr(i)/20))*noise33; %Addition of Noise
  
%--------------------------Receiver ---------------------------------------
data_received3 =y3;           %fadded data received with awgn noise
%---------------------Removing Cyclic Prefix-------------------------------
rx1=reshape(data_received3,(n+3),length(data_received3)/(n+3));
rx12=rx1';
rx13 = rx12(:,[(4:(n+3))]); 
rx14=rx13';
%-----------------Taking FFT ----------------------------------------------
fft_data_received =fft(rx14);
%----------------equilization of the channel-------------------------------
channel_response33=fft([gain1ss;gain2ss;gain3ss;gain4ss],n);
data_equilized22=fft_data_received.*conj(channel_response33);

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
recdata33=(data_equilized22'*code2')';
recdata33=real(recdata33)>0;
errors_user22(i) = size(find([data_user3- recdata33]),2); %Errors for User1
SBer22 = errors_user22/N;     

end

%% ------------------------Theoretical Result-------------------------------
snrlnr=10.^(snr/10);
TBerf = .5*erfc(sqrt(k1*snrlnr./(snrlnr+k1)));% theoretical BER fro Flat fadding

%% ------------------------Figures ------------------------

figure
grid on;
semilogy(snr,TBerf,'r-','LineWidth',2);
hold on;
semilogy(snr,SBer1,'gd-','LineWidth',2);%near 
hold on;
semilogy(snr,SBer2,'kd-','LineWidth',2);%far
hold on;
semilogy(snr,SBer11,'gd--','LineWidth',2);%near
hold on;
semilogy(snr,SBer12,'kd--','LineWidth',2);%far
hold on;
semilogy(snr,SBer22,'md-','LineWidth',2);%relay 1
axis([0 100 10^-8 1]);
set(gca,'fontsize',14) % say...default here is 10
grid on
legend('Theoratical BER for BPSK on Rician Channel (Urban)' ,'Simulated BER for near User - hop 1','Simulated BER for far User - hop 1','Simulated BER for near User - hop 2','Simulated BER for far User - hop 2','Simulated BER for relay User - hop 2','FontSize', 14);
xlabel('SNR, dB','FontSize', 18);
ylabel('Bit Error Rate','FontSize', 18);
title('BER Vs SNR on Rayleigh Channel','FontSize', 18)
