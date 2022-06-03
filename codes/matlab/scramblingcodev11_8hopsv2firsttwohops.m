%ref:https://in.mathworks.com/matlabcentral/fileexchange/23498-mc-cdma
clc;
close all;
clear all;
k1=5; %Rician factor %ref: https://www.researchgate.net/publication/263669548_Probability_Distribution_of_Rician_K-Factor_in_Urban_Suburban_and_Rural_Areas_Using_Real-World_Captured_Data
mean=sqrt(k1/(k1+1));% mean
sigma=sqrt(1/(2*(k1+1)));%variance
N=10^4;  % Number of Bits for data_user1
d1 = 0.8; d2 = 500;    %Distances of users from base station (BS)
d3 = 5;d4 = 5;d5 = 5;d6 = 5;d7 = 5;d8 = 5;d9 = 5;dA=0.8;
eta = 4;            %Path loss exponent

%% ------------------Hop 1------------------------------
users=1;            % Number of Users
%------------------Generation of Walsh code--------------------------------
n =48;                               %Number of  Data Sub-Carriers
walsh=hadamard(n);              
code1=walsh(2,:);                   %Taking 2nd row of walsh code for User1
code4=walsh(1,:);                   %Taking 2nd row of walsh code for User1
code2=walsh(4,:);   
code3=walsh(3,:); %Taking 3rd row of walsh code for User2
code5=walsh(5,:); 
code6=walsh(6,:); 
code7=walsh(7,:); 
code8=walsh(8,:); 
code9=walsh(9,:); 
codeA = walsh(10,:); 
codeB = walsh(11,:); 
codeC = walsh(12,:); 
codeD = walsh(13,:); 
codeE = walsh(14,:); 
codeF = walsh(15,:); 

%------------------Generating data for User1-------------------------------                     
data_user1 = rand(1,N)>0.5;          % Generation of data for user1
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
ch1=repmat(gain1,(n+3),1);     
data_channel1=sqrt(d1^-eta)*sqrt(p1/2)*x;  % Passing data through channel 

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
gain1s=sqrt(d2^-eta)*sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
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
data_channel = data_channel1+ data_channel2;
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
data_equilized1=fft_data_received2;
%/sqrt(d1^-eta)*sqrt(p1/2);
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
data_user21=data_user2bpsk';
spdata2=data_user21*code2;          % Spreading 
spdata22=(spdata2)';
ifftdata_user2=ifft(spdata22);      % Taking the IFFT
ifftdata22=ifftdata_user2';
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

%% ------user A using same code (code 2) but different channel

%original data for A                            % Number of Bits for  data_user2
data_userA= rand(1,M)>0.5;          % Generation of data for user2
data_userAbpsk = 2*data_userA-1;    % BPSK modulation 0 -> -1; 1 -> 0 
data_userAh2=data_userAbpsk';
spdataAh2=data_userAh2*code2;          % Spreading 
spdataAh2t=(spdataAh2)';
ifftdata_userAh2=ifft(spdataAh2t);      % Taking the IFFT
ifftdataAh2=ifftdata_userAh2';
yA=[ifftdataAh2(:,[(n-2):n]) ifftdataAh2];
transdataAh2=yA';
tx_userA=transdataAh2; 

%----------------------Adding data for Transmission of All User------------
xA=tx_userA;
%----------------------Creating a new Rayleigh Channel---------------------------
Taps=4;                                        % Number of Taps
data_channelA=sqrt(dA^-eta)*sqrt(p1/2)*xA;  % Passing data through channel 

data_channelhop2 = data_channel2+data_channelA;

%------------------------Addition of AWGN noise ---------------------------
data_noise2=data_channelhop2(:);
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
fft_data_receivedhop2 =fft(rx14);
%----------------equilization of the channel-------------------------------
channel_response22=fft([gain1ss;gain2ss;gain3ss;gain4ss],n);
data_equilized22=fft_data_receivedhop2.*conj(channel_response22);

data_equilizedA=fft_data_receivedhop2/sqrt(dA^-eta)*sqrt(p1/2);

%----------------BER of Data User1-----------------------------------------
recdata11=(data_equilized22'*code1')';
recdata120=real(recdata11)>0;
errors_user11(i) = size(find([data_user1- recdata120]),2); %Errors for User1
SBer11 = errors_user11/N;                              % simulated ber user1

%----------------BER of Data source User2-----------------------------------------
recdata12=(data_equilized22'*code3')';
recdata122=real(recdata12)>0;
errors_user12(i) = size(find([data_user2- recdata122]),2); %Errors for User1
SBer12 = errors_user12/N;                              % simulated ber user1

%----------------BER of Data relay User3-----------------------------------------
recdata3A=(data_equilizedA'*code2')';
xA_hat = ones(1,N);
xA_hat(recdata3A<0) = -1;
recdata3user=(data_equilized22'*code2'-xA_hat')';
estrecdata3=real(recdata3user)>0;
errors_user3(i) = size(find([data_user3- estrecdata3]),2); %Errors for User1
SBer3_user = errors_user3/M;     

%----------------BER of Data UserA-----------------------------------------

x3user_hat = ones(1,N);
x3user_hat(recdata3user<0) = -1;
recdataA=(data_equilizedA'*code2'-x3user_hat')';
estrecdataA=real(recdataA)>0;
errors_userA(i) = size(find([data_userA- estrecdataA]),2); %Errors for User1
SBerA = errors_userA/N;                              % simulated ber user1
                             
end

%% ------------------Hop 3------------------------------

users=4;            % Number of Users
%------------------Generation of Walsh code--------------------------------

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

gain1sss=sqrt(d4^-eta)*sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
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
%% ------------------Hop 4------------------------------
users=5;    

%user 1 estimated data %modulation 
estdata_user1bpskhop4 = 2*recdata1200-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user111=estdata_user1bpskhop4';
spdata1_user111=estdata_user111*code1;    % Spreading 
spdata121=(spdata1_user111)';
ifftdata_user111=ifft(spdata121);      % Taking the IFFT
ifftdata12=ifftdata_user111';
y11=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata11=y11';
esttx_user1hop4=transdata11;

%user 2 estimated data %modulation 
estdata_user2bpskhop4 = 2*recdata222-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user222=estdata_user2bpskhop4';
spdata1_user222=estdata_user222*code3;    % Spreading 
spdata222=(spdata1_user222)';
ifftdata_user22=ifft(spdata222);      % Taking the IFFT
ifftdata122=ifftdata_user22';
y222=[ifftdata122(:,[(n-2):n]) ifftdata122];
transdata222=y222';
esttx_user2hop4=transdata222;               

%user 3 estimated data %modulation 
estdata_user3bpskhop4 = 2*recdata333-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user333=estdata_user3bpskhop4';
spdata1_user333=estdata_user333*code2;    % Spreading 
spdata333=(spdata1_user333)';
ifftdata_user33=ifft(spdata333);      % Taking the IFFT
ifftdata33=ifftdata_user33';
y333=[ifftdata33(:,[(n-2):n]) ifftdata33];
transdata333=y333';
esttx_user3hop4=transdata333;                % Transmitting data for user1

%user 4 estimated data %modulation 
estdata_user4bpskhop4 = 2*recdata444-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user444=estdata_user4bpskhop4';
spdata1_user444=estdata_user444*code4;    % Spreading 
spdata444=(spdata1_user444)';
ifftdata_user44=ifft(spdata444);      % Taking the IFFT
ifftdata44=ifftdata_user44';
y444=[ifftdata44(:,[(n-2):n]) ifftdata44];
transdata444=y444';
esttx_user4hop4=transdata444;                % Transmitting data for user1

%user 5 original data %modulation 
data_user5= rand(1,N)>0.5;          % Generation of data for user2
data_user5bpsk = 2*data_user5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
data_user55=data_user5bpsk';
spdata5=data_user55*code5;          % Spreading 
spdata55=(spdata5)';
ifftdata_user5=ifft(spdata55);      % Taking the IFFT
ifftdata55=ifftdata_user5';
y5=[ifftdata55(:,[(n-2):n]) ifftdata55];
transdata5=y5';
tx_user5=transdata5;                % Transmitting data for user2

x5=esttx_user1hop4+esttx_user2hop4+esttx_user3hop4+esttx_user4hop4+tx_user5;

%rician fading channel 
%----------------------Creating a new Rayleigh Channel---------------------------
gain15=sqrt(d5^-eta)*sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
gain25=sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain35=sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3
gain45=sqrt(p4/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap4
x41=x5(:);
x42=reshape(x41,1,length(x41));
i=1:length(x42);        
delay1=1; 
for i=delay1+1:length(x42) % Producing one sample delay in Tap2 w.r.t. Tap1
   x43(i)=x5(i-delay1);
end
delay2=2;
for i=delay2+1:length(x42) % Producing two sample delay in Tap2 w.r.t. Tap1
   x44(i)=x5(i-delay2);
end
delay3=3;
for i=delay3+1:length(x42) % Producing three sample delay in Tap2 w.r.t. Tap1
   x45(i)=x5(i-delay3);
end
x1111=reshape(x43,(n+3),length(x43)/(n+3));
x2222=reshape(x44,(n+3),length(x44)/(n+3));
x3333=reshape(x45,(n+3),length(x45)/(n+3));
ch41=repmat(gain15,(n+3),1);     
ch42=repmat(gain25,(n+3),1);
ch43=repmat(gain35,(n+3),1);
ch44=repmat(gain45,(n+3),1);
data_channel5=x5.*ch41+x1111.*ch42+x2222.*ch43+x3333.*ch44;  % Passing data through channel 

%addition of noise 
data_noise5=data_channel5(:);
data_noise55=reshape(data_noise5,1,length(data_noise5));
noise55 = 1/sqrt(2)*[randn(1,length(data_noise55)) + j*randn(1,length(data_noise55))];

%ifft
%equalization 
%decoding
%demodulation
%ber rate calculation
for i = 1:length(snr)
y5 = data_noise55 + (sqrt(1)*10^(-snr(i)/20))*noise55; %Addition of Noise
  
%--------------------------Receiver ---------------------------------------
data_received5 = y5;           %fadded data received with awgn noise
%---------------------Removing Cyclic Prefix-------------------------------
rx1=reshape(data_received5,(n+3),length(data_received5)/(n+3));
rx12=rx1';
rx13 = rx12(:,[(4:(n+3))]); 
rx14=rx13';

%-----------------Taking FFT ----------------------------------------------
fft_data_received5 =fft(rx14);

%----------------equilization of the channel-------------------------------
channel_response5=fft([gain15;gain25;gain35;gain45],n);
data_equilized5=fft_data_received5.*conj(channel_response5);

%----------------BER of Data User1-----------------------------------------
recdata110=(data_equilized5'*code1')';
recdata1h4=real(recdata110)>0;
errors_user1h4(i) = size(find([data_user1- recdata1h4]),2); %Errors for User1
SBer1h4 = errors_user1h4/N;                              % simulated ber user1

%----------------BER of Data User2-----------------------------------------
recdata2221=(data_equilized5'*code3')';
recdata2h4=real(recdata2221)>0;
errors_user2h4(i) = size(find([data_user2- recdata2h4]),2); %Errors for User1
SBer2h4 = errors_user2h4/N;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata3=(data_equilized5'*code2')';
recdata3h4=real(recdata3)>0;
errors_user3h4(i) = size(find([data_user3- recdata3h4]),2); %Errors for User1
SBer3h5 = errors_user3h4/N;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata4=(data_equilized5'*code4')';
recdata4h4=real(recdata4)>0;
errors_user4h4(i) = size(find([data_user4- recdata4h4]),2); %Errors for User1
SBer4h5 = errors_user4h4/N;      

%----------------BER of Data User3-----------------------------------------
recdata5=(data_equilized5'*code5')';
recdata5h4=real(recdata5)>0;
errors_user5h4(i) = size(find([data_user5- recdata5h4]),2); %Errors for User1
SBer5h5= errors_user5h4/N;    
end

%% ------------------Hop 5------------------------------
users=6; 

%user 1 estimated data %modulation
estdata_user1bpskhop5 = 2*recdata1h4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user111=estdata_user1bpskhop5';
spdata1_user111=estdata_user111*code1;    % Spreading 
spdata121=(spdata1_user111)';
ifftdata_user111=ifft(spdata121);      % Taking the IFFT
ifftdata12=ifftdata_user111';
y11=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata11=y11';
esttx_user1hop5=transdata11;

%user 2 estimated data %modulation 
estdata_user2bpskhop5 = 2*recdata2h4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user222=estdata_user2bpskhop5';
spdata1_user222=estdata_user222*code3;    % Spreading 
spdata222=(spdata1_user222)';
ifftdata_user22=ifft(spdata222);      % Taking the IFFT
ifftdata122=ifftdata_user22';
y222=[ifftdata122(:,[(n-2):n]) ifftdata122];
transdata222=y222';
esttx_user2hop5=transdata222;          

%user 3 estimated data %modulation 
estdata_user3bpskhop5 = 2*recdata3h4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user3=estdata_user3bpskhop5';
spdata1_user3=estdata_user3*code2;    % Spreading 
spdata3=(spdata1_user3)';
ifftdata_user33=ifft(spdata3);      % Taking the IFFT
ifftdata3=ifftdata_user33';
y3h5=[ifftdata3(:,[(n-2):n]) ifftdata3];
transdata333=y3h5';
esttx_user3hop5=transdata333;          

%user 4 estimated data %modulation 
estdata_user4bpskhop5 = 2*recdata4h4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user4=estdata_user4bpskhop5';
spdata1_user4=estdata_user4*code4;    % Spreading 
spdata4=(spdata1_user4)';
ifftdata_user44=ifft(spdata4);      % Taking the IFFT
ifftdata4=ifftdata_user44';
y4h5=[ifftdata4(:,[(n-2):n]) ifftdata4];
transdata444=y4h5';
esttx_user4hop5=transdata444;      

%user 5 estimated data %modulation 
estdata_user5bpskhop5 = 2*recdata5h4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user5=estdata_user5bpskhop5';
spdata1_user5=estdata_user5*code5;    % Spreading 
spdata5=(spdata1_user5)';
ifftdata_user55=ifft(spdata5);      % Taking the IFFT
ifftdata5=ifftdata_user55';
y5h5=[ifftdata5(:,[(n-2):n]) ifftdata5];
transdata555=y5h5';
esttx_user5hop5=transdata555;  

%user 6 original data %modulation 
data_user6= rand(1,N)>0.5;          % Generation of data for user2
data_user6bpsk = 2*data_user6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
data_user66=data_user6bpsk';
spdata6=data_user66*code6;          % Spreading 
spdata66=(spdata6)';
ifftdata_user6=ifft(spdata66);      % Taking the IFFT
ifftdata66=ifftdata_user6';
y6=[ifftdata66(:,[(n-2):n]) ifftdata66];
transdata6=y6';
tx_user6=transdata6;                

% Transmitting data of all users
x6 = esttx_user1hop5+esttx_user2hop5+esttx_user3hop5+esttx_user4hop5+esttx_user5hop5+transdata6;

%rician fading channel 
%----------------------Creating a new Rayleigh Channel---------------------------
gain16=sqrt(d6^-eta)*sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
gain26=sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain36=sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3
gain46=sqrt(p4/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap4
x41=x6(:);
x42=reshape(x41,1,length(x41));
i=1:length(x42);        
delay1=1; 
for i=delay1+1:length(x42) % Producing one sample delay in Tap2 w.r.t. Tap1
   x43(i)=x6(i-delay1);
end
delay2=2;
for i=delay2+1:length(x42) % Producing two sample delay in Tap2 w.r.t. Tap1
   x44(i)=x6(i-delay2);
end
delay3=3;
for i=delay3+1:length(x42) % Producing three sample delay in Tap2 w.r.t. Tap1
   x45(i)=x6(i-delay3);
end
x1111=reshape(x43,(n+3),length(x43)/(n+3));
x2222=reshape(x44,(n+3),length(x44)/(n+3));
x3333=reshape(x45,(n+3),length(x45)/(n+3));
ch41=repmat(gain16,(n+3),1);     
ch42=repmat(gain26,(n+3),1);
ch43=repmat(gain36,(n+3),1);
ch44=repmat(gain46,(n+3),1);
data_channel6=x6.*ch41+x1111.*ch42+x2222.*ch43+x3333.*ch44;  % Passing data through channel 

%addition of noise 
data_noise6=data_channel6(:);
data_noise66=reshape(data_noise6,1,length(data_noise6));
noise66 = 1/sqrt(2)*[randn(1,length(data_noise66)) + j*randn(1,length(data_noise66))];

%ifft
%equalization 
%decoding
%demodulation
%ber rate calculation
for i = 1:length(snr)
y6 = data_noise66 + (sqrt(1)*10^(-snr(i)/20))*noise66; %Addition of Noise
  
%--------------------------Receiver ---------------------------------------
data_received6 = y6;           %fadded data received with awgn noise
%---------------------Removing Cyclic Prefix-------------------------------
rx1=reshape(data_received6,(n+3),length(data_received6)/(n+3));
rx12=rx1';
rx13 = rx12(:,[(4:(n+3))]); 
rx14=rx13';

%-----------------Taking FFT ----------------------------------------------
fft_data_received6 =fft(rx14);

%----------------equilization of the channel-------------------------------
channel_response6=fft([gain16;gain26;gain36;gain46],n);
data_equilized6=fft_data_received6.*conj(channel_response6);

%----------------BER of Data User1-----------------------------------------
recdata110=(data_equilized6'*code1')';
recdata1h5=real(recdata110)>0;
errors_user1h5(i) = size(find([data_user1- recdata1h5]),2); %Errors for User1
SBer1h5 = errors_user1h5/N;                              % simulated ber user1

%----------------BER of Data User2-----------------------------------------
recdata2221=(data_equilized6'*code3')';
recdata2h5=real(recdata2221)>0;
errors_user2h5(i) = size(find([data_user2- recdata2h5]),2); %Errors for User1
SBer2h5 = errors_user2h5/N;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata3=(data_equilized6'*code2')';
recdata3h5=real(recdata3)>0;
errors_user3h5(i) = size(find([data_user3- recdata3h5]),2); %Errors for User1
SBer3h5 = errors_user3h5/N;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata4=(data_equilized6'*code4')';
recdata4h5=real(recdata4)>0;
errors_user4h5(i) = size(find([data_user4- recdata4h5]),2); %Errors for User1
SBer4h5 = errors_user4h5/N;      

%----------------BER of Data User3-----------------------------------------
recdata5=(data_equilized6'*code5')';
recdata5h5=real(recdata5)>0;
errors_user5h5(i) = size(find([data_user5- recdata5h5]),2); %Errors for User1
SBer5h5= errors_user5h5/N; 

%----------------BER of Data User3-----------------------------------------
recdata6=(data_equilized6'*code6')';
recdata6h5=real(recdata6)>0;
errors_user6h5(i) = size(find([data_user6- recdata6h5]),2); %Errors for User1
SBer6h5= errors_user6h5/N;   
end

%% ------------------Hop 6------------------------------
users=7; 

%user 1 estimated data %modulation 
estdata_user1bpskhop6 = 2*recdata1h5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user111=estdata_user1bpskhop6';
spdata1_user111=estdata_user111*code1;    % Spreading 
spdata121=(spdata1_user111)';
ifftdata_user111=ifft(spdata121);      % Taking the IFFT
ifftdata12=ifftdata_user111';
y11=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata11=y11';
esttx_user1hop6=transdata11;

%user 2 estimated data %modulation 
estdata_user2bpskhop6 = 2*recdata2h5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user222=estdata_user2bpskhop6';
spdata1_user222=estdata_user222*code3;    % Spreading 
spdata222=(spdata1_user222)';
ifftdata_user22=ifft(spdata222);      % Taking the IFFT
ifftdata122=ifftdata_user22';
y222=[ifftdata122(:,[(n-2):n]) ifftdata122];
transdata222=y222';
esttx_user2hop6=transdata222;   

%user 3 estimated data %modulation 
estdata_user3bpskhop6 = 2*recdata3h5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user3=estdata_user3bpskhop6';
spdata1_user3=estdata_user3*code2;    % Spreading 
spdata3=(spdata1_user3)';
ifftdata_user33=ifft(spdata3);      % Taking the IFFT
ifftdata3=ifftdata_user33';
y3h5=[ifftdata3(:,[(n-2):n]) ifftdata3];
transdata333=y3h5';
esttx_user3hop6=transdata333;          

%user 4 estimated data %modulation 
estdata_user4bpskhop6 = 2*recdata4h5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user4=estdata_user4bpskhop6';
spdata1_user4=estdata_user4*code4;    % Spreading 
spdata4=(spdata1_user4)';
ifftdata_user44=ifft(spdata4);      % Taking the IFFT
ifftdata4=ifftdata_user44';
y4h5=[ifftdata4(:,[(n-2):n]) ifftdata4];
transdata444=y4h5';
esttx_user4hop6=transdata444;      

%user 5 estimated data %modulation 
estdata_user5bpskhop5 = 2*recdata5h5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user5=estdata_user5bpskhop5';
spdata1_user5=estdata_user5*code5;    % Spreading 
spdata5=(spdata1_user5)';
ifftdata_user55=ifft(spdata5);      % Taking the IFFT
ifftdata5=ifftdata_user55';
y5h5=[ifftdata5(:,[(n-2):n]) ifftdata5];
transdata555=y5h5';
esttx_user5hop6=transdata555;  

%user 6 estimated data %modulation 
estdata_user6bpskhop6 = 2*recdata6h5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user6=estdata_user6bpskhop6';
spdata1_user6=estdata_user6*code6;    % Spreading 
spdata6=(spdata1_user6)';
ifftdata_user66=ifft(spdata6);      % Taking the IFFT
ifftdata6=ifftdata_user66';
y5h6=[ifftdata6(:,[(n-2):n]) ifftdata6];
transdata666=y5h6';
esttx_user6hop6=transdata666;  

%user 7 original data %modulation 
data_user7= rand(1,N)>0.5;          % Generation of data for user2
data_user7bpsk = 2*data_user7-1;    % BPSK modulation 0 -> -1; 1 -> 0 
data_user77=data_user7bpsk';
spdata7=data_user77*code7;          % Spreading 
spdata77=(spdata7)';
ifftdata_user7=ifft(spdata77);      % Taking the IFFT
ifftdata77=ifftdata_user7';
y7=[ifftdata77(:,[(n-2):n]) ifftdata77];
transdata7=y7';
tx_user7=transdata7;                

% Transmitting data of all users
x7 = esttx_user1hop6+esttx_user2hop6+esttx_user3hop6+esttx_user4hop6+esttx_user5hop6+esttx_user6hop6+tx_user7;

%rician fading channel 
%----------------------Creating a new Rayleigh Channel---------------------------
gain17=sqrt(d6^-eta)*sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
gain27=sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain37=sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3
gain47=sqrt(p4/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap4
x41=x7(:);
x42=reshape(x41,1,length(x41));
i=1:length(x42);        
delay1=1; 
for i=delay1+1:length(x42) % Producing one sample delay in Tap2 w.r.t. Tap1
   x43(i)=x7(i-delay1);
end
delay2=2;
for i=delay2+1:length(x42) % Producing two sample delay in Tap2 w.r.t. Tap1
   x44(i)=x7(i-delay2);
end
delay3=3;
for i=delay3+1:length(x42) % Producing three sample delay in Tap2 w.r.t. Tap1
   x45(i)=x7(i-delay3);
end
x1111=reshape(x43,(n+3),length(x43)/(n+3));
x2222=reshape(x44,(n+3),length(x44)/(n+3));
x3333=reshape(x45,(n+3),length(x45)/(n+3));
ch41=repmat(gain17,(n+3),1);     
ch42=repmat(gain27,(n+3),1);
ch43=repmat(gain37,(n+3),1);
ch44=repmat(gain47,(n+3),1);
data_channel7=x7.*ch41+x1111.*ch42+x2222.*ch43+x3333.*ch44;% Passing data through channel 

%addition of noise 
data_noise7=data_channel7(:);
data_noise77=reshape(data_noise7,1,length(data_noise7));
noise77 = 1/sqrt(2)*[randn(1,length(data_noise77)) + j*randn(1,length(data_noise77))];

%ifft
%equalization 
%decoding
%demodulation
%ber rate calculation
for i = 1:length(snr)
y7 = data_noise77 + (sqrt(1)*10^(-snr(i)/20))*noise77; %Addition of Noise
  
%--------------------------Receiver ---------------------------------------
data_received7 = y7;           %fadded data received with awgn noise
%---------------------Removing Cyclic Prefix-------------------------------
rx1=reshape(data_received7,(n+3),length(data_received7)/(n+3));
rx12=rx1';
rx13 = rx12(:,[(4:(n+3))]); 
rx14=rx13';

%-----------------Taking FFT ----------------------------------------------
fft_data_received7 =fft(rx14);

%----------------equilization of the channel-------------------------------
channel_response7=fft([gain17;gain27;gain37;gain47],n);
data_equilized7=fft_data_received7.*conj(channel_response7);

%----------------BER of Data User1-----------------------------------------
recdata110=(data_equilized7'*code1')';
recdata1h6=real(recdata110)>0;
errors_user1h6(i) = size(find([data_user1- recdata1h6]),2); %Errors for User1
SBer1h6 = errors_user1h6/N;                              % simulated ber user1

%----------------BER of Data User2-----------------------------------------
recdata2221=(data_equilized7'*code3')';
recdata2h6=real(recdata2221)>0;
errors_user2h6(i) = size(find([data_user2- recdata2h6]),2); %Errors for User1
SBer2h6 = errors_user2h6/N;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata3=(data_equilized7'*code2')';
recdata3h6=real(recdata3)>0;
errors_user3h6(i) = size(find([data_user3- recdata3h6]),2); %Errors for User1
SBer3h6 = errors_user3h6/N;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata4=(data_equilized7'*code4')';
recdata4h6=real(recdata4)>0;
errors_user4h6(i) = size(find([data_user4- recdata4h6]),2); %Errors for User1
SBer4h6 = errors_user4h6/N;      

%----------------BER of Data User3-----------------------------------------
recdata5=(data_equilized7'*code5')';
recdata5h6=real(recdata5)>0;
errors_user5h6(i) = size(find([data_user5- recdata5h6]),2); %Errors for User1
SBer5h6= errors_user5h6/N; 

%----------------BER of Data User3-----------------------------------------
recdata6=(data_equilized7'*code6')';
recdata6h6=real(recdata6)>0;
errors_user6h6(i) = size(find([data_user6- recdata6h6]),2); %Errors for User1
SBer6h6= errors_user6h6/N;   

%----------------BER of Data User3-----------------------------------------
recdata7=(data_equilized7'*code7')';
recdata7h6=real(recdata7)>0;
errors_user7h6(i) = size(find([data_user6- recdata7h6]),2); %Errors for User1
SBer7h6= errors_user7h6/N;   
end

%% ------------------Hop 7------------------------------
users=8; 

%user 1 estimated data %modulation 
estdata_user1bpskhop7 = 2*recdata1h6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user111=estdata_user1bpskhop7';
spdata1_user111=estdata_user111*code1;    % Spreading 
spdata121=(spdata1_user111)';
ifftdata_user111=ifft(spdata121);      % Taking the IFFT
ifftdata12=ifftdata_user111';
y11=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata11=y11';
esttx_user1hop7=transdata11;

%user 2 estimated data %modulation 
estdata_user2bpskhop7 = 2*recdata2h6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user222=estdata_user2bpskhop7';
spdata1_user222=estdata_user222*code3;    % Spreading 
spdata222=(spdata1_user222)';
ifftdata_user22=ifft(spdata222);      % Taking the IFFT
ifftdata122=ifftdata_user22';
y222=[ifftdata122(:,[(n-2):n]) ifftdata122];
transdata222=y222';
esttx_user2hop7=transdata222;   

%user 3 estimated data %modulation 
estdata_user3bpskhop7 = 2*recdata3h6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user3=estdata_user3bpskhop7';
spdata1_user3=estdata_user3*code2;    % Spreading 
spdata3=(spdata1_user3)';
ifftdata_user33=ifft(spdata3);      % Taking the IFFT
ifftdata3=ifftdata_user33';
y3h7=[ifftdata3(:,[(n-2):n]) ifftdata3];
transdata333=y3h7';
esttx_user3hop7=transdata333;   

%user 4 estimated data %modulation 
estdata_user4bpskhop7 = 2*recdata4h6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user4=estdata_user4bpskhop7';
spdata1_user4=estdata_user4*code4;    % Spreading 
spdata4=(spdata1_user4)';
ifftdata_user44=ifft(spdata4);      % Taking the IFFT
ifftdata4=ifftdata_user44';
y4h7=[ifftdata4(:,[(n-2):n]) ifftdata4];
transdata444=y4h7';
esttx_user4hop7=transdata444;  

%user 5 estimated data %modulation 
estdata_user5bpskhop7 = 2*recdata5h6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user5=estdata_user5bpskhop7';
spdata1_user5=estdata_user5*code5;    % Spreading 
spdata5=(spdata1_user5)';
ifftdata_user55=ifft(spdata5);      % Taking the IFFT
ifftdata5=ifftdata_user55';
y5h7=[ifftdata5(:,[(n-2):n]) ifftdata5];
transdata555=y5h7';
esttx_user5hop7=transdata555;  

%user 6 estimated data %modulation 
estdata_user6bpskhop7 = 2*recdata6h6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user6=estdata_user6bpskhop7';
spdata1_user6=estdata_user6*code6;    % Spreading 
spdata6=(spdata1_user6)';
ifftdata_user66=ifft(spdata6);      % Taking the IFFT
ifftdata6=ifftdata_user66';
y6h7=[ifftdata6(:,[(n-2):n]) ifftdata6];
transdata666=y6h7';
esttx_user6hop7=transdata666; 

%user 7 estimated data %modulation 
estdata_user7bpskhop7 = 2*recdata7h6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user7=estdata_user7bpskhop7';
spdata1_user7=estdata_user7*code7;    % Spreading 
spdata7=(spdata1_user7)';
ifftdata_user77=ifft(spdata7);      % Taking the IFFT
ifftdata7=ifftdata_user77';
y7h7=[ifftdata7(:,[(n-2):n]) ifftdata7];
transdata777=y7h7';
esttx_user7hop7=transdata777; 

%user 8 original data %modulation 
%user 7 original data %modulation 
data_user8= rand(1,N)>0.5;          % Generation of data for user2
data_user8bpsk = 2*data_user8-1;    % BPSK modulation 0 -> -1; 1 -> 0 
data_user88=data_user8bpsk';
spdata8=data_user88*code8;          % Spreading 
spdata88=(spdata8)';
ifftdata_user8=ifft(spdata88);      % Taking the IFFT
ifftdata88=ifftdata_user8';
y8=[ifftdata88(:,[(n-2):n]) ifftdata88];
transdata8=y8';
tx_user8=transdata8;                

%rician fading channel 
% Transmitting data of all users
x8 = esttx_user1hop7+esttx_user2hop7+esttx_user3hop7+esttx_user4hop7+esttx_user5hop7+esttx_user6hop7+esttx_user7hop7+tx_user8;

%rician fading channel 
%----------------------Creating a new Rayleigh Channel---------------------------
gain18=sqrt(d7^-eta)*sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
gain28=sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain38=sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3
gain48=sqrt(p4/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap4
x41=x8(:);
x42=reshape(x41,1,length(x41));
i=1:length(x42);        
delay1=1; 
for i=delay1+1:length(x42) % Producing one sample delay in Tap2 w.r.t. Tap1
   x43(i)=x8(i-delay1);
end
delay2=2;
for i=delay2+1:length(x42) % Producing two sample delay in Tap2 w.r.t. Tap1
   x44(i)=x8(i-delay2);
end
delay3=3;
for i=delay3+1:length(x42) % Producing three sample delay in Tap2 w.r.t. Tap1
   x45(i)=x8(i-delay3);
end
x1111=reshape(x43,(n+3),length(x43)/(n+3));
x2222=reshape(x44,(n+3),length(x44)/(n+3));
x3333=reshape(x45,(n+3),length(x45)/(n+3));
ch41=repmat(gain18,(n+3),1);     
ch42=repmat(gain28,(n+3),1);
ch43=repmat(gain38,(n+3),1);
ch44=repmat(gain48,(n+3),1);
data_channel8=x8.*ch41+x1111.*ch42+x2222.*ch43+x3333.*ch44;% Passing data through channel 

%addition of noise 
%addition of noise 
data_noise8=data_channel8(:);
data_noise88=reshape(data_noise8,1,length(data_noise8));
noise88 = 1/sqrt(2)*[randn(1,length(data_noise88)) + j*randn(1,length(data_noise88))];

%ifft
%equalization 
%decoding
%demodulation
%ber rate calculation
for i = 1:length(snr)
y8 = data_noise88 + (sqrt(1)*10^(-snr(i)/20))*noise88; %Addition of Noise
  
%--------------------------Receiver ---------------------------------------
data_received8 = y8;           %fadded data received with awgn noise
%---------------------Removing Cyclic Prefix-------------------------------
rx1=reshape(data_received8,(n+3),length(data_received8)/(n+3));
rx12=rx1';
rx13 = rx12(:,[(4:(n+3))]); 
rx14=rx13';

%-----------------Taking FFT ----------------------------------------------
fft_data_received8 =fft(rx14);

%----------------equilization of the channel-------------------------------
channel_response8=fft([gain18;gain28;gain38;gain48],n);
data_equilized8=fft_data_received8.*conj(channel_response8);

%----------------BER of Data User1-----------------------------------------
recdata110=(data_equilized8'*code1')';
recdata1h7=real(recdata110)>0;
errors_user1h7(i) = size(find([data_user1- recdata1h7]),2); %Errors for User1
SBer1h7 = errors_user1h7/N;                              % simulated ber user1

%----------------BER of Data User2-----------------------------------------
recdata2221=(data_equilized8'*code3')';
recdata2h7=real(recdata2221)>0;
errors_user2h7(i) = size(find([data_user2- recdata2h7]),2); %Errors for User1
SBer2h7 = errors_user2h7/N;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata3=(data_equilized8'*code2')';
recdata3h7=real(recdata3)>0;
errors_user3h7(i) = size(find([data_user3- recdata3h7]),2); %Errors for User1
SBer3h7 = errors_user3h7/N;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata4=(data_equilized8'*code4')';
recdata4h7=real(recdata4)>0;
errors_user4h7(i) = size(find([data_user4- recdata4h7]),2); %Errors for User1
SBer4h7 = errors_user4h7/N;      

%----------------BER of Data User3-----------------------------------------
recdata5=(data_equilized8'*code5')';
recdata5h7=real(recdata5)>0;
errors_user5h7(i) = size(find([data_user5- recdata5h7]),2); %Errors for User1
SBer5h7 = errors_user5h7/N; 

%----------------BER of Data User3-----------------------------------------
recdata6=(data_equilized8'*code6')';
recdata6h7=real(recdata6)>0;
errors_user6h7(i) = size(find([data_user6- recdata6h7]),2); %Errors for User1
SBer6h7 = errors_user6h7/N;   

%----------------BER of Data User3-----------------------------------------
recdata7=(data_equilized8'*code7')';
recdata7h7=real(recdata7)>0;
errors_user7h7(i) = size(find([data_user7- recdata7h7]),2); %Errors for User1
SBer7h7 = errors_user7h7/N;   

%----------------BER of Data User3-----------------------------------------
recdata8=(data_equilized8'*code8')';
recdata8h7=real(recdata8)>0;
errors_user8h7(i) = size(find([data_user8- recdata8h7]),2); %Errors for User1
SBer8h7= errors_user8h7/N;   
end


%% ------------------Hop 8------------------------------
users=9;  
%user 1 estimated data %modulation 
estdata_user1bpskhop7 = 2*recdata1h7-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user111=estdata_user1bpskhop7';
spdata1_user111=estdata_user111*code1;    % Spreading 
spdata121=(spdata1_user111)';
ifftdata_user111=ifft(spdata121);      % Taking the IFFT
ifftdata12=ifftdata_user111';
y1h8=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata11=y1h8';
esttx_user1hop8=transdata11;

%user 2 estimated data %modulation 
estdata_user2bpskhop7 = 2*recdata2h7-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user222=estdata_user2bpskhop7';
spdata1_user222=estdata_user222*code3;    % Spreading 
spdata222=(spdata1_user222)';
ifftdata_user22=ifft(spdata222);      % Taking the IFFT
ifftdata122=ifftdata_user22';
y2h8=[ifftdata122(:,[(n-2):n]) ifftdata122];
transdata222=y2h8';
esttx_user2hop8=transdata222;   

%user 3 estimated data %modulation 
estdata_user3bpskhop7 = 2*recdata3h7-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user3=estdata_user3bpskhop7';
spdata1_user3=estdata_user3*code2;    % Spreading 
spdata3=(spdata1_user3)';
ifftdata_user33=ifft(spdata3);      % Taking the IFFT
ifftdata3=ifftdata_user33';
y3h8=[ifftdata3(:,[(n-2):n]) ifftdata3];
transdata333=y3h8';
esttx_user3hop8=transdata333;   

%user 4 estimated data %modulation 
estdata_user4bpskhop7 = 2*recdata4h7-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user4=estdata_user4bpskhop7';
spdata1_user4=estdata_user4*code4;    % Spreading 
spdata4=(spdata1_user4)';
ifftdata_user44=ifft(spdata4);      % Taking the IFFT
ifftdata4=ifftdata_user44';
y4h8=[ifftdata4(:,[(n-2):n]) ifftdata4];
transdata444=y4h8';
esttx_user4hop8=transdata444;  

%user 5 estimated data %modulation 
estdata_user5bpskhop7 = 2*recdata5h7-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user5=estdata_user5bpskhop7';
spdata1_user5=estdata_user5*code5;    % Spreading 
spdata5=(spdata1_user5)';
ifftdata_user55=ifft(spdata5);      % Taking the IFFT
ifftdata5=ifftdata_user55';
y5h8=[ifftdata5(:,[(n-2):n]) ifftdata5];
transdata555=y5h8';
esttx_user5hop8=transdata555;  

%user 6 estimated data %modulation 
estdata_user6bpskhop7 = 2*recdata6h7-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user6=estdata_user6bpskhop7';
spdata1_user6=estdata_user6*code6;    % Spreading 
spdata6=(spdata1_user6)';
ifftdata_user66=ifft(spdata6);      % Taking the IFFT
ifftdata6=ifftdata_user66';
y6h8=[ifftdata6(:,[(n-2):n]) ifftdata6];
transdata666=y6h8';
esttx_user6hop8=transdata666; 

%user 7 estimated data %modulation 
estdata_user7bpskhop7 = 2*recdata7h7-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user7=estdata_user7bpskhop7';
spdata1_user7=estdata_user7*code7;    % Spreading 
spdata7=(spdata1_user7)';
ifftdata_user77=ifft(spdata7);      % Taking the IFFT
ifftdata7=ifftdata_user77';
y7h8=[ifftdata7(:,[(n-2):n]) ifftdata7];
transdata777=y7h8';
esttx_user7hop8=transdata777; 

%user 8 estimated data %modulation 
estdata_user8bpskhop7 = 2*recdata8h7-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user8=estdata_user8bpskhop7';
spdata1_user8=estdata_user8*code8;    % Spreading 
spdata8=(spdata1_user8)';
ifftdata_user88=ifft(spdata8);      % Taking the IFFT
ifftdata8=ifftdata_user88';
y7h8=[ifftdata8(:,[(n-2):n]) ifftdata8];
transdata888=y7h8';
esttx_user8hop8=transdata888; 

%user 9 original data %modulation 
data_user9= rand(1,N)>0.5;          % Generation of data for user2
data_user9bpsk = 2*data_user9-1;    % BPSK modulation 0 -> -1; 1 -> 0 
data_user99=data_user9bpsk';
spdata9=data_user99*code9;          % Spreading 
spdata99=(spdata9)';
ifftdata_user9=ifft(spdata99);      % Taking the IFFT
ifftdata99=ifftdata_user9';
y9=[ifftdata99(:,[(n-2):n]) ifftdata99];
transdata9=y9';
tx_user9=transdata9;                

%transmitting signal of all users
x9 = esttx_user1hop8+esttx_user2hop8+esttx_user3hop8+esttx_user4hop8+esttx_user5hop8+esttx_user6hop8+esttx_user7hop8+esttx_user8hop8+tx_user9;

%rician fading channel 
%----------------------Creating a new Rayleigh Channel---------------------------
gain19=sqrt(d9^-eta)*sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
gain29=sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain39=sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3
gain49=sqrt(p4/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap4
x41=x9(:);
x42=reshape(x41,1,length(x41));
i=1:length(x42);        
delay1=1; 
for i=delay1+1:length(x42) % Producing one sample delay in Tap2 w.r.t. Tap1
   x43(i)=x9(i-delay1);
end
delay2=2;
for i=delay2+1:length(x42) % Producing two sample delay in Tap2 w.r.t. Tap1
   x44(i)=x9(i-delay2);
end
delay3=3;
for i=delay3+1:length(x42) % Producing three sample delay in Tap2 w.r.t. Tap1
   x45(i)=x9(i-delay3);
end
x1111=reshape(x43,(n+3),length(x43)/(n+3));
x2222=reshape(x44,(n+3),length(x44)/(n+3));
x3333=reshape(x45,(n+3),length(x45)/(n+3));
ch41=repmat(gain19,(n+3),1);     
ch42=repmat(gain29,(n+3),1);
ch43=repmat(gain39,(n+3),1);
ch44=repmat(gain49,(n+3),1);
data_channel9=x9.*ch41+x1111.*ch42+x2222.*ch43+x3333.*ch44;% Passing data through channel 

%addition of noise 
%addition of noise 
data_noise9=data_channel9(:);
data_noise99=reshape(data_noise9,1,length(data_noise9));
noise99 = 1/sqrt(2)*[randn(1,length(data_noise99)) + j*randn(1,length(data_noise99))];

%ifft
%equalization 
%decoding
%demodulation
%ber rate calculation
for i = 1:length(snr)
y9 = data_noise99 + (sqrt(1)*10^(-snr(i)/20))*noise99; %Addition of Noise
  
%--------------------------Receiver ---------------------------------------
data_received9 = y9;           %fadded data received with awgn noise
%---------------------Removing Cyclic Prefix-------------------------------
rx1=reshape(data_received9,(n+3),length(data_received9)/(n+3));
rx12=rx1';
rx13 = rx12(:,[(4:(n+3))]); 
rx14=rx13';

%-----------------Taking FFT ----------------------------------------------
fft_data_received9 =fft(rx14);

%----------------equilization of the channel-------------------------------
channel_response9=fft([gain19;gain29;gain39;gain49],n);
data_equilized9=fft_data_received9.*conj(channel_response9);

%----------------BER of Data User1-----------------------------------------
recdata110=(data_equilized9'*code1')';
recdata1h8=real(recdata110)>0;
errors_user1h8(i) = size(find([data_user1- recdata1h8]),2); %Errors for User1
SBer1h8 = errors_user1h8/N;                              % simulated ber user1

%----------------BER of Data User2-----------------------------------------
recdata2221=(data_equilized9'*code3')';
recdata2h8=real(recdata2221)>0;
errors_user2h8(i) = size(find([data_user2- recdata2h8]),2); %Errors for User1
SBer2h8 = errors_user2h8/N;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata3=(data_equilized9'*code2')';
recdata3h8=real(recdata3)>0;
errors_user3h8(i) = size(find([data_user3- recdata3h8]),2); %Errors for User1
SBer3h8 = errors_user3h8/N;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata4=(data_equilized9'*code4')';
recdata4h8=real(recdata4)>0;
errors_user4h8(i) = size(find([data_user4- recdata4h7]),2); %Errors for User1
SBer4h8 = errors_user4h8/N;      

%----------------BER of Data User3-----------------------------------------
recdata5=(data_equilized9'*code5')';
recdata5h8=real(recdata5)>0;
errors_user5h8(i) = size(find([data_user5- recdata5h8]),2); %Errors for User1
SBer5h8 = errors_user5h8/N; 

%----------------BER of Data User3-----------------------------------------
recdata6=(data_equilized9'*code6')';
recdata6h7=real(recdata6)>0;
errors_user6h8(i) = size(find([data_user6- recdata6h7]),2); %Errors for User1
SBer6h8 = errors_user6h8/N;   

%----------------BER of Data User3-----------------------------------------
recdata7=(data_equilized9'*code7')';
recdata7h7=real(recdata7)>0;
errors_user7h8(i) = size(find([data_user7- recdata7h7]),2); %Errors for User1
SBer7h8 = errors_user7h8/N;   

%----------------BER of Data User3-----------------------------------------
recdata8=(data_equilized9'*code8')';
recdata8h7=real(recdata8)>0;
errors_user8h8(i) = size(find([data_user8- recdata8h7]),2); %Errors for User1
SBer8h8= errors_user8h8/N;   

%----------------BER of Data User3-----------------------------------------
recdata9=(data_equilized9'*code9')';
recdata9h8=real(recdata9)>0;
errors_user9h8(i) = size(find([data_user9- recdata9h8]),2); %Errors for User1
SBer9h8= errors_user9h8/N; 
end

%% ------------------------Theoretical Result-------------------------------
snrlnr=10.^(snr/10);
TBerf = .5*erfc(sqrt(k1*snrlnr./(snrlnr+k1)));% theoretical BER fro Flat fadding
%% ------------------------Figures ------------------------
%SBer1

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
