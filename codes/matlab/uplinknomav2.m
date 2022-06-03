%ref:https://in.mathworks.com/matlabcentral/fileexchange/23498-mc-cdma
clc;
close all;
clear all;
k1=5; %Rician factor %ref: https://www.researchgate.net/publication/263669548_Probability_Distribution_of_Rician_K-Factor_in_Urban_Suburban_and_Rural_Areas_Using_Real-World_Captured_Data
mean=sqrt(k1/(k1+1));% mean
sigma=sqrt(1/(2*(k1+1)));%variance
N=10^5;  % Number of Bits for data_user1
d1 = 100; d2 = 500;    %Distances of users from base station (BS)
d3 = 5;
eta = 4;            %Path loss exponent

%% ------------------Hop 1------------------------------
users=1;            % Number of Users
%------------------Generation of Walsh code--------------------------------
n =16;                               %Number of  Data Sub-Carriers
walsh=hadamard(n);              
code1=walsh(2,:);                   %Taking 2nd row of walsh code for User1

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
data_channel =data_channel1+ data_channel2;
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