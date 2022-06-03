%ref:https://in.mathworks.com/matlabcentral/fileexchange/23498-mc-cdma
clc;close all;
clear all;
N=10^5; 
%% ------------------Hop 1------------------------------
users=1;            % Number of Users
%------------------Generation of Walsh code--------------------------------
n =4;                               %Number of  Data Sub-Carriers
walsh=hadamard(n);              
code1=walsh(2,:);                   %Taking 2nd row of walsh code for User1

%------------------Generating data for User1-------------------------------
                            % Number of Bits for  data_user1
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
p1=0.5/2.3;                                       % Power of Tap1
p2=0.9/2.3;                                     % Power of Tap2
p3=0.7/2.3;                                     % Power of Tap3
p4=0.2/2.3;
gain1=sqrt(p1/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap1
gain2=sqrt(p2/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap2
gain3=sqrt(p3/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap3
gain4=sqrt(p4/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap4
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
x1=reshape(x13,(n+3),length(x13)/(n+3));
x2=reshape(x14,(n+3),length(x14)/(n+3));
x3=reshape(x15,(n+3),length(x15)/(n+3));
ch1=repmat(gain1,(n+3),1);     
ch2=repmat(gain2,(n+3),1);
ch3=repmat(gain3,(n+3),1);
ch4=repmat(gain4,(n+3),1);
data_channel=x.*ch1+x1.*ch2+x2.*ch3+x3.*ch4;  % Passing data through channel 
%------------------------Addition of AWGN noise ---------------------------
data_noise1=data_channel(:);
data_noise2=reshape(data_noise1,1,length(data_noise1));
noise = 1/sqrt(2)*[randn(1,length(data_noise2)) + j*randn(1,length(data_noise2))]; 
snr = [0:20];                 % multiple Eb/N0 values

for i = 1:length(snr)
y = data_noise2 + (sqrt(1)*10^(-snr(i)/20))*noise; %Addition of Noise
  
%--------------------------Receiver ---------------------------------------
data_received =y;           %fadded data received with awgn noise
%---------------------Removing Cyclic Prefix-------------------------------
rx1=reshape(data_received,(n+3),length(data_received)/(n+3));
rx12=rx1';
rx13 = rx12(:,[(4:(n+3))]); 
rx14=rx13';
%-----------------Taking FFT ----------------------------------------------
fft_data_received =fft(rx14);
%----------------equilization of the channel-------------------------------
channel_response=fft([gain1;gain2;gain3;gain4],n);
data_equilized=fft_data_received.*conj(channel_response);
%----------------BER of Data User1-----------------------------------------
recdata11=(data_equilized'*code1')';
recdata12=real(recdata11)>0;
errors_user1(i) = size(find([data_user1- recdata12]),2); %Errors for User1
SBer1 = errors_user1/N;                              % simulated ber user1
                              % simulated ber user2
end

%% ------------------Hop 2------------------------------

users=2;            % Number of Users
%------------------Generation of Walsh code--------------------------------
n =4;                               %Number of  Data Sub-Carriers
walsh=hadamard(n);              
code1=walsh(2,:);                   %Taking 2nd row of walsh code for User1
code2=walsh(4,:);                   %Taking 3rd row of walsh code for User2
%------------------est for User1-------------------------------

%------------------Spreading & IFFT for User1------------------------------
estdata_user11=recdata12';
spdata1_user1=estdata_user11*code1;    % Spreading 
spdata12=(spdata1_user1)';
ifftdata_user1=ifft(spdata12);      % Taking the IFFT
ifftdata12=ifftdata_user1';
%------------------Append Cyclic Prefix1 for User1-------------------------
y1=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata1=y1';
esttx_user1=transdata1;                % Transmitting data for user1
%------------------Generating data for User2-------------------------------
M=N;                             % Number of Bits for  data_user2
data_user2= rand(1,M)>0.5;          % Generation of data for user2
data_user2bpsk = 2*data_user2-1;    % BPSK modulation 0 -> -1; 1 -> 0 
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
x2=esttx_user1+tx_user2;
%----------------------Creating a new Rayleigh Channel---------------------------
Taps=4;                                        % Number of Taps
p1=0.5/2.3;                                       % Power of Tap1
p2=0.9/2.3;                                     % Power of Tap2
p3=0.7/2.3;                                     % Power of Tap3
p4=0.2/2.3;
gain1=sqrt(p1/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap1
gain2=sqrt(p2/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap2
gain3=sqrt(p3/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap3
gain4=sqrt(p4/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap4
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
data_received2 =y2;           %fadded data received with awgn noise
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
recdata120=real(recdata110)>0;
errors_user11(i) = size(find([data_user1- recdata120]),2); %Errors for User1
SBer11 = errors_user11/N;                              % simulated ber user1
%----------------BER of Data User2-----------------------------------------
recdata21=(data_equilized2'*code2')';
recdata22=real(recdata21)>0;
errors_user22(i) = size(find([data_user2- recdata22]),2); %Errors for User1
SBer22 = errors_user22/M;                               % simulated ber user2
end

%% ------------------Hop 3------------------------------

users=3;            % Number of Users
%------------------Generation of Walsh code--------------------------------
n =4;                               %Number of  Data Sub-Carriers
walsh=hadamard(n);              
code1=walsh(2,:);                   %Taking 2nd row of walsh code for User1
code2=walsh(4,:);                   %Taking 3rd row of walsh code for User2
code3=walsh(3,:);                   %Taking 2nd row of walsh code for User1

%------------------est for User1-------------------------------

%------------------Spreading & IFFT for User1------------------------------
estdata_user111=recdata120';
spdata1_user1=estdata_user111*code1;    % Spreading 
spdata12=(spdata1_user1)';
ifftdata_user1=ifft(spdata12);      % Taking the IFFT
ifftdata12=ifftdata_user1';
%------------------Append Cyclic Prefix1 for User1-------------------------
y11=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata11=y11';
esttx_user11=transdata11;                % Transmitting data for user1

%------------------est for User2-------------------------------

%------------------Spreading & IFFT for User1------------------------------
estdata_user222=recdata22';
spdata1_user222=estdata_user222*code1;    % Spreading 
spdata222=(spdata1_user222)';
ifftdata_user22=ifft(spdata222);      % Taking the IFFT
ifftdata122=ifftdata_user22';
%------------------Append Cyclic Prefix1 for User1-------------------------
y222=[ifftdata122(:,[(n-2):n]) ifftdata122];
transdata222=y222';
esttx_user22=transdata222;                % Transmitting data for user1

%------------------Generating data for User3-------------------------------
M=N;                             % Number of Bits for  data_user2
data_user3= rand(1,M)>0.5;          % Generation of data for user2
data_user3bpsk = 2*data_user3-1;    % BPSK modulation 0 -> -1; 1 -> 0 
%-----------------Spreading & IFFT for User2-------------------------------
data_user31=data_user3bpsk';
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
Taps=4;                                        % Number of Taps
p1=0.5/2.3;                                       % Power of Tap1
p2=0.9/2.3;                                     % Power of Tap2
p3=0.7/2.3;                                     % Power of Tap3
p4=0.2/2.3;
gain1=sqrt(p1/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap1
gain2=sqrt(p2/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap2
gain3=sqrt(p3/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap3
gain4=sqrt(p4/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap4
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
data_received3 =y3;           %fadded data received with awgn noise
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
recdata1200=real(recdata110)>0;
errors_user111(i) = size(find([data_user1- recdata1200]),2); %Errors for User1
SBer111 = errors_user111/N;                              % simulated ber user1
%----------------BER of Data User2-----------------------------------------
recdata21=(data_equilized3'*code2')';
recdata222=real(recdata21)>0;
errors_user222(i) = size(find([data_user2- recdata222]),2); %Errors for User1
SBer222 = errors_user222/M;                               % simulated ber user2

recdata31=(data_equilized3'*code3')';
recdata333=real(recdata31)>0;
errors_user333(i) = size(find([data_user3- recdata333]),2); %Errors for User1
SBer333 = errors_user333/M;                               % simulated ber user2
end
% ------------------------Theoretical Result-------------------------------
snrlnr=10.^(snr/10);
TBer = 0.5*erfc(sqrt(snrlnr)); % Theoretical BER for AWGN
TBerf = 0.5.*(1-sqrt(snrlnr./(snrlnr+1)));% theoretical BER fro Flat fadding
%-------------------Displaying Result--------------------------------------       
figure
grid on;
semilogy(snr,TBer,'c*-','LineWidth',2);
hold on;
semilogy(snr,TBerf,'r-','LineWidth',2);
hold on;
semilogy(snr,SBer1,'kd-','LineWidth',2);
hold on;
semilogy(snr,SBer11,'kd--','LineWidth',2);
hold on;
semilogy(snr,SBer22,'go--','LineWidth',2);
hold on;
semilogy(snr,SBer111,'kd-.','LineWidth',2);
hold on;
semilogy(snr,SBer222,'go-.','LineWidth',2);
hold on;
semilogy(snr,SBer333,'mo-.','LineWidth',2);
axis([0 20 10^-5 1]);
grid on
legend('Theoratical BER for BPSK on AWGN ','Theoratical BER for BPSK on Rayleigh Channel ' ,'Simulated BER for User1 - hop 1','Simulated BER for User1 - hop 2','Simulated BER for User2 - hop 2','Simulated BER for User1 - hop 3','Simulated BER for User2 - hop 3','Simulated BER for User3 - hop 3','FontSize', 11);
xlabel('SNR, dB','FontSize', 18);
ylabel('Bit Error Rate','FontSize', 18);
title('BER Vs SNR on Rayleigh Channel','FontSize', 18)