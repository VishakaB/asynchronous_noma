%ref:https://in.mathworks.com/matlabcentral/fileexchange/23498-mc-cdma
clc;close all;
clear all;

%% ------------------Hop 1------------------------------

users=1;            % Number of Users of hop 1
  
%------------------Generation of Walsh code--------------------------------
n =4;                               %Number of  Data Sub-Carriers
walsh=hadamard(n);              
code1=walsh(1,:);                   %Taking 2nd row of walsh code for User1
                  %Taking 3rd row of walsh code for User2
N =10^4;  

%------------------ Function to generate data of users-------------------------------

[tx_user1,data_user1] = userdataGen(N,code1,n);

%----------------------Creating Rayleigh Channel for users---------------------------
p1=0.5/2.3;                                     % Power of Tap1
p2=0.9/2.3;                                     % Power of Tap2
p3=0.7/2.3;                                     % Power of Tap3
p4=0.2/2.3;
gain1=sqrt(p1/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap1
gain2=sqrt(p2/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap2
gain3=sqrt(p3/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap3
gain4=sqrt(p4/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap4
channel_response=fft([gain1;gain2;gain3;gain4],n);

[data_channel1] = rayleighFading(N,tx_user1,n);

%------------------------Addition of AWGN noise ---------------------------
[datawithnoise] = awgnNoiseAddition(data_channel1);

%--------------------------Receiver (3rd user) ---------------------------------------
snr = [0:40]; % multiple Eb/N0 values

for i = 1:length(snr)
    
    y = datawithnoise + (sqrt(1)*10^(-snr(i)/20))*datawithnoise; %Addition of Noise
    data_equilized1 = dataequilize(y,n,channel_response,gain1,gain2,gain3,gain4);

    [SBer1,est1] = berCompute(data_equilized1,data_user1,code1,N,i);

end

%% ------------------Hop 2------------------------------
users=1;            % Number of Users of hop 2
  
%------------------Generation of Walsh code--------------------------------
n =4;                               %Number of  Data Sub-Carriers
walsh=hadamard(n);              
code2=walsh(2,:);                   %Taking 2nd row of walsh code for User1

N =10^4;  

%------------------ Function to generate data of users-------------------------------
[tx_user2,data_user2] = userdataGen(N,code2,n);

[tx_user1_hop2] = estDataTransmit(N,code1,n,est1);
x2=tx_user1_hop2+tx_user2;

%----------------------Creating Rayleigh Channel for users---------------------------
p1=0.5/2.3;                                     % Power of Tap1
p2=0.9/2.3;                                     % Power of Tap2
p3=0.7/2.3;                                     % Power of Tap3
p4=0.2/2.3;
gain1=sqrt(p1/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap1
gain2=sqrt(p2/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap2
gain3=sqrt(p3/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap3
gain4=sqrt(p4/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap4
channel_response=fft([gain1;gain2;gain3;gain4],n);

[data_channel2] = rayleighFading(N,x2,n);

%------------------------Addition of AWGN noise ---------------------------
[datawithnoise2] = awgnNoiseAddition(data_channel2);

%--------------------------Receiver ---------------------------------------
snr = [0:40]; % multiple Eb/N0 values

for i = 1:length(snr)
    
 y = datawithnoise2 + (sqrt(1)*10^(-snr(i)/20))*datawithnoise2; %Addition of Noise
 data_equilized2 = dataequilize(y,n,channel_response,gain1,gain2,gain3,gain4);

 [SBer1,est1]  = berCompute(data_equilized2,data_user1,code1,N,i);
 [SBer2,est2]  = berCompute(data_equilized2,data_user2,code2,N,i);

end

%% ------------------Hop 3------------------------------

users=1;            % Number of Users of hop 1
  
%------------------Generation of Walsh code--------------------------------
n =4;                               %Number of  Data Sub-Carriers
walsh=hadamard(n);              
code3=walsh(3,:);                   %Taking 2nd row of walsh code for User1

N =10^4;  

%------------------ Function to generate data of users-------------------------------
[tx_user3,data_user3] = userdataGen(N,code3,n);

[tx_user1_hop3] = estDataTransmit(N,code1,n,est1);
[tx_user2_hop3] = estDataTransmit(N,code2,n,est2);

%----------------------Creating Rayleigh Channel for users---------------------------
p1=0.5/2.3;                                     % Power of Tap1
p2=0.9/2.3;                                     % Power of Tap2
p3=0.7/2.3;                                     % Power of Tap3
p4=0.2/2.3;
gain1=sqrt(p1/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap1
gain2=sqrt(p2/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap2
gain3=sqrt(p3/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap3
gain4=sqrt(p4/2)*[randn(1,N) + j*randn(1,N)];   % Gain for Tap4
channel_response=fft([gain1;gain2;gain3;gain4],n);

%----------------------superposition of the messages of different users---
x3=tx_user1_hop3+tx_user2_hop3+tx_user3;

[data_channel3] = rayleighFading(N,x3,n);


%transmit based on time slot
%time slot 1
%------------------------Addition of AWGN noise ---------------------------
[datawithnoise3] = awgnNoiseAddition(data_channel3);

%--------------------------Receiver ---------------------------------------
snr = [0:40]; % multiple Eb/N0 values

for i = 1:length(snr)
    
y =datawithnoise3 + (sqrt(1)*10^(-snr(i)/20))*datawithnoise3; %Addition of Noise
data_equilized3 = dataequilize(y,n,channel_response,gain1,gain2,gain3,gain4);

[SBer1,est1] = berCompute(data_equilized3,data_user1,code1,N,i);
[SBer2,est2] = berCompute(data_equilized3,data_user2,code2,N,i);
[SBer3,est3] = berCompute(data_equilized3,data_user3,code3,N,i);

end

%% ------------------Hop 4------------------------------




