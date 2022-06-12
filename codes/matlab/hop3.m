function [data_user4,data_userB,SBer1,SBer2,SBer3,SBerA,SBer4,SBerB,recdata1200,recdata222,estrecdata3,est4A,estrecdata4,estrecdataB] = ...
    hop3(est1,est2,est3,estA,code1,code2,code3,codeA,code4,n,N,d4,...
    eta,p1,sigma,mean,dA,snr,data_user1,data_user2,data_user3,data_userA)
%fprintf('okay\n')

%------------------est for User1-------------------------------
estdata_user1bpskhop3 = 2*est1-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user111=estdata_user1bpskhop3';
spdata1_user111=estdata_user111*code1;    % Spreading 
spdata121=(spdata1_user111)';
ifftdata_user111=ifft(spdata121);      % Taking the IFFT
ifftdata12=ifftdata_user111';
y11=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata11=y11';
esttx_user11=transdata11;                % Transmitting data for user1

%------------------est for User2-------------------------------
estdata_user2bpskhop3 = 2*est2-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user222=estdata_user2bpskhop3';
spdata1_user222=estdata_user222*code3;    % Spreading 
spdata222=(spdata1_user222)';
ifftdata_user22=ifft(spdata222);      % Taking the IFFT
ifftdata122=ifftdata_user22';
y222=[ifftdata122(:,[(n-2):n]) ifftdata122];
transdata222=y222';
esttx_user22=transdata222;                % Transmitting data for user1

%------------------est for User3-------------------------------
estdata_user3bpskhop3 = 2*est3-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user333=estdata_user3bpskhop3';
spdata1_user333=estdata_user333*code2;    % Spreading 
spdata333=(spdata1_user333)';
ifftdata_user33=ifft(spdata333);      % Taking the IFFT
ifftdata33=ifftdata_user33';
y333=[ifftdata33(:,[(n-2):n]) ifftdata33];
transdata333=y333';
esttx_user33=transdata333;                % Transmitting data for user1

%------------------est for UserA-------------------------------
estdata_userAbpskhop3 = 2*estA-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_userA=estdata_userAbpskhop3';
spdata1_userA=estdata_userA*codeA;    % Spreading 
spdataA=(spdata1_userA)';
ifftdata_userA=ifft(spdataA);      % Taking the IFFT
ifftdataA=ifftdata_userA';
yA=[ifftdataA(:,[(n-2):n]) ifftdataA];
transdataA=yA';
esttx_userA=transdataA;                % Transmitting data for user1

%------------------Generating original data for User4-------------------------------
data_user4= rand(1,N)>0.5;          % Generation of data for user2
data_user4bpsk = 2*data_user4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
data_user44=data_user4bpsk';
spdata4=data_user44*code4;          % Spreading 
spdata44=(spdata4)';
ifftdata_user4=ifft(spdata44);      % Taking the IFFT
ifftdata44=ifftdata_user4';
y4=[ifftdata44(:,[(n-2):n]) ifftdata44];
transdata4=y4';
tx_user4=transdata4;                % Transmitting data for user2

%----------------------Adding data for Transmission of All User------------
x4=esttx_user11+esttx_user22+esttx_user33+esttx_userA+tx_user4;
%----------------------Creating a new Rayleigh Channel---------------------------
p2=p1;
p3=p1;
p4=p1;
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

%% ------user b using same code (code 4) but different channel 
%original data for B                            % Number of Bits for  data_user2
data_userB= rand(1,N)>0.5;          % Generation of data for user2
data_userBbpsk = 2*data_userB-1;    % BPSK modulation 0 -> -1; 1 -> 0 
data_userBh2=data_userBbpsk';
spdataBh2=data_userBh2*code4;          % Spreading 
spdataBh2t=(spdataBh2)';
ifftdata_userBh2=ifft(spdataBh2t);      % Taking the IFFT
ifftdataBh2=ifftdata_userBh2';
yB=[ifftdataBh2(:,[(n-2):n]) ifftdataBh2];
transdataBh2=yB';
tx_userB=transdataBh2; 

%----------------------Adding data for Transmission of All User------------
xB=tx_userB;
%----------------------Creating a new Rayleigh Channel---------------------------
data_channelB=sqrt(dA^-eta)*sqrt(p1/2)*xB;  % Passing data through channel 

data_channelhop4 = data_channel4+data_channelB;

%------------------------Addition of AWGN noise ---------------------------
data_noise4=data_channelhop4(:);
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

%----------------equilization of the channel4-------------------------------
channel_response4=fft([gain1sss;gain2sss;gain3sss;gain4sss],n);
data_equilized4=fft_data_received4.*conj(channel_response4);

%----------------equilization of the channelB--------------------------
data_equilizedB=fft_data_received4/(sqrt(dA^-eta)*sqrt(p1/2));
%% 
%----------------BER of Data User1-----------------------------------------
recdata110=(data_equilized4'*code1')';
recdata1200=real(recdata110)>0;
errors_user111(i) = size(find([data_user1- recdata1200]),2); %Errors for User1
SBer1 = errors_user111/N;                              % simulated ber user1

%----------------BER of Data User2-----------------------------------------
recdata2221=(data_equilized4'*code3')';
recdata222=real(recdata2221)>0;
errors_user222(i) = size(find([data_user2- recdata222]),2); %Errors for User1
SBer2 = errors_user222/N;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata3=(data_equilized4'*code2')';
estrecdata3=real(recdata3)>0;
errors_user3(i) = size(find([data_user3- estrecdata3]),2); %Errors for User1
SBer3 = errors_user3/N;      
%----------------BER of Data UserA-----------------------------------------
recdata4A=(data_equilized4'*codeA')';
est4A=real(recdata4A)>0;
errors_user4A(i) = size(find([data_userA- est4A]),2); %Errors for User1
SBerA = errors_user4A/N;                               % simulated ber user2

%----------------BER of Data User 4 AND B SIC DECODING-----------------------------------------
recdata4B=(data_equilizedB'*code4')';
xB_hat = ones(1,N);
xB_hat(recdata4B<0) = -1;
recdata4user=(data_equilized4'*code4'-xB_hat')';
estrecdata4=real(recdata4user)>0;
errors_user4(i) = size(find([data_user4- estrecdata4]),2); %Errors for User1
SBer4 = errors_user4/N;   

x4user_hat = ones(1,N);
x4user_hat(recdata4user<0) = -1;
recdataB=(data_equilizedB'*code4'-x4user_hat')';
estrecdataB=real(recdataB)>0;
errors_userB(i) = size(find([data_userB- estrecdataB]),2); %Errors for User1
SBerB = errors_userB/N;                              % simulated ber user1

end

end
