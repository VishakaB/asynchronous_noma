function [data_user5,data_userC,SBer1h4,SBer2h4,SBer3h4,SBerAh4,SBer4h4,SBerBh4,...
    SBer5h4,SBerCh4,est1h4,est2h4,est3h4,estAh4,est4h4,estBh4,est5h4,estCh4]= ...
    hop4_nopower(est1,...
    est2,est3,estA,est4,estB,code1,code2,code3,codeA,...
    code4,codeB,code5,codeC,n,N,d5,eta,p1,sigma,mean,dA,snr,data_user1,data_user2,data_user3,data_userA,data_user4,data_userB)

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
esttx1=transdata11;                % Transmitting data for user1

%------------------est for User2-------------------------------
estdata_user2bpskhop3 = 2*est2-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user222=estdata_user2bpskhop3';
spdata1_user222=estdata_user222*code3;    % Spreading 
spdata222=(spdata1_user222)';
ifftdata_user22=ifft(spdata222);      % Taking the IFFT
ifftdata122=ifftdata_user22';
y222=[ifftdata122(:,[(n-2):n]) ifftdata122];
transdata222=y222';
esttx2=transdata222;                % Transmitting data for user1

%------------------est for User3-------------------------------
estdata_user3bpskhop3 = 2*est3-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user333=estdata_user3bpskhop3';
spdata1_user333=estdata_user333*code2;    % Spreading 
spdata333=(spdata1_user333)';
ifftdata_user33=ifft(spdata333);      % Taking the IFFT
ifftdata33=ifftdata_user33';
y333=[ifftdata33(:,[(n-2):n]) ifftdata33];
transdata333=y333';
esttx3=transdata333;                % Transmitting data for user1

%------------------est for UserA-------------------------------
estdata_userAbpskhop3 = 2*estA-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_userA=estdata_userAbpskhop3';
spdata1_userA=estdata_userA*codeA;    % Spreading 
spdataA=(spdata1_userA)';
ifftdata_userA=ifft(spdataA);      % Taking the IFFT
ifftdataA=ifftdata_userA';
yA=[ifftdataA(:,[(n-2):n]) ifftdataA];
transdataA=yA';
esttxA=transdataA;                % Transmitting data for user1

%user 4 estimated data %modulation 
estdata_user4bpskhop4 = 2*est4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user444=estdata_user4bpskhop4';
spdata1_user444=estdata_user444*code4;    % Spreading 
spdata444=(spdata1_user444)';
ifftdata_user44=ifft(spdata444);      % Taking the IFFT
ifftdata44=ifftdata_user44';
y444=[ifftdata44(:,[(n-2):n]) ifftdata44];
transdata444=y444';
esttx4=transdata444;                % Transmitting data for user1

%USER B ESTIMATED DATA
estdata_userBbpskhop4 = 2*estB-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user444=estdata_userBbpskhop4';
spdata1_user444=estdata_user444*codeB;    % Spreading 
spdata444=(spdata1_user444)';
ifftdata_user44=ifft(spdata444);      % Taking the IFFT
ifftdata44=ifftdata_user44';
y444=[ifftdata44(:,[(n-2):n]) ifftdata44];
transdata444=y444';
esttxB=transdata444;                % Transmitting data for user1

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
tx_user5=transdata5;   

%----------------------Adding data for Transmission of All User------------
x5=esttx1+esttx2+esttx3+esttxA+esttx4+esttxB+tx_user5;

%rician fading channel 
%----------------------Creating a new Rayleigh Channel---------------------------
p2=p1;p3=p1;p4=p1;
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

%% ------user C using same code (code 5) but different channel 
%original data for C                            % Number of Bits for  data_user2
data_userC = rand(1,N)>0.5;          % Generation of data for user2
data_userCbpsk = 2*data_userC-1;    % BPSK modulation 0 -> -1; 1 -> 0 
data_userCh5=data_userCbpsk';
spdataCh5=data_userCh5*code5;          % Spreading 
spdataCh5t=(spdataCh5)';
ifftdata_userCh5=ifft(spdataCh5t);      % Taking the IFFT
ifftdataCh5=ifftdata_userCh5';
yC=[ifftdataCh5(:,[(n-2):n]) ifftdataCh5];
transdataCh5=yC';
tx_userC=transdataCh5; 
xC=tx_userC;

%----------------------Creating a new Rayleigh Channel---------------------------
data_channelC=sqrt(dA^-eta)*sqrt(p1/2)*xC;  % Passing data through channel 

data_channelhop5 = data_channel5+data_channelC;

%addition of noise 
data_noise5=data_channelhop5(:);
data_noise55=reshape(data_noise5,1,length(data_noise5));
noise55 = 1/sqrt(2)*[randn(1,length(data_noise55)) + j*randn(1,length(data_noise55))];

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

%equalization over channel C
data_equilizedC=fft_data_received5/(sqrt(dA^-eta)*sqrt(p1/2));

%----------------BER of Data User1-----------------------------------------
recdata110=(data_equilized5'*code1')';
est1h4=real(recdata110)>0;
errors_user1h4(i) = size(find([data_user1- est1h4]),2); %Errors for User1
SBer1h4 = errors_user1h4/N;                              % simulated ber user1

%----------------BER of Data User2-----------------------------------------
recdata2221=(data_equilized5'*code3')';
est2h4=real(recdata2221)>0;
errors_user2h4(i) = size(find([data_user2- est2h4]),2); %Errors for User1
SBer2h4 = errors_user2h4/N;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata3=(data_equilized5'*code2')';
est3h4=real(recdata3)>0;
errors_user3h4(i) = size(find([data_user3- est3h4]),2); %Errors for User1
SBer3h4 = errors_user3h4/N;                               % simulated ber user2

%----------------BER of Data UserA-----------------------------------------
recdataAh4=(data_equilized5'*codeA')';
estAh4=real(recdataAh4)>0;
errors_userAh4(i) = size(find([data_userA- estAh4]),2); %Errors for User1
SBerAh4 = errors_userAh4/N;   

%----------------BER of Data User4-----------------------------------------
recdata4h4=(data_equilized5'*code4')';
est4h4=real(recdata4h4)>0;
errors_user4h4(i) = size(find([data_user4- est4h4]),2); %Errors for User1
SBer4h4 = errors_user4h4/N;      

%----------------BER of Data UserB-----------------------------------------
recdataBh4=(data_equilized5'*codeB')';
estBh4=real(recdataBh4)>0;
errors_userBh4(i) = size(find([data_userB- estBh4]),2); %Errors for User1
SBerBh4 = errors_userBh4/N;   

%----------------BER of Data User5 SIC DECODING 5 USER AND C -----------------------------------------
recdata5h4=(data_equilized5'*code5')';
est5h4=real(recdata5h4)>0;
errors_user5h4(i) = size(find([data_user5- est5h4]),2); %Errors for User1
SBer5h4 = errors_user5h4/N;  

recdataCh4=(data_equilizedC'*code5')';
estCh4=real(recdataCh4)>0;
errors_userCh4(i) = size(find([data_userC- estCh4]),2); %Errors for User1
SBerCh4 = errors_userCh4/N;  
end

end
