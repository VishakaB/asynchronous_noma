function [data_user7,data_userE,SBer1h6,SBer2h6,SBer3h6,SBerAh6,SBer4h6,SBerBh6,...
    SBer5h6,SBerCh6,SBer6h6,SBerDh6,SBer7h6,SBerEh6,est1h6,est2h6,est3h6,estAh6,est4h6,estBh6,est5h6,estCh6,...
    est6h6,estDh6,est7h6,estEh6]= ...
    hop6(est1h5,est2h5,est3h5,estAh5,est4h5,estBh5,est5h5,estCh5,...
    est6h5,estDh5,code1,code2,code3,codeA,...
    code4,codeB,code5,codeC,code6,codeD,code7,codeE,n,N,d6,eta,p1,sigma,mean,dA,snr,data_user1,data_user2,...
    data_user3,data_userA,data_user4,data_userB,data_user5,data_userC,data_user6,data_userD)

%fprintf('okay\n')
p2=p1;
p3=p1;
p4=p1;
M=N;
%------------------est for User1-------------------------------
estdata_user1bpskhop3 = 2*est1h5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user111=estdata_user1bpskhop3';
spdata1_user111=estdata_user111*code1;    % Spreading 
spdata121=(spdata1_user111)';
ifftdata_user111=ifft(spdata121);      % Taking the IFFT
ifftdata12=ifftdata_user111';
y11=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata11=y11';
esttx1=transdata11;                % Transmitting data for user1

%------------------est for User2-------------------------------
estdata_user2bpskhop3 = 2*est2h5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user222=estdata_user2bpskhop3';
spdata1_user222=estdata_user222*code3;    % Spreading 
spdata222=(spdata1_user222)';
ifftdata_user22=ifft(spdata222);      % Taking the IFFT
ifftdata122=ifftdata_user22';
y222=[ifftdata122(:,[(n-2):n]) ifftdata122];
transdata222=y222';
esttx2=transdata222;                % Transmitting data for user1

%------------------est for User3-------------------------------
estdata_user3bpskhop3 = 2*est3h5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user333=estdata_user3bpskhop3';
spdata1_user333=estdata_user333*code2;    % Spreading 
spdata333=(spdata1_user333)';
ifftdata_user33=ifft(spdata333);      % Taking the IFFT
ifftdata33=ifftdata_user33';
y333=[ifftdata33(:,[(n-2):n]) ifftdata33];
transdata333=y333';
esttx3=transdata333;                % Transmitting data for user1

%------------------est for UserA-------------------------------
estdata_userAbpskhop3 = 2*estAh5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_userA=estdata_userAbpskhop3';
spdata1_userA=estdata_userA*codeA;    % Spreading 
spdataA=(spdata1_userA)';
ifftdata_userA=ifft(spdataA);      % Taking the IFFT
ifftdataA=ifftdata_userA';
yA=[ifftdataA(:,[(n-2):n]) ifftdataA];
transdataA=yA';
esttxA=transdataA;                % Transmitting data for user1

%------------------est for User4-------------------------------
estdata_user4bpskhop4 = 2*est4h5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user444=estdata_user4bpskhop4';
spdata1_user444=estdata_user444*code4;    % Spreading 
spdata444=(spdata1_user444)';
ifftdata_user44=ifft(spdata444);      % Taking the IFFT
ifftdata44=ifftdata_user44';
y444=[ifftdata44(:,[(n-2):n]) ifftdata44];
transdata444=y444';
esttx4=transdata444;                % Transmitting data for user1

%------------------est for UserB-------------------------------
estdata_userBbpskhop4 = 2*estBh5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user444=estdata_userBbpskhop4';
spdata1_user444=estdata_user444*codeB;    % Spreading 
spdata444=(spdata1_user444)';
ifftdata_user44=ifft(spdata444);      % Taking the IFFT
ifftdata44=ifftdata_user44';
y444=[ifftdata44(:,[(n-2):n]) ifftdata44];
transdata444=y444';
esttxB=transdata444;             

%------------------est for User5-------------------------------
estdata_user5bpskhop5 = 2*est5h5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user5=estdata_user5bpskhop5';
spdata1_user5=estdata_user5*code5;    % Spreading 
spdata5=(spdata1_user5)';
ifftdata_user55=ifft(spdata5);      % Taking the IFFT
ifftdata5=ifftdata_user55';
y5h5=[ifftdata5(:,[(n-2):n]) ifftdata5];
transdata555=y5h5';
esttx5=transdata555;         

%------------------est for UserC-------------------------------
estdata_userCbpskhop5 = 2*estCh5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user5=estdata_userCbpskhop5';
spdata1_user5=estdata_user5*codeC;    % Spreading 
spdata5=(spdata1_user5)';
ifftdata_user55=ifft(spdata5);      % Taking the IFFT
ifftdata5=ifftdata_user55';
y5h5=[ifftdata5(:,[(n-2):n]) ifftdata5];
transdata555=y5h5';
esttxC=transdata555;  

%------------------est for User6-------------------------------
estdata_user6bpskhop6 = 2*est6h5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user6=estdata_user6bpskhop6';
spdata1_user6=estdata_user6*code6;    % Spreading 
spdata6=(spdata1_user6)';
ifftdata_user66=ifft(spdata6);      % Taking the IFFT
ifftdata6=ifftdata_user66';
y5h6=[ifftdata6(:,[(n-2):n]) ifftdata6];
transdata666=y5h6';
esttx6=transdata666;  

%------------------est for UserD-------------------------------
estdata_userDbpskhop6 = 2*estDh5-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user3=estdata_userDbpskhop6';
spdata1_user3=estdata_user3*codeD;    % Spreading 
spdata3=(spdata1_user3)';
ifftdata_user33=ifft(spdata3);      % Taking the IFFT
ifftdata3=ifftdata_user33';
y3h5=[ifftdata3(:,[(n-2):n]) ifftdata3];
transdata333=y3h5';
esttxD=transdata333; 

%------------------est for User7-------------------------------
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

%----------------------Adding data for Transmission of All User------------
x7=esttx1+esttx2+esttx3+esttxA+esttx4+esttxB+esttx5+esttxC+esttx6+esttxD+tx_user7;

%rician fading channel 
gain16=sqrt(d6^-eta)*sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
gain26=sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain36=sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3
gain46=sqrt(p4/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap4
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
ch41=repmat(gain16,(n+3),1);     
ch42=repmat(gain26,(n+3),1);
ch43=repmat(gain36,(n+3),1);
ch44=repmat(gain46,(n+3),1);
data_channel6=x7.*ch41+x1111.*ch42+x2222.*ch43+x3333.*ch44;  % Passing data through channel 

%% ------user E using same code (code 7) but different channel 
%user E original data
data_userE = rand(1,M)>0.5;          % Generation of data for user2
data_userEbpsk = 2*data_userE-1;    % BPSK modulation 0 -> -1; 1 -> 0 
data_userEh6=data_userEbpsk';
spdataEh6=data_userEh6*code7;          % Spreading 
spdataEh6t=(spdataEh6)';
ifftdata_userEh6=ifft(spdataEh6t);      % Taking the IFFT
ifftdataEh6=ifftdata_userEh6';
yE=[ifftdataEh6(:,[(n-2):n]) ifftdataEh6];
transdataEh6=yE';
tx_userE=transdataEh6; 

%----------------------Adding data for Transmission of All User------------
xE=tx_userE;
data_channelE=sqrt(dA^-eta)*sqrt(p1/2)*xE;  % Passing data through channel ;

%----------------------Adding data for Transmission of All User------------

%transmitted signal of data of all users
data_channelhop6 = data_channel6+data_channelE;

%addition of noise 
data_noise6=data_channelhop6 (:);
data_noise66=reshape(data_noise6,1,length(data_noise6));
noise66 = 1/sqrt(2)*[randn(1,length(data_noise66)) + j*randn(1,length(data_noise66))];

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

%equalization over channel C
data_equilizedE=fft_data_received6/(sqrt(dA^-eta)*sqrt(p1/2));

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

%equalization over channel D
data_equilizedE=fft_data_received6/(sqrt(dA^-eta)*sqrt(p1/2));

%----------------BER of Data User1-----------------------------------------
recdata110=(data_equilized6'*code1')';
est1h6=real(recdata110)>0;
errors_user1h6(i) = size(find([data_user1- est1h6]),2); %Errors for User1
SBer1h6 = errors_user1h6/N;                              % simulated ber user1

%----------------BER of Data User2-----------------------------------------
recdata2221=(data_equilized6'*code3')';
est2h6=real(recdata2221)>0;
errors_user2h6(i) = size(find([data_user2- est2h6]),2); %Errors for User1
SBer2h6 = errors_user2h6/N;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata3=(data_equilized6'*code2')';
est3h6=real(recdata3)>0;
errors_user3h6(i) = size(find([data_user3- est3h6]),2); %Errors for User1
SBer3h6 = errors_user3h6/N;                               % simulated ber user2

%ber of user A
recdataA=(data_equilized6'*codeA')';
estAh6=real(recdataA)>0;
errors_userAh6(i) = size(find([data_userA-estAh6]),2); %Errors for User1
SBerAh6 = errors_userAh6/N;   

%----------------BER of Data User4-----------------------------------------
recdata4=(data_equilized6'*code4')';
est4h6=real(recdata4)>0;
errors_user4h6(i) = size(find([data_user4- est4h6]),2); %Errors for User1
SBer4h6 = errors_user4h6/N;      

%ber of user B
recdataB=(data_equilized6'*codeB')';
estBh6=real(recdataB)>0;
errors_userBh6(i) = size(find([data_userB- estBh6]),2); %Errors for User1
SBerBh6 = errors_userBh6/N;   

%----------------BER of Data User5-----------------------------------------
recdata5=(data_equilized6'*code5')';
est5h6=real(recdata5)>0;
errors_user5h6(i) = size(find([data_user5- est5h6]),2); %Errors for User1
SBer5h6= errors_user5h6/N; 

%ber of user C
recdataC=(data_equilized6'*codeC')';
estCh6=real(recdataC)>0;
errors_userCh6(i) = size(find([data_userC- estCh6]),2); %Errors for User1
SBerCh6 = errors_userCh6/N;   

%----------------BER of Data User6-----------------------------------------
recdata6=(data_equilized6'*code6')';
est6h6=real(recdata6)>0;
errors_user6h6(i) = size(find([data_user6- est6h6]),2); %Errors for User1
SBer6h6= errors_user6h6/N;   

%BER of data user D
recdataDh6=(data_equilized6'*codeD')';
estDh6=real(recdataDh6)>0;
errors_userDh6(i) = size(find([data_userD- estDh6]),2); %Errors for User1
SBerDh6= errors_userDh6/N;                             % simulated ber user1

%----------------BER SIC coding User 7 and user E -----------------------------------------
%BER of data user E
recdata7E=(data_equilizedE'*code7')';
xE_hat = ones(1,N);
xE_hat(recdata7E<0) = -1;
recdata7user=(data_equilized6'*code7'-xE_hat')';
est7h6=real(recdata7user)>0;
errors_user7(i) = size(find([data_user7- est7h6]),2); %Errors for User1
SBer7h6= errors_user7/N;   

x7user_hat = ones(1,N);
x7user_hat(recdata7user<0) = -1;
recdataE=(data_equilizedE'*code7'-x7user_hat')';
estEh6=real(recdataE)>0;
errors_userE(i) = size(find([data_userE- estEh6]),2); %Errors for User1
SBerEh6= errors_userE/N;  

end

end
