function [data_user6,data_userD,SBer1h5,SBer2h5,SBer3h5,SBerAh5,SBer4h5,SBerBh5,...
    SBer5h5,SBerCh5,SBer6h5,SBerDh5,est1h5,est2h5,est3h5,estAh5,est4h5,estBh5,est5h5,estCh5,...
    est6h5,estDh5]= ...
    hop5(est1h4,est2h4,est3h4,estAh4,est4h4,estBh4,est5h4,estCh4,code1,code2,code3,codeA,...
    code4,codeB,code5,codeC,code6,codeD,n,N,d6,eta,p1,sigma,mean,dA,snr,data_user1,data_user2,...
    data_user3,data_userA,data_user4,data_userB,data_user5,data_userC)

%fprintf('okay\n')
p2=p1;
p3=p1;
p4=p1;
M=N;
%------------------est for User1-------------------------------
estdata_user1bpskhop3 = 2*est1h4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user111=estdata_user1bpskhop3';
spdata1_user111=estdata_user111*code1;    % Spreading 
spdata121=(spdata1_user111)';
ifftdata_user111=ifft(spdata121);      % Taking the IFFT
ifftdata12=ifftdata_user111';
y11=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata11=y11';
esttx1=transdata11;                % Transmitting data for user1

%------------------est for User2-------------------------------
estdata_user2bpskhop3 = 2*est2h4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user222=estdata_user2bpskhop3';
spdata1_user222=estdata_user222*code3;    % Spreading 
spdata222=(spdata1_user222)';
ifftdata_user22=ifft(spdata222);      % Taking the IFFT
ifftdata122=ifftdata_user22';
y222=[ifftdata122(:,[(n-2):n]) ifftdata122];
transdata222=y222';
esttx2=transdata222;                % Transmitting data for user1

%------------------est for User3-------------------------------
estdata_user3bpskhop3 = 2*est3h4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user333=estdata_user3bpskhop3';
spdata1_user333=estdata_user333*code2;    % Spreading 
spdata333=(spdata1_user333)';
ifftdata_user33=ifft(spdata333);      % Taking the IFFT
ifftdata33=ifftdata_user33';
y333=[ifftdata33(:,[(n-2):n]) ifftdata33];
transdata333=y333';
esttx3=transdata333;                % Transmitting data for user1

%------------------est for UserA-------------------------------
estdata_userAbpskhop3 = 2*estAh4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_userA=estdata_userAbpskhop3';
spdata1_userA=estdata_userA*codeA;    % Spreading 
spdataA=(spdata1_userA)';
ifftdata_userA=ifft(spdataA);      % Taking the IFFT
ifftdataA=ifftdata_userA';
yA=[ifftdataA(:,[(n-2):n]) ifftdataA];
transdataA=yA';
esttxA=transdataA;                % Transmitting data for user1

%user 4 estimated data %modulation 
estdata_user4bpskhop4 = 2*est4h4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user444=estdata_user4bpskhop4';
spdata1_user444=estdata_user444*code4;    % Spreading 
spdata444=(spdata1_user444)';
ifftdata_user44=ifft(spdata444);      % Taking the IFFT
ifftdata44=ifftdata_user44';
y444=[ifftdata44(:,[(n-2):n]) ifftdata44];
transdata444=y444';
esttx4=transdata444;                % Transmitting data for user1

%USER B ESTIMATED DATA
estdata_userBbpskhop4 = 2*estBh4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user444=estdata_userBbpskhop4';
spdata1_user444=estdata_user444*codeB;    % Spreading 
spdata444=(spdata1_user444)';
ifftdata_user44=ifft(spdata444);      % Taking the IFFT
ifftdata44=ifftdata_user44';
y444=[ifftdata44(:,[(n-2):n]) ifftdata44];
transdata444=y444';
esttxB=transdata444;             

%user 5 estimated data %modulation 
estdata_user5bpskhop5 = 2*est5h4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user5=estdata_user5bpskhop5';
spdata1_user5=estdata_user5*code5;    % Spreading 
spdata5=(spdata1_user5)';
ifftdata_user55=ifft(spdata5);      % Taking the IFFT
ifftdata5=ifftdata_user55';
y5h5=[ifftdata5(:,[(n-2):n]) ifftdata5];
transdata555=y5h5';
esttx5=transdata555;         

%user C estimated data modulation
estdata_userCbpskhop5 = 2*estCh4-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user5=estdata_userCbpskhop5';
spdata1_user5=estdata_user5*codeC;    % Spreading 
spdata5=(spdata1_user5)';
ifftdata_user55=ifft(spdata5);      % Taking the IFFT
ifftdata5=ifftdata_user55';
y5h5=[ifftdata5(:,[(n-2):n]) ifftdata5];
transdata555=y5h5';
esttxC=transdata555;  

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

%----------------------Adding data for Transmission of All User------------
x6=esttx1+esttx2+esttx3+esttxA+esttx4+esttxB+esttx5+esttxC+tx_user6;

%rician fading channel 
%rician fading channel 
%----------------------Creating a new Rayleigh Channel---------------------------
gain15=sqrt(d6^-eta)*sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
gain25=sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain35=sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3
gain45=sqrt(p4/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap4
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
ch41=repmat(gain15,(n+3),1);     
ch42=repmat(gain25,(n+3),1);
ch43=repmat(gain35,(n+3),1);
ch44=repmat(gain45,(n+3),1);
data_channel6=x6.*ch41+x1111.*ch42+x2222.*ch43+x3333.*ch44;  % Passing data through channel 

%% ------user D using same code (code 6) but different channel 
%user D original data
data_userD = rand(1,M)>0.5;          % Generation of data for user2
data_userDbpsk = 2*data_userD-1;    % BPSK modulation 0 -> -1; 1 -> 0 
data_userDh5=data_userDbpsk';
spdataDh5=data_userDh5*code6;          % Spreading 
spdataDh5t=(spdataDh5)';
ifftdata_userDh5=ifft(spdataDh5t);      % Taking the IFFT
ifftdataDh5=ifftdata_userDh5';
yD=[ifftdataDh5(:,[(n-2):n]) ifftdataDh5];
transdataDh5=yD';
tx_userD=transdataDh5; 

%----------------------Adding data for Transmission of All User------------
xD=tx_userD;
data_channelD=sqrt(dA^-eta)*sqrt(p1/2)*xD;  % Passing data through channel ;

%----------------------Adding data for Transmission of All User------------

%transmitted signal of data of all users
data_channelhop5 = data_channel6+data_channelD;

%addition of noise 
data_noise5=data_channelhop5 (:);
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

%equalization over channel D
data_equilizedD=fft_data_received5/(sqrt(dA^-eta)*sqrt(p1/2));

%----------------BER of Data User1-----------------------------------------
recdata110=(data_equilized5'*code1')';
est1h5=real(recdata110)>0;
errors_user1h5(i) = size(find([data_user1- est1h5]),2); %Errors for User1
SBer1h5 = errors_user1h5/N;                              % simulated ber user1

%----------------BER of Data User2-----------------------------------------
recdata2221=(data_equilized5'*code3')';
est2h5=real(recdata2221)>0;
errors_user2h5(i) = size(find([data_user2- est2h5]),2); %Errors for User1
SBer2h5 = errors_user2h5/N;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata3=(data_equilized5'*code2')';
est3h5=real(recdata3)>0;
errors_user3h5(i) = size(find([data_user3- est3h5]),2); %Errors for User1
SBer3h5 = errors_user3h5/N;                               % simulated ber user2

%ber of user A
recdataA=(data_equilized5'*codeA')';
estAh5=real(recdataA)>0;
errors_userAh5(i) = size(find([data_userA-estAh5]),2); %Errors for User1
SBerAh5 = errors_userAh5/N;   

%----------------BER of Data User4-----------------------------------------
recdata4=(data_equilized5'*code4')';
est4h5=real(recdata4)>0;
errors_user4h5(i) = size(find([data_user4- est4h5]),2); %Errors for User1
SBer4h5 = errors_user4h5/N;      

%ber of user B
recdataB=(data_equilized5'*codeB')';
estBh5=real(recdataB)>0;
errors_userBh5(i) = size(find([data_userB- estBh5]),2); %Errors for User1
SBerBh5 = errors_userBh5/N;   

%----------------BER of Data User5-----------------------------------------
recdata5=(data_equilized5'*code5')';
est5h5=real(recdata5)>0;
errors_user5h5(i) = size(find([data_user5- est5h5]),2); %Errors for User1
SBer5h5= errors_user5h5/N; 

%ber of user C
recdataC=(data_equilized5'*codeC')';
estCh5=real(recdataC)>0;
errors_userCh5(i) = size(find([data_userC- estCh5]),2); %Errors for User1
SBerCh5 = errors_userCh5/N;   

%----------------BER of Data SIC decoding user6 and user D-----------------------------------------
%ber of user D
recdata6D=(data_equilizedD'*code6')';
xD_hat = ones(1,N);
xD_hat(recdata6D<0) = -1;
recdata5user=(data_equilized5'*code6'-xD_hat')';
est6h5=real(recdata5user)>0;
errors_user6(i) = size(find([data_user6- est6h5]),2); %Errors for User1
SBer6h5= errors_user6/M;   

x5user_hat = ones(1,N);
x5user_hat(recdata5user<0) = -1;
recdataD=(data_equilizedD'*code6'-x5user_hat')';
estDh5=real(recdataD)>0;
errors_userD(i) = size(find([data_userD- estDh5]),2); %Errors for User1
SBerDh5 = errors_userD/N;                              % simulated ber user1

end

end
