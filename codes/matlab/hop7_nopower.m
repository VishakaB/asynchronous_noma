function [data_user8,data_userF,SBer1h7,SBer2h7,SBer3h7,SBerAh7,SBer4h7,SBerBh7,...
    SBer5h7,SBerCh7,SBer6h7,SBerDh7,SBer7h7,SBerEh7,SBer8h7,SBerFh7,est1h6,est2h6,est3h6,estAh6,est4h6,estBh6,est5h6,estCh6,...
    est6h6,estDh6,est7h6,estEh6,est8h7,estFh7]=hop7_nopower(est1h6,est2h6,est3h6,estAh6,est4h6,estBh6,est5h6,estCh6,...
    est6h6,estDh6,est7h6,estEh6,code1,code2,code3,codeA,...
    code4,codeB,code5,codeC,code6,codeD,code7,codeE,code8,codeF,n,N,d6,eta,p1,sigma,mean,dA,snr,data_user1,data_user2,...
    data_user3,data_userA,data_user4,data_userB,data_user5,data_userC,data_user6,data_userD,data_user7,data_userE)
%fprintf('okay\n')
p2=p1;
p3=p1;
p4=p1;
M=N;
%------------------est for User1-------------------------------
estdata_user1bpskhop3 = 2*est1h6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user111=estdata_user1bpskhop3';
spdata1_user111=estdata_user111*code1;    % Spreading 
spdata121=(spdata1_user111)';
ifftdata_user111=ifft(spdata121);      % Taking the IFFT
ifftdata12=ifftdata_user111';
y11=[ifftdata12(:,[(n-2):n]) ifftdata12];
transdata11=y11';
esttx1=transdata11;                % Transmitting data for user1

%------------------est for User2-------------------------------
estdata_user2bpskhop3 = 2*est2h6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user222=estdata_user2bpskhop3';
spdata1_user222=estdata_user222*code3;    % Spreading 
spdata222=(spdata1_user222)';
ifftdata_user22=ifft(spdata222);      % Taking the IFFT
ifftdata122=ifftdata_user22';
y222=[ifftdata122(:,[(n-2):n]) ifftdata122];
transdata222=y222';
esttx2=transdata222;                % Transmitting data for user1

%------------------est for User3-------------------------------
estdata_user3bpskhop3 = 2*est3h6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user333=estdata_user3bpskhop3';
spdata1_user333=estdata_user333*code2;    % Spreading 
spdata333=(spdata1_user333)';
ifftdata_user33=ifft(spdata333);      % Taking the IFFT
ifftdata33=ifftdata_user33';
y333=[ifftdata33(:,[(n-2):n]) ifftdata33];
transdata333=y333';
esttx3=transdata333;                % Transmitting data for user1

%------------------est for UserA-------------------------------
estdata_userAbpskhop3 = 2*estAh6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_userA=estdata_userAbpskhop3';
spdata1_userA=estdata_userA*codeA;    % Spreading 
spdataA=(spdata1_userA)';
ifftdata_userA=ifft(spdataA);      % Taking the IFFT
ifftdataA=ifftdata_userA';
yA=[ifftdataA(:,[(n-2):n]) ifftdataA];
transdataA=yA';
esttxA=transdataA;                % Transmitting data for user1

%------------------est for User4-------------------------------
estdata_user4bpskhop4 = 2*est4h6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user444=estdata_user4bpskhop4';
spdata1_user444=estdata_user444*code4;    % Spreading 
spdata444=(spdata1_user444)';
ifftdata_user44=ifft(spdata444);      % Taking the IFFT
ifftdata44=ifftdata_user44';
y444=[ifftdata44(:,[(n-2):n]) ifftdata44];
transdata444=y444';
esttx4=transdata444;                % Transmitting data for user1

%------------------est for UserB-------------------------------
estdata_userBbpskhop4 = 2*estBh6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user444=estdata_userBbpskhop4';
spdata1_user444=estdata_user444*codeB;    % Spreading 
spdata444=(spdata1_user444)';
ifftdata_user44=ifft(spdata444);      % Taking the IFFT
ifftdata44=ifftdata_user44';
y444=[ifftdata44(:,[(n-2):n]) ifftdata44];
transdata444=y444';
esttxB=transdata444;             

%------------------est for User5-------------------------------
estdata_user5bpskhop5 = 2*est5h6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user5=estdata_user5bpskhop5';
spdata1_user5=estdata_user5*code5;    % Spreading 
spdata5=(spdata1_user5)';
ifftdata_user55=ifft(spdata5);      % Taking the IFFT
ifftdata5=ifftdata_user55';
y5h5=[ifftdata5(:,[(n-2):n]) ifftdata5];
transdata555=y5h5';
esttx5=transdata555;         

%------------------est for UserC-------------------------------
estdata_userCbpskhop5 = 2*estCh6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user5=estdata_userCbpskhop5';
spdata1_user5=estdata_user5*codeC;    % Spreading 
spdata5=(spdata1_user5)';
ifftdata_user55=ifft(spdata5);      % Taking the IFFT
ifftdata5=ifftdata_user55';
y5h5=[ifftdata5(:,[(n-2):n]) ifftdata5];
transdata555=y5h5';
esttxC=transdata555;  

%------------------est for User6-------------------------------
estdata_user6bpskhop6 = 2*est6h6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user6=estdata_user6bpskhop6';
spdata1_user6=estdata_user6*code6;    % Spreading 
spdata6=(spdata1_user6)';
ifftdata_user66=ifft(spdata6);      % Taking the IFFT
ifftdata6=ifftdata_user66';
y5h6=[ifftdata6(:,[(n-2):n]) ifftdata6];
transdata666=y5h6';
esttx6=transdata666;  

%------------------est for UserD-------------------------------
estdata_userDbpskhop6 = 2*estDh6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user3=estdata_userDbpskhop6';
spdata1_user3=estdata_user3*codeD;    % Spreading 
spdata3=(spdata1_user3)';
ifftdata_user33=ifft(spdata3);      % Taking the IFFT
ifftdata3=ifftdata_user33';
y3h5=[ifftdata3(:,[(n-2):n]) ifftdata3];
transdata333=y3h5';
esttxD=transdata333; 

%------------------est for UserD-------------------------------
estdata_user7bpskhop6 = 2*est7h6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user3=estdata_user7bpskhop6';
spdata1_user3=estdata_user3*code7;    % Spreading 
spdata3=(spdata1_user3)';
ifftdata_user33=ifft(spdata3);      % Taking the IFFT
ifftdata3=ifftdata_user33';
y3h5=[ifftdata3(:,[(n-2):n]) ifftdata3];
transdata333=y3h5';
esttx7=transdata333; 

%------------------est for UserD-------------------------------
estdata_userEbpskhop6 = 2*estEh6-1;    % BPSK modulation 0 -> -1; 1 -> 0 
estdata_user3=estdata_userEbpskhop6';
spdata1_user3=estdata_user3*codeE;    % Spreading 
spdata3=(spdata1_user3)';
ifftdata_user33=ifft(spdata3);      % Taking the IFFT
ifftdata3=ifftdata_user33';
y3h5=[ifftdata3(:,[(n-2):n]) ifftdata3];
transdata333=y3h5';
esttxE=transdata333; 

%------------------est for User8-------------------------------
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

%----------------------Adding data for Transmission of All User------------
x8=esttx1+esttx2+esttx3+esttxA+esttx4+esttxB+esttx5+esttxC...
    +esttx6+esttxD+esttx7+esttxE+tx_user8;

%rician fading channel 
%rician fading channel 
%----------------------Creating a new Rayleigh Channel---------------------------
gain18=sqrt(d6^-eta)*sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
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
data_channel8=x8.*ch41+x1111.*ch42+x2222.*ch43+x3333.*ch44;  % Passing data through channel 


%% ------user E using same code (code 8) but different channel 
%user E original data
data_userF = rand(1,M)>0.5;          % Generation of data for user2
data_userFbpsk = 2*data_userF-1;    % BPSK modulation 0 -> -1; 1 -> 0 
data_userFh6=data_userFbpsk';
spdataFh6=data_userFh6*code8;          % Spreading 
spdataFh6t=(spdataFh6)';
ifftdata_userFh6=ifft(spdataFh6t);      % Taking the IFFT
ifftdataFh6=ifftdata_userFh6';
yF=[ifftdataFh6(:,[(n-2):n]) ifftdataFh6];
transdataFh6=yF';
tx_userF=transdataFh6; 

%----------------------Adding data for Transmission of All User------------
xF=tx_userF;
data_channelF=sqrt(dA^-eta)*sqrt(p1/2)*xF;  % Passing data through channel ;

%----------------------Adding data for Transmission of All User------------

%transmitted signal of data of all users
data_channelhop7 = data_channel8+data_channelF;

%addition of noise 
data_noise7=data_channelhop7 (:);
data_noise77=reshape(data_noise7,1,length(data_noise7));
noise77 = 1/sqrt(2)*[randn(1,length(data_noise77)) + j*randn(1,length(data_noise77))];

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
channel_response7=fft([gain18;gain28;gain38;gain48],n);
data_equilized7=fft_data_received7.*conj(channel_response7);

%equalization over channel C
data_equilizedF=fft_data_received7/(sqrt(dA^-eta)*sqrt(p1/2));

%----------------BER of Data User1-----------------------------------------
recdata110=(data_equilized7'*code1')';
est1h7=real(recdata110)>0;
errors_user1h6(i) = size(find([data_user1- est1h7]),2); %Errors for User1
SBer1h7 = errors_user1h6/N;                              % simulated ber user1

%----------------BER of Data User2-----------------------------------------
recdata2221=(data_equilized7'*code3')';
est2h7=real(recdata2221)>0;
errors_user2h6(i) = size(find([data_user2- est2h7]),2); %Errors for User1
SBer2h7 = errors_user2h6/N;                               % simulated ber user2

%----------------BER of Data User3-----------------------------------------
recdata3=(data_equilized7'*code2')';
est3h7=real(recdata3)>0;
errors_user3h6(i) = size(find([data_user3- est3h7]),2); %Errors for User1
SBer3h7 = errors_user3h6/N;                               % simulated ber user2

%ber of user A
recdataA=(data_equilized7'*codeA')';
estAh7=real(recdataA)>0;
errors_userAh6(i) = size(find([data_userA-estAh7]),2); %Errors for User1
SBerAh7 = errors_userAh6/N;   

%----------------BER of Data User4-----------------------------------------
recdata4=(data_equilized7'*code4')';
est4h7=real(recdata4)>0;
errors_user4h6(i) = size(find([data_user4- est4h7]),2); %Errors for User1
SBer4h7 = errors_user4h6/N;      

%ber of user B
recdataB=(data_equilized7'*codeB')';
estBh7=real(recdataB)>0;
errors_userBh6(i) = size(find([data_userB- estBh7]),2); %Errors for User1
SBerBh7 = errors_userBh6/N;   

%----------------BER of Data User5-----------------------------------------
recdata5=(data_equilized7'*code5')';
est5h7=real(recdata5)>0;
errors_user5h6(i) = size(find([data_user5- est5h7]),2); %Errors for User1
SBer5h7= errors_user5h6/N; 

%ber of user C
recdataC=(data_equilized7'*codeC')';
estCh7=real(recdataC)>0;
errors_userCh6(i) = size(find([data_userC- estCh7]),2); %Errors for User1
SBerCh7 = errors_userCh6/N;   

%----------------BER of Data User6-----------------------------------------
recdata6=(data_equilized7'*code6')';
est6h7=real(recdata6)>0;
errors_user6h6(i) = size(find([data_user6- est6h7]),2); %Errors for User1
SBer6h7= errors_user6h6/N;   

%BER of data user D
recdataDh6=(data_equilized7'*codeD')';
estDh7=real(recdataDh6)>0;
errors_userDh6(i) = size(find([data_userD- estDh7]),2); %Errors for User1
SBerDh7= errors_userDh6/N;                             % simulated ber user1

%----------------BER of Data User7-----------------------------------------
recdata7h7=(data_equilized7'*code7')';
est7h7=real(recdata7h7)>0;
errors_user7h7(i) = size(find([data_user7- est7h7]),2); %Errors for User1
SBer7h7= errors_user7h7/N;      

%BER of data user E
recdataEh7=(data_equilized7'*codeE')';
estEh7=real(recdataEh7)>0;
errors_userEh7(i) = size(find([data_userE- estEh7]),2); %Errors for User1
SBerEh7= errors_userEh7/N;

%----------------BER SIC coding User 8 and user F -----------------------------------------
%BER of data user E
recdata8user=(data_equilized7'*code8')';
est8h7=real(recdata8user)>0;
errors_user7(i) = size(find([data_user8- est8h7]),2); %Errors for User1
SBer8h7= errors_user7/N;   

recdataF=(data_equilizedF'*code8')';
estFh7=real(recdataF)>0;
errors_userF(i) = size(find([data_userF- estFh7]),2); %Errors for User1
SBerFh7= errors_userF/N;  

end

end
