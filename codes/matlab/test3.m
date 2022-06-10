clear all

N = 10^1;                             % Number of Bits for  data_user1
qam_order = 16;
data_user1 = randi([0 1],N,1);       % Generation of data for user1
data_user1qam = qammod(data_user1,qam_order);    % 4 QAM 0 -> -1; 1 -> 0 

x11=data_user1qam(:);
x12=reshape(x11,1,length(x11));

i = 1:length(x12);        
delay1=1; 

%modulated signal % delayed %time offset
for i=delay1+1:length(x12) % Producing one sample delay in Tap2 w.r.t. Tap1
   x13(i)=data_user1qam(i-delay1);
end

%original data delayed
for i=delay1+1:length(x12) % Producing one sample delay in Tap2 w.r.t. Tap1
   x_data1(i)=data_user1(i-delay1);
end

K = 1;
max_tx_power = 100;
max_dist     = 500;%meters
dist_k = max_dist*abs(randn(K,1));

%sorted transmit power vector %descending 
dist_vec = sort(dist_k,'ascend'); 

%rayleigh fading 
max_eta = 5;

%Path loss exponent
eta_k   = max_eta*abs(randn(K,1));
eta_vec = sort(eta_k,'ascend'); 

%unsorted transmit power vector
transmitpow_k = max_tx_power*abs(randn(K,1));

%sorted transmit power vector %descending 
power_vec = sort(transmitpow_k,'descend'); 

%tx power vec
%power_vec = tx_power_percentage_vec.*tx_pow_k;
pathloss_exp = sqrt(dist_vec.^-eta_vec);

%channel coefficients of each user vec
h_vec =  pathloss_exp.*sqrt(power_vec/2).*(randn(1,N)+1i*randn(1,N))/sqrt(2);

%superimposed data
super_signal = h_vec.*x13;

%AWGN noise 
noise= 1/sqrt(2)*[randn(K,length(super_signal)) + j*randn(K,length(super_signal))];

snr = [0:20];                 % multiple Eb/N0 values

for i = 1:length(snr)
%received signal 
y = super_signal +(sqrt(1)*10^(-snr(i)/20))*noise;

%equalization 
eq_vec = y./h_vec;

%ber analysis 
est_sym = qamdemod(eq_vec,qam_order);
data_user1;
ber = size(find([x_data1- est_sym]),2);

SBer(i) = ber/N;
end

SBer