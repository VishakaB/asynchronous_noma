
%asynchronous % imperfect time syncrhonization between different signals 
k1=5;
mean=sqrt(k1/(k1+1));% mean
sigma=sqrt(1/(2*(k1+1)));%variance
p1=1;                                       % Power of Tap1
p2=0.5;                                     % Power of Tap2
p3=0.3;                                     % Power of Tap3

%% data generation 
nbstreams = 3;%number of users
N= 3; %number of bits per user

time_offsets = [1,2,3];
for m = 1:nbstreams  
    if m>1
        time_offsetsarr(m) = time_offsetsarr(m-1)+randi(5);
    else
        time_offsetsarr(m) = randi(5);  
    end
end

for k =1:nbstreams
data_user(k,:) = rand(1,N)>0.5;%row=nbuser, coloumn= data bits 
data_userbpsk(k,:) = 2*data_user(k,:)-1;
end

%power transmitted
for l=1:nbstreams
 
 if l>1
    s0(l) = s0(l-1)+randi(1);
    P_transmitted(l) = p1/s0(l);
 else
    s0(l) = randi(5);
    P_transmitted(l) = p1;
 end
end

for kkk=1:nbstreams
 H_arr(kkk,:) = sqrt(P_transmitted(kkk)/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma]; 
end

for nn=1:nbstreams
x0=data_user(nn,:);%row=nbuser, coloumn= data bits 
x0reshape=reshape(x0,1,length(x0));      
if (nn>1)
  delay(nn)=1+delay(nn-1); 
else
  delay(nn)=0;  
end
for i=1:length(x0reshape) % Producing one sample delay in Tap2 w.r.t. Tap1
   x0new(nn,i+delay(nn))=x0(i);
   disp(x0new)
   gainnew(nn,i+delay(nn)) = H_arr(nn,i);
end   
end

%user 1 data
data_user1 = rand(1,N)>0.5;          % Generation of data for user1
data_user1bpsk = 2*data_user1-1;    % BPSK modulation 0 -> -1; 1 -> 0 

%user 2 data 
data_user2 = rand(1,N)>0.5;          % Generation of data for user1
data_user2bpsk = 2*data_user2-1;    % BPSK modulation 0 -> -1; 1 -> 0 

%user 3 data 
data_user3 = rand(1,N)>0.5;          % Generation of data for user1
data_user3bpsk = 2*data_user3-1;    % BPSK modulation 0 -> -1; 1 -> 0 

% user data generation, automate

%% time delay of receiving data 
gain1s=sqrt(p1/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap1
gain2s=sqrt(p2/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap2
gain3s=sqrt(p3/2)*[randn(1,N)*sigma+mean + j*randn(1,N)*sigma];   % Gain for Tap3

H = zeros(nbstreams,N);
P_received = zeros(1,nbstreams);
P_received =[sum(abs(gain1s)),sum(abs(gain2s)),sum(abs(gain3s))];
P_transmitted = [p1/2,p2/2,p3/2];
H = [sum(abs(gain1s)),sum(abs(gain2s)),sum(abs(gain3s))];

save P_transmitted.mat;
save P_received.mat;
save H.mat;

x2=data_user2(:);
x2reshape=reshape(x2,1,length(x2));
i=1:length(x2reshape);        
delay1=1; 
for i=1:length(x2reshape) % Producing one sample delay in Tap2 w.r.t. Tap1
   x2new(i+delay1)=x2(i);
   gain2new(i+delay1) = gain2s(i);
end

x3=data_user3(:);
x3reshape=reshape(x3,1,length(x3));
i=1:length(x3reshape);   
delay2=2;
for i=1:length(x3reshape) % Producing one sample delay in Tap2 w.r.t. Tap1
   x3new(i+delay2)=x3(i);
   gain3new(i+delay2) = gain3s(i);
end

%% transmission across channel 
data_user1(numel(x3new)) = 0;
gain1s(numel(x3new)) = 0;
x2new(numel(x3new)) = 0;
gain2new(numel(x3new)) = 0;
data_channel=data_user1.*gain1s+x2new.*gain2new+x3new.*gain3new;  % Passing data through channel 

%add AWGN
noise = 1/sqrt(2)*[randn(1,length(data_channel)) + j*randn(1,length(data_channel))]; 
