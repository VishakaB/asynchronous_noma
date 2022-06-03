%asynchronous % imperfect time syncrhonization between different signals 
k1=5;
mean=sqrt(k1/(k1+1));% mean
sigma=sqrt(1/(2*(k1+1)));%variance
p1=1;                                       % Power of Tap1
p2=0.5;                                     % Power of Tap2
p3=0.3;                                     % Power of Tap3

%% data generation 
nbstreams = 10;%number of users
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
   userdatawdelay(nn,i+delay(nn))=x0(i);
   %disp(x0new);
   Hwdelay(nn,i+delay(nn)) = H_arr(nn,i);
end   
end

data_channel=gainnew.*x0new;  % Passing data through channel 
superimposeddatachannel =sum(data_channel);

%add AWGN
noise = 1/sqrt(2)*[randn(1,length(superimposeddatachannel)) + j*randn(1,length(superimposeddatachannel))]; 

save data_user.mat;
save data_userbpsk.mat;
save time_offsetsarr.mat;
save delay.mat;
save P_transmitted.mat;
save userdatawdelay.mat;
save Hwdelay.mat;
save superimposeddatachannel.mat;
save noise.mat;

n_iter = 1;%number of iterations
e_th = 25;%energy threshold
A = ones(nbstreams);
max_n_sym = sum(sum(tril(A)));

alpha = 4;%pathloss exponent
Pt = 100;%maximum transmission power
d0 = 1500;%distance radius in meters
power_farthest = sqrt(Pt/d0)^alpha;%received power at the far end of connection range 
N0 = 1;
SNR = power_farthest/N0;
rate_thresh = log(1+SNR);%fixed threshold based on communication radius
bandwidth = 10^6;

%time offsets for each user
deltak = sort(randn(1,N));%time offsets in the increasing order
user_link_rate = rate_user_link(P_transmitted,H_arr,nbstreams,N0,delay,bandwidth,1);

%% optimization algorithm
if N ~=0
       
    cvx_begin quiet
    
    variable n_sym integer 
    variable uk(nbstreams,nbstreams) binary
    
    %objective
    maximize(n_sym)
    
    %constraints
    subject to
    
    0.00319*n_sym*n_iter <= log(e_th);%energy threshold
    n_sym >= 2;
    n_sym <= max_n_sym;
    
    sum(uk,1)*(abs(user_link_rate) - rate_thresh) >= 0;%sinr threshold
    
    uk <= tril(A);
    sum(sum(uk)) == n_sym;
    
    cvx_end

    n_sym;
    uk
       
end
