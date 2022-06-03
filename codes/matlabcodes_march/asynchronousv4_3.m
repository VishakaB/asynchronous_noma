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