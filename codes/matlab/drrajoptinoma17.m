%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% OPTIMAL POWER ALLOCATION NOMA %%%%%%%%%%%%%%
%%%%%%%%%%%% GRADIENT DECENT MINIMIZATION %%%%%%%%%%%%%%%%
clear all;
close all;
clc;
W=1;% System Bandwidth 
sig=1;

%%%%%%%%%%%%%SELECT MODEL %%%%%%%%%%%%%%%%
k=1;
%%%%% INITIALIZE LAGRANGIAN MULTIPLIERS %%%%%%%%%%%%%
lm_opt1(1)=0;
lm_opt2(1)=0;
mu_opt(1)=0;
p1(1)=0;p2(1)=0;

%%%%%%%%%%%% STOP CRITERIAN  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ep2=.01;

%%%%%%%%%%%% STEP SIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1(1)=.5;t2=.1;q(1)=.1;p11(1)=.1;p12(1)=.1;

for i1=1:1:10
%%%%%%%%%%%%%%%  Channel Coefficients %%%%%%%%%%%%%%%
h1=1.1+(0.5)*(abs(randn(1,1)+j*randn(1,1)))^2; %% 2.1412 2.0554 2.5872 1.8351 2.3413 Near User
h2=(0.5)*(abs(randn(1,1)+j*randn(1,1)))^2; %% 0.2935 0.6773 0.2439 0.6334 1.2899 Far User1
% h1=2.9293;h2=0.0819;
if h1>h2
    break
end
end
    
% h1=1.1;h2=1;  

%%%%%% ITERATION PROCESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
itr=1000;
for k=1:1:2
for snr=2:1:22
    ep1=snr;
    c=1;
    n=1;
for i=1:1:itr
    
%%%% Power Updation
%%% p1 far user (high power)
%p2 near user (low power)
p1_opt(n)= ((W*(1.4427)*(h2/((p1(n)*h2)+(p2(n)*h2)+sig)))*(1+0))...
    -mu_opt(n);
        
p2_opt(n)=(((W*(1.4427)*((h1/((p2(n)*h1)+sig))+(h2/((p1(n)*h2)...
    +(p2(n)*h2)+sig))-(h2/((p2(n)*h2)+sig)))))*(1+0))-mu_opt(n);

n=n+1;
p1(n)=max(0,(p1(n-1)+(t1*p1_opt(n-1))));
p2(n)=max(0,(p2(n-1)+(t1*p2_opt(n-1))));
P3=(p1(n)+p2(n));
if P3 < ep1
C11=(W*(1.4427)*log(1+((p1(n)*h1)/sig)));%% capacity of near user
C12=(W*(1.4427)*log(1+((p2(n)*h2)/((p1(n)*h2)+sig))));%% capacity of far user
C13=(p1(n)+p2(n));%% Updated power
C1=(C11-(c*W)); %% Constraint one
C2=(C12-(c*W));%% Constraint two
C3=(snr-C13);%% Constraint three
fprintf('Iteration: %d\n',i);
fprintf('-----------------------------------------------------\n');
fprintf(' Gradient1: %f\n ',p1_opt(n-1));
fprintf(' Gradient2: %f\n ',p2_opt(n-1));
fprintf(' Updated Power Far User:a= %f\n ',p1(n));
fprintf(' Updated Power Near User:b= %f\n ',p2(n));
fprintf(' Capacity of near user= C_sum1: %f\n ',C11);
fprintf(' Capacity of far user= C_sum2: %f\n ',C12);
fprintf(' Updatd Power= p1+p2: %f\n ',C13);
fprintf(' Contraint One near user = C_sum-C_min(approch to zero): %f\n ',C1);
fprintf(' Contraint Two far user =C_sum-C_min(approch to zero): %f\n ',C2);
fprintf(' Contraint Three =snr-(p1+p2)(approch to zero): %f\n ',C3);
lm_opt1(n)=max(lm_opt1(n-1)+p11(n-1)*C1);
lm_opt2(n)=max(lm_opt2(n-1)+p12(n-1)*C2);
mu_opt(n)=max(mu_opt(n-1)+q(n-1)*C3);
if lm_opt1(n-1)~=0
    p11(n)=-lm_opt1(n-1)/C1;
    p12(n)=-lm_opt2(n-1)/C2;
%         p11(n)=.1;
%     p12(n)=.1;
else
    p11(n)=.1;
    p12(n)=.1;
end
if mu_opt(n-1)~=0
   q(n)=-mu_opt(n-1)/(snr-(p1(n)+p2(n)));
%        q(n)=.1;
else
    q(n)=.1;
end
C4=abs((lm_opt1(n)-lm_opt1(n-1))+(lm_opt2(n)-lm_opt2(n-1)));
C5=abs(mu_opt(n)-mu_opt(n-1));
fprintf('-----------------------------------------------------\n');
fprintf(' Convergence One: %f\n ',C4);
fprintf(' Convergence Two: %f\n ',C5);
fprintf(' Convergence Three: %f\n ',C13);
fprintf(' Updated lamda One: %f\n ',lm_opt1(n));
fprintf(' Updated lamda Two: %f\n ',lm_opt2(n));
fprintf(' Updated Mu: %f\n ',mu_opt(n));
fprintf(' Lagrangian multiplier p11:= %f\n ',p11(n));
fprintf(' Lagrangian multiplier p11:= %f\n ',p12(n));
fprintf(' Lagrangian multiplier q:= %f\n ',q(n));
% if C4<=ep2 && C5<=ep2
%      break;
% end
else
    break
end
end
if k==1
    P4=(snr-P3);
    fprintf('Iteration: %d\n',i);
    fprintf('-----------------------------------------------------\n');
    fprintf(' Contraint Three =snr-(p1+p2)(approch to zero): %f\n ',P4);
    a=p2(n);%far user high power
    b=p1(n);%near user less power
    c_sum(snr)=(W*(1.4427)*(log(1+((b*h1)/sig))+log(1+((a*h2)/((b*h2)+sig)))));
else
    a=.8*snr;
    b=.2*snr;
    c_sum2(snr)=(W*(1.4427)*(log(1+((b*h1)/sig))+log(1+((a*h2)/((b*h2)+sig)))));
end
end
end

plot((c_sum));
axis([2 20 0 6]);
grid on;
hold on;
plot((c_sum2));
axis([2 20 0 6]);
xlabel(' Transmit Power Constraint (dB)')
ylabel ('Sum Capacity (b/s/Hz)')
