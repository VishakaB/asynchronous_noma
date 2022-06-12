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

else
    disp('yes')
    fprintf('P3: %f\n',P3);
    fprintf('iter: %i\n',itr);
    break
end
end

end
end


