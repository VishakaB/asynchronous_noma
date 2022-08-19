%{channel_gain  
%min_sinr  = %2Mbits
%signal = channel gain
%SINR = signal_power/(interference+noise)
%}
clear all
close all 
clc

for i = 1:10
    
dfar = 23;
cvx_begin quiet
   variable dnear
   dual variables var1 var2 var3
   minimize(-log(1+dnear)*dfar)%maximize secrecy
   subject to
      var1: dnear/dfar >= 0.1;
      var2: dnear <= 10;
      var3: dfar <= 20;
cvx_end
distancegap(i) = dfar - dnear;
capacity(i) = log(1+dnear)*dfar;
end

figure(1)
plot(distancegap,capacity,'k--')