%block error rate 
%blocks define with a fixed length : 3
N = 18;
transmitblocks = [100;010;100;111;010;101];

%transmission under channel fading


%addition of awgn 


%received blocks
receiveblocks = [101;011;100;110;010;101];
%check error blocks using same?
nbeblocks  = 6; 
totalblocks = 6;
avg_bler = nbeblocks/totalblocks;