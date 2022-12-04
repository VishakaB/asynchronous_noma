%bemaforming vector
N =3;%number of bits

h1 = (randn(1,N) +1i*randn(1,N))/sqrt(2);
h2 = (randn(1,N) +1i*randn(1,N))/sqrt(2);
h3 = (randn(1,N) +1i*randn(1,N))/sqrt(2);

Hm = [h2;h3]';
invp2 = inv(conj(Hm')*Hm);
invp1 = Hm ;
invp3 = conj((Hm)');

PI = ones(3,3) - invp1*invp2*invp3;

num1 = PI*h1';
den1 = norm(PI*h1');

wm = num1/den1
y1 =abs(mean(mean(conj(h1').*wm)))

Hm = [h1;h3]';
invp2 = inv(conj(Hm')*Hm);
invp1 = Hm ;
invp3 = conj((Hm)');
PI2 = ones(3,3) - invp1*invp2*invp3;
num2 = PI2*h2';
den2 = norm(PI2*h2');
wj = num2/den2

y2 = abs(mean(mean(conj(h1').*wj)))