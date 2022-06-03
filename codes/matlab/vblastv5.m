clc;
clear all;
close all;
N = 10^5;
M=4;
x = randi([0 3],1,N);
xmod=qammod(x,4);
xmod=reshape(xmod,8,N/8);
xmod=kron(xmod,[1,1,1,1,1,1,1,1]);
for i=1:64
    h(i,:)=1/sqrt(2)*(randn(1,N/8) + 1i*randn(1,N/8));
end
H=reshape(h,8,N);
y=reshape(sum(H.*xmod,1),8,N/8);
H=reshape(h,8,8,N/8);
snr=linspace(0,18,19);
ber=zeros(1,length(snr));
for ii=1:length(snr)
    N1=1/sqrt(2)*(randn(1,N)+1i*randn(1,N));
    N1=reshape(N1,8,N/8);
    ynoisy=y+10^(-(snr(ii)-10*log10(16))/20)*N1;
    ynoisy=reshape(ynoisy,8,1,N/8);
    finy=[];
    for kk=1:N/8
        c=zeros(8,1);
        r=ynoisy(:,:,kk);
        Heq=transpose(H(:,:,kk));
        for jj=1:8
            B=pinv(Heq);
            b=abs(B.*B);
            b=b(:,1)+b(:,2)+b(:,3)+b(:,4);
            b=b+c;
            [ss, dd] = min(b,[],1);
            c(dd)=100;
            temp=qamdemod(B(dd,:)*r,4);
            ydemod(dd)=temp;
            r=r-qammod(temp,4).*Heq(:,dd);
            Heq(:,dd)=0;
        end
        finy=[finy,ydemod];
    end
    [num ty]=symerr(x,finy);
    ber(ii)=ty/N;
    ty;
end
semilogy(snr,ber,'g-*');
grid on;hold on;
title('Plot of Bit error rate for Uncoded 8X8(4-QAM) System','FontSize',12);
legend('sim (nTx=8, nRx=8, Uncoded(4-QAM)) Using V-Blast/ZF','location','southwest');
xlabel('SNR(dB) ---->','Color','k','FontSize',11);
ylabel('Symbol Error rate ---->','Color','k','FontSize',11);