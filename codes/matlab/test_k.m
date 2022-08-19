
k_vec = 1:0.5:10
snr = 1:10

ber_mat = zeros(3,10)

for i= 1: length(k_vec)
    for j= 1: length(snr)
        %zero forcing
        avg_ber(j) = 0.05
        ber_mat(i,j)=avg_ber(j)%row= k, coloum = ber
    end
end

figure(1)
for h= 1:length(k_vec)
plot(snr,ber_mat(h,:))
hold on;
end