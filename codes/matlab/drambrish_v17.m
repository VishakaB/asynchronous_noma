%date: 13 july 2022
%goal: optimize the secrecy 
%for irc
close all
%clc
%clear all

%plotting 

for cases = 1:4

C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colros.
%% performance analysis
%secrecy capacity vs SINRd
hold on;
grid on;
plot(transmit_snrdb_vec,smooth(smooth(smooth(cs_smooth(cases,:)))),'color',C{cases},'marker','o');
hold on;
plot(transmit_snrdb_vec,smooth(smooth(smooth(wtcs_smooth(cases,:)))),'color',C{cases},'marker','*');
xlabel('SNR legitimate');
ylabel('Secrecy capacity');

end

legend({'opt \alpha_{{e}_1} = 0.001';'wt opt \alpha_{{e}_1} = 0.001';...
    'opt \alpha_{{e}_2} = 0.01';'wt opt \alpha_{{e}_2} = 0.01';...
    'opt \alpha_{{e}_3} = 0.1';'wt opt \alpha_{{e}_3} = 0.1';...
    'opt \alpha_{{e}_4} = 0.5';'wt opt \alpha_{{e}_4} = 0.5'},'Fontsize',9);

