
K = 3;%number of users
E_max =10;
total_throughput = 2;
decision_uk =[1;0;0];
K_vec= [3;2;1];
total_energ_consump = E_max - E_max^(exp(-log(2)/1000*decision_uk'*K_vec));
energy_eff = total_throughput/(total_energ_consump +0.01)


K = 3%number of users
E_max =10;
total_throughput = 2;
decision_uk =[1;1;0];
K_vec= [3;2;1];
total_energ_consump = E_max - E_max^(exp(-log(2)/1000*decision_uk'*K_vec));
energy_eff = total_throughput/(total_energ_consump +0.01)


K = 3%number of users
E_max =10;
total_throughput = 2;
decision_uk =[1;1;1];
K_vec= [3;2;1];
total_energ_consump = E_max - E_max^(exp(-log(2)/1000*decision_uk'*K_vec));
energy_eff = total_throughput/(total_energ_consump +0.01)

K = 5 %number of users
E_max =10;
total_throughput = 2;
decision_uk =[1;1;1;1;1];
K_vec= [5;4;3;2;1];
total_energ_consump = E_max - E_max^(exp(-log(2)/1000*decision_uk'*K_vec));
energy_eff = total_throughput/(total_energ_consump +0.01)

K = 10;%number of users
E_max =10;
total_throughput = 2;
decision_uk =ones(K,1);
K_vec= [10;9;8;7;6;5;4;3;2;1];
total_energ_consump = E_max - E_max^(exp(-log(2)/1000*decision_uk'*K_vec));
energy_eff = total_throughput/(total_energ_consump +0.01)