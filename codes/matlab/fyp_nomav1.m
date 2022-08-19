clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem parameter initialization
% L: number of users.
% G: channel gain.
% F: nonnegative matrix. F_{lj} = G_{lj} if l ~= j, and F_{lj} = 0 if l = j
% v: nonnegative vector. v_l = 1/G_{ll}
% beta: a vector that is a weight assigned to links to reflect priority.
% pmax: upper bound of the total power constraints.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = 4; 
G = zeros(L,L); 
F = zeros(L,L); 
v = zeros(L,1); 
beta = rand(L,1); 
pmax = 2;

for l = 1:1:L
    for j = 1:1:L
        if l~= j
            G(l,j) = 0.1+0.1*rand(1,1);
        else
            G(l,j) = 0.6+0.3*rand(1,1);
        end
    end
end

for l = 1:1:L
    for j = 1:1:L
        if l ~= j
            F(l,j) = G(l,j);
        else
            F(l,j) = 0;
        end
    end
    v(l) = 1/G(l,l);
end

itenum = 15;
power_evolution = [];
p = rand(L,1);
power_evolution = [power_evolution;p'];
pnew = zeros(L,1);

for i = 1:1:itenum 
    for l = 1:1:L
        tmp = v(l);
        for j = 1:1:L
    % Algorithm step 1: update power for each user in a distributed way 
           tmp = tmp+F(l,j)*v(l)*p(j);
        end
        pnew(l) = beta(l)*p(l)/(p(l)/tmp);
    end
    % Algorithm step 2: normalize power centrally at the base station
    pnew = pnew*(pmax/sum(pnew));
    power_evolution = [power_evolution;pnew'];
    p = pnew;
end


% plot the power evolution 
figure
plot1 = plot(1:1:(itenum+1),power_evolution(:,1),'-o',1:1:(itenum+1),power_evolution(:,2),...
        '-^',1:1:(itenum+1),power_evolution(:,3),'-+',1:1:(itenum+1),power_evolution(:,4),'-*','linewidth',1.5);
legend('User 1','User 2','User 3','User 4');