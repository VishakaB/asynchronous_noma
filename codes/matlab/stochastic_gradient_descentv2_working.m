startlam = 1.05;
%n_iter   = 30;
tolerance = 0.1;
lam = startlam;
grad_lam = 2*lam;
K_vec = [3,2,1];
learn_rate = 0.1;

[x,y,lam,a,decision_uk] = grad_descent_lam(grad_lam,...
    startlam,learn_rate,n_iter,tolerance,K_vec)

function [convergedukfin,nbiterationslam,lam,nbiterationsuk,decision_uk] = grad_descent_lam(grad_lam,...
    lam,learn_rate,n_iter,tolerance,K_vec,decision_uk)

    convergeduk = false;    
    convergedlam = false;   
    K = 3;
    nbiterationslam = 1;
    nbiterationsuk  = 1;
    %decision_uk = [0.1,0.1,0.1];%initial
    %decision_uk = decision_uk>0.8
    while(convergedlam==false)
        grad_uk = -K_vec-lam...
                            + 2*lam'*decision_uk;
                        
        fprintf('0 grad_uk %f\n',grad_uk)                
        diffuk  = -learn_rate*(grad_uk);
        fprintf('1 diffuk %f\n',diffuk)     
        
        if (abs(diffuk)<= tolerance)
            convergeduk = true;
            nbiterationsuk = nbiterationsuk+1;
            
            while(convergeduk == true)
                grad_lam = -sum(decision_uk)...
                + sum(decision_uk.^2);
                diff = -learn_rate*grad_lam;

                if (abs(diff)<= tolerance)
                    converged_lam = lam;
                    nbiterationslam = nbiterationslam+1;
                    disp('yea')
                    convergedlam = true;
                    convergeduk = false;   
                    convergedukfin = true;
                    fprintf('decision_uk %f\n',decision_uk)
                    decision_uk =decision_uk>0.8
                else
                    converged_lam = 0;
                    nbiterationslam = nbiterationslam+1;
                    lam = lam + diff;
                    fprintf('lam: %f\n',lam);
                end               
            end            
        else
            decision_uk    = decision_uk + diffuk;
            %decision_uk = decision_uk>0.8
            if ((decision_uk) <zeros(K,1))
                 decision_uk = zeros(K,1);
            end
            nbiterationsuk = nbiterationsuk+1;
        end
       
    end   
     
end