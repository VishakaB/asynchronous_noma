function [convergedukfin,nbiterationslam,lam,nbiterationsuk,decision_uk] =...
    grad_descent_uklam(grad_lam,...
    lam,learn_rate,n_iter,tolerance,K_vec,decision_uk,K,sym_dur_vec, energyth)
    if K>5
        startlam = 1.15;
        mew = 1.15;
    else 
        startlam = 0.9;
        mew = 0.9;
    end
    
    convergeduk = false;    
    convergedlam = false;
    convergedmew = false;
    nbiterationslam = 1;
    nbiterationsuk  = 1;
    
    while(convergedlam==false & convergedmew==false)
        grad_uk = -K_vec-lam...
                            + 2*lam'*decision_uk + ...
                            mew*(sum(sym_dur_vec) - energyth) ;
                                        
        diffuk  = -learn_rate*(grad_uk);       
        if (abs(diffuk)<= tolerance)
            convergeduk = true;
            nbiterationsuk = nbiterationsuk+1;
            
            while(convergeduk == true)
                grad_lam = -sum(decision_uk)...
                + sum(decision_uk.^2);
                grad_mew = (decision_uk'*sym_dur_vec) - energyth ;
                diff = -learn_rate*grad_lam;
                
                diffmew = -learn_rate*grad_mew;
                
                if (abs(diff)<= tolerance)
                    converged_lam = lam;
                    nbiterationslam = nbiterationslam+1;
                    convergedlam = true;
                    convergeduk = false;   
                    convergedukfin = true;
                    decision_uk =decision_uk>0.8;
                else
                    converged_lam = 0;
                    nbiterationslam = nbiterationslam+1;
                    lam = lam + diff;                   
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