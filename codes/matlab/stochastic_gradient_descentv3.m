startlam = 0.5;
n_iter   = 11;
learn_rate = 0.2;
tolerance  = 0.1;
lam = startlam;

K_vec = [3,2,1];

[x,y,z] = grad_descent_lam(grad_lam,...
    startlam,learn_rate,n_iter,tolerance,K_vec)

function [convergeduk,nbiterationslam,nbiterationsuk] = grad_descent_lam(grad_lam,...
    lam,learn_rate,n_iter,tolerance,K_vec)
    
    convergeduk = false;  
    convergedlam =false
    K = 3;
    nbiterationslam = 1;
    nbiterationsuk  = 1;
    decision_uk = [0.2,0.7,0.9];%initial
    nbiterationsuk = 0;
    while( convergeduk==false)
        grad_lam = -sum(decision_uk)...
        + sum(decision_uk.^2);
        difflam = -learn_rate*grad_lam;
        
        fprintf('grad_lam tolerance %f %f %f\n',difflam, grad_lam,tolerance);
        fprintf('grad_lam tolerance %i\n',abs(difflam)<= tolerance);
        if (abs(difflam)<= tolerance)
            if convergedlam ==true
                break;
            for j = 1:n_iter
                                
                grad_uk = -K_vec-lam...
                            +2*lam'*decision_uk;
                    
                diffuk  = -learn_rate*(grad_uk);
                if (abs(diffuk)<= tolerance)
                    convergeduk = true;
                    nbiterationsuk = nbiterationsuk+1;
                    disp('yea')
                    
                else
                    decision_uk    = decision_uk + diffuk;
                    if ((decision_uk) <zeros(K,1))
                         decision_uk = zeros(K,1);
                    end
                    nbiterationsuk = nbiterationsuk+1;
                end
            end   

        else
            converged_lam = 0;
            nbiterationslam = nbiterationslam+1;
            lam = lam + difflam;
            fprintf('lam: %f\n',lam);
        end
    end  
end