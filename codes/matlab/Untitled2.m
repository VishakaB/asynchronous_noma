Iz=eye(441); % Iz identity matrix refer to z part
BB=[B,Iz] ; % B is binary matrix with [441*263] dimension 
for i=1:704
      f_to_be_optimized(i)=1;
  end
for i=1:441      
      b_in_optimization_problem(i)=0;
  end
b_in_optimization_problem=b_in_optimization_problem';
t = bintprog(-f_to_be_optimized,BB,b_in_optimization_problem)