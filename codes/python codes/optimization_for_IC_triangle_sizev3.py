#optimization for determining the IC triangle size
#ref: https://www.cvxpy.org/examples/basic/linear_program.html
#Import packages.
import cvxpy as cp
import numpy as np
import math 

L = 1;
N= 5;
e_sym=1;
lsta = [1,3];
for i in range(3,N+1):
  #lsta = [1, 3, 6, 24, 120, 360]#update list
  lsta.append(math.factorial(i));
print(lsta)
# Generate a random non-trivial linear program.
np.random.seed(1)

# Define and solve the CVXPY problem.
n_sym = cp.Variable(1)
e_sic = cp.Variable(1)

#for conventional power noma
prob = cp.Problem(cp.Maximize(n_sym),
                 [e_sic <= 10000,
                  2 <= n_sym, n_sym <= math.factorial(N),
                  e_sic == N*(N+L-2)*math.factorial(N)*e_sym]) #n_sym should be in in lsta# can put it outside in a if condition
prob.solve()
#print("A solution nsym is",n_sym.value)
if (n_sym.value != None and math.ceil(n_sym.value) in lsta):
  # Print result.
  print("\nThe optimal value is", prob.value)
  print("A solution nsym is",math.ceil(n_sym.value))
