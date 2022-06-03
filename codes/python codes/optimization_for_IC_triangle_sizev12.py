#optimization for determining the IC triangle size
#ref: https://www.cvxpy.org/examples/basic/linear_program.html
#ref: https://stackoverflow.com/questions/58599059/why-am-i-getting-this-dcperror
#ref:https://math.stackexchange.com/questions/2816421/linear-exponential-and-interval-decay-of-distance-function-above-a-threshold

#Import packages.
import cvxpy as cp
import numpy as np
import math 

L = 1;
N = 20;

e_th = math.factorial(N-1);
lsta = [1,3];
e_th =25;
for i in range(3,N+1):
  lsta.append(math.factorial(i));
print(lsta)

# debug
n = N;
x  = 1*np.ones((n,n));
Pk = np.random.rand(n,n);

# Sort in ascending order
Pk = np.sort(Pk)

# Reverse the sorted array
Pk = Pk[::-1]

deltak = np.random.rand(n,n);

# Sort in ascending order
deltak = np.sort(deltak)

# Reverse the sorted array
deltak = deltak[::-1]

# Generate a random non-trivial linear program.
np.random.seed(1)

# Define and solve the CVXPY problem.
n_sym = cp.Variable(1);
e_residual = cp.Variable(1);
e_sic = cp.Variable((n,1),"e_sic");
p_signal = cp.Variable((1),"p_signal");
p_interference  = cp.Variable((1),"p_interference");
N_0  = cp.Variable(1)
xk = cp.Variable((n,n),"xk",boolean=True)
y = cp.Variable((n,n), "y", boolean=True)

e = np.ones((1,n))
array = 1*np.ones((n,n));
array[np.triu_indices(n)] = 0
np.fill_diagonal(array, 1)
n_iter = 1
k0 =-0.01;
prob = cp.Problem(cp.Maximize(n_sym),
                 [0.00319*n_sym*n_iter <= cp.log (e_th),
                  N_0 <= e_th,
                  (N_0 + L -2) <= e_th,
                  N_0 - 2 <= e_th,
                  2 <= n_sym, n_sym <= cp.sum(cp.sum(array, axis = 0)),
                  p_signal - cp.multiply(Pk[:,0],xk[:,0].T)>=0,
                  p_interference - cp.multiply(Pk[:,1],deltak[:,1])@xk[:,1].T>=0,
                  cp.multiply(Pk[:,[2,n-1]],deltak[:,[2,n-1]])@xk[:,[2,n-1]].T <= 0.1,
                  N_0 == cp.sum(xk[-1,:]),
                  2<=N_0,
                  N<=n_sym,
                  n_sym==cp.sum(cp.sum(xk, axis = 0)), 
                  #make all the values in indices above the diagonal zero, diagonal should be non zero   
                  #A @ x <= b
                  xk<=array
                 ])
prob.solve()

print("A solution xk 1",xk.value)
print("A solution nsym 1",n_sym.value)
if (n_sym.value != None and math.ceil(n_sym.value) in lsta):
  # Print result.
  print("\nThe optimal value is", prob.value)
  print("A solution nsym is",math.ceil(n_sym.value))

n_sym = cp.Variable(1);
a = cp.Variable(1);
b = cp.Variable(1);
e_sic = cp.Variable((n,1),"e_sic");
p_signal = cp.Variable((1),"p_signal");
p_interference  = cp.Variable((1),"p_interference");
p_residual   = cp.Variable((1),"p_interference");
N_0  = cp.Variable(1)
xk = cp.Variable((n,n),"xk",boolean=True)
y = cp.Variable((n,n), "y", boolean=True)
e = np.ones((1,n))
array = 1*np.ones((n,n));
array[np.triu_indices(n)] = 0
np.fill_diagonal(array, 1)
prob = cp.Problem(cp.Maximize(n_sym),
                 [0.00319*n_sym*n_iter <= cp.log (e_th),
                  N_0 <= e_th,
                  (N_0 + L -2) <= e_th,
                  N_0 - 2 <= e_th,
                  2 <= n_sym, n_sym <= cp.sum(cp.sum(array, axis = 0)),
                  p_signal - cp.multiply(Pk[:,0],xk[:,0].T)>=0.1,
                  p_interference - cp.multiply(Pk[:,1],deltak[:,1])@xk[:,1].T>=0.1,
                  p_residual - cp.multiply(Pk[:,[2,n-1]],deltak[:,[2,n-1]])@xk[:,[2,n-1]].T >= 0.1,
                  #- 0.1*(p_signal+p_interference)<=0,#if condition
                  N_0 == cp.sum(xk[-1,:]),
                  2<=N_0,
                  N<=n_sym,
                  n_sym==cp.sum(cp.sum(xk, axis = 0)), 
                  #make all the values in indices above the diagonal zero, diagonal should be non zero   
                  #A @ x <= b
                  xk<=array
                 ])
prob.solve()

print("A solution xk 2",xk.value)
print("A solution nsym 2",n_sym.value)
if (n_sym.value != None and math.ceil(n_sym.value) in lsta):
  # Print result.
  print("\nThe optimal value is", prob.value)
  print("A solution nsym is",math.ceil(n_sym.value))
