#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from numba import jit


@jit
def samurai_jacobi(A,tol):
    count =0
    k,l = np.unravel_index(np.argmax(A-np.diag(np.diag(A)), axis=None), A.shape)
    while(abs(A[k,l]) > tol):

        tau = (A[l,l]-A[k,k])/(2*A[k,l])

        t   = min(-tau+np.sqrt(1+tau**2),-tau-np.sqrt(1+tau**2))
        c   = 1/np.sqrt(1+t**2)
        s   = t*c
        
        S   = np.identity(A.shape[0])
        S[k,k] = c
        S[l,l] = c
        S[k,l] = s
        S[l,k] = -s
        A = S.T.dot(A.dot(S))
        #print(A)
        k,l = np.unravel_index(np.argmax((A-np.diag(np.diag(A)))**2, axis=None), A.shape)
        count += 1
    print(count)
    return A,count


N = 4
h = 1/N
from scipy.sparse import diags
diagonals = [np.ones(N)*(2), np.ones(N-1)*(-1), np.ones(N-1)*(-1)]
A = diags(diagonals, [0, -1, 1]).toarray()

print(A)
#A = np.array([[1,2,3],[2,1,2],[3,2,1]])#np.random.rand(n,n)



#print(2/(h**2) + 2*(-1/h**2)*np.cos(np.pi/(N+1)))
#print(np.sort(np.diag(samurai_jacobi(A,1e-6))))
#w, v = np.linalg.eig(A)
#print(np.sort(w)[0]/h**2)
print(samurai_jacobi(A,1e-6)[0])


N     = np.array(range(10,110,10))
count = []
err   = []
for i in N:
    print(i)
    h2 = 1./((i+1)**2)
    diagonals = [np.ones(i)*(2), np.ones(i-1)*(-1), np.ones(i-1)*(-1)]
    A = diags(diagonals, [0, -1, 1]).toarray()
    M, tell = samurai_jacobi(A,1e-6)
    w, v = np.linalg.eig(A/h2)
    Eig1 = np.sort(w)[0]
    Eig2 = np.sort(np.diag(M))[0]/h2
    err.append((abs(np.sort(w)[0]-np.sort(np.diag(M))[0]/h2)/np.sort(w)[0]))
    count.append(tell)
print(2/h2-(2/h2)*np.cos(np.pi/(i+1)))
print(np.sort(np.diag(M))/h2-np.sort(w))

print(err)

    
import numpy as np
import matplotlib.pyplot as plt


fig, ax1 = plt.subplots()


color = 'tab:red'
ax1.set_ylabel('Iterations', color=color)
ax1.plot(N, count, color=color)
ax1.plot(N,1.2*N**2,"--", color=color)
ax1.tick_params(axis='y', labelcolor=color)
plt.legend(["Iterations","$1.2N^2$"])
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Relative error', color=color)  # we already handled the x-label with ax1
ax2.plot(N, err, color=color)
plt.yscale("log")
ax2.tick_params(axis='y', labelcolor=color)

plt.legend(["Relative error of $\lambda_0$"])
ax1.set_xlabel("N in NxN matrix")
ax1.set_title("Development in iterations and eror for increasing N")
fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.show()