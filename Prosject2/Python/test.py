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

"""
N     = np.array(range(10,110,10))
count = []
err   = []
for i in N:
    print(i)
    h = 1/(i)
    diagonals = [np.ones(i)*(2), np.ones(i-1)*(-1), np.ones(i-1)*(-1)]
    A = diags(diagonals, [0, -1, 1]).toarray()
    M, tell = samurai_jacobi(A,1e-6)
    w, v = np.linalg.eig(A*h**2)
    err.append(max(abs(np.sort(w))-abs(np.sort(np.diag(M))*h**2)))
    count.append(tell)
print(np.sort(w)-np.sort(np.diag(M)*h**2))



    
import numpy as np
import matplotlib.pyplot as plt


fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_ylabel('iterations', color=color)
ax1.plot(N, count, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('error', color=color)  # we already handled the x-label with ax1
ax2.plot(N, err, color=color)
plt.yscale("log")
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()"""