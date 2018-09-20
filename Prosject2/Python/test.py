#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np



def samurai_jacobi(A,tol):
    
    k,l = np.unravel_index(np.argmax(A, axis=None), A.shape)
    while(A[k,l] > tol):
        tau = (A[l,l]-A[k,k])/(2*[k,l])
        t   = min(-tau+np.sqrt(1-tau**2),-tau-np.sqrt(1-tau**2))
        c   = 1/np.sqrt(1+t**2)
        s   = t*s
        
        S   = 
    
    return A


n = 3

A = np.array([[1,4,3],[1,2,0],[0,0,1]])#np.random.rand(n,n)

k,l = np.unravel_index(np.argmax(A, axis=None), A.shape)
print(k,l)



print(samurai_jacobi(A,1e-8))
#w, v = np.linalg.eig(A)
#print(w,"\n",v)