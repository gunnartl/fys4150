import numpy as np


a = 4 
b = 2
c = .5

steps = 3000

N = 400

S = np.zeros(steps)
I = np.zeros(steps)
R = np.zeros(steps)

S[0] = 300
I[0] = 100

dt = min(4/(a*N),1/(b*N),1/(c*N))
#%%
for j in range(steps-1):
    i = 0
    r = 0
    s = 0
    #S to I
    if np.random.random()< a*S[j]*I[j]*dt/N:
        i+=1
        s-=1
    
    # I to R    
    if np.random.random()< b*I[j]*dt:
        r += 1
        i -= 1

    if np.random.random()< c*R[j]*dt:
        s += 1
        r -= 1
    
    S[j+1] = S[j] + s
    R[j+1] = R[j] + r
    I[j+1] = I[j] + i
    
        
        
#%%        
import matplotlib.pyplot as plt 
plt.plot(S)
plt.plot(I)
plt.plot(R)
plt.plot(S+I+R)
### Legger inn et bedre programm er når jeg er ferdig
#def mc_sir(S0,I0,R0,n,)