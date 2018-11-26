import numpy as np

def S_d(S,I):
    return c*(N-S-I) - a*S*I/N

def I_d(S,I):
    return a*S*I/N -b*i

def rk4(S, I, f):

    k1 = f(S, I)
    k2 = f(S + dt/2, I + (dt/2)*k1)
    k3 = f(S + dt/2, I + (dt/2)*k2)
    k4 = f(S + dt,   I + dt * k3)
    

    return (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

dt = 1
a = 4 
b = 1
c = .5

steps = 100

N = 400

S = np.zeros(400)
I = np.zeros_like(S)

S[0] = 300
I[0] = 100

for i in range(1,steps): 
    S[i] = S[i] + rk4(S[i],I[i],S_d)
    I[i] = I[i] + rk4(S[i],I[i],I_d)
    
    
import matplotlib.pyplot as plt

plt.plot(S)
plt.plot(I)
plt.show()