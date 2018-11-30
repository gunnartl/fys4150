import numpy as np

def S_d(S,I,R):
    return c*(R) - a*S*I/N - d*S + e*N

def I_d(S,I):
    return a*S*I/N -b*I -d*I - d_i*I

def R_d(S,I,R):
    return b*I-c*R-d*R

def rk4(S, I, f):

    k1 = f(S, I)
    k2 = f(S + dt/2, I + (dt/2)*k1)
    k3 = f(S + dt/2, I + (dt/2)*k2)
    k4 = f(S + dt,   I + dt * k3)
    

    return (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

dt  = .01
a   = 4 
b   = 1
c   = .5
d   = 0.01
d_i = 0.1
e   = 0.03

steps = 1400

N = 400

S = np.zeros(steps)
I = np.zeros_like(S)
R = np.zeros_like(S)

S[0] = 300
I[0] = 100
R[0] = 0

for i in range(1,steps):
    
    sk1 = S_d(S[i-1], I[i-1],R[i-1])
    ik1 = I_d(S[i-1], I[i-1])
    rk1 = R_d(S[i-1], I[i-1],R[i-1])
    
    sk2 = S_d(S[i-1] + dt/2, I[i-1] + (dt/2)*sk1,R[i-1])
    ik2 = I_d(S[i-1] + dt/2, I[i-1] + (dt/2)*ik1)
    rk2 = R_d(S[i-1] + dt/2, I[i-1] + (dt/2)*sk1,R[i-1])
    
    sk3 = S_d(S[i-1] + dt/2, I[i-1] + (dt/2)*sk2,R[i-1])
    ik3 = I_d(S[i-1] + dt/2, I[i-1] + (dt/2)*ik2)
    rk3 = R_d(S[i-1] + dt/2, I[i-1] + (dt/2)*sk2,R[i-1])
    
    
    sk4 = S_d(S[i-1] + dt,   I[i-1] + dt * sk3,R[i-1])
    ik4 = I_d(S[i-1] + dt,   I[i-1] + dt * ik3)
    rk4 = S_d(S[i-1] + dt,   I[i-1] + dt * sk3,R[i-1])
    
    
    S[i] = S[i-1] + (dt/6)*(sk1 + 2*sk2 + 2*sk3 + sk4)
    I[i] = I[i-1] + (dt/6)*(ik1 + 2*ik2 + 2*ik3 + ik4)
    R[i] = R[i-1] + (dt/6)*(rk1 + 2*rk2 + 2*rk3 + rk4)
    


    
import matplotlib.pyplot as plt

plt.plot(S)
plt.plot(I)
plt.plot(R)
plt.plot(S+I+R)
plt.show()