import numpy as np

def S_d(S,I,R):
    return c*(R) - a*S*I/N - d*S + e*N 

def I_d(S,I):
    return a*S*I/N -b*I -d*I - d_i*I

def R_d(S,I,R):
    return b*I-c*R-d*R

dt  = .01
a   = 4 
b   = 1
c   = 1.5
d   = 0.0
d_i = 0.
e   = 0.0

steps = 700

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
    
    sk2 = S_d(S[i-1] + (dt/2)*sk1, I[i-1] + (dt/2)*ik1,R[i-1]+ (dt/2)*rk1)
    ik2 = I_d(S[i-1] + (dt/2)*sk1, I[i-1] + (dt/2)*ik1)
    rk2 = R_d(S[i-1] + (dt/2)*sk1, I[i-1] + (dt/2)*ik1,R[i-1]+ (dt/2)*rk1)
    
    sk3 = S_d(S[i-1] + (dt/2)*sk2, I[i-1] + (dt/2)*ik2,R[i-1]+ (dt/2)*rk2)
    ik3 = I_d(S[i-1] + (dt/2)*sk2, I[i-1] + (dt/2)*ik2)
    rk3 = R_d(S[i-1] + (dt/2)*sk2, I[i-1] + (dt/2)*ik2,R[i-1]+ (dt/2)*rk2)
    
    
    sk4 = S_d(S[i-1] + (dt)* sk3,   I[i-1] + dt * ik3,R[i-1]+ dt * rk3)
    ik4 = I_d(S[i-1] + (dt)* sk3,   I[i-1] + dt * ik3)
    rk4 = R_d(S[i-1] + (dt)* sk3,   I[i-1] + dt * ik3,R[i-1]+ dt * rk3)
    
    
    S[i] = S[i-1] + (dt/6)*(sk1 + 2*sk2 + 2*sk3 + sk4)
    I[i] = I[i-1] + (dt/6)*(ik1 + 2*ik2 + 2*ik3 + ik4)
    R[i] = R[i-1] + (dt/6)*(rk1 + 2*rk2 + 2*rk3 + rk4)
    


    
import matplotlib.pyplot as plt

plt.plot(S)
plt.plot(I)
plt.plot(R)
plt.plot(S+I+R)
plt.show()