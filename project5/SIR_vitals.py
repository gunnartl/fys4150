import numpy as np

class SIR_vitals:
    def __init__(self, population_size, S0, I0,a,b,c,d=0,e=0,d_i=0,R0=0):
        self.a, self.b,self.c, self.ps = a,b,c,population_size
        self.d,self.e,self.d_i   = d,e,d_i
        self.S0, self.I0,self.R0 = S0, I0,R0
        
    def S_d(self,S,I,R):
        return self.c*(R) - self.a*S*I/self.ps + self.e

    def I_d(self,S,I):
        return self.a*S*I/self.ps -self.b*I
    
    def R_d(self,S,I,R):
        return self.b*I-self.c*R-self.d*R
    
    
    def solve(self,steps,dt):
        S = np.zeros(steps)
        S[0] = self.S0
        I = np.zeros(steps)
        I[0] = self.I0
        R = np.zeros(steps)
        R[0] = self.R0
        
        for i in range(1,steps):
            
            sk1 = self.S_d(S[i-1], I[i-1],R[i-1])
            ik1 = self.I_d(S[i-1], I[i-1])
            rk1 = self.R_d(S[i-1], I[i-1],R[i-1])
            
            sk2 = self.S_d(S[i-1] + (dt/2)*sk1, I[i-1] + (dt/2)*ik1,R[i-1]+ (dt/2)*rk1)
            ik2 = self.I_d(S[i-1] + (dt/2)*sk1, I[i-1] + (dt/2)*ik1)
            rk2 = self.R_d(S[i-1] + (dt/2)*sk1, I[i-1] + (dt/2)*ik1,R[i-1]+ (dt/2)*rk1)
            
            sk3 = self.S_d(S[i-1] + (dt/2)*sk2, I[i-1] + (dt/2)*ik2,R[i-1]+ (dt/2)*rk2)
            ik3 = self.I_d(S[i-1] + (dt/2)*sk2, I[i-1] + (dt/2)*ik2)
            rk3 = self.R_d(S[i-1] + (dt/2)*sk2, I[i-1] + (dt/2)*ik2,R[i-1]+ (dt/2)*rk2)
            
            sk4 = self.S_d(S[i-1] + dt*sk3,   I[i-1] + dt * ik3,R[i-1]+ dt * rk3)
            ik4 = self.I_d(S[i-1] + dt*sk3,   I[i-1] + dt * ik3)
            rk4 = self.R_d(S[i-1] + dt*sk3,   I[i-1] + dt * ik3,R[i-1]+ dt * rk3)
            
            S[i] = S[i-1] + (dt/6)*(sk1 + 2*sk2 + 2*sk3 + sk4)
            I[i] = I[i-1] + (dt/6)*(ik1 + 2*ik2 + 2*ik3 + ik4)
            R[i] = R[i-1] + (dt/6)*(rk1 + 2*rk2 + 2*rk3 + rk4)
            
        return S,I,R

if __name__ == "__main__":
    dt  = .01
    a   = 4 
    b   = 1
    c   = 1.5
    d   = 0.0
    d_i = 0.
    e   = 0.0
    
    steps = 700
    N = 400
    
    S0 = 300
    I0 = 100
    R0 = 0
    
    start = SIR_vitals(N,S0,I0,a,b,c,d,e,d_i,R0)
    
    S,I,R = start.solve(steps,dt)
    
        
    import matplotlib.pyplot as plt
    
    plt.plot(S)
    plt.plot(I)
    plt.plot(R)
    plt.plot(S+I+R)
    plt.show()