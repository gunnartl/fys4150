import numpy as np

class SIRS:
    """
    SIRS-model with RK4 solver
    a = infection pressure
    b = recovery rate
    c = imunity loss rate
    d = death reate in population
    e = birth rate of population
    f = vacination rate of population
    d_i = increaseed death rate of infected individuals
    """
    def __init__(self, population_size, S0, I0,a,b,c,d=0,e=0,f=0,d_i=0,R0=0):
        self.a,self.b,self.c, self.ps = a,b,c,population_size
        self.d,self.e,self.f,self.d_i   = d,e,f,d_i
        self.S0, self.I0,self.R0 = S0, I0, R0    
        
    # S derivert    
    def S_d(self,S,I,R,a,f):
        return self.c*(R) - a*S*I/self.ps -f*self.ps + self.e*self.ps - self.d*S
    # I derivert
    def I_d(self,S,I,a):
        return a*S*I/self.ps -self.b*I - self.d*I - self.d_i*I
    #R derivert
    def R_d(self,S,I,R,f):
        return self.b*I-self.c*R +f*self.ps -self.d*R
    
    # RK4
    def solve(self,steps,dt):
        
        if type(self.a) == int or type(self.a) ==float: #if a is a constant sets it to an arrray, for generalized code
            self.a = np.ones(steps)*self.a
        
        if type(self.f) == int or type(self.f) == float: # same with f
            self.f = np.ones(steps)*self.f
        
        S = np.zeros(steps)

        S[0] = self.S0
        I = np.zeros(steps)
        I[0] = self.I0
        R = np.zeros(steps)
        R[0] = self.R0
        
        for i in range(1,steps):
            
            sk1 = self.S_d(S[i-1], I[i-1],R[i-1],self.a[i],self.f[i])
            ik1 = self.I_d(S[i-1], I[i-1],self.a[i])
            rk1 = self.R_d(S[i-1], I[i-1],R[i-1],self.f[i])
            
            sk2 = self.S_d(S[i-1] + (dt/2)*sk1, I[i-1] + (dt/2)*ik1,R[i-1]+ (dt/2)*rk1,self.a[i],self.f[i])
            ik2 = self.I_d(S[i-1] + (dt/2)*sk1, I[i-1] + (dt/2)*ik1,self.a[i])
            rk2 = self.R_d(S[i-1] + (dt/2)*sk1, I[i-1] + (dt/2)*ik1,R[i-1]+ (dt/2)*rk1,self.f[i])
            
            sk3 = self.S_d(S[i-1] + (dt/2)*sk2, I[i-1] + (dt/2)*ik2,R[i-1]+ (dt/2)*rk2,self.a[i],self.f[i])
            ik3 = self.I_d(S[i-1] + (dt/2)*sk2, I[i-1] + (dt/2)*ik2,self.a[i])
            rk3 = self.R_d(S[i-1] + (dt/2)*sk2, I[i-1] + (dt/2)*ik2,R[i-1]+ (dt/2)*rk2,self.f[i])
            
            sk4 = self.S_d(S[i-1] + dt*sk3,   I[i-1] + dt * ik3,R[i-1]+ dt * rk3,self.a[i],self.f[i])
            ik4 = self.I_d(S[i-1] + dt*sk3,   I[i-1] + dt * ik3,self.a[i])
            rk4 = self.R_d(S[i-1] + dt*sk3,   I[i-1] + dt * ik3,R[i-1]+ dt * rk3,self.f[i])
            
            S[i] = S[i-1] + (dt/6)*(sk1 + 2*sk2 + 2*sk3 + sk4)
            I[i] = I[i-1] + (dt/6)*(ik1 + 2*ik2 + 2*ik3 + ik4)
            R[i] = R[i-1] + (dt/6)*(rk1 + 2*rk2 + 2*rk3 + rk4)
            
            self.ps = S[i]+R[i]+I[i] # update size of population
            
        return S,I,R

if __name__ == "__main__":
    # steps 7000, og a = 8 + np.sin(np.linspace(0,np.pi,steps)*4)*8 it sjuk graf
    steps = 40000
    N = 400
    dt  = .001
    a   = 3 + np.sin(np.linspace(0,np.pi,steps)*6.5)*2
    b   = 1      #
    c   = .5     # 
    d   = 0.03
    d_i = .03
    e   = 0.04

    
    f   = 0.3 + np.sin(np.linspace(0,np.pi,steps)*6.5)*0.2
    
    time = np.linspace(0,steps*dt,steps)

    
    S0 = 300
    I0 = 100
    R0 = 0
    
    start = SIRS(N,S0,I0,a,b,c,d,e,f,d_i,R0)
    
    S,I,R = start.solve(steps,dt)
    
    #%%    
    import matplotlib.pyplot as plt
    
    plt.plot(time,S)
    plt.plot(time,I)
    plt.plot(time,R)
    plt.plot(time,S+I+R)

    plt.grid()
    plt.show()
