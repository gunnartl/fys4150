import numpy as np



class SIR_simple:
    def __init__(self, population_size, S0, I0,a,b,c):
        self.a, self.b,self.c, self.ps = a,b,c,population_size
        self.S0, self.I0 = S0, I0
        
    def S_d(self,S,I):
        return self.c*(self.ps-S-I) - self.a*S*I/self.ps

    def I_d(self,S,I):
        return self.a*S*I/self.ps -self.b*I
    
    
    def solve(self,steps,dt):
        S = np.zeros(steps)
        S[0] = self.S0
        I = np.zeros(steps)
        I[0] = self.I0
        
        for i in range(1,steps):
            
            sk1 = self.S_d(S[i-1], I[i-1])
            ik1 = self.I_d(S[i-1], I[i-1])
            
            sk2 = self.S_d(S[i-1] + (dt/2)*sk1, I[i-1] + (dt/2)*ik1)
            ik2 = self.I_d(S[i-1] + (dt/2)*sk1, I[i-1] + (dt/2)*ik1)
            
            sk3 = self.S_d(S[i-1] + (dt/2)*sk2, I[i-1] + (dt/2)*ik2)
            ik3 = self.I_d(S[i-1] + (dt/2)*sk2, I[i-1] + (dt/2)*ik2)
            
            sk4 = self.S_d(S[i-1] + dt*sk3,   I[i-1] + dt * ik3)
            ik4 = self.I_d(S[i-1] + dt*sk3,   I[i-1] + dt * ik3)    
            
            S[i] = S[i-1] + (dt/6)*(sk1 + 2*sk2 + 2*sk3 + sk4)
            I[i] = I[i-1] + (dt/6)*(ik1 + 2*ik2 + 2*ik3 + ik4)
            
        return S,I


import matplotlib.pyplot as plt
b = [1,2,3,4]
t = np.linspace(0,2000*0.01,2000)
for i in range(len(b)):
    test = SIR_simple(400,300,100,4,b[i],0.5)
    S, I = test.solve(2000,0.01)
    plt.subplot(2,2,i+1)
    plt.grid()
    plt.title("b = %i"%b[i])
    plt.plot(t,S)
    plt.plot(t,I)
    plt.plot(t,400-I-S)
   
    if i == 2 or i ==3:
         plt.xlabel("Time")
    if i == 0 or i==2:
        plt.ylabel("# of people")
    if i == 3:
        plt.legend(["S","I","R"])
        
plt.show()