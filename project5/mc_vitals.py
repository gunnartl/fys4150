import numpy as np
import time



class mc_SIR:
    def __init__(self, population_size, S0, I0,a,b,c,d=0,e=0,d_I=0,R0=0):
        self.a, self.b,self.c, self.ps = a,b,c,population_size
        self.e,self.d,self.d_I = e,d,d_I
        self.S0, self.I0,self.R0 = S0, I0,R0
        
        self.dt= min(4/(a*self.ps),1/(b*self.ps),1/(c*self.ps))
    

    def propagate(self, steps):
        S = np.zeros(steps)
        S[0] = self.S0
        I = np.zeros(steps)
        I[0] = self.I0
        R = np.zeros(steps)
        R[0] = self.R0
        
        
        for j in range(steps-1):
            i = 0
            r = 0
            s = 0
            #S to I
            if np.random.random()< self.a*S[j]*I[j]*self.dt/self.ps:
                i+=1
                s-=1
                
            # I to R    
            if np.random.random()< self.b*I[j]*self.dt:
                r += 1
                i -= 1
                    
            if np.random.random()< self.c*R[j]*self.dt:
                s += 1
                r -= 1
                        
            S[j+1] = S[j] + s - self.d*S[i-1]+ e*self.ps
            R[j+1] = R[j] + r - (self.d)*R[i-1]
            I[j+1] = I[j] + i - (self.d+self.d_I)*I[i-1]
        
        return S,I,R
    
    
if __name__ == "__main__":
    steps = 5600
    
    sav = np.zeros(steps)
    iav = np.zeros(steps)
    rav = np.zeros(steps)
    
    pop_size = 400
    S0 = 300
    I0 = 100
    R0 = 0
    
    a   = 8 #infecsiousness
    b   = 1
    c   = 0.5
    d   = 0.00001 #death rate
    e   = 0.000014 # birth rate
    d_I = 0.003 # death rate of infected
    
    test = mc_SIR(pop_size,S0,I0,a,b,c,d,e,d_I,R0)
    start = time.time()
    mcs = 100
    for i in range(mcs):
        S,I,R = test.propagate(steps) 
        sav += S
        iav += I
        rav += R
    print(time.time()-start)
            
            
    import matplotlib.pyplot as plt 
    plt.plot(sav/mcs)
    plt.plot(iav/mcs)
    plt.plot(rav/mcs)
    plt.plot(S+I+R)
