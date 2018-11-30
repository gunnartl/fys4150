import numpy as np
import time



class mc_SIR:
    def __init__(self, population_size, S0, I0,a,b,c,R0=0):
        self.a, self.b,self.c, self.ps = a,b,c,population_size
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
                        
            S[j+1] = S[j] + s
            R[j+1] = R[j] + r
            I[j+1] = I[j] + i
        
        return S,I,R

        
if __name__ == "__main__":
    steps = 5600
    
    sav = np.zeros(steps)
    iav = np.zeros(steps)
    rav = np.zeros(steps)
    
    test = mc_SIR(400,300,100,4,1,0.5)
    start = time.time()
    for i in range(1000):
        S,I,R = test.propagate(steps) 
        sav += S
        iav += I
        rav += R
        print(time.time()-start)
            
            
    import matplotlib.pyplot as plt 
    plt.plot(sav/1000)
    plt.plot(iav/1000)
    plt.plot(rav/1000)
    #plt.plot(S+I+R)
