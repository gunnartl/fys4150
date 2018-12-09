import numpy as np
import time as t



class mc_SIR:
    def __init__(self, population_size, S0, I0,a,b,c,d=0,e=0,f=0,d_I=0,R0=0):
        self.a,self.b,self.c, self.ps = a,b,c,population_size
        self.f,self.e,self.d,self.d_I = f,e,d,d_I
        self.S0, self.I0,self.R0 = S0, I0,R0
                
        
        
        
    def propagate(self, steps):
        
        if type(self.a) == int or type(self.a) ==float: #if a is a constant sets it to an arrray, for generalized code
            self.a = np.ones(steps)*self.a
        
        if type(self.f) == int or type(self.f) == float: # same with f
            self.f = np.ones(steps)*self.f
        
        self.dt= np.min([np.min(4/(self.a*self.ps)),1/(self.b*self.ps),1/(c*self.ps)])

        time = np.zeros(steps)
        S = np.zeros(steps)
        S[0] = self.S0
        I = np.zeros(steps)
        I[0] = self.I0
        R = np.zeros(steps)
        R[0] = self.R0
        rands = np.random.random((steps,8)) # eliminates multiple calls to random.random
        
        for j in range(steps-1):
            #self.dt= np.min([np.min(4/(self.a*self.ps)),1/(self.b*self.ps),1/(c*self.ps)])
            i = 0
            r = 0
            s = 0
            
            #S to I
            if rands[j,0]< self.a[j]*S[j]*I[j]*self.dt/self.ps:
                i+=1
                s-=1
                
            # I to R    
            if rands[j,1]< self.b*I[j]*self.dt:
                r += 1
                i -= 1
            
            # R to S
            if rands[j,2]< self.c*R[j]*self.dt:
                s += 1
                r -= 1
            
            # birth in whole population
            if rands[j,3] < (self.e*self.ps)*self.dt:
                s += 1
            
            # death in S
            if rands[j,4] < (self.d*S[j])*self.dt:
                s -= 1
            
            #death in R    
            if rands[j,5] < (self.d*R[j])*self.dt:
                r -= 1
            
            #death in I
            if rands[j,6] < ((self.d+self.d_I)*I[j])*self.dt:
                i -= 1
           
            ##vaccination 
            if rands[j,7] < self.f[j]*self.ps*self.dt:
              r += 1
              s -= 1
            
            
            S[j+1]    = S[j] + s #-self.f[j]*self.dt
            R[j+1]    = R[j] + r #+self.f[j]*self.dt
            I[j+1]    = I[j] + i 
            time[j+1] = time[j] + self.dt
            
            
            self.ps = S[j+1] +I[j+1]+R[j+1] # update population size
            
        return S,I,R, time
    
    
if __name__ == "__main__":
    steps = 50000
    

    
    pop_size = 400
    S0 = 300
    I0 = 100
    R0 = 0
    
    a   = 3 + np.sin(np.linspace(0,np.pi,steps)*6.5)*2#+ np.sin(np.linspace(0,np.pi,steps)*4)*4#infecsiousness
    b   = [1]
    c   = 0.5
    d   = 0.03 #death rate
    e   = 0.04 # birth rate
    
    #x   = #np.linspace(0,steps,steps)
    f   = 0.15 + np.sin(np.linspace(0,np.pi,steps)*6.5)*0.15#np.piecewise(x,[x<5000,x>=5000],[0,0.4])
    
    d_I = .03 # death rate of infected
    #file = open("outdata.txt","w")
    #file.write("b S ssd i isd r rsd \n")
    for j in b:
        sav = np.zeros(steps)
        iav = np.zeros(steps)
        rav = np.zeros(steps)
        tav = np.zeros(steps)
        population = mc_SIR(pop_size,S0,I0,a,j,c,d,e,f,d_I,R0)
        start = t.time()
        mcs = 10
        for i in range(mcs):
            S,I,R,time = population.propagate(steps) 
            sav += S
            iav += I
            rav += R
            tav += time
        print(t.time()-start)
        
        #start = 50000
        
        #sexp = np.mean(sav[start:])/mcs
        #ssd  = np.sqrt(np.var(sav[start:]))/np.sqrt(mcs)
    
        #iexp = np.mean(iav[start:])/mcs
        #isd = np.sqrt(np.var(iav[start:]))/np.sqrt(mcs)
              
        #rexp = np.mean(rav[start:])/mcs
        #rsd = np.sqrt(np.var(rav[start:]))/np.sqrt(mcs)
        #file.write("%i %.2f %.2f %.2f %.2f %.2f %.2f \n" %(j,sexp,ssd,iexp,isd,rexp,rsd))    
    #file.close()
    #%%
    import matplotlib.pyplot as plt 
    plt.plot(tav/mcs,sav/mcs)
    plt.plot(tav/mcs,iav/mcs)
    plt.plot(tav/mcs,rav/mcs)
    plt.plot(tav/mcs,(rav+sav+iav)/mcs)
    plt.grid()
    plt.ylabel("# of people")
    plt.xlabel("Time")
    #plt.title("b = 3")
    plt.legend(["S","I","R","Population Size"])
    plt.show()
    """
    b = [1,2,3,4]
    t = np.linspace(0,2000*0.01,2000)
    for i in range(len(b)):
            test = mc_SIR(400,300,100,4,b[i],0.5)
            S, I, R, time= test.propagate(8000)
            plt.subplot(2,2,i+1)
            plt.grid()
            plt.title("b = %i"%b[i])
            plt.plot(time,S)
            plt.plot(time,I)
            plt.plot(time,400-I-S)
            
            if i == 2 or i ==3:
                plt.xlabel("Time")
            if i == 0 or i==2:
                    plt.ylabel("# of people")
            if i == 3:
                    plt.legend(["S","I","R"])
    """
