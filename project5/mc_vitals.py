from mc_sir import *

from test import *

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
plt.plot(np.arange(5600)/3.9,sav/1000)
plt.plot(np.arange(5600)/3.9,iav/1000)
plt.plot(np.arange(5600)/3.9,rav/1000)