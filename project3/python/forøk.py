#2458402.500000000 = A.D. 2018-Oct-11 00:00:00.0000 TDB 
# X = 9.528047055398201E-01 Y = 3.053612869840809E-01 Z =-9.272902073041313E-05
import numpy as np


def acc(p):
    acc = G*(-p)/(np.linalg.norm(p)**3)
    return acc


N     = 1000 #per year
years = 1.5

VX=-5.428888690270241E-03 
VY= 1.636353485946535E-02 
VZ=-4.491683144318728E-07

E_p = np.zeros((int(N*years),2)) 
E_p[0] = [1,0]

SunM = 1
EarthM = 6e24/2e30
G = np.pi**2 * 4

E_v = np.zeros((int(N*years),2))
E_v[0] = [0,np.sqrt(G)]


#VERLET things
v_p = np.zeros((int(N*years),2))
v_p[0] = [1,0]
v_v = np.zeros((int(N*years),2))
v_v[0] = [0,np.sqrt(G)]

v_a=np.zeros_like(v_p)
v_a[0] = acc(v_p[0])

dt = 1/(N-1)

for i in range(int(N*years)-1):
    E_p[i+1] = E_p[i]+E_v[i]*dt # euelr
    E_v[i+1] = E_v[i]+ G*(-E_p[i])*dt/(np.linalg.norm(E_p[i])**3)
    
    
    v_p[i+1] = v_p[i] + dt*v_v[i] + 0.5*dt**2*v_a[i]
    v_a[i+1] = acc(v_p[i+1])
    v_v[i+1] = v_v[i] + 0.5*(v_a[i]+v_a[i+1])*dt

    

E_v = np.linalg.norm(E_v,axis=(1))
import matplotlib.pyplot as plt
plt.plot(E_p[:,0],E_p[:,1])
plt.plot(v_p[:,0],v_p[:,1])
plt.scatter(0,0,color="y")
plt.show()
#plt.plot(0.5*E_v**2+G/np.linalg.norm(E_p,axis=(1)))