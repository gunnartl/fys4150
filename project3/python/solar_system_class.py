import numpy as np
import time
from numba import jit

class body:
    def __init__(self,mass, pos0, vel0):
        self.mass = mass/2e30
        self.p0    = pos0
        self.v0    = vel0
        self.G = 4* np.pi**2
        
    def kinetic_energy(self,v):
        K = 0.5*self.mass*np.linalg.norm(v,axis=1)**2
        return K
    
    def potential_energy(self,r,M=1): #M = 1 for the sun
        U = -self.G*M*self.mass/np.linalg.norm(r,axis=1)
        return U
    
    def angular_momentum(self,v,r):
        angmom = np.cross(r,v*self.mass,axis=1)
        return angmom
    

def system(*args):
    system = []
    for i in args:
        system.append(i)
    
    mass_centre = np.zeros(3)
    for i in system:
        mass_centre += i.mass*i.p0
    
    for i in system:
        i.p0 = i.p0-mass_centre
    
    return system

class solve:
    def __init__(self,planet,N,dt):
        self.dt = dt
        self.N  = N
        
        self.pos = np.zeros((N,3))
        self.pos[0] = planet.p0
        
        self.vel = np.zeros((N,3))
        self.vel[0] = planet.v0

        self.G = 4*np.pi**2
    #@jit    
    def acceleration(self,pos):
        a = self.G*(-pos/(np.linalg.norm(pos)**3))
        return a
    #@jit
    def euler(self):
        pos = self.pos
        vel = self.vel
        dt  = self.dt
        start = time.time()
        for i in range(int(self.N)-1):
            pos[i+1] = pos[i]+vel[i]*dt # euelr
            vel[i+1] = vel[i]+ self.acceleration(pos[i])*dt
        stop = time.time()
        print("Euler : ",stop-start,"secs")
        
        return pos,vel
    #@jit
    def verlet(self):
        pos    = self.pos
        vel    = self.vel
        dt     = self.dt
        acc    = np.zeros_like(pos)
        acc[0] = self.acceleration(pos[0])
        start  = time.time()
        dt2    = dt**2
        for i in range(self.N-1):
            
            pos[i+1] = pos[i] + dt*vel[i] + 0.5*dt2*acc[i]
            acc[i+1] = self.acceleration(pos[i+1])
            vel[i+1] = vel[i] + 0.5*(acc[i]+acc[i+1])*dt
        stop = time.time()
        print("Verlet : ",stop-start,"secs")
        return pos,vel,acc


        
if __name__ == "__main__":
    
    VX=-5.428888690270241E-03*365
    VY= 1.636353485946535E-02*365
    VZ=-4.491683144318728E-07*365
    
    EarthM = 6e24
    G = np.pi**2 * 4

        
    Earth = body(EarthM, np.array([1.,0,0]),np.array([0,np.sqrt(G),0]))
    
    print(np.sqrt((2/Earth.mass)*abs(Earth.potential_energy(np.array([Earth.p0])))))
    
    print(np.sqrt(G)+2.602)
    
    året    = solve(Earth,2000000,1e-5)
    E_p,E_v = året.euler()
    året    = solve(Earth,int(365),1/365)
    v_p,v_v,acc = året.verlet()
    vpot     = Earth.potential_energy(v_p)
    vkin     = Earth.kinetic_energy(v_v)
    vangmom  = Earth.angular_momentum(v_v,v_p)
    #Epot     = Earth.potential_energy(E_p)
    #Ekin     = Earth.kinetic_energy(E_v)
    #Eangmom  = Earth.angular_momentum(E_v,E_p)
    
    print(np.mean(vpot+vkin))

    import matplotlib.pyplot as plt
    plt.plot(E_p[:,0],E_p[:,1])
    plt.plot(v_p[:,0],v_p[:,1])
    plt.scatter(0,0,color="y")
    #plt.plot(vpot-np.mean(vpot))
    #plt.plot(vkin-np.mean(vkin))
    #plt.plot(vpot+vkin - np.mean(vpot)-np.mean(vkin))
    #plt.plot(angmom[:,2])
    #cur_axes = plt.gca()
    #cur_axes.axes.get_yaxis().set_ticks([])
    #plt.legend(["Potential","Kinetic","Sum"],loc="lower right")
    #plt.xlabel("Time [Days]")
    #plt.ylabel("Energy")
    #plt.annotate('Verlet', xy=(0, 0), xytext=(0, 0.000000006))
    #plt.title("Variation from the mean in kinetic and potential\n energy of the earth over one year")
    #plt.subplot(2,1,2)
    #plt.plot(Epot-np.mean(Epot))
    #plt.plot(Ekin-np.mean(Ekin))
    #plt.plot(Epot+Ekin - np.mean(Epot)-np.mean(Ekin))
    #plt.plot(angmom[:,2])
    #cur_axes = plt.gca()
    #cur_axes.axes.get_yaxis().set_ticks([])
    #plt.legend(["Potential","Kinetic","Sum"],loc="lower right")
    #plt.xlabel("Time [Days]")
    #plt.ylabel("Energy")
    #plt.annotate('Euler', xy=(0, 0), xytext=(0, 0.000008))
    plt.show()