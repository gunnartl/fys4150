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
    
    def make_pos_vec(self,N):
        self.pos_vec = np.zeros((N,3))
        self.pos_vec[0] = self.p0

    def make_vel_vec(self,N):
        self.vel_vec = np.zeros((N,3))
        self.vel_vec[0] = self.v0

    def make_acc_vec(self,N,a0):
        self.acc_vec = np.zeros((N,3))
        self.acc_vec[0] = a0

        

def system(*args):
    system = []
    for i in args:
        system.append(i)
    
    return system

        
class solve:
    def __init__(self,system,N,dt):
        self.dt = dt
        self.N  = N
        self.system = system

        self.G = 4*np.pi**2
        
        for i in self.system:
            i.make_pos_vec(N)
            i.make_vel_vec(N)
        for i in self.system:
            i.make_acc_vec(N,self.acceleration(i,0))
            
            
    def acceleration(self,current,i):
        acc = -current.pos_vec[i] * (self.G/np.linalg.norm(current.pos_vec[i])**3)

        for j in self.system:
            if j != current:
                r = current.pos_vec[i]-j.pos_vec[i]
                acc -= r*self.G*j.mass/(np.linalg.norm(r)**3)
        return acc

    
    def verlet(self):
        dt     = self.dt
        dt2    = dt**2
        
        start  = time.time()
        for i in range(self.N-1):
            for j in self.system:
                j.pos_vec[i+1] = j.pos_vec[i] + dt*j.vel_vec[i] + 0.5*dt2*j.acc_vec[i]
            for j in self.system:
                j.acc_vec[i+1] = self.acceleration(j,i+1)
            for j in self.system:
                j.vel_vec[i+1] = j.vel_vec[i] + 0.5*(j.acc_vec[i]+j.acc_vec[i+1])*dt
        stop = time.time()
        print("Verlet : ",stop-start,"secs")
        


        
if __name__ == "__main__":
	a = 2
    
