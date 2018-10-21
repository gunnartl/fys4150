#import relatively_classy as rc
import multiple_body_class as rc
import numpy as np

mass_mercury= 3.3e23
p0_mercury =  np.array([0.3075,0,0])
v0_mercury =  np.array([0,12.44,0])

mass_sun     = 2e30
p0_sun       = np.array([0,0,0])
v0_sun       = -(mass_mercury*v0_mercury)/mass_sun

mercury  = rc.body(mass_mercury,p0_mercury,v0_mercury)
sun      = rc.body(mass_sun,p0_sun,v0_sun)

system = rc.system(sun,mercury)
N = int(1.5e6)
forsøk = rc.solve(system,N,1e-6)
forsøk.verlet()

import matplotlib.pyplot as plt

for i in range(len(system)):
    plt.plot(system[i].pos_vec[:,0],system[i].pos_vec[:,1])
plt.show()


periheels = []

print(system[1].pos_vec[i,:].shape)

for i in range(1,N-1):

    r_now  = np.linalg.norm(system[1].pos_vec[i,:]-system[0].pos_vec[i,:])
    r_next = np.linalg.norm(system[1].pos_vec[i+1,:]-system[0].pos_vec[i+1,:])
    r_prev = np.linalg.norm(system[1].pos_vec[i-1,:]-system[0].pos_vec[i-1,:])
    if r_now < r_next and r_now < r_prev:
        periheels.append(system[1].pos_vec[i,:])
#%% 
periheels = np.array(periheels)
x0 = periheels[1][0]
y0 = periheels[1][1]

xf = periheels[-2][0]
yf = periheels[-2][1]


print(periheels)
plt.plot(tans)
plt.show()

print(len(periheels))
print ((np.arctan(y0/x0)-np.arctan(yf/xf)) *180/np.pi)