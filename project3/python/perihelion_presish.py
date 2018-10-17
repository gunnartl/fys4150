import relatively_classy as rc

mass_mercury= 3.3e23
p0_mercury =  np.array([0.3075,0,0])
v0_mercury =  np.array([0,12.44,0])

mass_sun     = 2e30
p0_sun       = np.array([0,0,0])
v0_sun       = -(mass_mercury*v0_mercury)/mass_sun

mercury  = mpc.body(mass_mercury,p0_mercury,v0_mercury)
sun      = mpc.body(mass_sun,p0_sun,v0_sun)

system = rc.system(sun,mercury)

forsøk = rc.solve(system,2,1/365)

import matplotlib.pyplot as plt

for i in range(len(system)):
    plt.plot(system[i].pos_vec[:,0],system[i].pos_vec[:,1])
plt.show()
