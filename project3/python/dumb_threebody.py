import numpy as np
import dust_multiple_body_class as dmb

mass_earth  = 6e24
p0_earth    = np.array([9.528047055398201E-01,3.053612869840809E-01,-9.272902073041313E-05])
v0_earth    = np.array([-5.428888690270241E-03,1.636353485946535E-02,-4.491683144318728E-07])*365
 
#jupiter
mass_jupiter = 1.9e27*1000       
p0_jupiter   = np.array([-2.679859418467178E+00,-4.648870533678862E+00,7.922561878424454E-02])
v0_jupiter   = np.array([6.448353939132267E-03,-3.409518873130117E-03,-1.300875035507015E-04])*365

earth   = dmb.body(mass_earth,p0_earth,v0_earth)
jupiter = dmb.body(mass_jupiter,p0_jupiter,v0_jupiter)

system  = dmb.system(earth,jupiter)


N = int(20e5)
dt  = 1e-5

året = dmb.solve(system,N,dt)
året.verlet()


import matplotlib.pyplot as plt
for i in range(len(system)):
    plt.plot(system[i].pos_vec[:,0],system[i].pos_vec[:,1])

plt.legend(["Earth","Jupiter"])
plt.xlabel("X [AU]")
plt.ylabel("Y [AU]")
plt.show()
