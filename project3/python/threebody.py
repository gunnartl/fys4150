import numpy as np
import multiple_body_class as mpc

mass_mercur= 3.3e23
p0_mercur =  np.array([-2.138801978256535E-01,-4.028538544608614E-01,-1.397398026440086E-02])
v0_mercur =  np.array([1.927259979627735E-02,-1.164174161133437E-02,-2.720041450765680E-03])*365

mass_earth  = 6e24
p0_earth    = np.array([9.528047055398201E-01,3.053612869840809E-01,-9.272902073041313E-05])
v0_earth    = np.array([-5.428888690270241E-03,1.636353485946535E-02,-4.491683144318728E-07])*365
 
#jupiter
mass_jupiter = 1.9e27
p0_jupiter   = np.array([-2.679859418467178E+00,-4.648870533678862E+00,7.922561878424454E-02])
v0_jupiter   = np.array([6.448353939132267E-03,-3.409518873130117E-03,-1.300875035507015E-04])*365

mass_saturn  = 5.5e26
p0_saturn    = np.array([1.533542535074663E+00,-9.937711138425550E+00,1.117470354327505E-01])
v0_saturn    = np.array([5.206732001657008E-03,8.336898877431193E-04,-2.220848255804162E-04])*365

mass_sun     = 2e30
p0_sun       = np.array([0,0,0])
v0_sun       = -(mass_saturn*v0_saturn+mass_mercur*v0_mercur+mass_earth*v0_earth+mass_jupiter*v0_jupiter)/mass_sun


earth   = mpc.body(mass_earth,p0_earth,v0_earth)
jupiter = mpc.body(mass_jupiter,p0_jupiter,v0_jupiter)
sun     = mpc.body(mass_sun,p0_sun,v0_sun)
mercur  = mpc.body(mass_mercur,p0_mercur,v0_mercur)
saturn  = mpc.body(mass_saturn,p0_saturn,v0_saturn)

system  = mpc.system(sun,mercur,earth,jupiter,saturn)

#system   = mpc.system(sun,mercur)#,earth)




forsøk  = mpc.solve(system,1000,1/365)
#print(system[1].pos_vec[0])
forsøk.verlet()
#print(system[0].pos_vec[1]+system[0].vel_vec[1]/365 + 0.5*system[0].acc_vec[1]/(356**2))

import matplotlib.pyplot as plt

for i in range(len(system)):
    plt.plot(system[i].pos_vec[:,0],system[i].pos_vec[:,1])
plt.show()
