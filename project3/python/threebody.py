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

mass_venus   = 4.9E24
p0_venus     =  np.array([6.771183423096696E-01,2.635570892119448E-01,-3.564015770708658E-02])
v0_venus     =  np.array([-7.233924240439758E-03,1.883053303306194E-02,6.755080544414482E-04])*365

mass_mars  = 6.6E23
p0_mars    = np.array([1.386064640638050E+00,-7.195522876836125E-02,-3.574720875419098E-02])
v0_mars    = np.array([1.325375632613797E-03,1.516910673406354E-02,2.852845787253853E-04])*365

mass_uranus   = 8.8E25
p0_uranus     = np.array([1.716586975644875E+01,1.001373868288497E+01,-1.851949198554858E-01])
v0_uranus     = np.array([-2.010719471096281E-03,3.213990583409865E-03,3.805302445255308E-05])*365

mass_neptun   = 1.03E26
p0_neptun     = np.array([2.892424054587586E+01,-7.707281387885470E+00,-5.078715497631403E-01])
v0_neptun     = np.array([7.870585771863339E-04,3.052076356482153E-03,-8.060978991971341E-05])*365

mass_pluto   = 1.31E22
p0_pluto     = np.array([1.165895866372055E+01,-3.157299883744404E+01,6.054513403127080E-03])
v0_pluto     = np.array([3.023518143470085E-03,4.354001925477489E-04,-9.190982697153275E-04])*365


mass_sun     = 2e30
p0_sun       = np.array([0,0,0])
v0_sun       = -(mass_saturn*v0_saturn+mass_mercur*v0_mercur+mass_earth*\
                 v0_earth+mass_jupiter*v0_jupiter+v0_neptun*mass_neptun+\
                 v0_uranus*mass_uranus+v0_mars*mass_mars+v0_venus*mass_venus+\
                 p0_pluto*mass_pluto)/mass_sun


earth   = mpc.body(mass_earth,p0_earth,v0_earth)
jupiter = mpc.body(mass_jupiter,p0_jupiter,v0_jupiter)
sun     = mpc.body(mass_sun,p0_sun,v0_sun)
mercur  = mpc.body(mass_mercur,p0_mercur,v0_mercur)
saturn  = mpc.body(mass_saturn,p0_saturn,v0_saturn)
venus  = mpc.body(mass_venus,p0_venus,v0_venus)
mars    = mpc.body(mass_mars,p0_mars,v0_mars)
pluto  = mpc.body(mass_pluto,p0_pluto,v0_pluto)
neptun  = mpc.body(mass_neptun,p0_neptun,v0_neptun)
uranus  = mpc.body(mass_uranus,p0_uranus,v0_uranus)

system  = mpc.system(sun,mercur,earth,jupiter,saturn,venus,mars,pluto,uranus,neptun)

#system   = mpc.system(sun,mercur)#,earth)

system  = mpc.system(sun,earth,jupiter)



forsøk  = mpc.solve(system,30000,1/100)
#print(system[1].pos_vec[0])
forsøk.verlet()
#print(system[0].pos_vec[1]+system[0].vel_vec[1]/365 + 0.5*system[0].acc_vec[1]/(356**2))
"""
import matplotlib.pyplot as plt

for i in range(len(system)):
    plt.plot(system[i].pos_vec[:,0],system[i].pos_vec[:,1])
plt.show()
"""
#%%
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')

for i in range(len(system)):
    ax.plot(system[i].pos_vec[:,0],system[i].pos_vec[:,1],system[i].pos_vec[:,2])



plt.show()