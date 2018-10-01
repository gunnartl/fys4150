#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

r = np.linspace(0,5,350)
mat = np.loadtxt("eigenvectors0-2.txt")
import matplotlib.pyplot as plt
plt.plot(r,mat[:,0]**2)
plt.plot(r,mat[:,1]**2)
plt.plot(r,mat[:,2]**2)
plt.xlabel("r")
plt.legend(["1","2","3"])
plt.ylabel("Radial probability")
plt.title("Radial probability for the three lowest energy states")
plt.show()