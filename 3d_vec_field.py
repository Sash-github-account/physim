#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 20:27:43 2021

@author: sashwathnalinkanth
"""
#from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation


fig = plt.figure()
ax = fig.gca(projection='3d')

x, y, z = np.meshgrid(np.arange(-1, 1, 0.5),
                      np.arange(-1, 1, 0.5),
                      np.arange(-1, 1, 0.5))

for x_cor in range(len(x[0][0])):
    for y_cor in range(len(x[0])):
        for z_cor in range(len(x)):
            print(x[x_cor][y_cor][z_cor], y[x_cor][y_cor][z_cor], z[x_cor][y_cor][z_cor])
                       
                         

u = np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
v = -np.cos(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)
w = (np.sqrt(2.0 / 3.0) * np.cos(np.pi * x) * np.cos(np.pi * y) *
     np.sin(np.pi * z))

ax.quiver(x, y, z, u, v, w, length=0.1, color = 'black')

plt.show()



fig, ax = plt.subplots(subplot_kw=dict(projection="3d"))

def get_arrow(theta):
    x, y, z = np.meshgrid(np.arange(-1, 1, 0.5),
                      np.arange(-1, 1, 0.5),
                      np.arange(-1, 1, 0.5))
    u = np.sin(2*theta) * np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    v = np.sin(3*theta) * -np.cos(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)
    w = np.cos(3*theta) * (np.sqrt(2.0 / 3.0) * np.cos(np.pi * x) * np.cos(np.pi * y) *
     np.sin(np.pi * z))
    return x,y,z,u,v,w

quiver = ax.quiver(*get_arrow(0))

ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_zlim(-2, 2)

def update(theta):
    global quiver
    quiver.remove()
    quiver = ax.quiver(*get_arrow(theta))

ani = FuncAnimation(fig, update, frames=np.linspace(0,2*np.pi,200), interval=50)
plt.show()