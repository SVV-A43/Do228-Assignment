


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 21:36:40 2020

@author: Axel
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

N_z = 81
N_x = 41
C_a = 0.515
L_a = 2.691
th_z = []
th_x = []
z =[]
x = []

P = []

datastring = open('aero_loading_data/aerodynamicloaddo228.dat')

df= pd.DataFrame(datastring)
df.columns=['Data']
list_number= range(1,81)
result=[]
idx=0
for i in df['Data']:
    i = i.strip()
    row = i.split(',')
    
    for j in row:
            p_value = float(j)
            #print(p_value)
            P.append(p_value)


for i in range(1,N_z+2):
    th_z_i = (i-1)/N_z*np.pi
    th_z.append(th_z_i)

for i in range(1,N_z+1):
    z_i = -1/2*((C_a/2*(1-np.cos(th_z[i-1])))+(C_a/2*(1-np.cos(th_z[i]))))
    z.append(z_i)

for i in range(1,N_x+2):
    th_x_i = (i-1)/N_x*np.pi
    th_x.append(th_x_i)

for i in range(1,N_x+1):
   x_i = 1/2*((L_a/2*(1-np.cos(th_x[i-1])))+(L_a/2*(1-np.cos(th_x[i]))))
   x.append(x_i)


'''
for i in z:
    for j in x:
        zx_mesh = plt.plot(i,j,'r.', markersize=0.1)
plt.ylabel('x[m]')
plt.xlabel('z[m]')

'''
ax = plt.axes(projection='3d')


zlist = []
xlist = []

        
for q in range(1,len(z)+1):
    for r in range(1,len(x)+1):        
        zlist.append(z[q-1])
        
for k in range(1,len(z)+1):
    for l in range(1,len(x)+1):
            xlist.append(x[l-1])

P_arr=np.array(P)
x_arr=np.array(xlist)
z_arr=np.array(zlist)

ax = plt.axes(projection='3d')
ax.scatter3D(x_arr, z_arr, P_arr, c=P_arr, cmap='Greens')
#ax.plot_trisurf(x_arr, z_arr, P_arr, cmap='viridis', edgecolor='none')
ax.set_xlabel('Spanwise coordinate x [m]')
ax.set_ylabel('Chordwise coordinate z [m] ')
ax.set_zlabel('Aerodynamic Load p [kN/m^2]')
ax.view_init(10, 45)
ax.set_title('Aerodynamic loading')
plt.show()




